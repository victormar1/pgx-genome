#!/bin/bash
# Module pharmacogenetique pour genome entier. Un echantillon, une chaine d'etages.
#
#   pgx_genome.sh --cram X.cram --vcf X.vcf.gz --sortie DOSSIER [options]
#
# Le module ne modifie rien en amont : il consomme l'alignement et le fichier de
# variants deja produits par la plateforme, et rend un fichier structure, un
# compte rendu, et le perimetre reellement mesure.
#
# Regle de conception : aucun resultat faux ne doit pouvoir sortir sans etre
# signale. Chaque etage teste son code de retour, inscrit son etat, et refuse de
# se rabattre sur une entree de secours. Le code de retour du module ne vaut zero
# que si tous les etages requis ont abouti. Le typage complementaire est le seul
# etage optionnel : non arme, il s'inscrit IGNORE et ne pese pas sur le resultat.
#
# Options :
#   --echantillon ID   identifiant du prelevement    (defaut : nom de l'alignement)
#   --fasta CHEMIN     reference d'alignement, requise pour un CRAM
#   --identite FICHIER JSON {patient, naissance, contexte} pour le compte rendu ;
#                      fourni, il leve la mention prototype. Un fichier et non
#                      des arguments : la table des processus est lisible par tous
#   --contexte TEXTE   preindication portee sur le compte rendu
#   --gq N             seuil de qualite de genotype           (defaut 20)
#   --profondeur N     seuil de profondeur par position       (defaut 10)
#   --couverture-min N couverture mediane minimale du perimetre (defaut 18,
#                      0 pour accepter un genome hors du domaine valide)
#   --melange-max N    part d'heterozygotes desequilibres au-dela de laquelle
#                      l'echantillon est tenu pour melange     (defaut 35, 0 leve)
#   --ignorer-filtre   retient les sites que l'appelant a rejetes (FILTER)
#   --fils N           fils pour samtools, Cyrius et OptiType (defaut 4)
#   --reprise          reutilise les intermediaires si les entrees sont identiques
#   --forcer           ecrit dans un dossier non vide sans reprise
#   --etage LISTE      ne fait tourner que ces etages, separes par des virgules :
#                      recevabilite filtre tranche qc cyp2d6 hla pypgx appels
#                      pharmcat rendu provenance. La recevabilite ouvre le
#                      dossier ; les autres le reprennent la ou il en est.
#   --version          version du module
#
# Variables : PGX_FASTA PGX_CYRIUS PGX_GQ PGX_PROFONDEUR PGX_COUVERTURE_MIN
#             PGX_MELANGE_MAX PGX_IGNORER_FILTRE PGX_FILS PGX_MOTEUR
#             PGX_IMG_BCFTOOLS PGX_IMG_PHARMCAT PGX_IMG_OPTITYPE PGX_PYPGX

set -uo pipefail

RACINE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RES="$RACINE/ressources"

CRAM=""; VCF=""; SORTIE=""; ECH=""; REPRISE=0; FORCER=0; IDENTITE=""; CONTEXTE=""
ETAGES=""; RECEVABLE=0
FASTA="${PGX_FASTA:-}"
GQ="${PGX_GQ:-20}"
PROF="${PGX_PROFONDEUR:-10}"
# Couverture mediane sous laquelle le genome sort du domaine de validite.
# Dix-huit est le point le plus bas reellement mesure sur le perimetre clinique :
# le seuil n'extrapole pas sous lui. Les neuf executions refusees par la mesure
# plafonnent a treize, la marge est donc de cinq fois.
COUV_MIN="${PGX_COUVERTURE_MIN:-18}"
# Part d'heterozygotes desequilibres au-dela de laquelle l'echantillon n'en est
# pas un seul. Trente-cinq laisse passer les deux cent vingt-neuf echantillons
# purs mesures, dont le plus haut est a 31,9, et arrete un melange a trente pour
# cent, mesure a 41,5. Il ne voit pas un melange a vingt pour cent : voir
# doc/VALIDATION.md.
MELANGE_MAX="${PGX_MELANGE_MAX:-35}"
IGNORER_FILTRE="${PGX_IGNORER_FILTRE:+--ignorer-filtre}"
FILS="${PGX_FILS:-4}"
# Complement RNPGx (classes 1 et 2 hors definitions PharmCAT) : present, il est
# mesure ; absent, le module se comporte comme avant.
COMPLEMENT=""
[ -s "$RES/rnpgx_complement.vcf" ] && COMPLEMENT="$RES/rnpgx_complement.vcf"
# Typage complementaire (MT-RNR1, BCHE, MTHFR, POR) : quatre genes du core panel
# que l'interpreteur ne sait pas appeler. L'etage ne s'arme que si le catalogue
# et le typeur sont tous deux presents ; sinon il est ignore, sans consequence
# sur le reste. Le typeur est une image (PGX_IMG_PYPGX), lancee comme les trois
# autres par le moteur de conteneurs, ou une commande de l'hote (PGX_PYPGX).
PYPGX="${PGX_PYPGX:-}"
IMG_PYPGX="${PGX_IMG_PYPGX:-}"
CAT_PYPGX=""
{ [ -n "$PYPGX" ] || [ -n "$IMG_PYPGX" ]; } && [ -s "$RES/pypgx_genes.json" ] \
  && CAT_PYPGX="$RES/pypgx_genes.json"
IMG_BCF="${PGX_IMG_BCFTOOLS:-quay.io/biocontainers/bcftools:1.24--h118bc1c_2}"
IMG_PHARMCAT="${PGX_IMG_PHARMCAT:-pgkb/pharmcat:3.4.0}"
IMG_OPTITYPE="${PGX_IMG_OPTITYPE:-quay.io/biocontainers/optitype:1.3.5--hdfd78af_3}"
CYRIUS="${PGX_CYRIUS:-}"

while [ $# -gt 0 ]; do
  case "$1" in
    --cram) CRAM="$2"; shift 2;;
    --vcf) VCF="$2"; shift 2;;
    --sortie) SORTIE="$2"; shift 2;;
    --echantillon) ECH="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --gq) GQ="$2"; shift 2;;
    --profondeur) PROF="$2"; shift 2;;
    --couverture-min) COUV_MIN="$2"; shift 2;;
    --melange-max) MELANGE_MAX="$2"; shift 2;;
    --fils) FILS="$2"; shift 2;;
    --identite) IDENTITE="$2"; shift 2;;
    --contexte) CONTEXTE="$2"; shift 2;;
    --reprise) REPRISE=1; shift;;
    --forcer) FORCER=1; shift;;
    --ignorer-filtre) IGNORER_FILTRE="--ignorer-filtre"; shift;;
    --etage) ETAGES="${ETAGES:+$ETAGES,}$2"; shift 2;;
    -h|--help) sed -n '2,42p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    --version) cat "$RACINE/VERSION" 2>/dev/null || echo inconnue; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done

# Un etage mal nomme ne doit pas devenir un etage saute en silence : la liste
# se verifie avant tout le reste, et un nom inconnu est un refus.
CONNUS="recevabilite filtre tranche qc cyp2d6 hla pypgx appels pharmcat rendu provenance"
DECOUPE=0
if [ -n "$ETAGES" ]; then
  DECOUPE=1
  for e in ${ETAGES//,/ }; do
    case " $CONNUS " in
      *" $e "*) ;;
      *) echo "etage inconnu : $e" >&2
         echo "attendus : $CONNUS" >&2
         exit 2;;
    esac
  done
fi
a_faire() { [ "$DECOUPE" = 0 ] || [[ ",$ETAGES," == *",$1,"* ]]; }

[ -n "$CRAM" ] && [ -n "$VCF" ] && [ -n "$SORTIE" ] || {
  echo "usage : pgx_genome.sh --cram X.cram --vcf X.vcf.gz --sortie DOSSIER" >&2; exit 2; }
[ -s "$CRAM" ] || { echo "alignement introuvable : $CRAM" >&2; exit 2; }
[ -s "$VCF" ]  || { echo "fichier de variants introuvable : $VCF" >&2; exit 2; }
[ -n "$ECH" ] || ECH="$(basename "$CRAM" | sed 's/\.\(cram\|bam\)$//')"
# Un seuil non numerique serait lu zero par awk : la porte tomberait sans bruit.
case "$COUV_MIN" in
  ''|*[!0-9]*) echo "couverture minimale invalide : $COUV_MIN" >&2; exit 2;;
esac
case "$MELANGE_MAX" in
  ''|*[!0-9.]*) echo "part de melange invalide : $MELANGE_MAX" >&2; exit 2;;
esac

mkdir -p "$SORTIE"/{travail,sortie}
JOURNAL="$SORTIE/journal.txt"
ETATS="$SORTIE/travail/etats.tsv"

# --------------------------------------------------------------------- verrou
# Deux executions sur le meme dossier se detruisent mutuellement : elles purgent
# le dossier de travail l'une sous l'autre et melangent leurs lignes d'etat. Le
# verrou est pris par mkdir, qui est atomique.
#
# Une execution complete tient le verrou du dossier ; un etage lance seul tient
# le sien. Deux instances d'un meme etage s'excluent, des etages differents
# peuvent tourner ensemble, et rien ne tourne pendant que la recevabilite ou une
# execution complete purge le dossier de travail.
VERROUS=()
HOTE=$(hostname 2>/dev/null || echo inconnu)
vivant() {   # $1 dossier de verrou : tenu par un autre processus vivant ?
  local autre
  [ -d "$1" ] || return 1
  # Un verrou pose depuis un autre hote ne se juge pas d'ici : le numero de
  # processus n'y veut rien dire, et sur un dossier partage entre noeuds, le
  # reprendre ecraserait le travail d'une tache vivante ailleurs. Il est tenu
  # pour vivant ; apres une panne de noeud, l'exploitant le retire a la main.
  [ -s "$1/hote" ] && [ "$(cat "$1/hote")" != "$HOTE" ] && return 0
  autre=$(cat "$1/pid" 2>/dev/null)
  [ -n "$autre" ] && [ "$autre" != "$$" ] && kill -0 "$autre" 2>/dev/null
}
prendre() {  # $1 dossier de verrou
  if ! mkdir "$1" 2>/dev/null; then
    vivant "$1" && return 1
    # le detenteur precedent est mort : on reprend le verrou
    rm -rf "$1" && mkdir "$1" 2>/dev/null || return 1
  fi
  echo "$$" > "$1/pid"
  echo "$HOTE" > "$1/hote"
  date -Iseconds > "$1/depuis"
  VERROUS+=("$1")
}
nettoyer() { local v; for v in "${VERROUS[@]}"; do rm -rf "$v"; done; }
trap nettoyer EXIT INT TERM
refus_verrou() {  # $1 message
  echo "$1" >&2
  exit 3
}
if [ "$DECOUPE" = 0 ]; then
  prendre "$SORTIE/.verrou" || refus_verrou \
    "une autre execution travaille deja dans $SORTIE (processus $(cat "$SORTIE/.verrou/pid" 2>/dev/null) sur $(cat "$SORTIE/.verrou/hote" 2>/dev/null))"
  for v in "$SORTIE"/.verrou.*; do
    vivant "$v" && refus_verrou "l'etage ${v##*.verrou.} travaille deja dans $SORTIE"
  done
else
  vivant "$SORTIE/.verrou" && refus_verrou \
    "une execution complete travaille deja dans $SORTIE (processus $(cat "$SORTIE/.verrou/pid" 2>/dev/null))"
  for e in ${ETAGES//,/ }; do
    prendre "$SORTIE/.verrou.$e" || refus_verrou "l'etage $e travaille deja dans $SORTIE"
  done
  # Verrous pris d'abord, autres verifies ensuite : deux processus qui se
  # croisent se voient l'un l'autre et refusent tous deux, sans rien abimer.
  if a_faire recevabilite; then
    for v in "$SORTIE"/.verrou.*; do
      vivant "$v" && refus_verrou \
        "l'etage ${v##*.verrou.} travaille dans $SORTIE : la recevabilite purgerait son travail"
    done
  else
    vivant "$SORTIE/.verrou.recevabilite" && refus_verrou \
      "la recevabilite est en cours dans $SORTIE"
  fi
fi

# Des etages lances en parallele ecrivent dans la meme table d'etats. Un ajout
# n'est pas atomique sur un systeme de fichiers en reseau : deux lignes peuvent
# s'y entrelacer, et la trace lirait un etat qui n'existe pas. L'ajout passe
# donc par un verrou mkdir, atomique partout.
inscrire() {  # $1 fichier  $2 ligne
  local v="$1.ecriture" i=0 tenu=0
  while [ $i -lt 200 ]; do
    mkdir "$v" 2>/dev/null && { tenu=1; break; }
    i=$((i+1)); sleep 0.05
  done
  printf '%s\n' "$2" >> "$1"
  [ "$tenu" = 1 ] && rmdir "$v" 2>/dev/null
  return 0
}
dire() {
  local l
  l=$(printf '[%s] %s' "$(date +%H:%M:%S)" "$*")
  printf '%s\n' "$l"
  inscrire "$JOURNAL" "$l"
}
etat() { inscrire "$ETATS" "$(printf '%s\t%s\t%s\t%s' "$1" "$2" "$3" "$4")"; }

D_CRAM="$(cd "$(dirname "$CRAM")" && pwd)"
D_VCF="$(cd "$(dirname "$VCF")" && pwd)"
D_SORTIE="$(cd "$SORTIE" && pwd)"
DOSSIERS=("$D_CRAM" "$D_VCF" "$RES")
ECRIVABLES=("$D_SORTIE")
if [ -n "$FASTA" ]; then
  D_FASTA="$(cd "$(dirname "$FASTA")" && pwd)"
  DOSSIERS+=("$D_FASTA")
fi

# Le moteur de conteneurs est interchangeable. Docker reclame un demon
# privilegie, que beaucoup de plateformes de calcul n'autorisent pas ; elles
# disposent d'Apptainer, qui s'execute sans privilege. Les deux montent les
# memes dossiers, avec une syntaxe differente — et Apptainer monte le dossier
# courant de l'hote par defaut, ce qu'on desactive pour que l'execution ne
# depende pas du repertoire d'appel.
MOTEUR="${PGX_MOTEUR:-docker}"
# La commande qui lance une image, sans la commande interne : les etages
# l'executent, et le typeur complementaire la recoit comme prefixe. Un dossier
# en lecture seule et une variable d'environnement peuvent s'y ajouter.
CONT=()
args_conteneur() {   # $1 dossier de travail  $2 image  [$3 dossier en lecture seule]  [$4 VAR=valeur]
  local cwd="$1" image="$2" extra="${3:-}" env="${4:-}" d
  CONT=()
  case "$MOTEUR" in
    docker)
      if [ -f "$image" ]; then
        echo "pgx_genome : $image est un fichier ; Docker nomme une image, un fichier .sif demande PGX_MOTEUR=apptainer" >&2
        return 127
      fi
      # Sous l'identite de l'appelant, et non en root : un demon classique
      # ecrirait les sorties au nom de root, et le module ne pourrait plus
      # purger son propre dossier de travail. Le dossier courant tient lieu de
      # dossier personnel, seul endroit garanti en ecriture a cette identite.
      CONT=(docker run --rm --user "$(id -u):$(id -g)" -e "HOME=$cwd")
      for d in "${DOSSIERS[@]}" ${extra:+"$extra"}; do CONT+=(-v "$d":"$d":ro); done
      for d in "${ECRIVABLES[@]}"; do CONT+=(-v "$d":"$d"); done
      [ -n "$env" ] && CONT+=(-e "$env")
      CONT+=(-w "$cwd" "$image")
      ;;
    apptainer|singularity)
      # Un fichier .sif local est employe tel quel ; un nom d'image est tire
      # du depot, ce qui demande le reseau au moment de l'execution.
      CONT=("$MOTEUR" exec --cleanenv --no-home)
      for d in "${DOSSIERS[@]}" ${extra:+"$extra"}; do CONT+=(--bind "$d":"$d":ro); done
      for d in "${ECRIVABLES[@]}"; do CONT+=(--bind "$d":"$d"); done
      [ -n "$env" ] && CONT+=(--env "$env")
      if [ -f "$image" ]; then CONT+=(--pwd "$cwd" "$image")
      else CONT+=(--pwd "$cwd" "docker://$image"); fi
      ;;
    *)
      echo "moteur de conteneurs inconnu : $MOTEUR" >&2; return 127
      ;;
  esac
}
conteneur() {   # $1 dossier de travail  $2 image  $3... commande
  local cwd="$1" image="$2"; shift 2
  args_conteneur "$cwd" "$image" || return 127
  "${CONT[@]}" "$@"
}
BCF() { conteneur "$D_SORTIE" "$IMG_BCF" bcftools "$@"; }

# Le moteur se verifie ici, avant tout etage. Sans ce controle, un moteur
# inconnu ou absent faisait echouer la recevabilite sur « le fichier de
# variants porte 0 echantillons » : un message qui envoie chercher la faute
# dans les donnees du patient alors qu'elle est dans la configuration.
case "$MOTEUR" in
  docker|apptainer|singularity) ;;
  *) echo "pgx_genome : moteur de conteneurs inconnu : $MOTEUR" >&2
     echo "             attendus : docker, apptainer, singularity" >&2
     exit 2;;
esac
if ! command -v "$MOTEUR" >/dev/null 2>&1; then
  echo "pgx_genome : $MOTEUR demande par PGX_MOTEUR mais introuvable" >&2
  exit 2
fi

T="$SORTIE/travail"

# Le typeur complementaire en image : sa commande est celle de son conteneur,
# avec sa ressource de donnees montee. Elle n'est pas dans l'image, parce
# qu'elle pese plusieurs gigaoctets et evolue a part ; sans elle, le typeur
# echouerait au premier gene, et le motif serait moins clair qu'ici.
if [ -n "$IMG_PYPGX" ]; then
  if [ -z "${PYPGX_BUNDLE:-}" ] || [ ! -d "$PYPGX_BUNDLE" ]; then
    echo "pgx_genome : PGX_IMG_PYPGX demande PYPGX_BUNDLE, le dossier de donnees du typeur" >&2
    exit 2
  fi
  args_conteneur "$T" "$IMG_PYPGX" "$PYPGX_BUNDLE" "PYPGX_BUNDLE=$PYPGX_BUNDLE"
  PYPGX="${CONT[*]} pypgx"
fi
IMAGES="$IMG_BCF,$IMG_PHARMCAT,$IMG_OPTITYPE${IMG_PYPGX:+,$IMG_PYPGX}"

# Les chemins que les etages se transmettent sont fixes, et poses ici : un etage
# lance seul doit les connaitre sans que l'etage qui les produit ait tourne dans
# le meme processus.
FIL="$T/filtre.vcf.gz"
TRANCHE="$T/tranche.bam"
BEDT="$T/tranche.bed"
DEPTH="$T/profondeur.txt"
QUAL="$T/qualifie.vcf.gz"
D6="$T/cyp2d6.tsv"
PGV="$T/pypgx_entree.vcf.gz"
PO="$T/appels_externes.tsv"
PRE="$T/$(basename "$QUAL" .gz).preprocessed.vcf.bgz"
RAP="$SORTIE/sortie/$ECH.report.json"
ARRET="$T/arret.txt"

# Un alignement sur un autre assemblage rend des resultats aux mauvaises
# coordonnees sans lever d'erreur : les regions visees n'y existent pas, la
# tranche revient vide, et tout ce qui suit conclut a l'absence de couverture.
# La verification ne nomme aucun assemblage : elle confronte les contigs vises
# et leurs longueurs a ceux de la reference fournie.
ecarts_assemblage() {
  local entete="$T/entete_alignement.txt" fai=""
  [ -n "$FASTA" ] && [ -s "$FASTA.fai" ] && fai="$FASTA.fai"
  samtools view -H ${FASTA:+--reference "$FASTA"} "$CRAM" 2>/dev/null \
    | grep '^@SQ' > "$entete"
  [ -s "$entete" ] || { printf 'en-tete illisible'; return 0; }
  awk -v FS='\t' -v fai="$fai" '
    FILENAME == ARGV[1] {
      sn = ""; ln = ""
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^SN:/) sn = substr($i, 4)
        if ($i ~ /^LN:/) ln = substr($i, 4)
      }
      if (sn != "") L[sn] = ln
      next
    }
    fai != "" && FILENAME == fai { R[$1] = $2; next }
    { C[$1] = 1 }
    END {
      n = 0
      for (c in C) {
        if (n >= 3) break
        if (!(c in L)) { printf "%s absent de l alignement ; ", c; n++ }
        else if ((c in R) && L[c] + 0 != R[c] + 0) {
          printf "%s mesure %s dans l alignement et %s dans la reference ; ",
                 c, L[c], R[c]
          n++
        }
      }
    }' "$entete" ${fai:+"$fai"} "$RES/pharmcat_positions.bed" 2>/dev/null
}
# Le verdict de l'etage de controle se relit dans perimetre.json : un etage
# lance seul n'a pas en memoire ce que l'etage 2b a mesure. Quatre champs,
# jamais vides : un champ absent decalerait les suivants.
lire_verdict() {
  local v
  HORS_DOMAINE=""; MELANGE=""; COUV=""; MEL_PART=""
  [ -s "$T/perimetre.json" ] || return 0
  v=$(python3 -c "import json,sys;sys.path.insert(0,r'$RACINE/bin');import qc_perimetre as q;d=json.load(open(r'$T/perimetre.json',encoding='utf-8'));c=q.couverture_retenue(d);m=d.get('melange_part');print('-' if c is None else c, 1 if d.get('couverture_hors_domaine') else 0, '-' if m is None else m, 1 if d.get('melange_suspect') else 0)" 2>>"$JOURNAL")
  set -- $v
  COUV="${1:--}"; MEL_PART="${3:--}"
  [ "$COUV" = "-" ] && COUV=""
  [ "$MEL_PART" = "-" ] && MEL_PART=""
  [ "${2:-0}" = 1 ] && HORS_DOMAINE="couverture mediane ${COUV}x, domaine valide a partir de ${COUV_MIN}x"
  [ "${4:-0}" = 1 ] && MELANGE="part d'heterozygotes desequilibres ${MEL_PART}%, au-dela de ${MELANGE_MAX}% : echantillon tenu pour melange"
  return 0
}
REGION_TEST=$(awk 'NR==1{printf "%s:%d-%d", $1, $2+1, $3}' "$RES/pharmcat_positions.bed")

# --------------------------------------------------------------- 0. recevabilite
# Les intermediaires ne sont reutilisables que s'ils viennent des memes entrees.
# Sans cette empreinte, une reprise apres correction du manifeste melangerait
# deux echantillons sans que rien ne le signale.
empreinte() {
  # taille et date du fichier, plus le debut de son contenu : suffisant pour
  # distinguer deux entrees, sans lire trente gigaoctets
  local f="$1"
  [ -s "$f" ] || { echo "absent"; return; }
  printf '%s %s %s' "$(stat -c '%s' "$f")" "$(stat -c '%Y' "$f")" \
    "$(head -c 1048576 "$f" | sha256sum | cut -c1-32)"
}
# La reference fait partie des entrees : la tranche, la profondeur et tout ce
# qui en decoule changent avec elle. Son absence de la signature laissait une
# reprise reutiliser le travail de la veille apres un changement de reference,
# et un alignement decode avec la mauvaise reference ne leve rien. Les
# ressources y entrent aussi : une position ajoutee au perimetre change la
# tranche sans toucher aux entrees. Le seuil d'equilibre allelique n'y
# figure pas : il n'est pas reglable en ligne de commande, et l'etage de
# controle qualite qui le porte est rejoue a chaque execution.
# Le contenu, et non le nom ni la taille : une position corrigee a taille
# egale doit invalider la reprise. Cent soixante kilooctets, cent cinquante
# millisecondes.
EMPREINTE_RES=$(find "$RES" -type f -print0 2>/dev/null | sort -z \
  | xargs -0 sha256sum 2>/dev/null | sha256sum | cut -c1-16)
SIGNATURE="$(empreinte "$CRAM")|$(empreinte "$VCF")|$(empreinte "$FASTA")|gq=$GQ|prof=$PROF|couv=$COUV_MIN|mel=$MELANGE_MAX|filtre=${IGNORER_FILTRE:-respecte}|res=$EMPREINTE_RES"
SIG_FICHIER="$T/signature.txt"
# La recevabilite ouvre le dossier : elle seule le purge et remet la trace a
# zero. Un etage lance seul le reprend la ou il en est, a condition qu'il porte
# la recevabilite de ces memes entrees : sinon il travaillerait sur les
# intermediaires d'un autre echantillon sans que rien ne le signale.
if a_faire recevabilite; then
  if [ -s "$SIG_FICHIER" ]; then
    if [ "$(cat "$SIG_FICHIER")" = "$SIGNATURE" ] && [ "$REPRISE" = 1 ]; then
      : # meme entree, reprise demandee : on garde les intermediaires
    else
      rm -rf "$T"; mkdir -p "$T"
    fi
  elif [ -n "$(ls -A "$T" 2>/dev/null)" ] && [ "$FORCER" = 0 ]; then
    rm -rf "$T"; mkdir -p "$T"
  fi
  rm -f "$ARRET"
  printf '%s' "$SIGNATURE" > "$SIG_FICHIER"
  : > "$JOURNAL"
  printf 'etage\tetat\tsecondes\tdetail\n' > "$ETATS"

  dire "echantillon $ECH"
  dire "alignement  $CRAM"
  dire "variants    $VCF"
else
  if [ ! -s "$SIG_FICHIER" ] || [ "$(cat "$SIG_FICHIER")" != "$SIGNATURE" ]; then
    echo "le dossier $SORTIE ne porte pas la recevabilite de ces entrees :" >&2
    echo "lancer d'abord --etage recevabilite avec les memes entrees et options" >&2
    exit 2
  fi
  # Un refus, a l'entree ou apres le controle qualite, vaut pour tous les
  # etages qui suivent : aucun ne doit produire de resultat pour un genome refuse.
  if [ -s "$ARRET" ]; then
    dire "etage(s) $ETAGES : arret, $(cat "$ARRET")"
    exit 2
  fi
  if ! grep -q "^recevabilite	OK" "$ETATS" 2>/dev/null; then
    echo "la recevabilite n'a pas abouti dans $SORTIE" >&2
    exit 2
  fi
  RECEVABLE=1
  dire "etage(s) $ETAGES"
fi
# Les lignes d'etat ecrites par cette execution, et elles seules, decident de
# son code de retour quand elle ne fait tourner qu'une partie des etages.
DEBUT_ETATS=$(wc -l < "$ETATS")

if a_faire recevabilite; then
t0=$(date +%s)
NOMS_VCF=$(BCF query -l "$D_VCF/$(basename "$VCF")" 2>>"$JOURNAL")
NB_VCF=$(printf '%s\n' "$NOMS_VCF" | grep -c . || true)
SM_CRAM=$(samtools view -H ${FASTA:+--reference "$FASTA"} "$CRAM" 2>>"$JOURNAL" \
          | tr '\t' '\n' | grep -oE '^SM:.+' | sed 's/SM://' | sed 's/\\t.*//' \
          | sort -u | head -3 | tr '\n' ',' | sed 's/,$//')
if [ "$NB_VCF" != "1" ]; then
  dire "0. recevabilite : ECHEC, le fichier de variants porte $NB_VCF echantillons"
  etat recevabilite ECHEC $(( $(date +%s)-t0 )) "$NB_VCF echantillons dans le VCF"
  RECEVABLE=0
elif [ -n "$SM_CRAM" ] && [ "$SM_CRAM" != "$NOMS_VCF" ] && [ "$SM_CRAM" != "$ECH" ] && [ "$NOMS_VCF" != "$ECH" ]; then
  dire "0. recevabilite : ECHEC, identites divergentes (VCF $NOMS_VCF, alignement $SM_CRAM, demande $ECH)"
  etat recevabilite ECHEC $(( $(date +%s)-t0 )) "VCF=$NOMS_VCF alignement=$SM_CRAM demande=$ECH"
  RECEVABLE=0
elif ! samtools quickcheck "$CRAM" 2>>"$JOURNAL"; then
  dire "0. recevabilite : ECHEC, alignement tronque ou illisible"
  etat recevabilite ECHEC $(( $(date +%s)-t0 )) "samtools quickcheck : alignement tronque"
  RECEVABLE=0
elif ECARTS=$(ecarts_assemblage); [ -n "$ECARTS" ]; then
  dire "0. recevabilite : ECHEC, assemblage different de la reference : $ECARTS"
  etat recevabilite ECHEC $(( $(date +%s)-t0 )) "assemblage : $ECARTS"
  RECEVABLE=0
elif ! samtools view -c ${FASTA:+--reference "$FASTA"} "$CRAM" "$REGION_TEST" >/dev/null 2>>"$JOURNAL"; then
  dire "0. recevabilite : ECHEC, acces aleatoire impossible (index absent ou illisible)"
  etat recevabilite ECHEC $(( $(date +%s)-t0 )) "index de l alignement absent ou illisible"
  RECEVABLE=0
else
  dire "0. recevabilite : un echantillon, identites concordantes ($NOMS_VCF), assemblage conforme"
  etat recevabilite OK $(( $(date +%s)-t0 )) "VCF=$NOMS_VCF alignement=$SM_CRAM"
  RECEVABLE=1
fi

# Un refus a l'entree arrete tout. La purge du dossier de travail n'est pas une
# precaution de confort : avec --reprise, la tranche d'une execution anterieure
# serait retrouvee par l'etage 2a, qui se declarerait abouti, et le module
# produirait un compte rendu pour un genome qu'il vient de refuser.
if [ "$RECEVABLE" != 1 ]; then
  # On purge les intermediaires, pas la trace : l'etat de recevabilite porte le
  # motif du refus, et c'est la seule chose que le module ait a produire ici.
  find "$T" -mindepth 1 -maxdepth 1        ! -name "$(basename "$ETATS")" ! -name "$(basename "$SIG_FICHIER")"        -exec rm -rf {} + 2>/dev/null
  rm -rf "$SORTIE/sortie"; mkdir -p "$SORTIE/sortie"
  printf '%s' "$SIGNATURE" > "$SIG_FICHIER"
  printf '%s' "entrees non recevables" > "$ARRET"
  for e in filtre tranche qc cyp2d6 hla appels pharmcat rendu; do
    etat "$e" ECHEC 0 "entrees non recevables"
  done
  dire "arret : entrees non recevables, aucun etage n'a demarre"
  python3 "$RACINE/bin/provenance.py" --sortie "$SORTIE" --echantillon "$ECH" \
    --cram "$CRAM" --vcf "$VCF" --fasta "$FASTA" --ressources "$RES" \
    --gq "$GQ" --profondeur "$PROF" \
    --images "$IMAGES" --moteur "$MOTEUR" --cyrius "$CYRIUS" --pypgx "$PYPGX" \
  --filtre "$([ -n "$IGNORER_FILTRE" ] && echo ignore || echo respecte)" \
  --perimetre-clinique "$RES/perimetre_rnpgx.json" >>"$JOURNAL" 2>&1
  dire "termine, reussite complete : False"
  exit 2
fi
fi

# ------------------------------------------------------------------ 1. filtre
# Les regions visees reunissent les definitions PharmCAT et le complement RNPGx
# (classes 1 et 2 absentes de PharmCAT). Le complement est optionnel : absent, le
# module se comporte exactement comme avant.
if a_faire filtre; then
t0=$(date +%s)
BEDR="$T/regions.bed"
cat "$RES/pharmcat_positions.bed" > "$BEDR"
[ -s "$RES/rnpgx_complement.bed" ] && cat "$RES/rnpgx_complement.bed" >> "$BEDR"
sort -k1,1 -k2,2n "$BEDR" -o "$BEDR"
if [ "$RECEVABLE" = 1 ] \
   && BCF view -R "$D_SORTIE/travail/regions.bed" -Oz -o "$FIL" "$D_VCF/$(basename "$VCF")" >>"$JOURNAL" 2>&1 \
   && BCF index -f -t "$FIL" >>"$JOURNAL" 2>&1; then
  n=$(BCF view -H "$FIL" 2>/dev/null | wc -l)
  dire "1. filtre : $n enregistrements dans les regions visees"
  etat filtre OK $(( $(date +%s)-t0 )) "$n enregistrements"
else
  dire "1. filtre : ECHEC"; etat filtre ECHEC $(( $(date +%s)-t0 )) "bcftools view ou index"
fi
fi

# --------------------------------------------------------- 2a. tranche unique
# Une seule traversee de l'alignement sert au controle qualite et au typage HLA.
# Le BED de tranche est construit ici : les contigs HLA alternatifs dependent de
# la reference de la plateforme et ne peuvent pas etre figes dans les ressources.
if a_faire tranche; then
t0=$(date +%s)
if [ ! -s "$TRANCHE" ] && [ "$RECEVABLE" = 1 ]; then
  cat "$RES/pharmcat_positions.bed" > "$BEDT"
  [ -s "$RES/rnpgx_complement.bed" ] && cat "$RES/rnpgx_complement.bed" >> "$BEDT"
  # Les positions definissantes du typage complementaire, et elles seules : la
  # tranche ne porte pas les regions entieres des quatre genes, seulement de
  # quoi prouver que chaque position a bien ete lue.
  [ -n "$CAT_PYPGX" ] && [ -s "$RES/pypgx_tranche.bed" ] && cat "$RES/pypgx_tranche.bed" >> "$BEDT"
  printf 'chr6\t29600000\t33100000\n' >> "$BEDT"
  samtools view -H ${FASTA:+--reference "$FASTA"} "$CRAM" 2>/dev/null \
    | awk '$0 ~ /SN:HLA-/ {sn=""; ln=""; for(i=1;i<=NF;i++){if($i ~ /^SN:/) sn=substr($i,4); if($i ~ /^LN:/) ln=substr($i,4)} if(sn!="" && ln!="") printf "%s\t0\t%s\n", sn, ln}' \
    >> "$BEDT"
  # Une region sur un contig absent de l'alignement ferait echouer l'extraction.
  # On ne garde que les contigs reellement presents : un chrM absent, ou un
  # nommage different sur la plateforme, se traduit alors par une position
  # « profondeur inconnue », jamais par un etage en echec.
  samtools view -H ${FASTA:+--reference "$FASTA"} "$CRAM" 2>/dev/null \
    | awk '/^@SQ/ {for(i=1;i<=NF;i++) if($i ~ /^SN:/) print substr($i,4)}' > "$T/contigs.txt"
  awk 'NR==FNR {ok[$1]=1; next} ok[$1]' "$T/contigs.txt" "$BEDT" > "$BEDT.tmp"
  ecartes=$(( $(grep -c . "$BEDT") - $(grep -c . "$BEDT.tmp") ))
  [ "$ecartes" -gt 0 ] && dire "2a. tranche : $ecartes region(s) ecartee(s), contig absent de l'alignement"
  mv "$BEDT.tmp" "$BEDT"
  sort -k1,1 -k2,2n "$BEDT" -o "$BEDT"
  if samtools view -b -M -L "$BEDT" ${FASTA:+--reference "$FASTA"} -o "$TRANCHE.tmp" "$CRAM" 2>>"$JOURNAL" \
     && samtools sort -@ "$FILS" -o "$TRANCHE" "$TRANCHE.tmp" 2>>"$JOURNAL" \
     && samtools index "$TRANCHE" 2>>"$JOURNAL"; then
    rm -f "$TRANCHE.tmp"
  else
    rm -f "$TRANCHE.tmp" "$TRANCHE"
  fi
fi
if [ -s "$TRANCHE" ] && [ -s "$TRANCHE.bai" ]; then
  nl=$(samtools view -c "$TRANCHE" 2>/dev/null)
  dire "2a. tranche : $(du -h "$TRANCHE" | cut -f1), $nl lectures, $(grep -c . "$BEDT") regions"
  etat tranche OK $(( $(date +%s)-t0 )) "$nl lectures"
else
  dire "2a. tranche : ECHEC"; etat tranche ECHEC $(( $(date +%s)-t0 )) "samtools view, sort ou index"
fi
fi

# ------------------------------------------------------- 2b. controle qualite
if a_faire qc; then
t0=$(date +%s)
PROFOK=0
BEDP="$T/positions_mesurees.bed"
cat "$RES/positions_exactes.bed" > "$BEDP"
[ -s "$RES/rnpgx_complement_positions.bed" ] && cat "$RES/rnpgx_complement_positions.bed" >> "$BEDP"
# La profondeur sur les positions du typage complementaire est ce qui distingue
# un « reference » lu d'un « reference » suppose : sans elle, l'etage 4b refuse
# de rendre le gene.
[ -n "$CAT_PYPGX" ] && [ -s "$RES/pypgx_positions.bed" ] && cat "$RES/pypgx_positions.bed" >> "$BEDP"
sort -k1,1 -k2,2n "$BEDP" -o "$BEDP"
if [ -s "$TRANCHE.bai" ]; then
  if samtools depth -a -b "$BEDP" -Q 0 -q 0 "$TRANCHE" > "$DEPTH" 2>>"$JOURNAL" \
     && [ -s "$DEPTH" ]; then PROFOK=1; else rm -f "$DEPTH"; fi
fi
if [ "$PROFOK" = 0 ]; then
  dire "2b. controle qualite : ECHEC, profondeur non mesuree sur l'alignement"
  etat qc ECHEC $(( $(date +%s)-t0 )) "samtools depth"
elif python3 "$RACINE/bin/qc_perimetre.py" --vcf "$FIL" --profondeur "$DEPTH" \
       --positions "$RES/pharmcat_positions.vcf" --sortie "$T" \
       ${COMPLEMENT:+--positions-complement "$COMPLEMENT"} \
       --gq "$GQ" --profondeur-min "$PROF" --couverture-min "$COUV_MIN" \
       --melange-max "$MELANGE_MAX" \
       --echantillon "$ECH" \
       --perimetre-clinique "$RES/perimetre_rnpgx.json" $IGNORER_FILTRE >>"$JOURNAL" 2>&1 \
     && conteneur "$D_SORTIE" "$IMG_BCF" bgzip -f "$T/qualifie.vcf" >>"$JOURNAL" 2>&1 \
     && BCF index -f -t "$QUAL" >>"$JOURNAL" 2>&1 && [ -s "$QUAL" ]; then
  ret=$(python3 -c "import json;d=json.load(open(r'$T/perimetre.json',encoding='utf-8'));print(d['positions_retenues'],d['positions_visees'],sum(1 for v in d['genes'].values() if v['statut']=='complet'),len(d['genes']),d.get('positions_sans_GQ',0))")
  dire "2b. controle qualite : $(echo $ret | awk '{print $1" positions retenues sur "$2", "$3" genes complets sur "$4", "$5" sans qualite"}')"
  etat qc OK $(( $(date +%s)-t0 )) "$ret"
  # Le verdict sur le domaine est rendu par l'etage lui-meme, ou il se teste.
  lire_verdict
  dire "    couverture mediane du perimetre : ${COUV:-inconnue}x"
  dire "    heterozygotes desequilibres : ${MEL_PART:-part non etablie}%"
  true
else
  dire "2b. controle qualite : ECHEC"
  etat qc ECHEC $(( $(date +%s)-t0 )) "qc_perimetre.py, bgzip ou index"
fi

# Deux motifs d'arret apres le controle qualite, et la meme raison de ne pas se
# contenter d'un avertissement : rien dans le compte rendu ne distinguerait un
# appel fonde sur une position lue d'un appel fonde sur une position supposee,
# ni un heterozygote vrai d'un faux heterozygote venu d'un second individu.
#
# Sous le seuil de couverture, un allele variant est lu comme reference et
# l'interpreteur rend un diplotype faux sans reserve. Au-dela de la part
# d'heterozygotes desequilibres, l'echantillon n'en est pas un seul.
# --couverture-min 0 et --melange-max 0 levent chacun sa porte.
MOTIF="${HORS_DOMAINE:-}${HORS_DOMAINE:+${MELANGE:+ ; }}${MELANGE:-}"
if [ -n "$MOTIF" ]; then
  # On garde la mesure, qui est le motif du refus, et on purge ce qui pourrait
  # etre pris pour un resultat.
  rm -rf "$SORTIE/sortie"; mkdir -p "$SORTIE/sortie"
  printf '%s' "$MOTIF" > "$ARRET"
  m=$([ -n "${HORS_DOMAINE:-}" ] && echo "couverture hors domaine" || echo "echantillon tenu pour melange")
  for e in cyp2d6 hla appels pharmcat rendu; do
    etat "$e" ECHEC 0 "$m"
  done
  [ -n "$CAT_PYPGX" ] && etat pypgx ECHEC 0 "$m" \
    || etat pypgx IGNORE 0 "typeur non defini"
  dire "arret : $MOTIF"
  python3 "$RACINE/bin/provenance.py" --sortie "$SORTIE" --echantillon "$ECH" \
    --cram "$CRAM" --vcf "$VCF" --fasta "$FASTA" --ressources "$RES" \
    --gq "$GQ" --profondeur "$PROF" --couverture-min "$COUV_MIN" \
    --couverture "$COUV" --melange-max "$MELANGE_MAX" --melange "${MEL_PART:-}" \
    --images "$IMAGES" --moteur "$MOTEUR" --cyrius "$CYRIUS" --pypgx "$PYPGX" \
    --filtre "$([ -n "$IGNORER_FILTRE" ] && echo ignore || echo respecte)" \
    --perimetre-clinique "$RES/perimetre_rnpgx.json" >>"$JOURNAL" 2>&1
  dire "termine, reussite complete : False"
  exit 2
fi
fi

# --------------------------------------------------------------- 3. CYP2D6
if a_faire cyp2d6; then
t0=$(date +%s)
if [ -z "$CYRIUS" ]; then
  dire "3. CYP2D6 : ECHEC, chemin Cyrius non fourni (PGX_CYRIUS)"
  etat cyp2d6 ECHEC 0 "PGX_CYRIUS non defini"
elif [ "$RECEVABLE" = 1 ]; then
  # Dans le dossier de travail, pas dans /tmp de l'hote : un noeud de calcul
  # peut l'avoir petit, en lecture seule ou non partage.
  d=$(mktemp -d -p "$T" cyrius.XXXXXX); echo "$CRAM" > "$d/m.txt"
  (cd "$CYRIUS" && python3 star_caller.py -m "$d/m.txt" -g 38 -o "$d" -p x -t "$FILS" \
      ${FASTA:+-r "$FASTA"} >>"$JOURNAL" 2>&1)
  rc=$?
  if [ $rc -eq 0 ] && [ -s "$d/x.tsv" ] && [ "$(wc -l < "$d/x.tsv")" -ge 2 ]; then
    cp "$d/x.tsv" "$D6"; cp "$d/x.json" "$T/cyp2d6.json" 2>/dev/null
    dip=$(awk -F'\t' 'NR==2{print $2}' "$D6")
    flt=$(awk -F'\t' 'NR==2{print $3}' "$D6")
    # Un diplotype qui n'est pas une paire unique n'est pas transmis. Choisir
    # entre "*1/*46" et "*43/*45" serait une affirmation que l'outil n'a pas
    # faite, et les deux lectures ne donnent pas le meme phenotype ;
    # l'interpreteur refuse d'ailleurs la ligne et s'arrete. L'etage a tourne
    # jusqu'au bout sans conclure : ce n'est pas une panne du module, et le
    # compte rendu rend alors le gene non analyse.
    motif=""
    case "$dip" in
      ""|None)  motif="aucun diplotype rendu";;
      *";"*)    motif="plusieurs diplotypes possibles : $dip";;
      */*)      : ;;
      *)        motif="diplotype illisible : $dip";;
    esac
    if [ -n "$motif" ]; then
      dire "3. CYP2D6 : sans resultat, $motif"
      etat cyp2d6 SANS_RESULTAT $(( $(date +%s)-t0 )) "$motif"
    else
      dire "3. CYP2D6 : $dip ($flt)"
      etat cyp2d6 OK $(( $(date +%s)-t0 )) "$dip ($flt)"
    fi
  else
    dire "3. CYP2D6 : ECHEC (code $rc)"; etat cyp2d6 ECHEC $(( $(date +%s)-t0 )) "star_caller.py code $rc"
  fi
  rm -rf "$d"
else
  dire "3. CYP2D6 : non lance, entrees non recevables"; etat cyp2d6 ECHEC 0 "entrees non recevables"
fi
fi

# ------------------------------------------------------------------ 4. HLA
if a_faire hla; then
t0=$(date +%s)
FQ1="$T/mhc_1.fq"; FQ2="$T/mhc_2.fq"
# Le typage HLA coute a lui seul la moitie du temps d'un genome. Quand la
# reprise est legitime — meme signature, donc memes entrees, meme reference et
# memes ressources — son resultat est reutilise tel quel. Sans cette condition,
# tout genome repris repayait cent soixante-dix secondes pour un resultat
# identique.
# L'empreinte immuable de l'image, et non son etiquette : une etiquette peut
# etre repointee sur une autre image sans changer de nom.
empreinte_image() {
  case "$MOTEUR" in
    docker) docker image inspect "$1" \
              --format '{{index .RepoDigests 0}}' 2>/dev/null \
            || docker image inspect "$1" --format '{{.Id}}' 2>/dev/null;;
    *) if [ -f "$1" ]; then printf 'sha256:%s\n' "$(sha256sum "$1" | cut -d' ' -f1)"
       else echo "$1"; fi;;   # une image de depot : Apptainer la resout a chaque appel
  esac
}
HLA_REPRIS=0
if [ "$REPRISE" = 1 ] && [ -s "$T/hla.tsv" ] && [ "$(wc -l < "$T/hla.tsv")" -ge 2 ] \
   && [ -s "$T/hla.image" ] \
   && [ "$(cat "$T/hla.image")" = "$(empreinte_image "$IMG_OPTITYPE")" ]; then
  HLA_REPRIS=1
else
  rm -rf "$T/hla"
fi
if [ "$HLA_REPRIS" = 0 ] && [ ! -s "$FQ1" ] && [ -s "$TRANCHE.bai" ]; then
  CONTIGS=$(samtools view -H "$TRANCHE" 2>/dev/null | grep -oE 'SN:HLA-[^[:space:]]+' | sed 's/SN://' | tr '\n' ' ')
  samtools view -u "$TRANCHE" chr6:29600000-33100000 $CONTIGS 2>>"$JOURNAL" \
    | samtools collate -u -O - 2>>"$JOURNAL" \
    | samtools fastq -1 "$FQ1" -2 "$FQ2" -0 /dev/null -s /dev/null -n - 2>>"$JOURNAL"
fi
NPAIRES=0
[ -s "$FQ1" ] && NPAIRES=$(( $(wc -l < "$FQ1") / 4 ))
if [ "$HLA_REPRIS" = 1 ]; then
  r=$(awk -F'\t' 'NR==2{print $2"/"$3" "$4"/"$5}' "$T/hla.tsv")
  dire "4. HLA : $r, repris de l'execution precedente"
  etat hla OK $(( $(date +%s)-t0 )) "$r (repris)"
elif [ "$NPAIRES" -lt 200 ]; then
  dire "4. HLA : ECHEC, seulement $NPAIRES paires extraites du complexe majeur"
  etat hla ECHEC $(( $(date +%s)-t0 )) "$NPAIRES paires"
else
  # La configuration de l'image fixe seize fils au placement des lectures. Le
  # typage ne les occupe pas : a quatre, il rend le meme resultat dans le meme
  # temps, et seize fils par genome effondrent un lot. La configuration de
  # l'image est reprise au seul nombre de fils pres.
  CONF_OT=()
  conteneur "$T" "$IMG_OPTITYPE" sh -c \
    'cat "$(dirname "$(readlink -f "$(command -v OptiTypePipeline.py)")")/config.ini"' \
    2>>"$JOURNAL" | awk -v n="$FILS" '
      /^\[/ { s = $0 }
      s == "[mapping]" && /^threads=/ { $0 = "threads=" n; f = 1 }
      { print }
      END { exit !f }' > "$T/optitype.ini" \
    && CONF_OT=(-c "$T/optitype.ini") \
    || dire "4. HLA : configuration d'OptiType illisible, nombre de fils de l'image"
  conteneur "$T" "$IMG_OPTITYPE" \
    OptiTypePipeline.py -i "$T/mhc_1.fq" "$T/mhc_2.fq" --dna -v \
    ${CONF_OT[@]+"${CONF_OT[@]}"} -o "$T/hla" >>"$JOURNAL" 2>&1
  rc=$?
  HLA=$(find "$T/hla" -name "*_result.tsv" 2>/dev/null | sort | tail -1)
  if [ $rc -eq 0 ] && [ -n "$HLA" ] && [ "$(wc -l < "$HLA")" -ge 2 ]; then
    cp "$HLA" "$T/hla.tsv"
    # L'image qui a produit ce resultat, pour que la reprise puisse verifier
    # qu'elle n'a pas change.
    empreinte_image "$IMG_OPTITYPE" > "$T/hla.image"
    r=$(awk -F'\t' 'NR==2{print $2"/"$3" "$4"/"$5}' "$T/hla.tsv")
    dire "4. HLA : $r, sur $NPAIRES paires"; etat hla OK $(( $(date +%s)-t0 )) "$r"
  else
    dire "4. HLA : ECHEC (code $rc)"; etat hla ECHEC $(( $(date +%s)-t0 )) "OptiType code $rc"
  fi
fi
fi

# --------------------------------------- 4b. typage complementaire (PyPGx)
# Quatre genes du core panel RNPGx que l'interpreteur ne sait pas appeler. Le
# typeur lit la region entiere du gene, pas les seules positions retenues : on
# lui prepare donc une extraction distincte, qui ne touche ni au VCF filtre ni
# au VCF qualifie. Les contigs du fichier de variants sont releves ici : un
# contig que l'appelant de la plateforme a ignore rendrait un « reference »
# indistinguable d'une absence de lecture.
if a_faire pypgx; then
t0=$(date +%s)
if [ -z "$CAT_PYPGX" ]; then
  dire "4b. typage complementaire : ignore (ni PGX_IMG_PYPGX ni PGX_PYPGX)"
  etat pypgx IGNORE 0 "typeur non defini"
elif [ "$RECEVABLE" != 1 ]; then
  dire "4b. typage complementaire : ECHEC, genome non recevable"
  etat pypgx ECHEC 0 "genome non recevable"
elif BCF view -R "$RES/pypgx_regions.bed" -Oz -o "$PGV" "$D_VCF/$(basename "$VCF")" >>"$JOURNAL" 2>&1 \
     && BCF index -f -t "$PGV" >>"$JOURNAL" 2>&1; then
  BCF index -s "$D_VCF/$(basename "$VCF")" 2>/dev/null | cut -f1 > "$T/contigs_vcf.txt"
  # Le meme critere de qualite que l'etage 2b, sur les seules positions qui
  # definissent les alleles : un genotype que le controle qualite ecarte ne doit
  # pas servir ici a nommer un allele.
  BCF query -R "$RES/pypgx_positions.bed" \
    -f '%CHROM\t%POS\t%FILTER[\t%GQ\t%DP]\n' "$PGV" > "$T/pypgx_qualites.tsv" 2>>"$JOURNAL"
  if python3 "$RACINE/bin/typage_pypgx.py" --catalogue "$CAT_PYPGX" \
       --variants "$PGV" --profondeur "$DEPTH" --contigs-vcf "$T/contigs_vcf.txt" \
       --qualites "$T/pypgx_qualites.tsv" --gq "$GQ" $IGNORER_FILTRE \
       --profondeur-min "$PROF" --travail "$T" --sortie "$T/pypgx.tsv" \
       --json "$T/pypgx.json" --pypgx "$PYPGX" >>"$JOURNAL" 2>&1 \
     && [ -s "$T/pypgx.tsv" ]; then
    r=$(awk -F'\t' 'NR>1 && $3=="rendu"{n++} NR>1{t++} END{printf "%d rendu(s) sur %d", n, t}' "$T/pypgx.tsv")
    dire "4b. typage complementaire : $r"
    etat pypgx OK $(( $(date +%s)-t0 )) "$r"
  else
    dire "4b. typage complementaire : ECHEC"
    etat pypgx ECHEC $(( $(date +%s)-t0 )) "typage_pypgx.py"
  fi
else
  dire "4b. typage complementaire : ECHEC a l'extraction"
  etat pypgx ECHEC $(( $(date +%s)-t0 )) "bcftools view ou index"
fi
fi

# ------------------------------------------------------- 5. appels externes
if a_faire appels; then
t0=$(date +%s)
if python3 "$RACINE/bin/appels_externes.py" --travail "$T" --sortie "$PO" >>"$JOURNAL" 2>&1; then
  nlig=$(grep -c . "$PO" 2>/dev/null || echo 0)
  # coherence : un etage marque OK doit avoir depose sa ligne
  attendu=0
  grep -q '^cyp2d6	OK' "$ETATS" && attendu=$((attendu+1))
  grep -q '^hla	OK' "$ETATS" && attendu=$((attendu+2))
  # MT-RNR1 n'est attendu que si l'etage 4b l'a effectivement rendu : un gene
  # non conclusif ne doit pas faire echouer l'etage des appels externes.
  [ -s "$T/pypgx.tsv" ] && awk -F'\t' '$1=="MT-RNR1" && $3=="rendu"{f=1} END{exit !f}' "$T/pypgx.tsv" \
    && attendu=$((attendu+1))
  if [ "$nlig" -ge "$attendu" ]; then
    dire "5. appels externes : $nlig ligne(s)"; etat appels OK $(( $(date +%s)-t0 )) "$nlig lignes"
  else
    dire "5. appels externes : ECHEC, $nlig ligne(s) pour $attendu attendue(s)"
    etat appels ECHEC $(( $(date +%s)-t0 )) "$nlig lignes pour $attendu attendues"
  fi
else
  dire "5. appels externes : ECHEC"; etat appels ECHEC $(( $(date +%s)-t0 )) "appels_externes.py"
fi
fi

# ------------------------------------------------------------ 6. PharmCAT
# Aucun repli : sans le VCF qualifie, l'etage echoue. Se rabattre sur le VCF
# filtre remettrait en service les genotypes ecartes sur leur qualite.
if a_faire pharmcat; then
t0=$(date +%s)
rm -f "$PRE" "$RAP"
if [ ! -s "$QUAL" ]; then
  dire "6. interpretation : ECHEC, VCF qualifie absent"
  etat pharmcat ECHEC 0 "qualifie.vcf.gz absent"
elif conteneur "$T" "$IMG_PHARMCAT" \
       /pharmcat/pharmcat_vcf_preprocessor -vcf "$QUAL" -o "$T" --absent-to-ref >>"$JOURNAL" 2>&1 \
     && [ -s "$PRE" ]; then
  ARG=(); [ -s "$PO" ] && ARG=(-po "$PO")
  if conteneur "$D_SORTIE" "$IMG_PHARMCAT" \
       java -jar /pharmcat/pharmcat.jar -vcf "$PRE" "${ARG[@]}" \
       -o "$D_SORTIE/sortie" -bf "$ECH" -reporterJson -del >>"$JOURNAL" 2>&1 \
     && [ -s "$RAP" ]; then
    dire "6. interpretation : $RAP"
    etat pharmcat OK $(( $(date +%s)-t0 )) "$ECH.report.json"
  else
    dire "6. interpretation : ECHEC de l'interpreteur"
    etat pharmcat ECHEC $(( $(date +%s)-t0 )) "pharmcat.jar"
  fi
else
  dire "6. interpretation : ECHEC au pretraitement"
  etat pharmcat ECHEC $(( $(date +%s)-t0 )) "pretraitement, $PRE absent"
fi
fi

# --------------------------------------------------------------- 7. rendu
if a_faire rendu; then
t0=$(date +%s)
if [ -s "$RAP" ] && [ -s "$T/perimetre.json" ]; then
  # L'identite, si elle a ete fournie, passe par un chemin de fichier et non
  # par sa valeur : un nom en argument serait lisible dans la table des
  # processus et dans les journaux de l'ordonnanceur. Fournie, elle leve la
  # mention prototype du document.
  if python3 "$RACINE/bin/compte_rendu.py" --rapport "$RAP" --perimetre "$T/perimetre.json" \
       --sortie "$SORTIE/sortie/CR_$ECH.pdf" --echantillon "$ECH" \
       --couverture-min "$COUV_MIN" \
       ${IDENTITE:+--identite "$IDENTITE"} ${CONTEXTE:+--contexte "$CONTEXTE"} >>"$JOURNAL" 2>&1 \
     && [ -s "$SORTIE/sortie/CR_$ECH.pdf" ]; then
    dire "7. compte rendu : $SORTIE/sortie/CR_$ECH.pdf"
    etat rendu OK $(( $(date +%s)-t0 )) "CR_$ECH.pdf"
  else
    dire "7. compte rendu : ECHEC"; etat rendu ECHEC $(( $(date +%s)-t0 )) "compte_rendu.py"
  fi
else
  dire "7. compte rendu : ECHEC, rapport ou perimetre absent"
  etat rendu ECHEC 0 "rapport ou perimetre absent"
fi
fi

# ------------------------------------------------------------- provenance
if a_faire provenance; then
# Lancee seule, la trace ne connait pas la couverture mesuree par l'etage 2b :
# elle la relit.
[ "$DECOUPE" = 1 ] && [ -z "${COUV:-}" ] && lire_verdict
python3 "$RACINE/bin/provenance.py" --sortie "$SORTIE" --echantillon "$ECH" \
  --cram "$CRAM" --vcf "$VCF" --fasta "$FASTA" --ressources "$RES" \
  --gq "$GQ" --profondeur "$PROF" --couverture-min "$COUV_MIN" \
  --couverture "${COUV:-}" --melange-max "$MELANGE_MAX" --melange "${MEL_PART:-}" \
  --images "$IMAGES" --moteur "$MOTEUR" --cyrius "$CYRIUS" --pypgx "$PYPGX" \
  --filtre "$([ -n "$IGNORER_FILTRE" ] && echo ignore || echo respecte)" \
  --perimetre-clinique "$RES/perimetre_rnpgx.json" >>"$JOURNAL" 2>&1

REUSSI=$(python3 -c "import json;print(json.load(open(r'$SORTIE/sortie/provenance.json',encoding='utf-8'))['reussite_complete'])" 2>/dev/null || echo False)
dire "termine, reussite complete : $REUSSI"
[ "$REUSSI" = "True" ] && exit 0 || exit 1
fi

# Sans la trace, le code de retour ne dit que les etages de cette execution :
# un ordonnanceur doit savoir si celui qu'il a lance a abouti, pas si le
# genome entier l'a fait. Un etage sans resultat ou ignore n'est pas un echec.
# Ses propres etages seulement : un etage lance en parallele ecrit dans la meme
# table, et son echec n'est pas celui de cette execution.
ECHECS=$(tail -n +"$((DEBUT_ETATS + 1))" "$ETATS" 2>/dev/null \
  | awk -F'\t' -v e=",$ETAGES," 'index(e, "," $1 ",") && $2 == "ECHEC" {print $1}' \
  | sort -u | tr '\n' ' ')
if [ -n "$ECHECS" ]; then
  dire "etage(s) en echec : $ECHECS"
  exit 1
fi
exit 0
