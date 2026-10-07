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
#   --fasta CHEMIN     reference d'alignement, requise pour un CRAM
#   --gq N             seuil de qualite de genotype           (defaut 20)
#   --profondeur N     seuil de profondeur par position       (defaut 10)
#   --couverture-min N couverture mediane minimale du perimetre (defaut 18,
#                      0 pour accepter un genome hors du domaine valide)
#   --fils N           fils pour samtools et Cyrius           (defaut 4)
#   --reprise          reutilise les intermediaires si les entrees sont identiques
#   --forcer           ecrit dans un dossier non vide sans reprise
#
# Variables : PGX_FASTA PGX_CYRIUS PGX_GQ PGX_PROFONDEUR PGX_COUVERTURE_MIN
#             PGX_FILS
#             PGX_IMG_BCFTOOLS PGX_IMG_PHARMCAT PGX_IMG_OPTITYPE PGX_PYPGX

set -uo pipefail

RACINE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RES="$RACINE/ressources"

CRAM=""; VCF=""; SORTIE=""; ECH=""; REPRISE=0; FORCER=0; IDENTITE=""
FASTA="${PGX_FASTA:-}"
GQ="${PGX_GQ:-20}"
PROF="${PGX_PROFONDEUR:-10}"
# Couverture mediane sous laquelle le genome sort du domaine de validite.
# Dix-huit est le point le plus bas reellement mesure sur le perimetre clinique :
# le seuil n'extrapole pas sous lui. Les neuf executions refusees par la mesure
# plafonnent a treize, la marge est donc de cinq fois.
COUV_MIN="${PGX_COUVERTURE_MIN:-18}"
IGNORER_FILTRE="${PGX_IGNORER_FILTRE:+--ignorer-filtre}"
FILS="${PGX_FILS:-4}"
# Complement RNPGx (classes 1 et 2 hors definitions PharmCAT) : present, il est
# mesure ; absent, le module se comporte comme avant.
COMPLEMENT=""
[ -s "$RES/rnpgx_complement.vcf" ] && COMPLEMENT="$RES/rnpgx_complement.vcf"
# Typage complementaire (MT-RNR1, BCHE, MTHFR, POR) : quatre genes du core panel
# que l'interpreteur ne sait pas appeler. L'etage ne s'arme que si le catalogue
# et la commande sont tous deux presents ; sinon il est ignore, sans consequence
# sur le reste.
PYPGX="${PGX_PYPGX:-}"
CAT_PYPGX=""
[ -n "$PYPGX" ] && [ -s "$RES/pypgx_genes.json" ] && CAT_PYPGX="$RES/pypgx_genes.json"
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
    --fils) FILS="$2"; shift 2;;
    --identite) IDENTITE="$2"; shift 2;;
    --reprise) REPRISE=1; shift;;
    --forcer) FORCER=1; shift;;
    --ignorer-filtre) IGNORER_FILTRE="--ignorer-filtre"; shift;;
    -h|--help) sed -n '2,29p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    --version) cat "$RACINE/VERSION" 2>/dev/null || echo inconnue; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done

[ -n "$CRAM" ] && [ -n "$VCF" ] && [ -n "$SORTIE" ] || {
  echo "usage : pgx_genome.sh --cram X.cram --vcf X.vcf.gz --sortie DOSSIER" >&2; exit 2; }
[ -s "$CRAM" ] || { echo "alignement introuvable : $CRAM" >&2; exit 2; }
[ -s "$VCF" ]  || { echo "fichier de variants introuvable : $VCF" >&2; exit 2; }
[ -n "$ECH" ] || ECH="$(basename "$CRAM" | sed 's/\.\(cram\|bam\)$//')"
# Un seuil non numerique serait lu zero par awk : la porte tomberait sans bruit.
case "$COUV_MIN" in
  ''|*[!0-9]*) echo "couverture minimale invalide : $COUV_MIN" >&2; exit 2;;
esac

mkdir -p "$SORTIE"/{travail,sortie}
JOURNAL="$SORTIE/journal.txt"
ETATS="$SORTIE/travail/etats.tsv"

# --------------------------------------------------------------------- verrou
# Deux executions sur le meme dossier se detruisent mutuellement : elles purgent
# le dossier de travail l'une sous l'autre et melangent leurs lignes d'etat. Le
# verrou est pris par mkdir, qui est atomique.
VERROU="$SORTIE/.verrou"
if ! mkdir "$VERROU" 2>/dev/null; then
  autre=$(cat "$VERROU/pid" 2>/dev/null)
  if [ -n "$autre" ] && kill -0 "$autre" 2>/dev/null; then
    echo "une autre execution travaille deja dans $SORTIE (processus $autre)" >&2
    exit 3
  fi
  # le detenteur precedent est mort : on reprend le verrou
  rm -rf "$VERROU" && mkdir "$VERROU" 2>/dev/null || { echo "verrou illisible" >&2; exit 3; }
fi
echo "$$" > "$VERROU/pid"
date -Iseconds > "$VERROU/depuis"
nettoyer() { rm -rf "$VERROU"; }
trap nettoyer EXIT INT TERM

dire() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$JOURNAL"; }
etat() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "$ETATS"; }

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
conteneur() {   # $1 dossier de travail  $2 image  $3... commande
  local cwd="$1" image="$2"; shift 2
  local m=() d
  case "$MOTEUR" in
    docker)
      for d in "${DOSSIERS[@]}"; do m+=(-v "$d":"$d":ro); done
      for d in "${ECRIVABLES[@]}"; do m+=(-v "$d":"$d"); done
      docker run --rm "${m[@]}" -w "$cwd" "$image" "$@"
      ;;
    apptainer|singularity)
      for d in "${DOSSIERS[@]}"; do m+=(--bind "$d":"$d":ro); done
      for d in "${ECRIVABLES[@]}"; do m+=(--bind "$d":"$d"); done
      "$MOTEUR" exec --cleanenv --no-home --pwd "$cwd" "${m[@]}" \
        "docker://$image" "$@"
      ;;
    *)
      echo "moteur de conteneurs inconnu : $MOTEUR" >&2; return 127
      ;;
  esac
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
SIGNATURE="$(empreinte "$CRAM")|$(empreinte "$VCF")|$(empreinte "$FASTA")|gq=$GQ|prof=$PROF|couv=$COUV_MIN|filtre=${IGNORER_FILTRE:-respecte}|res=$EMPREINTE_RES"
SIG_FICHIER="$T/signature.txt"
if [ -s "$SIG_FICHIER" ]; then
  if [ "$(cat "$SIG_FICHIER")" = "$SIGNATURE" ] && [ "$REPRISE" = 1 ]; then
    : # meme entree, reprise demandee : on garde les intermediaires
  else
    rm -rf "$T"; mkdir -p "$T"
  fi
elif [ -n "$(ls -A "$T" 2>/dev/null)" ] && [ "$FORCER" = 0 ]; then
  rm -rf "$T"; mkdir -p "$T"
fi
printf '%s' "$SIGNATURE" > "$SIG_FICHIER"
: > "$JOURNAL"
printf 'etage\tetat\tsecondes\tdetail\n' > "$ETATS"

dire "echantillon $ECH"
dire "alignement  $CRAM"
dire "variants    $VCF"

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
  for e in filtre tranche qc cyp2d6 hla appels pharmcat rendu; do
    etat "$e" ECHEC 0 "entrees non recevables"
  done
  dire "arret : entrees non recevables, aucun etage n'a demarre"
  python3 "$RACINE/bin/provenance.py" --sortie "$SORTIE" --echantillon "$ECH" \
    --cram "$CRAM" --vcf "$VCF" --fasta "$FASTA" --ressources "$RES" \
    --gq "$GQ" --profondeur "$PROF" \
    --images "$IMG_BCF,$IMG_PHARMCAT,$IMG_OPTITYPE" --cyrius "$CYRIUS" --pypgx "$PYPGX" \
  --filtre "$([ -n "$IGNORER_FILTRE" ] && echo ignore || echo respecte)" \
  --perimetre-clinique "$RES/perimetre_rnpgx.json" >>"$JOURNAL" 2>&1
  dire "termine, reussite complete : False"
  exit 2
fi

# ------------------------------------------------------------------ 1. filtre
# Les regions visees reunissent les definitions PharmCAT et le complement RNPGx
# (classes 1 et 2 absentes de PharmCAT). Le complement est optionnel : absent, le
# module se comporte exactement comme avant.
t0=$(date +%s)
FIL="$T/filtre.vcf.gz"
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

# --------------------------------------------------------- 2a. tranche unique
# Une seule traversee de l'alignement sert au controle qualite et au typage HLA.
# Le BED de tranche est construit ici : les contigs HLA alternatifs dependent de
# la reference de la plateforme et ne peuvent pas etre figes dans les ressources.
t0=$(date +%s)
TRANCHE="$T/tranche.bam"
BEDT="$T/tranche.bed"
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

# ------------------------------------------------------- 2b. controle qualite
t0=$(date +%s)
DEPTH="$T/profondeur.txt"
QUAL="$T/qualifie.vcf.gz"
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
       --echantillon "$ECH" \
       --perimetre-clinique "$RES/perimetre_rnpgx.json" $IGNORER_FILTRE >>"$JOURNAL" 2>&1 \
     && conteneur "$D_SORTIE" "$IMG_BCF" bgzip -f "$T/qualifie.vcf" >>"$JOURNAL" 2>&1 \
     && BCF index -f -t "$QUAL" >>"$JOURNAL" 2>&1 && [ -s "$QUAL" ]; then
  ret=$(python3 -c "import json;d=json.load(open(r'$T/perimetre.json',encoding='utf-8'));print(d['positions_retenues'],d['positions_visees'],sum(1 for v in d['genes'].values() if v['statut']=='complet'),len(d['genes']),d.get('positions_sans_GQ',0))")
  dire "2b. controle qualite : $(echo $ret | awk '{print $1" positions retenues sur "$2", "$3" genes complets sur "$4", "$5" sans qualite"}')"
  etat qc OK $(( $(date +%s)-t0 )) "$ret"
  # Le verdict sur le domaine est rendu par l'etage lui-meme, ou il se teste.
  VERDICT=$(python3 -c "import json,sys;sys.path.insert(0,r'$RACINE/bin');import qc_perimetre as q;d=json.load(open(r'$T/perimetre.json',encoding='utf-8'));c=q.couverture_retenue(d);print('' if c is None else c, 1 if d.get('couverture_hors_domaine') else 0)" 2>>"$JOURNAL")
  COUV="${VERDICT% *}"
  dire "    couverture mediane du perimetre : ${COUV:-inconnue}x"
  if [ "${VERDICT##* }" = 1 ]; then
    HORS_DOMAINE="couverture mediane ${COUV}x, domaine valide a partir de ${COUV_MIN}x"
  fi
else
  dire "2b. controle qualite : ECHEC"
  etat qc ECHEC $(( $(date +%s)-t0 )) "qc_perimetre.py, bgzip ou index"
fi

# Sous le seuil, un allele variant est lu comme reference et l interpreteur
# rend un diplotype faux, sans reserve : le phenotype change, et la conduite a
# tenir avec lui. Le genome est donc refuse et non rendu avec un avertissement,
# parce que rien dans le document ne distinguerait un appel sur une position lue
# d'un appel sur une position supposee. --couverture-min 0 leve la porte.
if [ -n "${HORS_DOMAINE:-}" ]; then
  # On garde la mesure, qui est le motif du refus, et on purge ce qui pourrait
  # etre pris pour un resultat.
  rm -rf "$SORTIE/sortie"; mkdir -p "$SORTIE/sortie"
  for e in cyp2d6 hla appels pharmcat rendu; do
    etat "$e" ECHEC 0 "couverture hors domaine"
  done
  [ -n "$CAT_PYPGX" ] && etat pypgx ECHEC 0 "couverture hors domaine" \
    || etat pypgx IGNORE 0 "PGX_PYPGX non defini"
  dire "arret : $HORS_DOMAINE"
  python3 "$RACINE/bin/provenance.py" --sortie "$SORTIE" --echantillon "$ECH" \
    --cram "$CRAM" --vcf "$VCF" --fasta "$FASTA" --ressources "$RES" \
    --gq "$GQ" --profondeur "$PROF" --couverture-min "$COUV_MIN" \
    --couverture "$COUV" \
    --images "$IMG_BCF,$IMG_PHARMCAT,$IMG_OPTITYPE" --cyrius "$CYRIUS" --pypgx "$PYPGX" \
    --filtre "$([ -n "$IGNORER_FILTRE" ] && echo ignore || echo respecte)" \
    --perimetre-clinique "$RES/perimetre_rnpgx.json" >>"$JOURNAL" 2>&1
  dire "termine, reussite complete : False"
  exit 2
fi

# --------------------------------------------------------------- 3. CYP2D6
t0=$(date +%s)
D6="$T/cyp2d6.tsv"
if [ -z "$CYRIUS" ]; then
  dire "3. CYP2D6 : ECHEC, chemin Cyrius non fourni (PGX_CYRIUS)"
  etat cyp2d6 ECHEC 0 "PGX_CYRIUS non defini"
elif [ "$RECEVABLE" = 1 ]; then
  d=$(mktemp -d); echo "$CRAM" > "$d/m.txt"
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

# ------------------------------------------------------------------ 4. HLA
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
    *) echo "$1";;   # Apptainer resout l'image a chaque appel
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
  conteneur "$T" "$IMG_OPTITYPE" \
    OptiTypePipeline.py -i "$T/mhc_1.fq" "$T/mhc_2.fq" --dna -v -o "$T/hla" >>"$JOURNAL" 2>&1
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

# --------------------------------------- 4b. typage complementaire (PyPGx)
# Quatre genes du core panel RNPGx que l'interpreteur ne sait pas appeler. Le
# typeur lit la region entiere du gene, pas les seules positions retenues : on
# lui prepare donc une extraction distincte, qui ne touche ni au VCF filtre ni
# au VCF qualifie. Les contigs du fichier de variants sont releves ici : un
# contig que l'appelant de la plateforme a ignore rendrait un « reference »
# indistinguable d'une absence de lecture.
t0=$(date +%s)
PGV="$T/pypgx_entree.vcf.gz"
if [ -z "$CAT_PYPGX" ]; then
  dire "4b. typage complementaire : ignore (PGX_PYPGX non defini)"
  etat pypgx IGNORE 0 "PGX_PYPGX non defini"
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

# ------------------------------------------------------- 5. appels externes
t0=$(date +%s)
PO="$T/appels_externes.tsv"
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

# ------------------------------------------------------------ 6. PharmCAT
# Aucun repli : sans le VCF qualifie, l'etage echoue. Se rabattre sur le VCF
# filtre remettrait en service les genotypes ecartes sur leur qualite.
t0=$(date +%s)
PRE="$T/$(basename "$QUAL" .gz).preprocessed.vcf.bgz"
RAP="$SORTIE/sortie/$ECH.report.json"
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

# --------------------------------------------------------------- 7. rendu
t0=$(date +%s)
if [ -s "$RAP" ] && [ -s "$T/perimetre.json" ]; then
  # L'identite, si elle a ete fournie, passe par un chemin de fichier et non
  # par sa valeur : un nom en argument serait lisible dans la table des
  # processus et dans les journaux de l'ordonnanceur.
  if python3 "$RACINE/bin/compte_rendu.py" --rapport "$RAP" --perimetre "$T/perimetre.json" \
       --sortie "$SORTIE/sortie/CR_$ECH.pdf" --echantillon "$ECH" \
       --couverture-min "$COUV_MIN" \
       ${IDENTITE:+--identite "$IDENTITE"} >>"$JOURNAL" 2>&1 \
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

# ------------------------------------------------------------- provenance
python3 "$RACINE/bin/provenance.py" --sortie "$SORTIE" --echantillon "$ECH" \
  --cram "$CRAM" --vcf "$VCF" --fasta "$FASTA" --ressources "$RES" \
  --gq "$GQ" --profondeur "$PROF" --couverture-min "$COUV_MIN" \
  --couverture "${COUV:-}" \
  --images "$IMG_BCF,$IMG_PHARMCAT,$IMG_OPTITYPE" --cyrius "$CYRIUS" --pypgx "$PYPGX" \
  --filtre "$([ -n "$IGNORER_FILTRE" ] && echo ignore || echo respecte)" \
  --perimetre-clinique "$RES/perimetre_rnpgx.json" >>"$JOURNAL" 2>&1

REUSSI=$(python3 -c "import json;print(json.load(open(r'$SORTIE/sortie/provenance.json',encoding='utf-8'))['reussite_complete'])" 2>/dev/null || echo False)
dire "termine, reussite complete : $REUSSI"
[ "$REUSSI" = "True" ] && exit 0 || exit 1
