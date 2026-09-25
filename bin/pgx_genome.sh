#!/bin/bash
# Module pharmacogenetique pour genome entier. Un echantillon, sept etages.
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
# que si les neuf etages ont abouti.
#
# Options :
#   --fasta CHEMIN     reference d'alignement, requise pour un CRAM
#   --gq N             seuil de qualite de genotype           (defaut 20)
#   --profondeur N     seuil de profondeur par position       (defaut 10)
#   --fils N           fils pour samtools et Cyrius           (defaut 4)
#   --reprise          reutilise les intermediaires si les entrees sont identiques
#   --forcer           ecrit dans un dossier non vide sans reprise
#
# Variables : PGX_FASTA PGX_CYRIUS PGX_GQ PGX_PROFONDEUR PGX_FILS
#             PGX_IMG_BCFTOOLS PGX_IMG_PHARMCAT PGX_IMG_OPTITYPE

set -uo pipefail

RACINE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RES="$RACINE/ressources"

CRAM=""; VCF=""; SORTIE=""; ECH=""; REPRISE=0; FORCER=0
FASTA="${PGX_FASTA:-}"
GQ="${PGX_GQ:-20}"
PROF="${PGX_PROFONDEUR:-10}"
IGNORER_FILTRE="${PGX_IGNORER_FILTRE:+--ignorer-filtre}"
FILS="${PGX_FILS:-4}"
# Complement RNPGx (classes 1 et 2 hors definitions PharmCAT) : present, il est
# mesure ; absent, le module se comporte comme avant.
COMPLEMENT=""
[ -s "$RES/rnpgx_complement.vcf" ] && COMPLEMENT="$RES/rnpgx_complement.vcf"
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
    --fils) FILS="$2"; shift 2;;
    --reprise) REPRISE=1; shift;;
    --forcer) FORCER=1; shift;;
    --ignorer-filtre) IGNORER_FILTRE="--ignorer-filtre"; shift;;
    -h|--help) sed -n '2,26p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done

[ -n "$CRAM" ] && [ -n "$VCF" ] && [ -n "$SORTIE" ] || {
  echo "usage : pgx_genome.sh --cram X.cram --vcf X.vcf.gz --sortie DOSSIER" >&2; exit 2; }
[ -s "$CRAM" ] || { echo "alignement introuvable : $CRAM" >&2; exit 2; }
[ -s "$VCF" ]  || { echo "fichier de variants introuvable : $VCF" >&2; exit 2; }
[ -n "$ECH" ] || ECH="$(basename "$CRAM" | sed 's/\.\(cram\|bam\)$//')"

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
MONTE=(-v "$D_CRAM":"$D_CRAM":ro -v "$D_VCF":"$D_VCF":ro -v "$D_SORTIE":"$D_SORTIE" -v "$RES":"$RES":ro)
if [ -n "$FASTA" ]; then
  D_FASTA="$(cd "$(dirname "$FASTA")" && pwd)"
  MONTE+=(-v "$D_FASTA":"$D_FASTA":ro)
fi
BCF() { docker run --rm "${MONTE[@]}" -w "$D_SORTIE" "$IMG_BCF" bcftools "$@"; }

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
SIGNATURE="$(empreinte "$CRAM")|$(empreinte "$VCF")|gq=$GQ|prof=$PROF|filtre=${IGNORER_FILTRE:-respecte}"
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
    --images "$IMG_BCF,$IMG_PHARMCAT,$IMG_OPTITYPE" --cyrius "$CYRIUS" \
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
       --gq "$GQ" --profondeur-min "$PROF" --echantillon "$ECH" \
       --perimetre-clinique "$RES/perimetre_rnpgx.json" $IGNORER_FILTRE >>"$JOURNAL" 2>&1 \
     && docker run --rm "${MONTE[@]}" "$IMG_BCF" bgzip -f "$T/qualifie.vcf" >>"$JOURNAL" 2>&1 \
     && BCF index -f -t "$QUAL" >>"$JOURNAL" 2>&1 && [ -s "$QUAL" ]; then
  ret=$(python3 -c "import json;d=json.load(open(r'$T/perimetre.json',encoding='utf-8'));print(d['positions_retenues'],d['positions_visees'],sum(1 for v in d['genes'].values() if v['statut']=='complet'),len(d['genes']),d.get('positions_sans_GQ',0))")
  dire "2b. controle qualite : $(echo $ret | awk '{print $1" positions retenues sur "$2", "$3" genes complets sur "$4", "$5" sans qualite"}')"
  etat qc OK $(( $(date +%s)-t0 )) "$ret"
else
  dire "2b. controle qualite : ECHEC"
  etat qc ECHEC $(( $(date +%s)-t0 )) "qc_perimetre.py, bgzip ou index"
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
rm -rf "$T/hla"
if [ ! -s "$FQ1" ] && [ -s "$TRANCHE.bai" ]; then
  CONTIGS=$(samtools view -H "$TRANCHE" 2>/dev/null | grep -oE 'SN:HLA-[^[:space:]]+' | sed 's/SN://' | tr '\n' ' ')
  samtools view -u "$TRANCHE" chr6:29600000-33100000 $CONTIGS 2>>"$JOURNAL" \
    | samtools collate -u -O - 2>>"$JOURNAL" \
    | samtools fastq -1 "$FQ1" -2 "$FQ2" -0 /dev/null -s /dev/null -n - 2>>"$JOURNAL"
fi
NPAIRES=0
[ -s "$FQ1" ] && NPAIRES=$(( $(wc -l < "$FQ1") / 4 ))
if [ "$NPAIRES" -lt 200 ]; then
  dire "4. HLA : ECHEC, seulement $NPAIRES paires extraites du complexe majeur"
  etat hla ECHEC $(( $(date +%s)-t0 )) "$NPAIRES paires"
else
  docker run --rm "${MONTE[@]}" -w "$T" "$IMG_OPTITYPE" \
    OptiTypePipeline.py -i "$T/mhc_1.fq" "$T/mhc_2.fq" --dna -v -o "$T/hla" >>"$JOURNAL" 2>&1
  rc=$?
  HLA=$(find "$T/hla" -name "*_result.tsv" 2>/dev/null | sort | tail -1)
  if [ $rc -eq 0 ] && [ -n "$HLA" ] && [ "$(wc -l < "$HLA")" -ge 2 ]; then
    cp "$HLA" "$T/hla.tsv"
    r=$(awk -F'\t' 'NR==2{print $2"/"$3" "$4"/"$5}' "$T/hla.tsv")
    dire "4. HLA : $r, sur $NPAIRES paires"; etat hla OK $(( $(date +%s)-t0 )) "$r"
  else
    dire "4. HLA : ECHEC (code $rc)"; etat hla ECHEC $(( $(date +%s)-t0 )) "OptiType code $rc"
  fi
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
elif docker run --rm "${MONTE[@]}" -w "$T" "$IMG_PHARMCAT" \
       /pharmcat/pharmcat_vcf_preprocessor -vcf "$QUAL" -o "$T" --absent-to-ref >>"$JOURNAL" 2>&1 \
     && [ -s "$PRE" ]; then
  ARG=(); [ -s "$PO" ] && ARG=(-po "$PO")
  if docker run --rm "${MONTE[@]}" -w "$D_SORTIE" "$IMG_PHARMCAT" \
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
  if python3 "$RACINE/bin/compte_rendu.py" --rapport "$RAP" --perimetre "$T/perimetre.json" \
       --sortie "$SORTIE/sortie/CR_$ECH.pdf" --echantillon "$ECH" >>"$JOURNAL" 2>&1 \
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
  --gq "$GQ" --profondeur "$PROF" \
  --images "$IMG_BCF,$IMG_PHARMCAT,$IMG_OPTITYPE" --cyrius "$CYRIUS" \
  --filtre "$([ -n "$IGNORER_FILTRE" ] && echo ignore || echo respecte)" \
  --perimetre-clinique "$RES/perimetre_rnpgx.json" >>"$JOURNAL" 2>&1

REUSSI=$(python3 -c "import json;print(json.load(open(r'$SORTIE/sortie/provenance.json',encoding='utf-8'))['reussite_complete'])" 2>/dev/null || echo False)
dire "termine, reussite complete : $REUSSI"
[ "$REUSSI" = "True" ] && exit 0 || exit 1
