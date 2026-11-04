#!/usr/bin/env bash
# Reconstruit les entrees du banc de validation depuis les donnees publiques du
# 1000 Genomes, puis fait tourner le module et calcule la fiche de performance.
#
#   bench/preparer_banc.sh --donnees DOSSIER --fasta GRCh38.fa [--genomes liste.txt]
#                          [--etapes alignements,variants,lot,fiche] [--parallele 6]
#
# DOSSIER recoit les alignements (cram/), les fichiers de variants (vcf/), le lot
# (lot/) et la fiche (fiche/). La liste des genomes vaut par defaut celle de
# bench/echantillons.csv. Chaque etape reprend la ou la precedente s'est arretee.
#
#   alignements  alignements 30x NYGC et leurs index, taille et MD5 verifies
#                contre l'index public du projet ; une passe a la fois, sans
#                relance automatique, qui reappend des octets en double ;
#   variants     appels de la plateforme NYGC sur les regions de l'interpreteur,
#                tires de ses fichiers par chromosome, et appels bcftools sur les
#                regions du typage complementaire et du complement RNPGx, chrM
#                compris — les fichiers NYGC ne couvrent pas ces regions ;
#   lot          bin/pgx_lot.sh sur le manifeste construit ;
#   fiche        bench/fiche_performance.py, puis bench/confronter_puce.sh.
#
# Variables du module attendues pour l'etape lot : PGX_CYRIUS, et pour le typage
# complementaire PGX_IMG_PYPGX et PYPGX_BUNDLE.
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RACINE="$(dirname "$ICI")"
DONNEES=""; FASTA=""; LISTE=""; ETAPES="alignements,variants,lot,fiche"; PAR=6
while [ $# -gt 0 ]; do
  case "$1" in
    --donnees) DONNEES="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --genomes) LISTE="$2"; shift 2;;
    --etapes) ETAPES="$2"; shift 2;;
    --parallele) PAR="$2"; shift 2;;
    -h|--help) sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$DONNEES" ] && [ -n "$FASTA" ] || { echo "usage : preparer_banc.sh --donnees DOSSIER --fasta GRCh38.fa" >&2; exit 2; }
[ -s "$FASTA.fai" ] || { echo "reference non indexee : $FASTA" >&2; exit 2; }
etape() { [[ ",$ETAPES," == *",$1,"* ]]; }
mkdir -p "$DONNEES"/{cram,vcf,travail}
W=$DONNEES/travail
IMG=${PGX_IMG_BCFTOOLS:-quay.io/biocontainers/bcftools:1.24--h118bc1c_2}
D_FASTA="$(cd "$(dirname "$FASTA")" && pwd)"
BCF() { docker run --rm --user "$(id -u):$(id -g)" -v "$DONNEES":"$DONNEES" -v "$RACINE":"$RACINE":ro \
          -v "$D_FASTA":"$D_FASTA":ro -w "$W" "$IMG" bcftools "$@"; }

if [ -n "$LISTE" ]; then grep -v '^#' "$LISTE" | cut -f1 | sed '/^$/d' > "$W/genomes.txt"
else grep -v '^#' "$ICI/echantillons.csv" | tail -n +2 | cut -d, -f1 > "$W/genomes.txt"; fi
echo "genomes : $(wc -l < "$W/genomes.txt")"

IDX_URL=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
for f in 1000G_2504_high_coverage.sequence.index 1000G_698_related_high_coverage.sequence.index; do
  [ -s "$W/$f" ] || curl -fsS -o "$W/$f" "$IDX_URL/$f" || { echo "index introuvable : $f" >&2; exit 1; }
done

# ------------------------------------------------------------ alignements
if etape alignements; then
  while read -r E; do
    f=$DONNEES/cram/$E.cram
    [ -s "$f" ] && [ -s "$f.crai" ] && continue
    l=$(grep -h "/$E.final.cram	" "$W"/*.sequence.index | head -1)
    url=$(printf '%s' "$l" | cut -f1 | sed 's#^ftp://#https://#'); md5=$(printf '%s' "$l" | cut -f2)
    [ -n "$url" ] || { echo "$E : absent de l'index du projet"; continue; }
    taille=$(curl -sI "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
    for passe in $(seq 1 40); do
      have=$(stat -c %s "$f.part" 2>/dev/null || echo 0)
      [ "$have" -eq "$taille" ] && break
      [ "$have" -gt "$taille" ] && rm -f "$f.part"
      curl -s -C - -o "$f.part" "$url"
    done
    if [ "$(stat -c %s "$f.part" 2>/dev/null || echo 0)" -ne "$taille" ]; then echo "$E : taille incomplete"; continue; fi
    if [ "$(md5sum "$f.part" | cut -d' ' -f1)" != "$md5" ]; then echo "$E : MD5 faux, fichier retire"; rm -f "$f.part"; continue; fi
    mv "$f.part" "$f" && curl -fsS -o "$f.crai" "$url.crai"
    samtools quickcheck "$f" && echo "$E : alignement conforme" || echo "$E : alignement illisible"
  done < "$W/genomes.txt"
fi

# ------------------------------------------------------------ variants
if etape variants; then
  NYGC=$IDX_URL/working/20201028_3202_raw_GT_with_annot
  PRE=20201028_CCDG_14151_B01_GRM_WGS_2020-08-05
  SUF=recalibrated_variants.vcf.gz
  # Appels de la plateforme, tous les genomes ensemble, chromosome par chromosome.
  # Le filtrage de cohorte est retire : la decision d'exclure appartient au
  # module, et un filtre de cohorte supprimait une position de CYP2B6.
  PARTS=""
  for c in $(cut -f1 "$RACINE/ressources/pharmcat_positions.bed" | sort -uV); do
    if [ ! -s "$W/br_$c.vcf.gz.tbi" ]; then
      grep "^$c	" "$RACINE/ressources/pharmcat_positions.bed" > "$W/reg_$c.bed"
      BCF view -R "$W/reg_$c.bed" -S "$W/genomes.txt" --force-samples -Oz -o "$W/br_$c.vcf.gz" \
        "$NYGC/${PRE}_${c}.${SUF}" && BCF index -f -t "$W/br_$c.vcf.gz" || echo "$c : extraction en echec"
    fi
    PARTS="$PARTS $W/br_$c.vcf.gz"
  done
  BCF concat -a -Oz -o "$W/tout.vcf.gz" $PARTS && BCF index -f -t "$W/tout.vcf.gz"
  BCF annotate -x FILTER -Oz -o "$W/plateforme.vcf.gz" "$W/tout.vcf.gz" && BCF index -f -t "$W/plateforme.vcf.gz"
  cut -f1-3 "$RACINE/ressources/pypgx_regions.bed" "$RACINE/ressources/rnpgx_complement.bed" | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"} {if (c==$1 && $2<=e) { if($3>e) e=$3 } else { if(c!="") print c,s,e; c=$1; s=$2; e=$3 }} END{if(c!="") print c,s,e}' \
    > "$W/regions_complement.bed"
  while read -r E; do
    o=$DONNEES/vcf/$E.vcf.gz
    [ -s "$o.tbi" ] && continue
    [ -s "$DONNEES/cram/$E.cram" ] || continue
    BCF view -s "$E" -Oz -o "$W/$E.p.vcf.gz" "$W/plateforme.vcf.gz" && BCF index -f -t "$W/$E.p.vcf.gz"
    # L'alignement peut etre un lien : le conteneur ne voit que ce qui est monte,
    # il recoit donc le chemin reel et son dossier.
    al=$(readlink -f "$DONNEES/cram/$E.cram"); d_al=$(dirname "$al")
    docker run --rm --user "$(id -u):$(id -g)" -v "$DONNEES":"$DONNEES" -v "$D_FASTA":"$D_FASTA":ro \
      -v "$d_al":"$d_al":ro -w "$W" "$IMG" sh -c "
      bcftools mpileup -f '$FASTA' -R '$W/regions_complement.bed' -a AD,DP -q 1 -Q 13 -d 500 '$al' 2>/dev/null \
        | bcftools call -mv -f GQ -Oz -o '$W/$E.c0.vcf.gz' 2>/dev/null
      bcftools view -e 'FMT/GQ=\".\"' '$W/$E.c0.vcf.gz' | bcftools sort -Oz -o '$W/$E.c1.vcf.gz' 2>/dev/null
      printf '%s\n' '$E' > '$W/$E.nom'
      bcftools reheader -s '$W/$E.nom' '$W/$E.c1.vcf.gz' -o '$W/$E.c2.vcf.gz' && bcftools index -f -t '$W/$E.c2.vcf.gz'
      bcftools concat -a -Oz -o '$W/$E.c3.vcf.gz' '$W/$E.p.vcf.gz' '$W/$E.c2.vcf.gz' 2>/dev/null
      bcftools norm -d exact -Oz -o '$o' '$W/$E.c3.vcf.gz' 2>/dev/null && bcftools index -f -t '$o'"
    rm -f "$W/$E".p.vcf.gz* "$W/$E".c[0-3].vcf.gz* "$W/$E.nom"
    [ -s "$o.tbi" ] && echo "$E : fichier de variants pret" || echo "$E : fichier de variants en echec"
  done < "$W/genomes.txt"
fi

# ------------------------------------------------------------ lot
if etape lot; then
  : > "$W/manifeste.tsv"
  while read -r E; do
    [ -s "$DONNEES/cram/$E.cram" ] || continue
    # Un alignement sans fichier de variants pret est soumis quand meme : un
    # alignement tronque doit etre refuse par le module, et il compte au banc.
    v=$DONNEES/vcf/$E.vcf.gz; [ -s "$v" ] || v=$W/plateforme.vcf.gz
    printf '%s\t%s\t%s\n' "$E" "$DONNEES/cram/$E.cram" "$v" >> "$W/manifeste.tsv"
  done < "$W/genomes.txt"
  PGX_FASTA="$FASTA" bash "$RACINE/bin/pgx_lot.sh" --manifeste "$W/manifeste.tsv" \
    --sortie "$DONNEES/lot" --parallele "$PAR" --reprise
fi

# ------------------------------------------------------------ fiche
if etape fiche; then
  python3 "$ICI/fiche_performance.py" --lot "$DONNEES/lot" --sortie "$DONNEES/fiche"
  bash "$ICI/confronter_puce.sh" --lot "$DONNEES/lot" --travail "$W"
fi
