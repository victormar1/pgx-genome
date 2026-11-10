#!/usr/bin/env bash
# Extrait les haplotypes d'une liste de genomes HPRC R2 : pour chaque region du
# perimetre, le statut de couverture et les variants portes par l'haplotype, et
# pour chaque locus sa sequence.
#
# Usage : lancer_extraction.sh --travail DOSSIER --fasta GRCh38.fa \
#                              --echantillons LISTE [--parallele 4]
#
# Les URL sont lues a distance : rien n'est rapatrie en entier. Un haplotype
# deja extrait est saute, ce qui rend la reprise sure.
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="${PGX_VALIDATION:-}"; REF="${PGX_FASTA:-}"; LISTE=""; PAR=4
while [ $# -gt 0 ]; do
  case "$1" in
    --travail) W="$2"; shift 2;;
    --fasta) REF="$2"; shift 2;;
    --echantillons) LISTE="$2"; shift 2;;
    --parallele) PAR="$2"; shift 2;;
    -h|--help) sed -n '2,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$W" ] && [ -n "$REF" ] && [ -n "$LISTE" ] || {
  echo "usage : lancer_extraction.sh --travail DOSSIER --fasta GRCh38.fa --echantillons LISTE" >&2; exit 2; }
[ -s "$REF.fai" ] || { echo "reference non indexee : $REF" >&2; exit 2; }

mkdir -p "$W/haplotypes" "$W/tmp"
MANIFESTE=$W/haplotypes/manifeste_$(basename "$LISTE" .txt).tsv
python3 "$ICI/manifeste.py" "$LISTE" "$ICI/meta/assemblies_r2.csv" \
        "$ICI/meta/aln_grch38_r2.csv" > "$MANIFESTE" || exit 1
echo "$(wc -l < "$MANIFESTE") haplotypes dans $MANIFESTE"

# Le perimetre des positions se refait a partir des deux BED du module, pour
# qu'il ne puisse pas s'en ecarter en silence.
cut -f1-3 "$ICI/../../ressources/pharmcat_positions.bed" \
          "$ICI/../../ressources/rnpgx_complement.bed" \
  | sort -k1,1 -k2,2n > "$W/tmp/positions.bed"

cd "$W/tmp" || exit 1
une() {
  local p=$1 bam=$2 fa=$3
  [ -s "$W/haplotypes/$p.statuts.tsv" ] && { echo "$p deja fait"; return 0; }
  python3 "$ICI/extraire_haplotype.py" --bam "$bam" --fasta "$fa" \
    --positions "$W/tmp/positions.bed" --loci "$ICI/ressources/loci.bed" \
    --ref "$REF" --prefixe "$W/haplotypes/$p" 2>> "$W/haplotypes/$p.log" \
    || echo "$p ECHEC (voir $p.log)"
  tail -1 "$W/haplotypes/$p.log"
}
while IFS=$'\t' read -r p bam fa; do
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 2; done
  une "$p" "$bam" "$fa" &
done < "$MANIFESTE"
wait
echo "termine"
