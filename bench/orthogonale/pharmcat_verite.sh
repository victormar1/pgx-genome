#!/usr/bin/env bash
# Vérité des genes a variants ponctuels : un VCF phase par genome, tire des deux
# haplotypes d'assemblage, puis PharmCAT avec le pretraitement et les options du
# module. Aucun appel exterieur n'est fait sur ce VCF.
#
# Usage : pharmcat_verite.sh --travail DOSSIER --fasta GRCh38.fa --echantillons LISTE
#
# Image : PGX_IMG_PHARMCAT, par defaut celle du module.
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="${PGX_VALIDATION:-}"; REF="${PGX_FASTA:-}"; LISTE=""
IMG="${PGX_IMG_PHARMCAT:-pgkb/pharmcat:3.4.0}"
while [ $# -gt 0 ]; do
  case "$1" in
    --travail) W="$2"; shift 2;;
    --fasta) REF="$2"; shift 2;;
    --echantillons) LISTE="$2"; shift 2;;
    -h|--help) sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$W" ] && [ -n "$REF" ] && [ -n "$LISTE" ] || {
  echo "usage : pharmcat_verite.sh --travail DOSSIER --fasta GRCh38.fa --echantillons LISTE" >&2; exit 2; }
[ -s "$REF.fai" ] || { echo "reference non indexee : $REF" >&2; exit 2; }
RES=$ICI/../../ressources

mkdir -p "$W/verite"
cd "$W/verite" || exit 1
for S in $(tr -d '\r' < "$LISTE"); do
  H=$W/haplotypes/$S
  if [ ! -s "$H.h1.statuts.tsv" ] || [ ! -s "$H.h2.statuts.tsv" ]; then echo "$S : haplotype manquant"; continue; fi
  python3 "$ICI/vcf_verite.py" "$S" "$H.h1" "$H.h2" "$REF.fai" "$W/verite/$S.vcf" \
    "$RES/pharmcat_positions.vcf" "$RES/rnpgx_complement.vcf" || { echo "$S : VCF en echec"; continue; }
  D=$W/verite/pc_$S; rm -rf "$D"; mkdir -p "$D"
  docker run --rm -u "$(id -u):$(id -g)" -v "$W/verite:$W/verite" -w "$D" "$IMG" \
    /pharmcat/pharmcat_vcf_preprocessor -vcf "$W/verite/$S.vcf" -o "$D" --absent-to-ref \
    > "$D/pretraitement.log" 2>&1
  if [ ! -s "$D/$S.preprocessed.vcf.bgz" ]; then echo "$S : pretraitement en echec"; continue; fi
  docker run --rm -u "$(id -u):$(id -g)" -v "$W/verite:$W/verite" -w "$D" "$IMG" \
    java -jar /pharmcat/pharmcat.jar -vcf "$D/$S.preprocessed.vcf.bgz" -o "$D" -bf "$S" -reporterJson \
    > "$D/pharmcat.log" 2>&1
  [ -s "$D/$S.report.json" ] && echo "$S : OK" || echo "$S : PharmCAT en echec"
done
