#!/usr/bin/env bash
# ABCG2 rs2231142 et POR rs1057868 contre la puce Illumina Omni 2.5 du 1000
# Genomes : une technologie independante du sequencage et de tout appeleur.
#
#   bench/confronter_puce.sh --lot DOSSIER --travail DOSSIER
#
# La puce est en GRCh37 : les deux positions y sont lues par leurs coordonnees
# GRCh37, sur le meme brin que GRCh38.
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOT=""; W=""
while [ $# -gt 0 ]; do
  case "$1" in
    --lot) LOT="$2"; shift 2;;
    --travail) W="$2"; shift 2;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$LOT" ] && [ -n "$W" ] || { echo "usage : confronter_puce.sh --lot DOSSIER --travail DOSSIER" >&2; exit 2; }
mkdir -p "$W"
IMG=${PGX_IMG_BCFTOOLS:-quay.io/biocontainers/bcftools:1.24--h118bc1c_2}
O=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
docker run --rm --user "$(id -u):$(id -g)" -v "$W":"$W" -w "$W" "$IMG" sh -c "
  bcftools query -r 4:89052323,7:75615006 -f '%CHROM:%POS[\t%GT]\n' '$O' > '$W/puce.tsv' &&
  bcftools query -l '$O' > '$W/puce_echantillons.txt'" || { echo "lecture de la puce en echec" >&2; exit 1; }
python3 "$ICI/confronter_puce.py" --lot "$LOT" --puce "$W/puce.tsv" --echantillons "$W/puce_echantillons.txt"
