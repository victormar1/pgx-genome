#!/usr/bin/env bash
# Verite d'une liste de genomes HPRC R2, de bout en bout : extraction des
# haplotypes, VCF de verite et PharmCAT, puis deCYPher et Immuannot.
#
# Usage : construire_verite.sh --travail DOSSIER --fasta GRCh38.fa --echantillons LISTE
#
# Chaque etape saute ce qui est deja fait : une execution interrompue se reprend
# en relancant la meme commande. Le journal s'ouvre ici meme, et non par une
# redirection de l'appelant, qui perd des lignes sous Windows.
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="${PGX_VALIDATION:-}"; REF="${PGX_FASTA:-}"; LISTE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --travail) W="$2"; shift 2;;
    --fasta) REF="$2"; shift 2;;
    --echantillons) LISTE="$2"; shift 2;;
    -h|--help) sed -n '2,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$W" ] && [ -n "$REF" ] && [ -n "$LISTE" ] || {
  echo "usage : construire_verite.sh --travail DOSSIER --fasta GRCh38.fa --echantillons LISTE" >&2; exit 2; }

mkdir -p "$W"
exec >> "$W/construction_$(basename "$LISTE" .txt).log" 2>&1
M=$W/haplotypes/manifeste_$(basename "$LISTE" .txt).tsv
echo "debut $(date '+%F %T')"
bash "$ICI/lancer_extraction.sh" --travail "$W" --fasta "$REF" --echantillons "$LISTE"
echo "extraction finie $(date '+%F %T') : $(find "$W/haplotypes" -name '*.statuts.tsv' | wc -l) haplotypes"
bash "$ICI/pharmcat_verite.sh" --travail "$W" --fasta "$REF" --echantillons "$LISTE" | grep -v ': OK$'
echo "PharmCAT fini $(date '+%F %T') : $(find "$W/verite" -maxdepth 1 -type d -name 'pc_*' | wc -l) genomes"
bash "$ICI/typer_haplotypes.sh" --travail "$W" --manifeste "$M" | grep -v ': fini$'
echo "typage fini $(date '+%F %T') : $(find "$W/decypher" -name '*.json' | wc -l) deCYPher, $(find "$W/immuannot" -name '*.gtf.gz' | wc -l) Immuannot"
