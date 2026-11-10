#!/usr/bin/env bash
# Fichier de variants du bras B : la reanalyse DRAGEN 3.7.6 du 1000 Genomes,
# lue a distance et restreinte au perimetre du module.
#
# Le VCF de la reanalyse est un fichier par genome, indexe : seules les regions
# utiles descendent, quelques mega-octets au lieu de 545. Le FILTER est garde tel
# quel, la decision d'exclure appartient au module.
#
# Usage : vcf_dragen.sh --donnees DOSSIER --echantillon ECH [--echantillons FICHIER]
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RACINE="$(cd "$ICI/../.." && pwd)"
DONNEES=""; ECH=""; LISTE=""
IMG=${PGX_IMG_BCFTOOLS:-quay.io/biocontainers/bcftools:1.24--h118bc1c_2}
SOURCE=https://1000genomes-dragen.s3.amazonaws.com/data/dragen-3.7.6/hg38-graph-based
while [ $# -gt 0 ]; do
  case "$1" in
    --donnees) DONNEES="$2"; shift 2;;
    --echantillon) ECH="$2"; shift 2;;
    --echantillons) LISTE="$2"; shift 2;;
    -h|--help) sed -n '2,10p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$DONNEES" ] || { echo "usage : vcf_dragen.sh --donnees DOSSIER --echantillon ECH" >&2; exit 2; }
[ -n "$ECH$LISTE" ] || { echo "ni --echantillon ni --echantillons" >&2; exit 2; }

W=$DONNEES/travail
mkdir -p "$DONNEES/vcf_dragen" "$W"
BCF() { docker run --rm --user "$(id -u):$(id -g)" -v "$DONNEES":"$DONNEES" \
          -v "$RACINE":"$RACINE":ro -w "$W" "$IMG" bcftools "$@"; }

# Perimetre : les trois BED du module, fusionnes.
REG=$W/regions_dragen.bed
if [ ! -s "$REG" ]; then
  cut -f1-3 "$RACINE/ressources/pharmcat_positions.bed" "$RACINE/ressources/pypgx_regions.bed" \
            "$RACINE/ressources/rnpgx_complement.bed" | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"} {if (c==$1 && $2<=e) { if($3>e) e=$3 } else { if(c!="") print c,s,e; c=$1; s=$2; e=$3 }} END{if(c!="") print c,s,e}' \
    > "$REG"
fi

un() {
  local E=$1 o=$DONNEES/vcf_dragen/$1.vcf.gz
  [ -s "$o.tbi" ] && { echo "$E : deja la"; return 0; }
  if ! BCF view -R "$REG" -Oz -o "$o" "$SOURCE/$E/$E.hard-filtered.vcf.gz" 2> "$W/$E.dragen.log"; then
    echo "$E : extraction DRAGEN en echec (voir travail/$E.dragen.log)"; rm -f "$o"; return 1
  fi
  BCF index -f -t "$o" || return 1
  echo "$E : $(BCF view -H "$o" | wc -l) variants dans le perimetre"
}

code=0
for E in ${ECH:-} $([ -n "$LISTE" ] && tr -d '\r' < "$LISTE"); do
  un "$E" || code=1
done
exit $code
