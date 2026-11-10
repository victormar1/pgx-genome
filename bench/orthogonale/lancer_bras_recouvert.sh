#!/usr/bin/env bash
# Les deux bras, genome par genome, le rapatriement du suivant recouvrant le
# calcul du courant.
#
# Un genome prend 8 min de reseau et 15 de calcul. Le genome suivant arrive
# pendant le calcul, et le pic disque passe a deux alignements.
#
# Les deux bras sont tentes separement.
# Un genome qui n'a qu'un bras sort de la comparaison entre bras, pas de la
# mesure de concordance de son bras.
#
# Usage : lancer_bras_recouvert.sh --donnees DOSSIER --travail DOSSIER \
#            --fasta GRCh38.fa --echantillons LISTE --lot-a DOSSIER --lot-b DOSSIER \
#            [--plages N] [--prepare-seulement] [--garder-cram]
#
#   --plages             plages simultanees par alignement (4 par defaut)
#   --prepare-seulement  ne fait que la tranche commune de l'appel conjoint
#   --garder-cram        n'efface pas les alignements
#
# Chaque etape saute ce qui est deja fait : la meme commande reprend.
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="$(dirname "$ICI")"
RACINE="$(dirname "$BENCH")"
DONNEES=""; W="${PGX_VALIDATION:-}"; REF="${PGX_FASTA:-}"; LISTE=""
LOT_A=""; LOT_B=""; PREPARE=0; GARDER=0; PLAGES=4
while [ $# -gt 0 ]; do
  case "$1" in
    --donnees) DONNEES="$2"; shift 2;;
    --travail) W="$2"; shift 2;;
    --fasta) REF="$2"; shift 2;;
    --echantillons) LISTE="$2"; shift 2;;
    --lot-a) LOT_A="$2"; shift 2;;
    --lot-b) LOT_B="$2"; shift 2;;
    --plages) PLAGES="$2"; shift 2;;
    --prepare-seulement) PREPARE=1; shift;;
    --garder-cram) GARDER=1; shift;;
    -h|--help) sed -n '2,23p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$DONNEES" ] && [ -n "$REF" ] && [ -n "$LISTE" ] || {
  echo "usage : lancer_bras_recouvert.sh --donnees D --fasta GRCh38.fa --echantillons LISTE --lot-a A --lot-b B" >&2; exit 2; }
[ "$PREPARE" = 1 ] || { [ -n "$LOT_A" ] && [ -n "$LOT_B" ] && [ -n "$W" ]; } || {
  echo "lancer_bras_recouvert.sh : --lot-a, --lot-b et --travail sont requis hors --prepare-seulement" >&2; exit 2; }

T=$DONNEES/travail
mkdir -p "$T" "$DONNEES/cram" "$DONNEES/vcf" "$DONNEES/vcf_dragen"
IMG_BCF=${PGX_IMG_BCFTOOLS:-quay.io/biocontainers/bcftools:1.24--h118bc1c_2}
BCF() { docker run --rm --user "$(id -u):$(id -g)" -v "$DONNEES":"$DONNEES" -w "$T" \
          "$IMG_BCF" bcftools "$@"; }

# Tranche commune : l'appel conjoint de la plateforme se lit une seule fois pour
# toute la liste. Elle est refaite si elle ne porte pas tous les genomes de la
# liste : une tranche laissee par une liste plus courte donnerait une plateforme
# amputee, et l'extraction par genome echouerait ensuite un genome a la fois.
voulus=$(tr -d '\r' < "$LISTE" | sed '/^$/d' | sort -u | wc -l)
presents=0
if [ -s "$T/plateforme.vcf.gz.tbi" ]; then
  presents=$(BCF query -l "$T/plateforme.vcf.gz" 2> /dev/null | wc -l)
fi
if [ "$presents" -lt "$voulus" ]; then
  echo "[$(date '+%F %T')] tranche commune : $presents genome(s) sur $voulus, a refaire"
  rm -f "$T"/br_chr*.vcf.gz "$T"/br_chr*.vcf.gz.tbi "$T"/reg_chr*.bed \
        "$T"/tout.vcf.gz "$T"/tout.vcf.gz.tbi "$T"/plateforme.vcf.gz "$T"/plateforme.vcf.gz.tbi
  bash "$BENCH/preparer_banc.sh" --donnees "$DONNEES" --fasta "$REF" \
       --genomes "$LISTE" --etapes variants || exit 1
  presents=$(BCF query -l "$T/plateforme.vcf.gz" 2> /dev/null | wc -l)
  [ "$presents" -ge "$voulus" ] || {
    echo "lancer_bras_recouvert.sh : tranche commune incomplete, $presents sur $voulus" >&2; exit 1; }
fi
echo "[$(date '+%F %T')] tranche commune : $presents genomes"
[ "$PREPARE" = 1 ] && { echo "tranche commune prete"; exit 0; }

rendu() { [ -s "$1/$2/sortie/$2.report.json" ]; }

PRE_PID=""; PRE_ECH=""
precharger() {
  local E=${1:-}
  [ -n "$E" ] || return 0
  rendu "$LOT_A" "$E" && rendu "$LOT_B" "$E" && return 0
  bash "$ICI/recuperer_cram_plages.sh" --donnees "$DONNEES" --echantillon "$E" \
       --plages "$PLAGES" > "$T/$E.rapatriement.log" 2>&1 &
  PRE_PID=$!; PRE_ECH=$E
}
attendre() {
  [ -n "$PRE_PID" ] && [ "$PRE_ECH" = "${1:-}" ] || return 0
  wait "$PRE_PID"
  local c=$?
  cat "$T/$PRE_ECH.rapatriement.log" 2> /dev/null
  PRE_PID=""; PRE_ECH=""
  return $c
}

# Le fichier de variants doit porter ce genome et au moins autant de lignes que
# la tranche dont il descend : un VCF ampute fait rendre « reference » des
# positions jamais lues.
variants_bras_a() {
  local E=$1 un_seul=$T/un_$1.txt noms lignes attendu
  printf '%s\n' "$E" > "$un_seul"
  if [ ! -s "$DONNEES/vcf/$E.vcf.gz.tbi" ]; then
    bash "$BENCH/preparer_banc.sh" --donnees "$DONNEES" --fasta "$REF" \
         --genomes "$un_seul" --etapes variants || return 1
  fi
  noms=$(BCF query -l "$DONNEES/vcf/$E.vcf.gz" 2> /dev/null | tr -d '\r' | paste -sd,)
  lignes=$(BCF view -H "$DONNEES/vcf/$E.vcf.gz" 2> /dev/null | wc -l)
  attendu=$(BCF view -H "$T/plateforme.vcf.gz" 2> /dev/null | wc -l)
  if [ "$noms" != "$E" ] || [ "$lignes" -lt "$attendu" ]; then
    echo "$E : fichier de variants suspect, echantillon(s) [$noms], $lignes lignes pour $attendu attendues"
    return 1
  fi
}

module() {
  local E=$1 lot=$2 vcf=$3
  rendu "$lot" "$E" && return 0
  mkdir -p "$lot"
  bash "$RACINE/bin/pgx_genome.sh" --cram "$DONNEES/cram/$E.cram" --vcf "$vcf" \
       --sortie "$lot/$E" --echantillon "$E" || echo "$E bras $lot : code $?"
  rendu "$lot" "$E"
}

un() {
  local E=$1 a=1 b=1
  if rendu "$LOT_A" "$E" && rendu "$LOT_B" "$E"; then
    echo "[$(date '+%F %T')] $E : deja rendu"; return 0
  fi
  if variants_bras_a "$E"; then
    module "$E" "$LOT_A" "$DONNEES/vcf/$E.vcf.gz" && a=0
  else
    echo "$E : bras A sans entree"
  fi
  if bash "$ICI/vcf_dragen.sh" --donnees "$DONNEES" --echantillon "$E"; then
    module "$E" "$LOT_B" "$DONNEES/vcf_dragen/$E.vcf.gz" && b=0
  else
    echo "$E : bras B sans entree"
  fi
  [ "$GARDER" = 0 ] && rm -f "$DONNEES/cram/$E.cram" "$DONNEES/cram/$E.cram.crai"
  echo "[$(date '+%F %T')] $E : bras A $( [ $a = 0 ] && echo rendu || echo manquant ), bras B $( [ $b = 0 ] && echo rendu || echo manquant )"
  [ "$a" = 0 ] || [ "$b" = 0 ]
}

mapfile -t ECHS < <(tr -d '\r' < "$LISTE" | sed '/^$/d')
nb=${#ECHS[@]}
echo "[$(date '+%F %T')] $nb genomes, $PLAGES plages par alignement"
for k in $(seq 0 $(( nb - 1 ))); do
  E=${ECHS[$k]}
  attendre "$E" || echo "$E : rapatriement en echec au prechargement"
  if ! rendu "$LOT_A" "$E" || ! rendu "$LOT_B" "$E"; then
    bash "$ICI/recuperer_cram_plages.sh" --donnees "$DONNEES" --echantillon "$E" \
         --plages "$PLAGES" || { echo "[$(date '+%F %T')] $E : ECHEC de rapatriement"; continue; }
  fi
  [ $(( k + 1 )) -lt "$nb" ] && precharger "${ECHS[$(( k + 1 ))]}"
  un "$E" || echo "[$(date '+%F %T')] $E : ECHEC"
done
attendre "${PRE_ECH:-}" || true
echo "[$(date '+%F %T')] $nb genomes parcourus"
