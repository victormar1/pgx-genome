#!/usr/bin/env bash
# Les deux bras, genome par genome, en flux : rapatrier l'alignement, preparer
# les deux fichiers de variants, lancer le module deux fois, puis effacer
# l'alignement avant de passer au suivant.
#
# Un alignement pese 15 Go et 211 genomes en pesent 3,3 To : les garder tous
# demanderait un disque que nous n'avons pas. Le pic reste ici a deux ou trois
# alignements.
#
# Usage : lancer_bras.sh --donnees DOSSIER --travail DOSSIER --fasta GRCh38.fa \
#                        --echantillons LISTE --lot-a DOSSIER --lot-b DOSSIER \
#                        [--prepare-seulement] [--garder-cram]
#
#   --prepare-seulement  ne fait que la tranche commune de l'appel conjoint,
#                        celle qui sert a tous les genomes de la liste
#   --garder-cram        n'efface pas les alignements
#
# Chaque etape saute ce qui est deja fait : la meme commande reprend.
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="$(dirname "$ICI")"
RACINE="$(dirname "$BENCH")"
DONNEES=""; W="${PGX_VALIDATION:-}"; REF="${PGX_FASTA:-}"; LISTE=""
LOT_A=""; LOT_B=""; PREPARE=0; GARDER=0
while [ $# -gt 0 ]; do
  case "$1" in
    --donnees) DONNEES="$2"; shift 2;;
    --travail) W="$2"; shift 2;;
    --fasta) REF="$2"; shift 2;;
    --echantillons) LISTE="$2"; shift 2;;
    --lot-a) LOT_A="$2"; shift 2;;
    --lot-b) LOT_B="$2"; shift 2;;
    --prepare-seulement) PREPARE=1; shift;;
    --garder-cram) GARDER=1; shift;;
    -h|--help) sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$DONNEES" ] && [ -n "$REF" ] && [ -n "$LISTE" ] || {
  echo "usage : lancer_bras.sh --donnees D --fasta GRCh38.fa --echantillons LISTE --lot-a A --lot-b B" >&2; exit 2; }
[ "$PREPARE" = 1 ] || { [ -n "$LOT_A" ] && [ -n "$LOT_B" ] && [ -n "$W" ]; } || {
  echo "lancer_bras.sh : --lot-a, --lot-b et --travail sont requis hors --prepare-seulement" >&2; exit 2; }

T=$DONNEES/travail
mkdir -p "$T" "$DONNEES/cram" "$DONNEES/vcf" "$DONNEES/vcf_dragen"
IMG_BCF=${PGX_IMG_BCFTOOLS:-quay.io/biocontainers/bcftools:1.24--h118bc1c_2}
BCF() { docker run --rm --user "$(id -u):$(id -g)" -v "$DONNEES":"$DONNEES" -w "$T" \
          "$IMG_BCF" bcftools "$@"; }

# Tranche commune : l'appel conjoint de la plateforme se lit une seule fois pour
# toute la liste. Sans alignement sur le disque, preparer_banc.sh s'arrete juste
# apres, avant la partie par genome.
#
# La tranche est refaite si elle ne porte pas tous les genomes de la liste. Les
# tranches par chromosome, elles, sont reutilisees telles quelles par le banc :
# une tranche laissee par une liste plus courte donnerait une plateforme amputee,
# et l'extraction par genome echouerait ensuite, un genome a la fois.
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
    echo "lancer_bras.sh : tranche commune incomplete, $presents sur $voulus" >&2; exit 1; }
fi
echo "[$(date '+%F %T')] tranche commune : $presents genomes"
[ "$PREPARE" = 1 ] && { echo "tranche commune prete"; exit 0; }

un() {
  local E=$1 un_seul=$T/un_genome.txt
  if [ -s "$LOT_B/$E/sortie/$E.report.json" ] && [ -s "$LOT_A/$E/sortie/$E.report.json" ]; then
    echo "[$(date '+%F %T')] $E : deja rendu"; return 0
  fi
  bash "$ICI/recuperer_cram.sh" --donnees "$DONNEES" --echantillon "$E" || return 1
  printf '%s\n' "$E" > "$un_seul"
  # Bras A : le fichier de variants de la cohorte, prepare par le banc du depot.
  if [ ! -s "$DONNEES/vcf/$E.vcf.gz.tbi" ]; then
    bash "$BENCH/preparer_banc.sh" --donnees "$DONNEES" --fasta "$REF" \
         --genomes "$un_seul" --etapes variants || return 1
  fi
  # Le fichier doit porter ce genome et au moins autant de lignes que la tranche
  # commune : un VCF ampute ferait rendre « reference » des positions jamais lues.
  local noms lignes attendu
  noms=$(BCF query -l "$DONNEES/vcf/$E.vcf.gz" 2> /dev/null | tr -d '\r' | paste -sd,)
  lignes=$(BCF view -H "$DONNEES/vcf/$E.vcf.gz" 2> /dev/null | wc -l)
  attendu=$(BCF view -H "$T/plateforme.vcf.gz" 2> /dev/null | wc -l)
  if [ "$noms" != "$E" ] || [ "$lignes" -lt "$attendu" ]; then
    echo "$E : fichier de variants suspect, echantillon(s) [$noms], $lignes lignes pour $attendu attendues"
    return 1
  fi
  # Bras B : la reanalyse DRAGEN, lue a distance.
  bash "$ICI/vcf_dragen.sh" --donnees "$DONNEES" --echantillon "$E" || return 1
  for bras in a b; do
    local lot vcf
    if [ "$bras" = a ]; then lot=$LOT_A; vcf=$DONNEES/vcf/$E.vcf.gz
    else lot=$LOT_B; vcf=$DONNEES/vcf_dragen/$E.vcf.gz; fi
    [ -s "$lot/$E/sortie/$E.report.json" ] && continue
    mkdir -p "$lot"
    bash "$RACINE/bin/pgx_genome.sh" --cram "$DONNEES/cram/$E.cram" --vcf "$vcf" \
         --sortie "$lot/$E" --echantillon "$E" || echo "$E bras $bras : code $?"
  done
  if [ "$GARDER" = 0 ]; then
    rm -f "$DONNEES/cram/$E.cram" "$DONNEES/cram/$E.cram.crai"
  fi
  echo "[$(date '+%F %T')] $E : termine"
}

n=0
for E in $(tr -d '\r' < "$LISTE"); do
  n=$((n + 1))
  un "$E" || echo "[$(date '+%F %T')] $E : ECHEC"
done
echo "[$(date '+%F %T')] $n genomes parcourus"
