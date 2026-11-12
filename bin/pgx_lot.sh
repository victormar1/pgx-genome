#!/bin/bash
# Traitement d'un lot d'echantillons, avec tableau de bord.
#
#   pgx_lot.sh --manifeste liste.tsv --sortie DOSSIER [--parallele N] [--reprise]
#
# Le manifeste porte trois colonnes separees par des tabulations :
#   identifiant   chemin de l'alignement   chemin du fichier de variants
#
# Chaque echantillon est traite independamment. Un echec n'arrete pas le lot ;
# il est consigne. Le tableau de bord porte une ligne par echantillon du
# manifeste, y compris ceux qui n'ont produit aucune sortie.

set -uo pipefail
RACINE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MAN=""; SORTIE=""; PAR="${PGX_PARALLELE:-4}"; REPRISE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --manifeste) MAN="$2"; shift 2;;
    --sortie) SORTIE="$2"; shift 2;;
    --parallele) PAR="$2"; shift 2;;
    --reprise) REPRISE="--reprise"; shift;;
    -h|--help) sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$MAN" ] && [ -n "$SORTIE" ] || {
  echo "usage : pgx_lot.sh --manifeste liste.tsv --sortie DOSSIER [--parallele N]" >&2; exit 2; }
[ -s "$MAN" ] || { echo "manifeste introuvable ou vide : $MAN" >&2; exit 2; }
mkdir -p "$SORTIE" || exit 2
JOURNAL="$SORTIE/lot.log"

# --------------------------------------------------------------- recevabilite
# Un identifiant sert de nom de dossier : il ne doit contenir ni separateur ni
# espace, et ne doit pas apparaitre deux fois.
python3 "$RACINE/bin/lot.py" manifeste "$MAN" \
  > "$SORTIE/manifeste_verifie.tsv" 2>"$SORTIE/manifeste_refuse.txt"
NREF=$(grep -c . "$SORTIE/manifeste_refuse.txt" 2>/dev/null); NREF=${NREF:-0}
NOK=$(grep -c . "$SORTIE/manifeste_verifie.tsv" 2>/dev/null); NOK=${NOK:-0}
echo "[$(date +%H:%M)] manifeste : $NOK recevables, $NREF refuses, $PAR en parallele" | tee -a "$JOURNAL"
[ "$NREF" -gt 0 ] && sed 's/^/  refuse : /' "$SORTIE/manifeste_refuse.txt" | tee -a "$JOURNAL"
[ "$NOK" -gt 0 ] || { echo "aucun echantillon recevable" >&2; exit 1; }

un() {
  local ligne="$1" ech cram vcf d rc
  ech=$(printf '%s' "$ligne" | cut -f1)
  cram=$(printf '%s' "$ligne" | cut -f2)
  vcf=$(printf '%s' "$ligne" | cut -f3)
  d="$SORTIE/$ech"
  # chaque echantillon ecrit dans son propre journal : en parallele, un journal
  # commun entrelace les lignes et aucune n'est attribuable
  #
  # Le module decide seul si l'echantillon est deja complet : il compare la
  # signature de ces entrees a celle du resultat present. Juger ici sur la seule
  # reussite garderait, apres correction du manifeste, le resultat d'un autre
  # fichier sous le meme identifiant.
  bash "$RACINE/bin/pgx_genome.sh" --cram "$cram" --vcf "$vcf" --sortie "$d" \
       --echantillon "$ech" --si-incomplet $REPRISE > "$SORTIE/$ech.lancement.tmp" 2>&1
  rc=$?
  if [ "$rc" = 0 ] && [ "$(head -c 12 "$SORTIE/$ech.lancement.tmp")" = "deja complet" ]; then
    # le journal du lancement qui a produit le resultat reste en place
    rm -f "$SORTIE/$ech.lancement.tmp"
    echo "[$(date +%H:%M)] [$ech] deja complet, saute" >> "$JOURNAL"
    return 0
  fi
  mv -f "$SORTIE/$ech.lancement.tmp" "$SORTIE/$ech.lancement.log"
  echo "[$(date +%H:%M)] [$ech] code $rc" >> "$JOURNAL"
}
export -f un; export SORTIE RACINE JOURNAL REPRISE

xargs -d '\n' -P "$PAR" -I{} bash -c 'un "$@"' _ {} < "$SORTIE/manifeste_verifie.tsv"

# --------------------------------------------------------------- tableau de bord
python3 "$RACINE/bin/lot.py" tableau "$SORTIE"
