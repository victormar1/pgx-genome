#!/usr/bin/env bash
# Rapatrie l'alignement NYGC 30x d'un genome du 1000 Genomes.
#
# Le miroir AWS sert environ huit fois plus vite que le serveur du projet, mais
# il n'heberge que le collectif de 2504 genomes : les apparentes passent par
# l'URL de l'index. Dans les deux cas le fichier n'est garde que si sa taille
# puis son MD5 repondent a ce que l'index annonce.
#
# Usage : recuperer_cram.sh --donnees DOSSIER --echantillon ECH [--echantillons FICHIER]
#   --donnees      dossier du banc ; l'alignement arrive dans cram/
#   --echantillon  un identifiant du 1000 Genomes
#   --echantillons un fichier d'identifiants, un par ligne, traites en serie
set -uo pipefail
DONNEES=""; ECH=""; LISTE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --donnees) DONNEES="$2"; shift 2;;
    --echantillon) ECH="$2"; shift 2;;
    --echantillons) LISTE="$2"; shift 2;;
    -h|--help) sed -n '2,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$DONNEES" ] || { echo "usage : recuperer_cram.sh --donnees DOSSIER --echantillon ECH" >&2; exit 2; }
[ -n "$ECH$LISTE" ] || { echo "ni --echantillon ni --echantillons" >&2; exit 2; }

W=$DONNEES/travail
mkdir -p "$DONNEES/cram" "$W"
IDX_URL=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
AWS=https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data
for f in 1000G_2504_high_coverage.sequence.index 1000G_698_related_high_coverage.sequence.index; do
  [ -s "$W/$f" ] || curl -fsS -o "$W/$f" "$IDX_URL/$f" || { echo "index introuvable : $f" >&2; exit 1; }
done

un() {
  local E=$1 f=$DONNEES/cram/$1.cram l url md5 taille have passe
  if [ -s "$f" ] && [ -s "$f.crai" ]; then echo "$E : deja la"; return 0; fi
  l=$(grep -h "/$E.final.cram	" "$W"/*.sequence.index | head -1)
  [ -n "$l" ] || { echo "$E : absent des index du projet"; return 1; }
  url=$(printf '%s' "$l" | cut -f1 | sed 's#^ftp://#https://#')
  md5=$(printf '%s' "$l" | cut -f2)
  # Le collectif de 2504 est sur AWS, sous le numero de passage de l'index.
  if grep -q "/$E.final.cram	" "$W/1000G_2504_high_coverage.sequence.index"; then
    url="$AWS/$(basename "$(dirname "$url")")/$E.final.cram"
  fi
  taille=$(curl -sI "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}' | tail -1)
  [ -n "$taille" ] || { echo "$E : taille inconnue a la source"; return 1; }
  # Une passe de reprise par tour, sans --retry : la source coupe, et un retry
  # interne de curl reappend ce qui est deja ecrit.
  for passe in $(seq 1 40); do
    have=$(stat -c %s "$f.part" 2>/dev/null || echo 0)
    [ "$have" -eq "$taille" ] && break
    [ "$have" -gt "$taille" ] && rm -f "$f.part"
    curl -s -C - -o "$f.part" "$url"
  done
  have=$(stat -c %s "$f.part" 2>/dev/null || echo 0)
  [ "$have" -eq "$taille" ] || { echo "$E : taille incomplete, $have sur $taille"; return 1; }
  if [ "$(md5sum "$f.part" | cut -d' ' -f1)" != "$md5" ]; then
    echo "$E : MD5 faux, fichier retire"; rm -f "$f.part"; return 1
  fi
  mv "$f.part" "$f" || return 1
  curl -fsS -o "$f.crai" "$url.crai" || { echo "$E : index de l'alignement absent"; return 1; }
  if command -v samtools > /dev/null; then
    samtools quickcheck "$f" || { echo "$E : alignement illisible"; return 1; }
  fi
  echo "$E : alignement conforme, $((taille / 1000000000)) Go"
}

code=0
for E in ${ECH:-} $([ -n "$LISTE" ] && tr -d '\r' < "$LISTE"); do
  un "$E" || code=1
done
exit $code
