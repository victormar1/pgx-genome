#!/usr/bin/env bash
# Rapatrie un alignement NYGC 30x par morceaux tires simultanement.
#
# Le serveur du projet bride la connexion, pas le debit : quatre flux rendent
# 6 Mo/s, huit en rendent 16 et seize en rendent 26, ou le poste plafonne. Le
# nombre de flux se choisit donc a l'appel. Mais des plages egales ne tiennent
# pas la distance,
# les connexions n'allant pas a la meme vitesse : les deux plus lentes finissent
# seules et le debit moyen retombe a 3 Mo/s. Le fichier est donc decoupe en
# morceaux courts qu'une file distribue aux travailleurs, de sorte qu'ils soient
# tous occupes jusqu'au bout.
#
# Un morceau interrompu se reprend par ce qui lui manque, jamais par -C - : la
# source coupe souvent, et une reprise interne de curl reappend ce qui est deja
# ecrit. Le fichier recolle n'est garde que si sa taille puis son MD5 repondent
# a ce que l'index annonce ; sinon il est retire et rien n'est rendu.
#
# Usage : recuperer_cram_plages.sh --donnees DOSSIER --echantillon ECH \
#                                  [--echantillons FICHIER] [--plages N] [--morceau Mo]
#   --donnees      dossier du banc ; l'alignement arrive dans cram/
#   --echantillon  un identifiant du 1000 Genomes
#   --echantillons un fichier d'identifiants, un par ligne, traites en serie
#   --plages       morceaux tires en meme temps (4 par defaut)
#   --morceau      taille d'un morceau en Mo (256 par defaut)
set -uo pipefail

# Un morceau, appele par la file. Il rend toujours 0 : un code d'echec ferait
# abandonner le lot entier a xargs, sans message. Ce sont les tailles qui disent
# ce qui manque, et la passe suivante le reprend.
#
# Le fragment arrive dans un fichier a part. Une reponse d'erreur HTTP a un
# corps : l'appendre decalerait tout le reste du morceau sans jamais changer sa
# taille finale, et c'est au MD5 du fichier entier qu'on s'en apercevrait, une
# heure plus tard.
#
# Le code et la longueur ne suffisent pas : une plage servie autre que celle
# demandee, de meme longueur, passait ces deux controles. Le fragment n'est
# recolle que si la source confirme mot pour mot la plage ET la taille totale,
# ce qui ecarte du meme coup une taille de source erronee. Tout fragment non
# retenu est nomme dans le journal, avec ce qui a ete servi.
if [ "${1:-}" = "--un-morceau" ]; then
  base=$2; pas=$3; taille=$4; url=$5; i=$6
  deb=$(( i * pas ))
  fin=$(( deb + pas - 1 )); [ "$fin" -ge "$taille" ] && fin=$(( taille - 1 ))
  p=$(printf '%s.c%05d' "$base" "$i")
  veut=$(( fin - deb + 1 ))
  have=$(stat -c %s "$p" 2> /dev/null || echo 0)
  [ "$have" -eq "$veut" ] && exit 0
  [ "$have" -gt "$veut" ] && { rm -f "$p"; have=0; }
  reste=$(( veut - have ))
  rm -f "$p.n" "$p.h"
  code=$(curl -s --max-time 1800 -r "$(( deb + have ))-$fin" \
              -D "$p.h" -o "$p.n" -w '%{http_code}' "$url")
  recu=$(stat -c %s "$p.n" 2> /dev/null || echo 0)
  servi=$(tr -d '\r' < "$p.h" 2> /dev/null \
          | awk 'tolower($1)=="content-range:"{print $2, $3}' | tail -1)
  if [ "${code:-000}" = "206" ] && [ "$recu" -gt 0 ] && [ "$recu" -le "$reste" ] \
     && [ "$servi" = "bytes $(( deb + have ))-$fin/$taille" ]; then
    cat "$p.n" >> "$p"
  elif [ "$recu" -gt 0 ]; then
    echo "morceau $i ecarte : code ${code:-000}, $recu octet(s), servi" \
         "[${servi:-aucune plage}] pour [bytes $(( deb + have ))-$fin/$taille]" >&2
  fi
  rm -f "$p.n" "$p.h"
  exit 0
fi

DONNEES=""; ECH=""; LISTE=""; N=4; MO=256
while [ $# -gt 0 ]; do
  case "$1" in
    --donnees) DONNEES="$2"; shift 2;;
    --echantillon) ECH="$2"; shift 2;;
    --echantillons) LISTE="$2"; shift 2;;
    --plages) N="$2"; shift 2;;
    --morceau) MO="$2"; shift 2;;
    -h|--help) sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$DONNEES" ] || { echo "usage : recuperer_cram_plages.sh --donnees DOSSIER --echantillon ECH" >&2; exit 2; }
[ -n "$ECH$LISTE" ] || { echo "ni --echantillon ni --echantillons" >&2; exit 2; }
case "$N" in ''|*[!0-9]*) echo "--plages attend un entier" >&2; exit 2;; esac
case "$MO" in ''|*[!0-9]*) echo "--morceau attend un entier" >&2; exit 2;; esac
[ "$N" -ge 1 ] && [ "$N" -le 16 ] || { echo "--plages hors de 1 a 16" >&2; exit 2; }
[ "$MO" -ge 16 ] && [ "$MO" -le 4096 ] || { echo "--morceau hors de 16 a 4096 Mo" >&2; exit 2; }

MOI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
W=$DONNEES/travail
mkdir -p "$DONNEES/cram" "$W"
IDX_URL=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
AWS=https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data
for f in 1000G_2504_high_coverage.sequence.index 1000G_698_related_high_coverage.sequence.index; do
  [ -s "$W/$f" ] || curl -fsS -o "$W/$f" "$IDX_URL/$f" || { echo "index introuvable : $f" >&2; exit 1; }
done

# La taille annoncee par la source. L'en-tete la donne d'ordinaire, mais une
# requete peut echouer sans raison durable : un seul echec faisait abandonner le
# genome. Le repli demande deux octets et lit la taille totale de la plage.
taille_source() {
  local u=$1 t
  t=$(curl -sI --max-time 60 "$u" | tr -d '\r' \
      | awk 'tolower($1)=="content-length:"{print $2}' | tail -1)
  case "${t:-}" in ''|*[!0-9]*) t=0;; esac
  if [ "$t" -le 0 ]; then
    t=$(curl -s --max-time 60 -r 0-1 -D - -o /dev/null "$u" | tr -d '\r' \
        | awk 'tolower($1)=="content-range:"{n=split($NF,a,"/"); print a[n]}' | tail -1)
    case "${t:-}" in ''|*[!0-9]*) t=0;; esac
  fi
  echo "$t"
}

un() {
  local E=$1 f=$DONNEES/cram/$1.cram l url md5 taille pas nb passe i have rejet
  local -a morceaux
  if [ -s "$f" ] && [ -s "$f.crai" ]; then echo "$E : deja la"; return 0; fi
  l=$(grep -h "/$E.final.cram	" "$W"/*.sequence.index | head -1)
  [ -n "$l" ] || { echo "$E : absent des index du projet"; return 1; }
  url=$(printf '%s' "$l" | cut -f1 | sed 's#^ftp://#https://#')
  md5=$(printf '%s' "$l" | cut -f2)
  if grep -q "/$E.final.cram	" "$W/1000G_2504_high_coverage.sequence.index"; then
    url="$AWS/$(basename "$(dirname "$url")")/$E.final.cram"
  fi
  taille=0
  for essai in 1 2 3; do
    taille=$(taille_source "$url")
    [ "$taille" -gt 0 ] && break
  done
  [ "$taille" -gt 0 ] || { echo "$E : taille inconnue a la source, trois essais"; return 1; }

  pas=$(( MO * 1048576 ))
  nb=$(( (taille + pas - 1) / pas ))
  morceaux=()
  for i in $(seq 0 $(( nb - 1 ))); do morceaux+=("$(printf '%s.c%05d' "$f" "$i")"); done

  for passe in $(seq 1 40); do
    seq 0 $(( nb - 1 )) | xargs -P "$N" -n 1 \
      bash "$MOI" --un-morceau "$f" "$pas" "$taille" "$url"
    have=0
    for i in "${!morceaux[@]}"; do
      have=$(( have + $(stat -c %s "${morceaux[$i]}" 2> /dev/null || echo 0) ))
    done
    [ "$have" -eq "$taille" ] && break
    echo "$E : passe $passe, $(( have / 1048576 )) Mo sur $(( taille / 1048576 ))"
  done
  if [ "$have" -ne "$taille" ]; then
    echo "$E : morceaux incomplets, $have sur $taille"; return 1
  fi

  cat "${morceaux[@]}" > "$f.part" || return 1
  have=$(stat -c %s "$f.part" 2> /dev/null || echo 0)
  [ "$have" -eq "$taille" ] || { echo "$E : recollage incomplet, $have sur $taille"; return 1; }
  rejet=$W/rejete/$E
  if [ "$(md5sum "$f.part" | cut -d' ' -f1)" != "$md5" ]; then
    echo "$E : MD5 faux, fichier retire"
    mkdir -p "$rejet"
    ( cd "$DONNEES/cram" && md5sum "$E".cram.c[0-9]* ) > "$rejet/ecarte.txt" 2> /dev/null
    rm -f "$f.part" "$f".c*
    return 1
  fi
  if [ -s "$rejet/ecarte.txt" ]; then
    ( cd "$DONNEES/cram" && md5sum "$E".cram.c[0-9]* ) > "$rejet/retenu.txt" 2> /dev/null
    echo "$E : $(diff <(sort -k2 "$rejet/ecarte.txt") <(sort -k2 "$rejet/retenu.txt") \
                 | grep -c '^>' || true) morceau(x) differaient du rapatriement ecarte"
    diff <(sort -k2 "$rejet/ecarte.txt") <(sort -k2 "$rejet/retenu.txt") \
      | sed -n 's/^> /    refait : /p'
    rm -rf "$rejet"
  fi
  rm -f "$f".c*
  mv "$f.part" "$f" || return 1
  curl -fsS -o "$f.crai" "$url.crai" || { echo "$E : index de l'alignement absent"; return 1; }
  if command -v samtools > /dev/null; then
    samtools quickcheck "$f" || { echo "$E : alignement illisible"; return 1; }
  fi
  echo "$E : alignement conforme, $((taille / 1000000000)) Go, $nb morceaux a $N de front"
}

code=0
for E in ${ECH:-} $([ -n "$LISTE" ] && tr -d '\r' < "$LISTE"); do
  un "$E" || code=1
done
exit $code
