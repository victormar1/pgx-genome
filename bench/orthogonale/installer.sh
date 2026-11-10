#!/usr/bin/env bash
# Installe, dans DOSSIER/outils, les deux outils tiers dont la verite par
# assemblages a besoin, a la version exacte qui a servi, et leurs donnees :
#   deCYPher   CYP2D6 et huit autres genes, depuis la sequence d'un haplotype
#   Immuannot  HLA, depuis la sequence du locus
# Plus k8 et paftools.js, que deCYPher appelle.
#
# Rien n'est ecrit hors de DOSSIER. Compter 400 Mo pour les references de
# deCYPher et 11 Mo pour les donnees IPD d'Immuannot.
#
# Usage : installer.sh --travail DOSSIER
set -euo pipefail
W="${PGX_VALIDATION:-}"
while [ $# -gt 0 ]; do
  case "$1" in
    --travail) W="$2"; shift 2;;
    -h|--help) sed -n '2,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$W" ] || { echo "usage : installer.sh --travail DOSSIER" >&2; exit 2; }
O=$W/outils
mkdir -p "$O/bin"
cd "$O"

# Les deux depots sont epingles sur le commit qui a servi. core.autocrlf est
# desactive a la source : un script bash arrive avec des fins de ligne Windows
# echoue sur un \r invisible.
depot() {
  local nom=$1 url=$2 commit=$3
  if [ ! -d "$nom/.git" ]; then
    git -c core.autocrlf=false clone -q "$url" "$nom"
  fi
  git -C "$nom" -c core.autocrlf=false fetch -q origin "$commit" 2> /dev/null || git -C "$nom" fetch -q
  git -C "$nom" config core.autocrlf false
  git -C "$nom" -c advice.detachedHead=false checkout -q "$commit"
  git -C "$nom" checkout-index -f -a
  echo "$nom : $(git -C "$nom" rev-parse --short HEAD)"
}
depot deCYPher https://github.com/HolidayChang/deCYPher.git 54a0f9d0150518d7a50c6fdd553dd5b6e3dff257
depot Immuannot https://github.com/YingZhou-Lab/Immuannot.git 31362c370123f6706a04ffdda6a3bf0c108328e0

if [ ! -x bin/k8 ]; then
  curl -fsSL https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 -o k8.tar.bz2
  tar -jxf k8.tar.bz2 && cp k8-0.2.4/k8-Linux bin/k8 && chmod +x bin/k8 && rm -rf k8-0.2.4 k8.tar.bz2
fi
# paftools.js de la version de minimap2 installee (2.24), pas la branche master.
[ -s bin/paftools.js ] || curl -fsSL https://raw.githubusercontent.com/lh3/minimap2/v2.24/misc/paftools.js -o bin/paftools.js

if [ ! -d deCYPher/ref ]; then
  curl -fsSL https://github.com/HolidayChang/deCYPher/releases/download/v0.1.0/ref.tar.gz -o deCYPher/ref.tar.gz
  tar -zxf deCYPher/ref.tar.gz -C deCYPher && rm deCYPher/ref.tar.gz
fi

if [ ! -d Immuannot/Data-2024Feb02 ]; then
  curl -fsSL https://zenodo.org/api/records/10948964/files/Data-2024Feb02.tar.gz/content -o Immuannot/data.tar.gz
  tar -zxf Immuannot/data.tar.gz -C Immuannot && rm Immuannot/data.tar.gz
fi

echo "k8 : $(bin/k8 -v 2>&1 | head -1)"
echo "paftools : $(bin/k8 bin/paftools.js version 2>&1 | head -1)"
echo "minimap2 : $(minimap2 --version 2>&1 | head -1)"
echo "deCYPher ref : $(ls deCYPher/ref | tr '\n' ' ')"
echo "Immuannot : $(ls Immuannot | tr '\n' ' ')"
du -sh deCYPher/ref Immuannot/Data-2024Feb02 2>/dev/null || true
