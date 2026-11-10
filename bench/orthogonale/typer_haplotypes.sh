#!/usr/bin/env bash
# Verite des genes que les variants ponctuels ne suffisent pas a decider :
# deCYPher pour CYP2D6 et huit autres genes, Immuannot pour le HLA, sur la
# sequence de chaque locus d'assemblage.
#
# Usage : typer_haplotypes.sh --travail DOSSIER --manifeste M.tsv [--parallele 5]
#
# Les deux outils sont installes par installer.sh dans DOSSIER/outils.
set -uo pipefail
ICI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="${PGX_VALIDATION:-}"; MANIFESTE=""; PAR=5
while [ $# -gt 0 ]; do
  case "$1" in
    --travail) W="$2"; shift 2;;
    --manifeste) MANIFESTE="$2"; shift 2;;
    --parallele) PAR="$2"; shift 2;;
    -h|--help) sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "option inconnue : $1" >&2; exit 2;;
  esac
done
[ -n "$W" ] && [ -n "$MANIFESTE" ] || {
  echo "usage : typer_haplotypes.sh --travail DOSSIER --manifeste M.tsv" >&2; exit 2; }
O=$W/outils
[ -x "$O/bin/k8" ] || { echo "outils absents : lancer installer.sh" >&2; exit 2; }
# deCYPher cherche ses dependances sous CONDA_PREFIX.
export CONDA_PREFIX=$O PATH=$O/bin:$PATH
mkdir -p "$W/decypher" "$W/immuannot"

une() {
  local p=$1 fa=$W/haplotypes/$1.loci.fa
  if [ ! -s "$W/decypher/$p/$p.json" ]; then
    python3 "$O/deCYPher/deCYPher.py" --input "$fa" --name "$p" --outdir "$W/decypher" --threads 4 \
      > "$W/decypher/$p.log" 2>&1 || echo "$p : deCYPher en echec"
  fi
  if [ ! -s "$W/immuannot/$p.gtf.gz" ]; then
    python3 "$ICI/extraire_hla.py" "$fa" "$W/immuannot/$p.hla.fa" || { echo "$p : HLA non extrait"; return 1; }
    ( cd "$W/immuannot" && bash "$O/Immuannot/scripts.pub.v3/immuannot.sh" -c "$p.hla.fa" \
        -r "$O/Immuannot/Data-2024Feb02" -o "$p" -t 4 > "$p.log" 2>&1 ) || echo "$p : Immuannot en echec"
  fi
  echo "$p : fini"
}

for p in $(cut -f1 "$MANIFESTE"); do
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 2; done
  une "$p" &
done
wait
echo "termine"
