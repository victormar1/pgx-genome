#!/bin/bash
# Traitement d'un lot d'echantillons, avec tableau de bord.
#
#   pgx_lot.sh --manifeste liste.tsv --sortie DOSSIER [--parallele 4] [--reprise]
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
python3 - "$MAN" > "$SORTIE/manifeste_verifie.tsv" 2>"$SORTIE/manifeste_refuse.txt" <<'PY'
import csv, os, re, sys
src = sys.argv[1]
vus, sortie, refus = set(), [], []
for i, l in enumerate(open(src, encoding="utf-8"), 1):
    l = l.rstrip("\n")
    if not l.strip() or l.startswith("#"):
        continue
    c = l.split("\t")
    if len(c) < 3:
        refus.append("ligne %d : moins de trois colonnes" % i); continue
    ech, cram, vcf = c[0].strip(), c[1].strip(), c[2].strip()
    if not re.fullmatch(r"[A-Za-z0-9._-]{1,64}", ech):
        refus.append("ligne %d : identifiant inutilisable comme nom de dossier : %r" % (i, ech)); continue
    if ech in vus:
        refus.append("ligne %d : identifiant en double : %s" % (i, ech)); continue
    if not os.path.exists(cram):
        refus.append("ligne %d : alignement introuvable : %s" % (i, cram)); continue
    if not os.path.exists(vcf):
        refus.append("ligne %d : fichier de variants introuvable : %s" % (i, vcf)); continue
    vus.add(ech)
    sortie.append((ech, cram, vcf))
for r in sortie:
    print("\t".join(r))
for r in refus:
    print(r, file=sys.stderr)
PY
NREF=$(grep -c . "$SORTIE/manifeste_refuse.txt" 2>/dev/null); NREF=${NREF:-0}
NOK=$(grep -c . "$SORTIE/manifeste_verifie.tsv" 2>/dev/null); NOK=${NOK:-0}
echo "[$(date +%H:%M)] manifeste : $NOK recevables, $NREF refuses, $PAR en parallele" | tee -a "$JOURNAL"
[ "$NREF" -gt 0 ] && sed 's/^/  refuse : /' "$SORTIE/manifeste_refuse.txt" | tee -a "$JOURNAL"
[ "$NOK" -gt 0 ] || { echo "aucun echantillon recevable" >&2; exit 1; }

un() {
  local ligne="$1" ech cram vcf d
  ech=$(printf '%s' "$ligne" | cut -f1)
  cram=$(printf '%s' "$ligne" | cut -f2)
  vcf=$(printf '%s' "$ligne" | cut -f3)
  d="$SORTIE/$ech"
  if [ -s "$d/sortie/provenance.json" ] \
     && python3 -c "import json,sys;sys.exit(0 if json.load(open(r'$d/sortie/provenance.json',encoding='utf-8'))['reussite_complete'] else 1)"; then
    echo "[$(date +%H:%M)] [$ech] deja complet, saute" >> "$JOURNAL"
    return 0
  fi
  # chaque echantillon ecrit dans son propre journal : en parallele, un journal
  # commun entrelace les lignes et aucune n'est attribuable
  bash "$RACINE/bin/pgx_genome.sh" --cram "$cram" --vcf "$vcf" --sortie "$d" \
       --echantillon "$ech" $REPRISE > "$SORTIE/$ech.lancement.log" 2>&1
  echo "[$(date +%H:%M)] [$ech] code $?" >> "$JOURNAL"
}
export -f un; export SORTIE RACINE JOURNAL REPRISE

xargs -d '\n' -P "$PAR" -I{} bash -c 'un "$@"' _ {} < "$SORTIE/manifeste_verifie.tsv"

# --------------------------------------------------------------- tableau de bord
python3 - "$SORTIE" <<'PY'
import csv, json, os, sys
sortie = sys.argv[1]
ATTENDUS = ["recevabilite", "filtre", "tranche", "qc", "cyp2d6", "hla",
            "appels", "pharmcat", "rendu"]
lignes = []
for l in open(os.path.join(sortie, "manifeste_verifie.tsv"), encoding="utf-8"):
    ech = l.split("\t")[0].strip()
    if not ech:
        continue
    p = os.path.join(sortie, ech, "sortie", "provenance.json")
    ligne = {"echantillon": ech}
    if not os.path.exists(p):
        ligne["reussite_complete"] = False
        ligne["date"] = ""
        for e in ATTENDUS:
            ligne[e] = "AUCUNE SORTIE"
            ligne[e + "_s"] = ""
    else:
        d = json.load(open(p, encoding="utf-8"))
        ligne["reussite_complete"] = d["reussite_complete"]
        ligne["date"] = d["date"]
        for e in ATTENDUS:
            v = d.get("etages", {}).get(e)
            ligne[e] = v["etat"] if v else "ABSENT"
            ligne[e + "_s"] = v["secondes"] if v else ""
    lignes.append(ligne)

cols = ["echantillon", "reussite_complete", "date"] + [x for e in ATTENDUS for x in (e, e + "_s")]
p = os.path.join(sortie, "tableau_de_bord.tsv")
with open(p, "w", newline="\n", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(lignes)

ok = sum(1 for l in lignes if l["reussite_complete"])
print("echantillons du manifeste : %d, reussite complete : %d" % (len(lignes), ok))
for e in ATTENDUS:
    ko = [l["echantillon"] for l in lignes if l.get(e) not in ("OK", "SANS_RESULTAT")]
    sr = [l["echantillon"] for l in lignes if l.get(e) == "SANS_RESULTAT"]
    if ko:
        print("  etage %-13s en echec chez %d : %s%s"
              % (e, len(ko), " ".join(ko[:8]), " …" if len(ko) > 8 else ""))
    if sr:
        # l etage a tourne jusqu au bout sans rien conclure : le gene est rendu
        # non analyse, le module n a pas failli
        print("  etage %-13s sans resultat chez %d : %s%s"
              % (e, len(sr), " ".join(sr[:8]), " …" if len(sr) > 8 else ""))
print("tableau de bord :", p)
PY
