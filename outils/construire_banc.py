# -*- coding: utf-8 -*-
"""Reconstruit la liste des prelevements du banc depuis le lot de validation.

Une liste tenue a la main finit par ne plus correspondre au lot : elle comptait
244 prelevements quand le lot en portait 256. Ici, chaque dossier du lot donne
une ligne, et son statut vient de sa trace d'execution ; les colonnes de verite
disent si le prelevement figure dans la table GeT-RM et dans le typage HLA de
reference.

Usage :
    construire_banc.py --lot DOSSIER --getrm verite.tsv --hla verite_hla.tsv \
        --populations metadonnees.tsv --sortie bench/echantillons.csv

Le lot est un dossier par prelevement, chacun avec sortie/provenance.json. Les
trois tables sont au format tabule, la premiere colonne portant l'identifiant.
"""
import argparse
import csv
import io
import json
import os
import sys

SOURCE = "1000G 30x (NYGC)"


def identifiants(chemin):
    with io.open(chemin, encoding="utf-8") as fh:
        return {l.split("\t")[0].strip() for l in fh if l.strip()} - {"sample", "fid"}


def populations(chemin):
    pop = {}
    with io.open(chemin, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pop[r["id"]] = r.get("pop", "")
    return pop


def statut(provenance):
    """Ce que la trace dit du prelevement : analyse, en echec, ou refuse."""
    with io.open(provenance, encoding="utf-8") as fh:
        d = json.load(fh)
    if d.get("reussite_complete"):
        return "analysé"
    etages = d.get("etages") or {}
    r = etages.get("recevabilite") or {}
    if r.get("etat") == "ECHEC":
        return "refusé : " + (r.get("detail") or "entrées non recevables").strip(" ;")
    echecs = d.get("etages_en_echec") or [e for e, v in etages.items() if v.get("etat") == "ECHEC"]
    return "échec : " + ", ".join(echecs)


def lignes(lot, getrm, hla, pop):
    out = []
    for e in sorted(os.listdir(lot)):
        p = os.path.join(lot, e, "sortie", "provenance.json")
        if not os.path.isfile(p):
            continue
        out.append([e, SOURCE, pop.get(e, ""), statut(p),
                    "oui" if e in getrm else "non", "oui" if e in hla else "non"])
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--lot", required=True)
    p.add_argument("--getrm", required=True)
    p.add_argument("--hla", required=True)
    p.add_argument("--populations", required=True)
    p.add_argument("--sortie", required=True)
    a = p.parse_args()

    rows = lignes(a.lot, identifiants(a.getrm), identifiants(a.hla), populations(a.populations))
    if not rows:
        print("aucun prelevement dans le lot", file=sys.stderr)
        return 1
    n = {"analysé": 0, "refusé": 0, "échec": 0}
    for r in rows:
        n[r[3].split(" ")[0]] += 1
    with io.open(a.sortie, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# Prélèvements soumis au banc — 1000 Genomes Project, séquençage 30x haute "
                 "couverture NYGC, aligné GRCh38. %d soumis : %d analysés, %d refusés à "
                 "l'entrée, %d en échec d'étage. Produit par outils/construire_banc.py.\n"
                 % (len(rows), n["analysé"], n["refusé"], n["échec"]))
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["Prélèvement", "Source", "Population", "Statut", "Vérité GeT-RM", "Vérité HLA Sanger"])
        w.writerows(rows)
    print("%d prelevements : %d analyses, %d refuses, %d en echec ; GeT-RM %d, HLA %d"
          % (len(rows), n["analysé"], n["refusé"], n["échec"],
             sum(r[4] == "oui" for r in rows), sum(r[5] == "oui" for r in rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
