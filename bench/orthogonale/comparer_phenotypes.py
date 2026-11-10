# -*- coding: utf-8 -*-
"""Concordance au niveau qui decide : le phenotype, et non le diplotype.

Le banc mesure d'abord la concordance de diplotype. Celle-ci penalise le module
sur de la nomenclature : un *36 de trop dans un tandem CYP2D6, un sous-allele
HLA rendu par le representant de son groupe. Or un compte rendu ne rend pas un
diplotype, il rend une conduite a tenir, et celle-ci depend du phenotype.

Les deux cotes ont un rapport PharmCAT complet : la verite dans
<verite>/pc_<E>/<E>.report.json, le module dans <lot>/<E>/sortie/<E>.report.json.

    python comparer_phenotypes.py --verite DOSSIER --lot DOSSIER [--bras a]

Limite a dire avec le chiffre : la verite est batie sans appel exterieur. Elle
rend « No Result » pour CYP2D6 et ne porte ni HLA ni POR, qui sortent donc du
denominateur au lieu de compter un desaccord qui n'en est pas un.
"""
import argparse
import glob
import io
import json
import os
import sys

CLINIQUE = ("ABCG2", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A4", "CYP3A5",
            "HLA-A", "HLA-B", "NUDT15", "POR", "SLCO1B1", "TPMT")


def charger(chemin):
    try:
        with io.open(chemin, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def phenotypes(rapport, gene):
    """L'ensemble des phenotypes proposes pour ce gene, ou None si absent."""
    g = (rapport or {}).get("genes") or {}
    if isinstance(g.get("CPIC"), dict):
        g = g["CPIC"]
    e = g.get(gene)
    if not isinstance(e, dict):
        return None
    vus = set()
    for d in e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []:
        for p in d.get("phenotypes") or []:
            p = (p or "").strip()
            if p:
                vus.add(p)
    return vus or None


def comparer(verite, lot):
    """Rend (par_gene, desaccords) pour un bras."""
    genomes = sorted(
        os.path.basename(os.path.dirname(os.path.dirname(p)))
        for p in glob.glob(os.path.join(lot, "*", "sortie", "*.report.json"))
        if os.path.getsize(p) > 1024
    )
    par_gene, desaccords = {}, []
    for e in genomes:
        rv = charger(os.path.join(verite, "pc_" + e, e + ".report.json"))
        rm = charger(os.path.join(lot, e, "sortie", e + ".report.json"))
        if rv is None or rm is None:
            continue
        for gene in CLINIQUE:
            pv, pm = phenotypes(rv, gene), phenotypes(rm, gene)
            c = par_gene.setdefault(gene, {"vus": 0, "accord": 0, "sans": 0})
            if pv is None or pm is None or pv == {"No Result"}:
                c["sans"] += 1
                continue
            c["vus"] += 1
            if pv == pm:
                c["accord"] += 1
            else:
                desaccords.append((e, gene, " | ".join(sorted(pv)),
                                   " | ".join(sorted(pm))))
    return len(genomes), par_gene, desaccords


def main():
    a = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    a.add_argument("--verite", required=True,
                   help="dossier des rapports PharmCAT de la verite")
    a.add_argument("--lot", required=True,
                   help="dossier des sorties du module pour un bras")
    a.add_argument("--bras", default="", help="etiquette affichee")
    o = a.parse_args()

    n, par_gene, desaccords = comparer(o.verite, o.lot)
    print("=== bras %s : %d genome(s) ===" % (o.bras.upper() or "?", n))
    tot = acc = 0
    for gene in CLINIQUE:
        c = par_gene.get(gene)
        if not c or not c["vus"]:
            print("  %-10s  aucun couple comparable" % gene)
            continue
        tot += c["vus"]
        acc += c["accord"]
        print("  %-10s %4d / %4d  %6.2f %%%s"
              % (gene, c["accord"], c["vus"], 100.0 * c["accord"] / c["vus"],
                 "   (%d sans phenotype des deux cotes)" % c["sans"]
                 if c["sans"] else ""))
    if tot:
        print("  %-10s %4d / %4d  %6.2f %%"
              % ("ensemble", acc, tot, 100.0 * acc / tot))
    print("  desaccords de phenotype : %d" % len(desaccords))
    for e, gene, pv, pm in desaccords:
        print("    %-10s %-9s verite [%s]  module [%s]" % (e, gene, pv, pm))
    return 0


if __name__ == "__main__":
    sys.exit(main())
