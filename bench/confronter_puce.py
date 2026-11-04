# -*- coding: utf-8 -*-
"""ABCG2 rs2231142 et POR rs1057868 : le module contre la puce Omni 2.5.

Pour chaque genome du lot genotype sur puce, on compare le nombre d'alleles
alternatifs :
  - ABCG2 : le genotype du controle qualite (travail/qc_positions.tsv), celui
    que l'interpreteur recoit ;
  - POR : le diplotype rendu par le typage complementaire (travail/pypgx.tsv),
    dont les alleles *28 et *36 portent la position.

Usage :
    confronter_puce.py --lot DOSSIER --puce puce.tsv --echantillons liste.txt
Le code de retour vaut 0 sans ecart, 1 sinon.
"""
import argparse, collections, csv, io, os, sys

POSITIONS = {"ABCG2": "4:89052323", "POR": "7:75615006"}
PORTEURS_POR = ("*28", "*36")


def nalt(g):
    return {"0/0": 0, "0|0": 0, "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1,
            "1/1": 2, "1|1": 2}.get((g or "").strip())


def lire_puce(puce, echantillons):
    ech = [l.strip() for l in io.open(echantillons, encoding="utf-8") if l.strip()]
    gt = {}
    for l in io.open(puce, encoding="utf-8"):
        c = l.rstrip("\n").split("\t")
        for gene, pos in POSITIONS.items():
            if c[0] == pos:
                gt[gene] = dict(zip(ech, c[1:]))
    return gt


def module_abcg2(travail):
    p = os.path.join(travail, "qc_positions.tsv")
    if not os.path.isfile(p):
        return None
    for r in csv.DictReader(io.open(p, encoding="utf-8"), delimiter="\t"):
        if r["rsid"] == "rs2231142" and r.get("retenue") == "oui":
            return nalt(r["genotype"])
    return None


def module_por(travail):
    p = os.path.join(travail, "pypgx.tsv")
    if not os.path.isfile(p):
        return None
    for r in csv.DictReader(io.open(p, encoding="utf-8"), delimiter="\t"):
        if r["gene"] == "POR" and r["statut"] == "rendu":
            return sum(r["genotype"].split("/").count(a) for a in PORTEURS_POR)
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--lot", required=True)
    p.add_argument("--puce", required=True)
    p.add_argument("--echantillons", required=True)
    a = p.parse_args()
    puce = lire_puce(a.puce, a.echantillons)
    lecteurs = {"ABCG2": module_abcg2, "POR": module_por}
    ecarts = 0
    for gene in POSITIONS:
        n, ecart, rep = 0, [], collections.Counter()
        for e in sorted(os.listdir(a.lot)):
            t = os.path.join(a.lot, e, "travail")
            v = nalt(puce.get(gene, {}).get(e))
            m = lecteurs[gene](t) if os.path.isdir(t) else None
            if v is None or m is None:
                continue
            n += 1
            rep[v] += 1
            if m != v:
                ecart.append("%s (puce %d, module %d)" % (e, v, m))
        print("%-6s %s : %d / %d concordants ; 0/1/2 alleles alternatifs : %d / %d / %d"
              % (gene, POSITIONS[gene], n - len(ecart), n, rep[0], rep[1], rep[2]))
        for x in ecart:
            print("   ecart", x)
        ecarts += len(ecart)
    return 1 if ecarts else 0


if __name__ == "__main__":
    sys.exit(main())
