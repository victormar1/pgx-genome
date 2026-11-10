# -*- coding: utf-8 -*-
"""Etat de la verite construite, avant toute mesure du module.

Dit ce que la verite contient et ce qui lui manque : fichiers absents, statuts
de couverture, genes ecartes, et pour CYP2D6 et le HLA les formes rencontrees.
Aucune sortie du module n'est lue ici : ce bilan doit pouvoir se faire avant que
le module ait tourne.

Usage : bilan_construction.py --travail DOSSIER --echantillons LISTE
"""
import argparse
import collections
import gzip
import json
import os
import re

import chemins
import verite_cyp2d6 as v6

PHARMCAT = ["CYP2C9", "CYP2C19", "CYP3A5", "CYP3A4", "TPMT", "NUDT15", "SLCO1B1", "ABCG2"]
RISQUE = {"HLA-A": ("31:01",), "HLA-B": ("15:02", "58:01", "57:01")}


def lignes_tsv(f):
    with open(f, encoding="utf-8") as fh:
        tete = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(tete, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def immuannot(trav, p):
    f = os.path.join(trav, "immuannot", p + ".gtf.gz")
    sortie = {}
    if os.path.exists(f):
        with gzip.open(f, "rt") as fh:
            for ligne in fh:
                c = ligne.rstrip("\n").split("\t")
                if len(c) > 8 and c[2] == "transcript":
                    g = re.search(r'gene_name "([^"]+)"', c[8])
                    a = re.search(r'consensus "([^"]+)"', c[8])
                    if g and a:
                        sortie.setdefault(g.group(1), []).append(a.group(1))
    return sortie


def deux_champs(a):
    return ":".join(a.split("*")[-1].split(":")[:2])


def etiquettes(trav, s):
    f = os.path.join(trav, "verite", f"pc_{s}", f"{s}.report.json")
    if not os.path.exists(f):
        return {}
    with open(f, encoding="utf-8") as fh:
        j = json.load(fh)
    sortie = {}
    for g, e in j["genes"].items():
        dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
        sortie[g] = sorted({x.get("label") for x in dips if x.get("label")})
    return sortie


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--travail")
    p.add_argument("--echantillons", required=True)
    a = p.parse_args()
    trav = chemins.travail(a.travail)
    ech = [l.strip() for l in open(a.echantillons) if l.strip()]
    print(f"genomes : {len(ech)}")

    manque = collections.Counter()
    st_p, st_l = collections.Counter(), collections.defaultdict(collections.Counter)
    exclus = collections.Counter()
    for s in ech:
        for h in (1, 2):
            q = f"{s}.h{h}"
            for nom, f in (("statuts", f"haplotypes/{q}.statuts.tsv"),
                           ("deCYPher", f"decypher/{q}/{q}.json"),
                           ("Immuannot", f"immuannot/{q}.gtf.gz")):
                if not os.path.exists(os.path.join(trav, f)):
                    manque[nom] += 1
            f = os.path.join(trav, f"haplotypes/{q}.statuts.tsv")
            if os.path.exists(f):
                for r in lignes_tsv(f):
                    (st_p if r["genre"] == "P" else st_l[r["nom"]])[r["statut"]] += 1
        for nom, f in (("etat", f"verite/{s}.vcf.etat.tsv"),
                       ("rapport", f"verite/pc_{s}/{s}.report.json")):
            if not os.path.exists(os.path.join(trav, f)):
                manque[nom] += 1
        f = os.path.join(trav, f"verite/{s}.vcf.etat.tsv")
        if os.path.exists(f):
            for r in lignes_tsv(f):
                exclus[r["gene"]] += 1

    print("fichiers manquants :", dict(manque) or "aucun")
    print("statuts des regions de positions (haplotypes x regions) :", dict(st_p))
    print("statuts des loci :")
    for nom in sorted(st_l):
        print(f"  {nom:14} {dict(st_l[nom])}")
    print("genes ecartes de la verite PharmCAT (genomes) :", dict(exclus.most_common()))

    print("\nappels PharmCAT par gene a variants ponctuels :")
    for g in PHARMCAT:
        n_ok = n_amb = n_inc = 0
        for s in ech:
            e = etiquettes(trav, s).get(g, [])
            if not e or any("Unknown" in x for x in e):
                n_inc += 1
            elif len(e) > 1 or any(" or " in x for x in e):
                n_amb += 1
            else:
                n_ok += 1
        print(f"  {g:8} appel unique {n_ok}, plusieurs {n_amb}, sans appel {n_inc} "
              f"(dont ecartes {exclus.get(g, 0)})")

    print("\nCYP2D6 (deCYPher, apres la regle de choix) :")
    concl, corrs, formes, non_concl = 0, [], collections.Counter(), []
    for s in ech:
        haps = []
        for h in (1, 2):
            cps, c = v6.copies(f"{s}.h{h}", trav)
            corrs += [(f"{s}.h{h}", x, y) for x, y in c]
            haps.append(v6.haplotype(cps))
            if cps:
                if cps == ["*5"]:
                    formes["deletion *5"] += 1
                elif len(cps) == 1:
                    formes["une copie"] += 1
                elif len(set(cps)) == 1:
                    formes["duplication"] += 1
                else:
                    formes["tandem ou hybride"] += 1
        if all(haps):
            concl += 1
        else:
            non_concl.append(s)
    print(f"  genomes conclusifs {concl} sur {len(ech)} ; formes d'haplotypes {dict(formes)}")
    print(f"  non conclusifs : {' '.join(non_concl) or 'aucun'}")
    print(f"  choix differents de l'appel de deCYPher : {len(corrs)} "
          + ", ".join(f"{q} {x}->{y}" for q, x, y in corrs))

    print("\nHLA (Immuannot, deux champs, une copie par haplotype) :")
    for gene in ("HLA-A", "HLA-B"):
        concl, porteurs, non = 0, collections.Counter(), []
        for s in ech:
            a1, a2 = immuannot(trav, f"{s}.h1"), immuannot(trav, f"{s}.h2")
            n1, n2 = len(a1.get(gene, [])), len(a2.get(gene, []))
            if n1 != 1 or n2 != 1:
                non.append(f"{s}({n1}+{n2})")
                continue
            concl += 1
            al = {deux_champs(a1[gene][0]), deux_champs(a2[gene][0])}
            for r in RISQUE[gene]:
                if r in al:
                    porteurs[r] += 1
        print(f"  {gene} : conclusifs {concl} sur {len(ech)} ; porteurs {dict(porteurs)} ; "
              f"non conclusifs {' '.join(non) or 'aucun'}")


if __name__ == "__main__":
    main()
