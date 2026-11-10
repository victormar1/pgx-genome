# -*- coding: utf-8 -*-
"""Mesures du protocole, a partir d'une ou deux tables de confrontation.

Sur le perimetre clinique : concordance du module avec la verite, par gene et en
tout, avec l'intervalle de confiance a 95 % de Wilson. Puis les allèles HLA a
risque, les genomes sans compte rendu, et, si deux bras sont donnes, les couples
ou ils ne disent pas la meme chose.

Un couple gene x genome ne compte que si la verite est conclusive et que le
module a rendu quelque chose : les deux cas ecartes sont comptes a part, car un
denominateur bricole rend des taux faux.

Usage : mesurer.py T.tsv [T_autre_bras.tsv] [--sensibilite]
  --sensibilite  refait les mesures sans les loci a chaine douteuse
"""
import collections
import csv
import math
import sys

Z = 1.959963984540054          # loi normale, 97,5 %
BONS = ("concordant", "memes copies, phase differente")
RISQUE = {"HLA-A": ("31:01",), "HLA-B": ("15:02", "58:01", "57:01")}


def wilson(k, n):
    """Intervalle de Wilson a 95 % pour k succes sur n."""
    if not n:
        return (0.0, 0.0, 0.0)
    p = k / n
    d = 1 + Z * Z / n
    centre = (p + Z * Z / (2 * n)) / d
    demi = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return (p, max(0.0, centre - demi), min(1.0, centre + demi))


def retenus(lignes, douteux=None):
    """Couples mesurables, et les motifs d'ecart."""
    gardes, ecartes = [], collections.Counter()
    for l in lignes:
        if l["perimetre"] != "oui":
            continue
        if douteux and (l["echantillon"], l["gene"]) in douteux:
            ecartes["locus a chaine douteuse"] += 1
        elif not l["assemblage"] or l["module_vs_assemblage"] == "pas de verite":
            # Une verite peut etre presente et pourtant ne rien conclure : un
            # « Unknown » rendu par l'outil de verite remplit la colonne sans
            # trancher. Juger la vacuite de la colonne comptait ce couple comme
            # un echec, et affichait son gene a 0 %.
            ecartes[l["note"] or "verite non conclusive"] += 1
        elif l["module"] == "sans compte rendu":
            ecartes["genome sans compte rendu"] += 1
        elif not l["module_vs_assemblage"]:
            ecartes["module sans appel"] += 1
        else:
            gardes.append(l)
    return gardes, ecartes


def taux(lignes, titre):
    print(f"\n{titre}")
    par_gene = collections.defaultdict(lambda: [0, 0])
    for l in lignes:
        d = par_gene[l["gene"]]
        d[1] += 1
        d[0] += l["module_vs_assemblage"] in BONS
    for g in sorted(par_gene):
        k, n = par_gene[g]
        p, bas, haut = wilson(k, n)
        print(f"  {g:9} {k:4} / {n:4}   {100 * p:6.2f} %   [{100 * bas:6.2f} ; {100 * haut:6.2f}]")
    k = sum(d[0] for d in par_gene.values())
    n = sum(d[1] for d in par_gene.values())
    p, bas, haut = wilson(k, n)
    print(f"  {'ensemble':9} {k:4} / {n:4}   {100 * p:6.2f} %   [{100 * bas:6.2f} ; {100 * haut:6.2f}]")
    return bas


def risque_hla(lignes):
    print("\nAlleles HLA a risque")
    for gene, alleles in sorted(RISQUE.items()):
        for a in alleles:
            vp = fn = fp = vn = 0
            for l in lignes:
                if l["gene"] != gene or not l["assemblage"] or not l["module"]:
                    continue
                dans_v = a in l["assemblage"].split("/")
                dans_m = a in [x.strip().lstrip("*") for x in l["module"].split("/")]
                vp += dans_v and dans_m
                fn += dans_v and not dans_m
                fp += not dans_v and dans_m
                vn += not dans_v and not dans_m
            if vp + fn + fp == 0:
                print(f"  {gene}*{a} : aucun porteur, {vn} genomes sans l'allele des deux cotes")
                continue
            sens = f"{100 * vp / (vp + fn):.1f} %" if vp + fn else "sans porteur"
            spec = f"{100 * vn / (vn + fp):.1f} %" if vn + fp else "-"
            print(f"  {gene}*{a} : porteurs {vp + fn}, retrouves {vp}, manques {fn}, "
                  f"a tort {fp} ; sensibilite {sens}, specificite {spec}")


def sans_fonction_en_normal(lignes):
    """Un allele sans fonction de la verite rendu en fonction normale est un
    echec qualitatif, independamment du taux."""
    cas = [l for l in lignes if l["module_vs_assemblage"] == "DISCORDANT"
           and set(x.strip().lstrip("*") for x in l["module"].replace(" ", "").split("/")) == {"1"}
           and l["assemblage"].replace(" ", "") not in ("*1/*1",)]
    print(f"\nVerite non-reference rendue *1/*1 : {len(cas)}")
    for l in cas:
        print(f"  {l['echantillon']} {l['gene']:9} verite {l['assemblage']:24} module {l['module']}")


def charger(f):
    return list(csv.DictReader(open(f, encoding="utf-8"), delimiter="\t"))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    options = {a for a in sys.argv[1:] if a.startswith("--")}
    if not args:
        raise SystemExit(__doc__)
    bras = {}
    for f in args:
        lignes = charger(f)
        nom = (lignes[0]["bras"] if lignes and lignes[0]["bras"] else f)
        bras[nom] = lignes

    for nom, lignes in bras.items():
        gardes, ecartes = retenus(lignes)
        genomes = len({l["echantillon"] for l in lignes})
        print(f"\n===== bras {nom} : {genomes} genomes")
        print("couples ecartes :", dict(ecartes) or "aucun")
        taux(gardes, "Concordance du module avec la verite, perimetre clinique")
        risque_hla([l for l in lignes if l["gene"] in RISQUE])
        sans_fonction_en_normal(gardes)
        if "--sensibilite" in options:
            douteux = {(l["echantillon"], l["gene"]) for l in lignes if "douteu" in l["note"]}
            if douteux:
                g2, e2 = retenus(lignes, douteux)
                print(f"\n-- sans les {len(douteux)} loci a chaine douteuse")
                taux(g2, "Concordance, perimetre clinique")

    if len(bras) == 2:
        (n1, l1), (n2, l2) = bras.items()
        cle = lambda l: (l["echantillon"], l["gene"])
        d1 = {cle(l): l for l in l1 if l["perimetre"] == "oui"}
        d2 = {cle(l): l for l in l2 if l["perimetre"] == "oui"}
        communs = sorted(set(d1) & set(d2))
        ecarts = [k for k in communs if d1[k]["module"] != d2[k]["module"]]
        print(f"\n===== bras {n1} contre bras {n2} : {len(communs)} couples communs, "
              f"{len(ecarts)} ou le rendu differe")
        for k in ecarts:
            print(f"  {k[0]} {k[1]:9} {n1} {d1[k]['module']:24} {n2} {d2[k]['module']:24} "
                  f"verite {d1[k]['assemblage']}")


if __name__ == "__main__":
    main()
