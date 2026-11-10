# -*- coding: utf-8 -*-
"""Confronte, genome par genome et gene par gene, la verite tiree des
assemblages HPRC, les vérités de reference quand elles existent (GeT-RM, typage
Sanger), et le rendu du module.

La verite ne vient pas de la meme source pour tous les genes :
- genes a variants ponctuels : PharmCAT sur le VCF phase des deux haplotypes ;
- CYP2D6 : deCYPher sur la sequence du locus, car le gene se decide par sa
  structure, que le VCF ne porte pas ;
- HLA-A et HLA-B : Immuannot, a deux champs.
Le rendu de PharmCAT pour CYP2D6 et le HLA est donc ignore : ce ne sont pas des
appels de verite.

Usage : confronter.py --travail DOSSIER --lot DOSSIER --echantillons LISTE
                      --sortie T.tsv [--bras NOM]
"""
import argparse
import csv
import gzip
import json
import os
import re
import sys

import chemins
import verite_cyp2d6 as v6

sys.path.insert(0, os.path.join(chemins.RACINE, "bench"))
import comparer  # noqa: E402
from comparer import compare  # noqa: E402

# Memes tables d'equivalence de nomenclature que le banc de la cohorte : sans
# elles, deux noms d'un meme diplotype comptent comme un desaccord.
for _f, _d in (("alleles_pharmcat.json", comparer.ALLELES_OUTIL),
               ("definitions_alleles.json", comparer.DEFINITIONS)):
    with open(os.path.join(chemins.VERITES, _f), encoding="utf-8") as _fh:
        _d.update(json.load(_fh))

AUTRE_SOURCE = ("CYP2D6", "HLA-A", "HLA-B", "POR")
RISQUE = {"HLA-A": ("31:01",), "HLA-B": ("15:02", "58:01", "57:01")}
# POR n'a pas de table PharmCAT : sa verite se lit sur la position elle-meme, et
# le module le rend par son typage complementaire, hors du rapport PharmCAT.
POR = ("chr7", "75985688", "C", "T", "*28")
PYPGX_SEUL = ("POR",)


def perimetre():
    f = os.path.join(chemins.RESSOURCES, "perimetre_rnpgx.json")
    return set(json.load(open(f, encoding="utf-8"))["genes"])


def etiquettes(rapport):
    """Diplotypes rendus par un rapport PharmCAT, un par gene."""
    if not os.path.exists(rapport):
        return None
    j = json.load(open(rapport, encoding="utf-8"))
    sortie = {}
    for g, e in j["genes"].items():
        dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
        sortie[g] = "; ".join(sorted({x.get("label") for x in dips if x.get("label")}))
    return sortie


def incomplets(trav, ech):
    """Genes dont une region n'est pas couverte par les deux haplotypes."""
    p = os.path.join(trav, "verite", ech + ".vcf.etat.tsv")
    return {l.split("\t")[0] for l in open(p).read().splitlines()[1:] if l.strip()} \
        if os.path.exists(p) else set()


def hla_assemblage(trav, ech, gene):
    """Allele a deux champs de chaque haplotype, ou None si Immuannot n'en
    trouve pas exactement un par haplotype."""
    trouves = []
    for h in (1, 2):
        f = os.path.join(trav, "immuannot", f"{ech}.h{h}.gtf.gz")
        if not os.path.exists(f):
            return None
        copies = []
        for ligne in gzip.open(f, "rt"):
            c = ligne.rstrip("\n").split("\t")
            if len(c) > 8 and c[2] == "transcript" and f'gene_name "{gene}"' in c[8]:
                a = re.search(r'consensus "([^"]+)"', c[8])
                if a:
                    copies.append(":".join(a.group(1).split("*")[-1].split(":")[:2]))
        if len(copies) != 1:
            return None
        trouves.append(copies[0])
    return "/".join(sorted(trouves))


def deux_champs(x):
    """« HLA-A*02:01:01:02N », « *02:01 » ou « 02:01 » -> « 02:01 »."""
    return ":".join(x.strip().split("*")[-1].split(":")[:2])


def locus_douteux(trav, ech, gene):
    """Vrai si, chez l'un des deux haplotypes, le locus de ce gene a du etre
    reconstruit le long du contig par un chainage que l'extraction juge
    douteux. L'analyse de sensibilite retire ces couples."""
    nom = {"CYP2D6": "CYP2D6_locus"}.get(gene, gene)
    for h in (1, 2):
        f = os.path.join(trav, "haplotypes", f"{ech}.h{h}.statuts.tsv")
        if not os.path.exists(f):
            continue
        with open(f, encoding="utf-8") as fh:
            for l in fh.read().splitlines()[1:]:
                c = l.split("\t")
                if len(c) > 5 and c[0] == "L" and c[1] == nom and c[5] == "chaine_douteuse":
                    return True
    return False


def pypgx_module(lot, ech, gene):
    """Genotype rendu par le typage complementaire du module. PharmCAT n'a pas de
    table POR : sans cette lecture, un gene du perimetre que le module rend
    vraiment serait compte comme sans appel."""
    f = os.path.join(lot, ech, "travail", "pypgx.tsv")
    if not os.path.exists(f):
        return ""
    with open(f, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("gene") == gene and r.get("statut") == "rendu":
                return (r.get("genotype") or "").strip()
    return ""


def por_assemblage(trav, ech):
    """Genotype de rs1057868, lu sur chaque haplotype. Rend None si une des deux
    regions n'est pas couverte."""
    chrom, pos, ref, alt, etoile = POR
    copies = []
    for h in (1, 2):
        st = os.path.join(trav, "haplotypes", f"{ech}.h{h}.statuts.tsv")
        if not os.path.exists(st):
            return None
        couvert = any(c[2] == chrom and int(c[3]) < int(pos) <= int(c[4]) and c[5] == "complet"
                      for c in (l.split("\t") for l in open(st).read().splitlines()[1:] if l.strip()))
        if not couvert:
            return None
        v = os.path.join(trav, "haplotypes", f"{ech}.h{h}.variants.tsv")
        porte = any(l.split("\t")[:4] == [chrom, pos, ref, alt]
                    for l in open(v).read().splitlines() if l.strip())
        copies.append(etoile if porte else "*1")
    return "/".join(sorted(copies))


def hla_reference(sanger, ech, gene):
    r = sanger.get(ech)
    if not r:
        return ""
    k1, k2 = ("A", "A.1") if gene == "HLA-A" else ("B", "B.1")
    deux = deux_champs
    a = sorted({deux(x) for x in (r.get(k1) or "").split("/") if x.strip()})
    b = sorted({deux(x) for x in (r.get(k2) or "").split("/") if x.strip()})
    # Le typage Sanger peut laisser plusieurs possibilites par copie : on garde
    # les deux listes, la comparaison teste l'appariement.
    return "|".join(["+".join(a), "+".join(b)]) if a and b else ""


def hla_verdict(asm, ref):
    if not asm or not ref:
        return ""
    x, y = asm.split("/")
    a, b = [c.split("+") for c in ref.split("|")]
    ok = (x in a and y in b) or (x in b and y in a)
    return "concordant" if ok else "DISCORDANT"


def egal_hors_phase(un, autre):
    """Memes copies de CYP2D6 des deux cotes, dans un autre ordre."""
    a, b = v6.multiensemble(un), v6.multiensemble(autre.replace(" ", ""))
    return a is not None and a == b


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--travail")
    p.add_argument("--lot")
    p.add_argument("--echantillons", required=True)
    p.add_argument("--sortie", required=True)
    p.add_argument("--bras", default="")
    a = p.parse_args()
    trav, lot = chemins.travail(a.travail), chemins.lot(a.lot)
    peri = perimetre()
    ech = [l.strip() for l in open(a.echantillons) if l.strip()]
    tsv = lambda f: {r["sample"]: r for r in csv.DictReader(open(f, encoding="utf-8"), delimiter="\t")}
    getrm = tsv(os.path.join(chemins.VERITES, "verite_getrm.tsv"))
    sanger = tsv(os.path.join(chemins.VERITES, "hla_sanger.tsv"))

    lignes = []
    for s in ech:
        asm = etiquettes(os.path.join(trav, "verite", f"pc_{s}", f"{s}.report.json"))
        if asm is None:
            print(f"{s} : pas de rapport de verite", file=sys.stderr)
            continue
        mod = etiquettes(os.path.join(lot, s, "sortie", f"{s}.report.json"))
        trous = incomplets(trav, s)
        # Tout gene du perimetre a sa ligne, meme sans verite ni rendu : un
        # gene absent de la table serait un denominateur ampute en silence.
        genes = sorted(set(asm) | set(AUTRE_SOURCE) | peri)
        for g in genes:
            source, note = "pharmcat_assemblage", ""
            if g in ("CYP2D6", "HLA-A", "HLA-B") and locus_douteux(trav, s, g):
                note = "locus a chaine douteuse"
            if g == "CYP2D6":
                source = "decypher"
                cps = [v6.copies(f"{s}.h{h}", trav) for h in (1, 2)]
                haps = [v6.haplotype(c[0]) for c in cps]
                v_asm = "/".join(sorted(haps)) if all(haps) else ""
                if not v_asm:
                    note = "deCYPher non conclusif"
            elif g in ("HLA-A", "HLA-B"):
                source = "immuannot"
                v_asm = hla_assemblage(trav, s, g) or ""
                if not v_asm:
                    note = "Immuannot : pas une copie par haplotype"
            elif g == "POR":
                source = "position"
                v_asm = por_assemblage(trav, s) or ""
                note = "" if v_asm else "region non couverte par les deux haplotypes"
            elif g in trous:
                v_asm, note = "", "region non couverte par les deux haplotypes"
            elif g not in asm:
                v_asm, note = "", "gene absent du rendu de PharmCAT sur la verite"
            else:
                v_asm = asm[g]
            v_mod = (mod or {}).get(g, "")
            if g in PYPGX_SEUL and mod is not None:
                v_mod = pypgx_module(lot, s, g)
            v_ref = hla_reference(sanger, s, g) if g in ("HLA-A", "HLA-B") \
                else (getrm.get(s, {}).get(g) or "").strip()
            src_ref = ("sanger" if g in ("HLA-A", "HLA-B") else "getrm") if v_ref else ""

            if g in ("HLA-A", "HLA-B"):
                vs_ref = hla_verdict(v_asm, v_ref)
                m = sorted(deux_champs(x) for x in v_mod.split("/") if x.strip())
                vs_asm = ("concordant" if m == sorted(v_asm.split("/")) else "DISCORDANT") \
                    if v_asm and len(m) == 2 else ""
            else:
                vs_ref = compare(g, v_asm, v_ref)[0] if v_asm and v_ref else ""
                vs_asm = compare(g, v_mod, v_asm)[0] if v_asm and v_mod else ""
                if g == "CYP2D6":
                    if vs_ref == "DISCORDANT" and egal_hors_phase(v_asm, v_ref):
                        vs_ref = "memes copies, phase differente"
                    if vs_asm == "DISCORDANT" and egal_hors_phase(v_asm, v_mod):
                        vs_asm = "memes copies, phase differente"
            lignes.append({
                "bras": a.bras, "echantillon": s, "gene": g,
                "perimetre": "oui" if g in peri else "non", "source_verite": source,
                "assemblage": v_asm, "reference": v_ref, "source_reference": src_ref,
                "module": v_mod if mod is not None else "sans compte rendu",
                "assemblage_vs_reference": vs_ref, "module_vs_assemblage": vs_asm,
                "risque_hla": ",".join(r for r in RISQUE.get(g, ()) if r in v_asm.split("/")),
                "note": note,
            })

    if not lignes:
        raise SystemExit("aucune ligne : verite ou liste vide")
    with open(a.sortie, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lignes[0]), delimiter="\t")
        w.writeheader()
        w.writerows(lignes)
    print(f"{len(ech)} genomes, {len(lignes)} couples gene x genome -> {a.sortie}")


if __name__ == "__main__":
    main()
