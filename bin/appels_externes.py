# -*- coding: utf-8 -*-
"""Etage 5 : met les appels dedies au format attendu par l interpreteur.

CYP2D6, HLA et MT-RNR1 ne sortent pas du fichier de variants tel que
l interpreteur le lit. Les outils dedies rendent chacun leur format ; cet etage
les traduit, et refuse de produire une ligne qu il ne peut pas justifier.

Le plus d un tandem CYP2D6 doit porter des espaces : sans eux l interpreteur
rend un resultat indetermine et ne leve aucune erreur.

MT-RNR1 est haploide : l interpreteur attend un allele seul, pas une paire. Ses
vingt-quatre noms d alleles sont ceux de PyPGx, a l identique.
"""
import argparse, csv, os, re, sys


def cyp2d6(chemin):
    if not os.path.exists(chemin):
        return None, "fichier absent"
    with open(chemin, encoding="utf-8") as fh:
        lignes = list(csv.DictReader(fh, delimiter="\t"))
    if not lignes:
        return None, "fichier vide"
    g = (lignes[0].get("Genotype") or "").strip()
    if not g or g == "None":
        return None, "aucun diplotype"
    # Une paire, et une seule. Cyrius separe par un point-virgule les lectures
    # entre lesquelles il n a pas tranche ; l interpreteur s arrete sur une
    # telle ligne, et en choisir une serait affirmer ce que l outil n a pas dit.
    if ";" in g or "," in g or g.count("/") != 1:
        return None, "plusieurs diplotypes possibles : " + g
    return re.sub(r"\s*\+\s*", " + ", g), None


def hla(chemin):
    if not os.path.exists(chemin):
        return {}, "fichier absent"
    with open(chemin, encoding="utf-8") as fh:
        lignes = list(csv.DictReader(fh, delimiter="\t"))
    if not lignes:
        return {}, "fichier vide"
    r = lignes[0]
    quatre = lambda x: (lambda m: "*" + m.group(1) if m else None)(re.search(r"\*(\d+:\d+)", x or ""))
    out, manque = {}, []
    for gene, k1, k2 in (("HLA-A", "A1", "A2"), ("HLA-B", "B1", "B2")):
        a, b = quatre(r.get(k1, "")), quatre(r.get(k2, ""))
        if a and b:
            out[gene] = "%s/%s" % (a, b)
        else:
            manque.append(gene)
    return out, ("alleles illisibles pour " + ", ".join(manque)) if manque else None


def mtrnr1(chemin):
    """L etage de typage complementaire a deja tranche la recevabilite : on ne
    reprend que les lignes qu il a marquees rendues."""
    if not os.path.exists(chemin):
        return None, "fichier absent"
    with open(chemin, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("gene") != "MT-RNR1":
                continue
            if r.get("statut") != "rendu":
                return None, r.get("motif") or "non conclusif"
            a = (r.get("allele_externe") or "").strip()
            return (a, None) if a else (None, "aucun allele")
    return None, "gene absent du typage complementaire"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--travail", required=True)
    p.add_argument("--sortie", required=True)
    a = p.parse_args()

    lignes, notes = [], []
    d6, err = cyp2d6(os.path.join(a.travail, "cyp2d6.tsv"))
    if d6:
        lignes.append("CYP2D6\t" + d6)
    else:
        notes.append("CYP2D6 : " + err)
    h, err = hla(os.path.join(a.travail, "hla.tsv"))
    for gene in ("HLA-A", "HLA-B"):
        if gene in h:
            lignes.append("%s\t%s" % (gene, h[gene]))
    if err:
        notes.append("HLA : " + err)
    # Le typage complementaire est optionnel : absent, le module se comporte
    # exactement comme avant.
    tsv = os.path.join(a.travail, "pypgx.tsv")
    if os.path.exists(tsv):
        m, err = mtrnr1(tsv)
        if m:
            lignes.append("MT-RNR1\t" + m)
        else:
            notes.append("MT-RNR1 : " + err)

    with open(a.sortie, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\n".join(lignes) + ("\n" if lignes else ""))
    for n in notes:
        print("  note :", n)
    print("appels externes ecrits : %d" % len(lignes))
    return 0


if __name__ == "__main__":
    sys.exit(main())
