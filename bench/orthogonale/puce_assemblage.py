# -*- coding: utf-8 -*-
"""ABCG2 rs2231142 et POR rs1057868 : la verite par assemblage contre la puce
Omni 2.5 du 1000 Genomes.

Ces deux positions n'ont pas de reference GeT-RM. La puce est une technologie
independante du sequencage comme de l'assemblage : elle dit si la verite tiree
des haplotypes lit la bonne base, en dehors de toute chaine de lecture courte.

La comparaison porte sur le nombre d'alleles alternatifs, 0, 1 ou 2 : la puce ne
phase pas. Les deux positions sont lues chez la puce par leurs coordonnees
GRCh37, sur le meme brin qu'en GRCh38.

Le fichier de la puce et sa liste d'echantillons se produisent par
`bench/confronter_puce.sh`, qui les lit a distance.

Usage : puce_assemblage.py --travail DOSSIER --puce puce.tsv
                           --echantillons-puce liste.txt --genomes LISTE
Le code de retour vaut 0 sans ecart, 1 sinon.
"""
import argparse
import collections
import io
import os
import sys

import chemins

sys.path.insert(0, os.path.join(chemins.RACINE, "bench"))
from confronter_puce import POSITIONS, lire_puce, nalt  # noqa: E402

# Memes rsID, en GRCh38, tels que les porte le VCF des positions du module.
EN_38 = {"ABCG2": ("chr4", "88131171", "G", "T"), "POR": ("chr7", "75985688", "C", "T")}


def nalt_assemblage(trav, ech, gene):
    """Nombre d'haplotypes qui portent l'allele alternatif, ou None si une des
    deux regions n'est pas couverte d'un seul alignement."""
    chrom, pos, ref, alt = EN_38[gene]
    porteurs = 0
    for h in (1, 2):
        base = os.path.join(trav, "haplotypes", f"{ech}.h{h}")
        if not os.path.exists(base + ".statuts.tsv"):
            return None
        with io.open(base + ".statuts.tsv", encoding="utf-8") as f:
            lignes = [l.rstrip("\n").split("\t") for l in f.read().splitlines()[1:] if l.strip()]
        if not any(c[2] == chrom and int(c[3]) < int(pos) <= int(c[4]) and c[5] == "complet"
                   for c in lignes):
            return None
        with io.open(base + ".variants.tsv", encoding="utf-8") as f:
            porteurs += any(l.split("\t")[:4] == [chrom, pos, ref, alt]
                            for l in f.read().splitlines() if l.strip())
    return porteurs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--travail")
    p.add_argument("--puce", required=True)
    p.add_argument("--echantillons-puce", required=True)
    p.add_argument("--genomes", required=True)
    a = p.parse_args()
    trav = chemins.travail(a.travail)
    puce = lire_puce(a.puce, a.echantillons_puce)
    ech = [l.strip() for l in io.open(a.genomes, encoding="utf-8") if l.strip()]

    bilan = collections.Counter()
    ecarts = []
    for gene in sorted(POSITIONS):
        for s in ech:
            attendu = nalt(puce.get(gene, {}).get(s))
            obtenu = nalt_assemblage(trav, s, gene)
            if attendu is None:
                bilan[f"{gene} : absent de la puce"] += 1
            elif obtenu is None:
                bilan[f"{gene} : region non couverte"] += 1
            elif attendu == obtenu:
                bilan[f"{gene} : accord"] += 1
            else:
                bilan[f"{gene} : ECART"] += 1
                ecarts.append((s, gene, attendu, obtenu))
    for cle in sorted(bilan):
        print(f"  {cle} : {bilan[cle]}")
    for s, gene, attendu, obtenu in ecarts:
        print(f"  {s} {gene} : puce {attendu} allele(s) alternatif(s), assemblage {obtenu}")
    return 1 if ecarts else 0


if __name__ == "__main__":
    sys.exit(main())
