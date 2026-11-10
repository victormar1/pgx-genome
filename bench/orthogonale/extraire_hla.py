# -*- coding: utf-8 -*-
"""Garde les seules sequences HLA d'un fichier de loci, pour Immuannot.

Usage : extraire_hla.py <loci.fa> <sortie.fa>
"""
import sys


def main():
    entree, sortie = sys.argv[1:3]
    garder = False
    with open(sortie, "w") as f:
        for ligne in open(entree):
            if ligne.startswith(">"):
                garder = ligne[1:].split()[0].startswith("HLA-")
            if garder:
                f.write(ligne)


if __name__ == "__main__":
    main()
