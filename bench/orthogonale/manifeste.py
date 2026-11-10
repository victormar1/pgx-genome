# -*- coding: utf-8 -*-
"""Manifeste des haplotypes a extraire, une ligne par haplotype :

    ECH.hN<tab>URL de l'alignement winnowmap<tab>URL de l'assemblage

Les deux tables viennent de l'index HPRC ; les adresses s3:// sont reecrites en
https, seule forme que lit le samtools du poste.

Usage : manifeste.py <liste_echantillons.txt> <assemblies.csv> <aln_grch38.csv>
"""
import csv
import sys


def https(u):
    return u.replace("s3://human-pangenomics", "https://human-pangenomics.s3.amazonaws.com")


def main():
    liste, assemblages, alignements = sys.argv[1:4]
    voulus = {l.strip() for l in open(liste) if l.strip()}
    fasta = {r["assembly_name"]: https(r["assembly"]) for r in csv.DictReader(open(assemblages))}
    for r in csv.DictReader(open(alignements)):
        if r["sample_id"] in voulus:
            print(f"{r['sample_id']}.h{r['haplotype']}\t{https(r['location'])}\t{fasta[r['assembly_name']]}")


if __name__ == "__main__":
    main()
