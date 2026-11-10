#!/usr/bin/env python3
"""VCF de verite d'un echantillon HPRC : variants des deux haplotypes
d'assemblage, phases h1|h2, sur les regions des positions PharmCAT et du
complement RNPGx.

Une region que l'un des deux haplotypes ne couvre pas d'un seul alignement de
bout en bout n'apporte aucun variant : ses positions sont ecrites ./. pour que
PharmCAT ne les lise pas comme reference (--absent-to-ref laisse les ./.). Les
genes touches sont declares incomplets dans le tableau d'etat.

Usage : vcf_verite.py <echantillon> <prefixe_h1> <prefixe_h2> <ref.fa.fai>
        <sortie.vcf> <positions.vcf> [<positions.vcf> ...]
"""
import collections
import sys

PRIMAIRES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def lire_statuts(p):
    st = {}
    for ligne in open(p):
        c = ligne.rstrip("\n").split("\t")
        if c[0] == "P":
            st[(c[2], int(c[3]), int(c[4]))] = c[5]
    return st


def lire_variants(p):
    v = set()
    for ligne in open(p):
        c = ligne.rstrip("\n").split("\t")
        if len(c) == 4:
            v.add((c[0], int(c[1]), c[2], c[3]))
    return v


def lire_positions(chemins):
    for p in chemins:
        for ligne in open(p, encoding="utf-8"):
            if ligne.startswith("#"):
                continue
            c = ligne.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
            yield c[0], int(c[1]), c[3], info.get("PX", "?")


def main():
    ech, h1, h2, fai, sortie = sys.argv[1:6]
    fichiers_positions = sys.argv[6:]
    st1, st2 = lire_statuts(h1 + ".statuts.tsv"), lire_statuts(h2 + ".statuts.tsv")
    regions = sorted(set(st1) | set(st2))
    bonnes = {r for r in regions if st1.get(r) == "complet" and st2.get(r) == "complet"}
    par_chrom = collections.defaultdict(list)
    for r in regions:
        par_chrom[r[0]].append(r)

    def chevauchent(chrom, s0, e0):
        return [r for r in par_chrom[chrom] if r[1] < e0 and r[2] > s0]

    def retenu(v):
        chev = chevauchent(v[0], v[1] - 1, v[1] - 1 + len(v[2]))
        return bool(chev) and all(r in bonnes for r in chev)

    v1 = {v for v in lire_variants(h1 + ".variants.tsv") if retenu(v)}
    v2 = {v for v in lire_variants(h2 + ".variants.tsv") if retenu(v)}
    groupes = collections.defaultdict(set)
    for v in v1 | v2:
        groupes[v[:3]].add(v[3])
    lignes = []
    for (chrom, pos, ref), alts in groupes.items():
        alts = sorted(alts)
        gt = lambda hv: str(next((i + 1 for i, a in enumerate(alts) if (chrom, pos, ref, a) in hv), 0))
        lignes.append((chrom, pos, ref, ",".join(alts), f"{gt(v1)}|{gt(v2)}"))

    occupees = {(l[0], l[1]) for l in lignes}
    incomplets = collections.Counter()
    vus = set()
    for chrom, pos, ref, gene in lire_positions(fichiers_positions):
        chev = chevauchent(chrom, pos - 1, pos)
        if chev and all(r in bonnes for r in chev):
            continue
        incomplets[gene] += 1
        if (chrom, pos) not in occupees and (chrom, pos) not in vus:
            vus.add((chrom, pos))
            lignes.append((chrom, pos, ref, ".", "./."))

    longueurs = dict(l.split("\t")[:2] for l in open(fai))
    ordre = {c: i for i, c in enumerate(PRIMAIRES)}
    lignes.sort(key=lambda l: (ordre.get(l[0], 99), l[1]))
    with open(sortie, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=verite_assemblage_HPRC_R2\n")
        for c in PRIMAIRES:
            f.write(f"##contig=<ID={c},length={longueurs[c]}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, haplotype 1 | haplotype 2">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + ech + "\n")
        for chrom, pos, ref, alt, gt in lignes:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    with open(sortie + ".etat.tsv", "w") as f:
        f.write("gene\tpositions_non_couvertes\n")
        for g, n in sorted(incomplets.items()):
            f.write(f"{g}\t{n}\n")
    print(f"{ech} : {len(lignes)} lignes, {len(v1)}+{len(v2)} variants, "
          f"genes incomplets : {dict(incomplets) or 'aucun'}", file=sys.stderr)


if __name__ == "__main__":
    main()
