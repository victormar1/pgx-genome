# -*- coding: utf-8 -*-
"""Construit ressources/rnpgx_complement.{vcf,bed} : les positions de classe 1 et 2
du core panel RNPGx 2026 absentes des definitions PharmCAT.

Coordonnees GRCh38 depuis Ensembl, chaque allele de reference verifie contre le
FASTA GRCh38 avant publication (bcftools norm --check-ref). Provenance distincte
de pharmcat_positions.* (MPL, derive de PharmCAT).

Usage : python3 outils/construire_complement_rnpgx.py [dossier_ressources]
"""
import os, sys, collections

RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ressources")
PAD = 50   # meme marge que pharmcat_positions.bed

# (gene, rsid, chrom, pos, ref, alt, classe, type, interpreteur, note)
P = [
 ("BCHE", "rs1799807", "chr3", 165830741, "T", "C", 1, "snv", "aucun", "p.Asp70Gly, allèle atypique"),
 ("BCHE", "rs1803274", "chr3", 165773492, "C", "T", 1, "snv", "aucun", "p.Ala539Thr, variante K"),
 ("MT-RNR1", "rs267606618", "chrM", 1095, "T", "C", 1, "snv", "aucun", "m.1095T>C, surdité sous aminoside"),
 ("MT-RNR1", "rs267606619", "chrM", 1494, "C", "T", 1, "snv", "aucun", "m.1494C>T, surdité sous aminoside"),
 ("MT-RNR1", "rs267606617", "chrM", 1555, "A", "G", 1, "snv", "aucun", "m.1555A>G, surdité sous aminoside"),
 ("NAT2", "rs1041983", "chr8", 18400285, "C", "T", 1, "snv", "pharmcat", "c.282C>T, hors définitions PharmCAT"),
 ("NAT2", "rs72554616", "chr8", 18400437, "A", "C", 1, "snv", "pharmcat", "c.434A>C, hors définitions PharmCAT"),
 ("NAT2", "rs1799929", "chr8", 18400484, "C", "T", 1, "snv", "pharmcat", "c.481C>T, hors définitions PharmCAT"),
 ("VKORC1", "rs9934438", "chr16", 31093557, "G", "A,C", 1, "snv", "pharmcat", "intron 1 (1173C>T), hors définitions PharmCAT"),
 ("CYP2C_cluster", "rs11188059", "chr10", 94709142, "G", "A", 2, "snv", "aucun", "haplotype CYP2C:TG, phasage nécessaire si utilisé"),
 ("CYP2C_cluster", "rs2860840", "chr10", 94735475, "C", "T", 2, "snv", "aucun", "haplotype CYP2C:TG, phasage nécessaire si utilisé"),
 ("GGCX", "rs11676382", "chr2", 85550510, "C", "G", 2, "snv", "aucun", "modulateur de la dose d'AVK"),
 ("MTHFR", "rs1801131", "chr1", 11794419, "T", "G", 2, "snv", "aucun", "p.Glu429Ala (A1298C)"),
 ("MTHFR", "rs1801133", "chr1", 11796321, "G", "A,C", 2, "snv", "aucun", "p.Ala222Val (C677T)"),
 ("POR", "rs1057868", "chr7", 75985688, "C", "T", 2, "snv", "aucun", "POR*28, gène absent des tables de l'interpréteur"),
 ("TYMS", "rs45445694", "chr18", 657646, "C", "<VNTR>", 2, "vntr", "aucun",
  "VNTR 28 pb du promoteur (2R à 9R), chr18:657646-657712, couverture seule"),
 ("TYMS", "rs2853542", "chr18", 657685, "G", "A,C,T", 2, "snv", "aucun",
  "situé dans le VNTR, appel peu fiable en lectures courtes"),
 ("TYMS", "rs11280056", "chr18", 673444, "TTAAAGTTA", "TTA,TTAAAGTTAAAGTTA", 2, "indel", "aucun",
  "indel 6 pb du 3'UTR, la représentation peut différer selon l'appelant"),
]
FIN = {"rs45445694": 657712, "rs11280056": 673452}   # borne droite des evenements multi-bases

P.sort(key=lambda x: (x[2].replace("chrM", "chrZZ"), x[3]))
contigs = sorted({x[2] for x in P}, key=lambda c: (c == "chrM", c))

with open(os.path.join(RES, "rnpgx_complement.vcf"), "w", newline="\n", encoding="utf-8") as fh:
    fh.write("##fileformat=VCFv4.2\n")
    fh.write("##source=Core panel RNPGx 2026 (Picard N. et coll., PMID 41491663, Table 2), "
             "positions de classe 1 et 2 absentes des definitions PharmCAT\n")
    fh.write("##reference=GRCh38\n")
    fh.write("##note=Coordonnees GRCh38 via Ensembl ; allele de reference verifie contre "
             "GRCh38_full_analysis_set_plus_decoy_hla.fa\n")
    for c in contigs:
        fh.write(f'##contig=<ID={c},assembly=GRCh38>\n')
    fh.write('##INFO=<ID=PX,Number=.,Type=String,Description="Gene">\n')
    fh.write('##INFO=<ID=PXCLASSE,Number=1,Type=Integer,Description="Classe d actionnabilite RNPGx : 1 ou 2">\n')
    fh.write('##INFO=<ID=PXTYPE,Number=1,Type=String,Description="snv, indel ou vntr">\n')
    fh.write('##INFO=<ID=PXINTERP,Number=1,Type=String,Description="pharmcat si le gene est interprete mais pas cette position ; aucun sinon">\n')
    fh.write('##INFO=<ID=PXNOTE,Number=1,Type=String,Description="Commentaire">\n')
    fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Borne droite des evenements multi-bases">\n')
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRNPGx\n")
    for g, rs, ch, pos, ref, alt, cl, ty, it, note in P:
        info = f"PX={g};PXCLASSE={cl};PXTYPE={ty};PXINTERP={it}"
        if rs in FIN:
            info += f";END={FIN[rs]}"
        # INFO n'admet ni espace, ni ';', ni '=' : on encode
        propre = note.replace(";", ",").replace("=", "-").replace(" ", "_")
        info += f";PXNOTE={propre}"
        fh.write(f"{ch}\t{pos}\t{rs}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0/0\n")

# BED des positions exactes (span complet pour indel et VNTR) et BED des regions +-50 pb
exact, reg = [], []
for g, rs, ch, pos, ref, alt, cl, ty, it, note in P:
    fin = FIN.get(rs, pos + len(ref) - 1)
    exact.append((ch, pos - 1, fin))
    reg.append((ch, max(0, pos - 1 - PAD), fin + PAD))


def fusionne(iv):
    out = []
    for ch, s, e in sorted(iv):
        if out and out[-1][0] == ch and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([ch, s, e])
    return out


for nom, iv in (("rnpgx_complement_positions.bed", exact), ("rnpgx_complement.bed", fusionne(reg))):
    with open(os.path.join(RES, nom), "w", newline="\n", encoding="utf-8") as fh:
        for ch, s, e in (iv if nom.endswith("_positions.bed") else iv):
            fh.write(f"{ch}\t{s}\t{e}\n")

par_gene = collections.Counter(x[0] for x in P)
print("positions :", len(P), "| genes :", len(par_gene), dict(par_gene))
print("classe 1 :", sum(1 for x in P if x[6] == 1), "| classe 2 :", sum(1 for x in P if x[6] == 2))
print("regions BED :", len(fusionne(reg)), "| pb couvertes :", sum(e - s for _, s, e in fusionne(reg)))
print("contigs :", contigs)
