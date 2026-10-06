# -*- coding: utf-8 -*-
"""Construit ressources/pypgx_*.{json,bed} : le perimetre de l'etage 4b.

Les quatre genes traites ici sont ceux du core panel RNPGx que l'interpreteur
PharmCAT ne sait pas appeler mais que PyPGx sait typer. Les positions et les
regions ne sont pas saisies a la main : elles sont lues dans PyPGx lui-meme, de
sorte qu'une montee de version du typeur met a jour le perimetre mesure.

MT-RNR1 est le seul a porter un phenotype. PharmCAT embarque sa table de
recommandations sans aucun appeleur : le resultat lui est rendu en appel
externe, comme HLA et CYP2D6.

Prerequis : pypgx installe et PYPGX_BUNDLE defini.
Usage : python3 outils/construire_pypgx.py [dossier_ressources]
"""
import json, os, sys, datetime

import pypgx

RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ressources")
PAD = 50   # meme marge que pharmcat_positions.bed

# gene -> (classe RNPGx, appel externe rendu a l'interpreteur, motif de presence)
GENES = {
    "MT-RNR1": (1, True,  "surdite sous aminoside ; PharmCAT a la table, pas l'appeleur"),
    "BCHE":    (1, False, "allele atypique et variante K ; gene absent de PharmCAT"),
    "MTHFR":   (2, False, "C677T et A1298C ; gene absent de PharmCAT"),
    "POR":     (2, False, "POR*28 ; gene absent de PharmCAT"),
}


def region(gene):
    """PyPGx nomme les contigs sans prefixe : M, 3, 1, 7. On revient au nommage
    GRCh38 de reference, seul present dans les alignements de la plateforme."""
    c, bornes = pypgx.get_region(gene, assembly="GRCh38").split(":")
    debut, fin = (int(x) for x in bornes.split("-"))
    return "chr" + c, debut, fin


def positions_panel(gene):
    """Les positions que le panel RNPGx nomme pour ce gene, lues dans le
    complement. Le typeur en connait bien davantage : sur MTHFR il nomme neuf
    alleles quand le panel n'en retient que deux, C677T et A1298C. Le compte
    rendu s'en tient au panel, et c'est cette liste qui le lui dit."""
    out = set()
    chemin = os.path.join(RES, "rnpgx_complement.vcf")
    if not os.path.exists(chemin):
        return out
    with open(chemin, encoding="utf-8") as fh:
        for ligne in fh:
            if ligne.startswith("#"):
                continue
            c = ligne.split("\t")
            if len(c) > 7 and ("PX=%s;" % gene) in c[7]:
                out.add("%s:%s" % (c[0], c[1]))
    return out


def positions_des_alleles(gene):
    """Quelle position definit quel allele. C'est ce qui permet de ne retenir
    que les alleles du panel sans jamais en ecarter un en silence."""
    out = {}
    for a in pypgx.list_alleles(gene):
        p = set()
        for mode in ("core", "tag"):
            try:
                v = pypgx.list_variants(gene, alleles=[a], mode=mode,
                                        assembly="GRCh38")
            except Exception:
                continue
            for x in v:
                c, pos = x.split("-")[:2]
                p.add("chr%s:%s" % (c, pos))
        if p:
            out[a] = sorted(p)
    return out


def positions(gene):
    """Les positions definissantes des alleles, au format PyPGx <contig>-<pos>-<ref>-<alt>.
    Ce sont elles, et non la region entiere, qui doivent etre lues pour qu'un
    appel « reference » soit une lecture et non une absence."""
    out = []
    for v in pypgx.list_variants(gene, mode="all", assembly="GRCh38"):
        c, pos, ref, alt = v.split("-", 3)
        out.append(("chr" + c, int(pos), ref, alt))
    # Le tri doit etre total : deux alleles a la meme coordonnee seraient
    # sinon departages par l'ordre d'iteration de l'ensemble, qui varie d'une
    # execution a l'autre. La coordonnee reste la cle principale, pour que le
    # fichier se lise dans l'ordre du genome.
    return sorted(set(out), key=lambda x: (x[1], x[2], x[3]))


catalogue, exact, larges = {}, [], []
for gene, (classe, externe, motif) in GENES.items():
    contig, debut, fin = region(gene)
    pos = positions(gene)
    alleles = pypgx.list_alleles(gene)
    phen = str(pypgx.load_gene_table().set_index("Gene").loc[gene].get("PhenotypeMethod"))
    catalogue[gene] = {
        "classe_rnpgx": classe,
        "motif": motif,
        "contig": contig,
        "debut": debut,
        "fin": fin,
        "alleles": len(alleles),
        # Le nom que le gene donne a son allele de reference : « *1 » pour
        # un gene a etoiles, « Reference » pour ceux nommes par rsID.
        "allele_reference": "*1" if "*1" in alleles else "Reference",
        "positions": ["%s:%d:%s:%s" % p for p in pos],
        # Le perimetre du compte rendu : les positions que le panel nomme,
        # et les alleles du typeur qui tombent dessus.
        "positions_panel": sorted(positions_panel(gene)),
        "alleles_positions": positions_des_alleles(gene),
        # PyPGx ne rend un phenotype que pour MT-RNR1 ; ailleurs il nomme
        # l'allele sans lui attacher de sens clinique, et le compte rendu doit
        # le dire plutot que de laisser croire a une interpretation.
        "phenotype": phen if phen and phen != "nan" else None,
        "appel_externe_pharmcat": externe,
        # MT-RNR1 est haploide : l'interpreteur attend un allele, pas une paire.
        "haploide": gene == "MT-RNR1",
    }
    for contig_p, p, ref, alt in pos:
        exact.append((contig_p, p - 1, p - 1 + len(ref)))
        larges.append((contig_p, max(0, p - 1 - PAD), p - 1 + len(ref) + PAD))


def fusionne(iv):
    out = []
    for ch, s, e in sorted(iv):
        if out and out[-1][0] == ch and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([ch, s, e])
    return out


with open(os.path.join(RES, "pypgx_genes.json"), "w", newline="\n", encoding="utf-8") as fh:
    json.dump({
        "source": "PyPGx %s" % getattr(pypgx, "__version__", "?"),
        "genere_le": datetime.date.today().isoformat(),
        "genes": catalogue,
    }, fh, ensure_ascii=False, indent=1)
    fh.write("\n")

# Positions exactes : ce qui doit etre lu pour qu'un appel soit recevable.
with open(os.path.join(RES, "pypgx_positions.bed"), "w", newline="\n", encoding="utf-8") as fh:
    for ch, s, e in fusionne(exact):
        fh.write("%s\t%d\t%d\n" % (ch, s, e))

# Regions de tranche : les lectures a extraire de l'alignement, +-50 pb.
with open(os.path.join(RES, "pypgx_tranche.bed"), "w", newline="\n", encoding="utf-8") as fh:
    for ch, s, e in fusionne(larges):
        fh.write("%s\t%d\t%d\n" % (ch, s, e))

# Regions completes des genes : l'extraction du fichier de variants. PyPGx lit
# la region entiere et non les seules positions definissantes.
with open(os.path.join(RES, "pypgx_regions.bed"), "w", newline="\n", encoding="utf-8") as fh:
    for gene in sorted(catalogue):
        d = catalogue[gene]
        fh.write("%s\t%d\t%d\t%s\n" % (d["contig"], d["debut"] - 1, d["fin"], gene))

tot = sum(len(d["positions"]) for d in catalogue.values())
print("genes :", len(catalogue), "| positions definissantes :", tot)
for g in sorted(catalogue):
    d = catalogue[g]
    print("  %-9s classe %d  %3d alleles  %2d positions  %6d pb  phenotype: %s%s"
          % (g, d["classe_rnpgx"], d["alleles"], len(d["positions"]),
             d["fin"] - d["debut"], d["phenotype"] or "aucun",
             "  -> appel externe" if d["appel_externe_pharmcat"] else ""))
print("BED positions :", len(fusionne(exact)), "| BED tranche :", len(fusionne(larges)),
      "pour", sum(e - s for _, s, e in fusionne(larges)), "pb")
