#!/usr/bin/env python3
"""Variants et sequences d'un haplotype d'assemblage HPRC, lus dans son
alignement winnowmap sur GRCh38 (BAM distant, lu par plages).

Pour chaque region du BED des positions, statut de couverture et variants de
l'haplotype contre la reference. Pour chaque locus nomme, la sequence de
l'haplotype entre les deux bornes du locus, a l'usage d'un typage sur
assemblage (CYP2D6, HLA).

Statuts : complet (un seul alignement non secondaire couvre toute la region,
sans N), partiel, absent, ambigu (plusieurs alignements couvrent), N.
"""
import argparse
import os
import subprocess
import sys
import tempfile

import pysam

GRAPPE = 1_000_000  # regions a moins d'un megabase partagent une lecture distante


def lire_plage(url, chrom, debut, fin, dossier):
    """Alignements non secondaires d'une plage du BAM distant. Le pysam de ce
    poste ne lit pas les URL : samtools rapatrie la tranche, pysam la lit."""
    tranche = os.path.join(dossier, "tranche.bam")
    subprocess.run(["samtools", "view", "-b", "-o", tranche, url,
                    f"{chrom}:{debut + 1}-{fin}"], check=True)
    with pysam.AlignmentFile(tranche, "rb", check_sq=False) as f:
        recs = [r for r in f if not r.is_secondary and not r.is_unmapped]
    os.remove(tranche)
    return recs


def lire_bed(chemin):
    regions = []
    with open(chemin) as f:
        for ligne in f:
            if not ligne.strip() or ligne.startswith(("#", "track")):
                continue
            c = ligne.rstrip("\n").split("\t")
            nom = c[3] if len(c) > 3 else f"{c[0]}:{int(c[1]) + 1}-{c[2]}"
            regions.append((c[0], int(c[1]), int(c[2]), nom))
    return regions


def grappes(regions):
    """Regroupe les regions proches pour ne lire chaque plage distante qu'une fois."""
    par_chrom = {}
    for r in regions:
        par_chrom.setdefault(r[0], []).append(r)
    for chrom, rs in par_chrom.items():
        rs.sort(key=lambda r: r[1])
        courant = [rs[0]]
        for r in rs[1:]:
            if r[1] - courant[-1][2] <= GRAPPE:
                courant.append(r)
            else:
                yield chrom, courant
                courant = [r]
        yield chrom, courant


def ref_vers_requete(rec, pos):
    """Indice dans query_sequence de la base alignee sur la position pos (0-based)."""
    rp, qp = rec.reference_start, 0
    for op, n in rec.cigartuples:
        if op in (0, 7, 8):
            if rp <= pos < rp + n:
                return qp + (pos - rp)
            rp += n
            qp += n
        elif op == 1 or op == 4:
            qp += n
        elif op in (2, 3):
            if rp <= pos < rp + n:
                return qp
            rp += n
    return None


def coord_contig(rec, q):
    """Coordonnee, dans l'orientation du FASTA de l'assemblage, de l'indice q
    de query_sequence (les coupures dures ne figurent pas dans la sequence)."""
    ct = rec.cigartuples
    h_debut = ct[0][1] if ct[0][0] == 5 else 0
    pos = h_debut + q
    return rec.infer_read_length() - 1 - pos if rec.is_reverse else pos


def chainer(chevauchent, s, e):
    """Locus decoupe en plusieurs alignements d'un meme contig (cas du CMH) :
    bornes du locus sur le contig, en suivant la chaine colineaire la plus longue."""
    groupes = {}
    for r in chevauchent:
        groupes.setdefault((r.query_name, r.is_reverse), []).append(r)
    recouvrement = lambda g: sum(min(r.reference_end, e) - max(r.reference_start, s) for r in g)
    g = sorted(max(groupes.values(), key=recouvrement), key=lambda r: r.reference_start)
    gauche, droite = g[0], g[-1]
    ps, pe = max(s, gauche.reference_start), min(e, droite.reference_end) - 1
    ca = coord_contig(gauche, ref_vers_requete(gauche, ps))
    cb = coord_contig(droite, ref_vers_requete(droite, pe))
    debuts = [coord_contig(r, 0) for r in g]
    colineaire = debuts == sorted(debuts, reverse=g[0].is_reverse)
    return {"contig": g[0].query_name, "inverse": g[0].is_reverse,
            "a": min(ca, cb), "b": max(ca, cb) + 1, "colineaire": colineaire,
            "tronque": (ps - s) + (e - 1 - pe), "morceaux": len(g),
            "couvert": recouvrement(g) / (e - s)}


def sequence_contig(fasta, contig, a, b, inverse):
    cmd = ["samtools", "faidx"] + (["-i"] if inverse else []) + [fasta, f"{contig}:{a + 1}-{b}"]
    sortie = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    return "".join(sortie.splitlines()[1:])


def variants(rec, seq, ref, chrom, s, e):
    """Variants de l'alignement dans [s, e) : (pos 1-based, ref, alt), et nombre de N."""
    sortie, n_inconnues = [], 0
    rp, qp = rec.reference_start, 0
    for op, n in rec.cigartuples:
        if rp >= e and op in (0, 2, 3, 7, 8):
            break
        if op in (0, 7, 8):
            a, b = max(rp, s), min(rp + n, e)
            if a < b:
                r_seq = ref.fetch(chrom, a, b).upper()
                q_seq = seq[qp + (a - rp): qp + (b - rp)].upper()
                for i, (x, y) in enumerate(zip(r_seq, q_seq)):
                    if y not in "ACGT":
                        n_inconnues += 1
                    elif x != y and x in "ACGT":
                        sortie.append((a + i + 1, x, y))
            rp += n
            qp += n
        elif op == 1:
            if s <= rp - 1 < e:
                ancre = ref.fetch(chrom, rp - 1, rp).upper()
                insere = seq[qp: qp + n].upper()
                n_inconnues += sum(1 for c in insere if c not in "ACGT")
                sortie.append((rp, ancre, ancre + insere))
            qp += n
        elif op == 2:
            if rp < e and rp + n > s:
                bloc = ref.fetch(chrom, rp - 1, rp + n).upper()
                sortie.append((rp, bloc, bloc[0]))
            rp += n
        elif op == 3:
            rp += n
        elif op == 4:
            qp += n
    return sortie, n_inconnues


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--positions", required=True, help="BED des regions a genotyper")
    ap.add_argument("--loci", required=True, help="BED des loci a extraire en sequence")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--fasta", required=True, help="FASTA bgzippe de l'assemblage, avec .fai et .gzi")
    ap.add_argument("--prefixe", required=True, help="prefixe des fichiers de sortie")
    a = ap.parse_args()

    ref = pysam.FastaFile(a.ref)
    dossier = tempfile.mkdtemp(prefix="tranche_", dir=os.path.dirname(os.path.abspath(a.prefixe)))
    positions = lire_bed(a.positions)
    loci = lire_bed(a.loci)
    toutes = [(c, s, e, ("P", n)) for c, s, e, n in positions] + \
             [(c, s, e, ("L", n)) for c, s, e, n in loci]

    statuts, vcf, fasta = [], [], []
    for chrom, groupe in grappes(toutes):
        debut = min(r[1] for r in groupe)
        fin = max(r[2] for r in groupe)
        recs = lire_plage(a.bam, chrom, debut, fin, dossier)
        sequences = {}
        for _, s, e, (genre, nom) in groupe:
            chevauchent = [r for r in recs if r.reference_start < e and r.reference_end > s]
            couvrent = [r for r in chevauchent if r.reference_start <= s and r.reference_end >= e]
            if not chevauchent:
                statut, rec = "absent", None
            elif not couvrent:
                statut, rec = "partiel", None
            elif len(couvrent) > 1:
                statut, rec = "ambigu", None
            else:
                statut, rec = "complet", couvrent[0]
            contigs = ",".join(sorted({f"{r.query_name}:{r.mapping_quality}" for r in chevauchent}))
            if genre == "L" and rec is None and chevauchent:
                c = chainer(chevauchent, s, e)
                douteux = not c["colineaire"] or c["tronque"] > 10_000 or not 0.5 < (c["b"] - c["a"]) / (e - s) < 2
                statut = "chaine_douteuse" if douteux else "chaine"
                seq = sequence_contig(a.fasta, c["contig"], c["a"], c["b"], c["inverse"])
                fasta.append((nom, f"{c['contig']}:{c['a']}-{c['b']} chaine morceaux={c['morceaux']} "
                                   f"couvert={c['couvert']:.2f} tronque={c['tronque']} statut={statut}", seq))
                contigs = f"{contigs};morceaux={c['morceaux']};couvert={c['couvert']:.2f};tronque={c['tronque']}"
            if rec is not None:
                cle = (rec.query_name, rec.reference_start, rec.flag)
                if cle not in sequences:
                    sequences[cle] = rec.query_sequence
                seq = sequences[cle]
                if seq is None:
                    statut = "sans_sequence"
                elif genre == "P":
                    vs, n_inc = variants(rec, seq, ref, chrom, s, e)
                    if n_inc:
                        statut = "N"
                    else:
                        vcf.extend((chrom, p, r_, alt) for p, r_, alt in vs)
                else:
                    qa, qb = ref_vers_requete(rec, s), ref_vers_requete(rec, e - 1)
                    morceau = seq[min(qa, qb): max(qa, qb) + 1]
                    brin = "-" if rec.is_reverse else "+"
                    fasta.append((nom, f"{rec.query_name}:{min(qa, qb)}-{max(qa, qb) + 1} brin_bam{brin}", morceau))
            statuts.append((genre, nom, chrom, s, e, statut, contigs))

    with open(a.prefixe + ".statuts.tsv", "w") as f:
        f.write("genre\tnom\tchrom\tdebut\tfin\tstatut\talignements\n")
        for l in statuts:
            f.write("\t".join(map(str, l)) + "\n")
    with open(a.prefixe + ".variants.tsv", "w") as f:
        for chrom, p, r_, alt in sorted(set(vcf), key=lambda v: (v[0], v[1])):
            f.write(f"{chrom}\t{p}\t{r_}\t{alt}\n")
    with open(a.prefixe + ".loci.fa", "w") as f:
        for nom, desc, seq in fasta:
            f.write(f">{nom} {desc}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
    os.rmdir(dossier)
    comptes = {}
    for l in statuts:
        comptes[l[5]] = comptes.get(l[5], 0) + 1
    print(os.path.basename(a.prefixe), comptes, f"{len(set(vcf))} variants", file=sys.stderr)


if __name__ == "__main__":
    main()
