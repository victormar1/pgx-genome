# -*- coding: utf-8 -*-
"""Doublures des outils externes, pour executer la chaine sans eux.

Le module appelle samtools, un moteur de conteneurs (bcftools, bgzip, OptiType,
PharmCAT) et Cyrius. Aucun n'est disponible en integration continue, et la
chaine n'etait donc verifiee qu'a la lecture de son texte. Ces doublures
rendent des sorties de la bonne forme, sans rien calculer : ce qui se teste est
l'orchestration — l'ordre des etages, leurs etats, les refus, la reprise, le
code de retour — et non la biologie, que mesure le banc.

Le nom sous lequel le fichier est appele choisit la doublure : samtools,
docker ou star_caller.py. Les variables FAUX_* reglent une execution :

  FAUX_SM          identite portee par l'alignement        (defaut HG00699)
  FAUX_PROFONDEUR  profondeur rendue a chaque position     (defaut 30)
  FAUX_PAIRES      paires extraites du complexe majeur     (defaut 500)
  FAUX_D6          diplotype rendu par Cyrius              (defaut *1/*2)
  FAUX_D6_FILTRE   filtre que Cyrius lui attache           (defaut PASS)
  FAUX_PANNE       defaillance a simuler : sort_tue, index_tue, collate,
                   optitype, vcf_illisible, quickcheck, d6_illisible
  FAUX_VCF         defaut du fichier de variants : sans_ad, deux_echantillons
  FAUX_JOURNAL     fichier ou chaque appel inscrit sa ligne de commande
"""
import gzip
import json
import os
import shutil
import signal
import sys

ICI = os.path.dirname(os.path.abspath(__file__))
RACINE = os.path.dirname(os.path.dirname(ICI))
RES = os.path.join(RACINE, "ressources")
RAPPORT = os.path.join(RACINE, "bench", "qualification", "HG00699", "sortie",
                       "HG00699.report.json")
PANNE = os.environ.get("FAUX_PANNE", "")


def journal(outil, args):
    chemin = os.environ.get("FAUX_JOURNAL")
    if chemin:
        with open(chemin, "a", encoding="utf-8") as fh:
            fh.write(json.dumps([outil] + list(args), ensure_ascii=False) + "\n")


def ecrire(chemin, contenu):
    with open(chemin, "wb") as fh:
        fh.write(contenu)


def apres(args, option):
    return args[args.index(option) + 1]


def contigs():
    vus = []
    for nom in sorted(os.listdir(RES)):
        if not nom.endswith(".bed"):
            continue
        with open(os.path.join(RES, nom), encoding="utf-8") as fh:
            for l in fh:
                c = l.split("\t")[0].strip()
                if c and c not in vus:
                    vus.append(c)
    return vus


def entete():
    lignes = ["@HD\tVN:1.6\tSO:coordinate"]
    for c in contigs():
        lignes.append("@SQ\tSN:%s\tLN:300000000" % c)
    lignes.append("@SQ\tSN:HLA-A*01:01:01:01\tLN:3503")
    lignes.append("@RG\tID:1\tSM:%s" % os.environ.get("FAUX_SM", "HG00699"))
    return "\n".join(lignes) + "\n"


def tuer_l_appelant():
    """L'execution est interrompue pendant que l'outil ecrit : c'est ce que
    laisse un noeud de calcul qui retire la tache."""
    os.kill(os.getppid(), signal.SIGTERM)


def samtools(a):
    cmd, a = a[0], a[1:]
    if cmd == "view":
        if "-H" in a:
            sys.stdout.write(entete())
        elif "-c" in a:
            print(1000)
        elif "-b" in a:
            ecrire(apres(a, "-o"), b"alignement extrait")
        else:
            sys.stdout.write("flux")
        return 0
    if cmd == "quickcheck":
        return 1 if PANNE == "quickcheck" else 0
    if cmd == "sort":
        sortie = apres(a, "-o")
        if PANNE == "sort_tue":
            ecrire(sortie, b"tri interrompu")
            tuer_l_appelant()
            return 1
        ecrire(sortie, b"alignement trie")
        return 0
    if cmd == "index":
        if PANNE == "index_tue":
            ecrire(a[-1] + ".bai", b"index interrompu")
            tuer_l_appelant()
            return 1
        ecrire(a[-1] + ".bai", b"index")
        return 0
    if cmd == "depth":
        p = os.environ.get("FAUX_PROFONDEUR", "30")
        with open(apres(a, "-b"), encoding="utf-8") as fh:
            for l in fh:
                c = l.split()
                if len(c) >= 3:
                    for pos in range(int(c[1]) + 1, int(c[2]) + 1):
                        print("%s\t%d\t%s" % (c[0], pos, p))
        return 0
    if cmd == "collate":
        flux = sys.stdin.read()
        if PANNE == "collate":
            return 1
        sys.stdout.write(flux)
        return 0
    if cmd == "fastq":
        sys.stdin.read()
        n = int(os.environ.get("FAUX_PAIRES", "500"))
        if PANNE == "collate":
            n = n // 2        # le flux s'est arrete en route : fichiers tronques
        for option in ("-1", "-2"):
            with open(apres(a, option), "w", encoding="ascii") as fh:
                for i in range(n):
                    fh.write("@l%d\nACGT\n+\nIIII\n" % i)
        return 0
    sys.stderr.write("faux samtools : commande non doublee : %s\n" % cmd)
    return 64


def lignes_vcf(chemin):
    ouvrir = gzip.open if chemin.endswith(("gz", "bgz")) else open
    with ouvrir(chemin, "rt", encoding="utf-8") as fh:
        for l in fh:
            yield l


def bcftools(a):
    cmd, a = a[0], a[1:]
    defaut = os.environ.get("FAUX_VCF", "")
    if cmd == "query" and "-l" in a:
        if PANNE == "vcf_illisible":
            sys.stderr.write("faux bcftools : fichier illisible\n")
            return 1
        for l in lignes_vcf(a[-1]):
            if l.startswith("#CHROM"):
                for nom in l.rstrip("\n").split("\t")[9:]:
                    print(nom)
                if defaut == "deux_echantillons":
                    print("SECOND")
                return 0
        return 1
    if cmd == "view" and "-h" in a:
        for l in lignes_vcf(a[-1]):
            if not l.startswith("#"):
                break
            if defaut == "sans_ad" and l.startswith("##FORMAT=<ID=AD,"):
                continue
            sys.stdout.write(l)
        return 0
    if cmd == "view" and "-H" in a:
        for l in lignes_vcf(a[-1]):
            if not l.startswith("#"):
                sys.stdout.write(l)
        return 0
    if cmd == "view" and "-R" in a:
        # Le fichier de qualification est deja restreint aux regions utiles.
        shutil.copyfile(a[-1], apres(a, "-o"))
        return 0
    if cmd == "index":
        ecrire(a[-1] + ".tbi", b"index")
        return 0
    sys.stderr.write("faux bcftools : commande non doublee : %s\n" % cmd)
    return 64


def conteneur(cmd):
    outil, a = cmd[0], cmd[1:]
    if outil == "bcftools":
        return bcftools(a)
    if outil == "bgzip":
        with open(a[-1], "rb") as src, gzip.open(a[-1] + ".gz", "wb") as dst:
            shutil.copyfileobj(src, dst)
        os.remove(a[-1])
        return 0
    if outil == "sh":
        sys.stdout.write("[mapping]\nrazers3=razers3\nthreads=16\n\n"
                         "[ilp]\nsolver=glpk\nthreads=1\n")
        return 0
    if outil == "OptiTypePipeline.py":
        if PANNE == "optitype":
            return 1
        d = os.path.join(apres(a, "-o"), "2026_01_01_00_00_00")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "2026_01_01_00_00_00_result.tsv"), "w",
                  encoding="utf-8") as fh:
            fh.write("\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n")
            fh.write("0\tA*02:01\tA*11:01\tB*40:01\tB*58:01\tC*03:02\tC*07:02"
                     "\t1000.0\t950.0\n")
        return 0
    if outil == "/pharmcat/pharmcat_vcf_preprocessor":
        nom = os.path.basename(apres(a, "-vcf"))
        if nom.endswith(".gz"):
            nom = nom[:-3]
        ecrire(os.path.join(apres(a, "-o"), nom + ".preprocessed.vcf.bgz"),
               b"vcf pretraite")
        return 0
    if outil == "java":
        shutil.copyfile(RAPPORT, os.path.join(apres(a, "-o"),
                                              apres(a, "-bf") + ".report.json"))
        return 0
    sys.stderr.write("faux conteneur : commande non doublee : %s\n" % outil)
    return 64


def docker(a):
    if a[:2] == ["image", "inspect"]:
        print("faux@sha256:" + "0" * 64)
        return 0
    if a[0] != "run":
        sys.stderr.write("faux docker : commande non doublee : %s\n" % a[0])
        return 64
    i, cwd = 1, None
    while i < len(a) and a[i].startswith("-"):
        if a[i] == "--rm":
            i += 1
            continue
        if a[i] == "-w":
            cwd = a[i + 1]
        i += 2
    # Comme le vrai : le dossier de travail d'un conteneur est un chemin absolu,
    # et la commande s'y execute, pas dans le dossier d'appel de l'hote.
    if not cwd or not os.path.isabs(cwd):
        sys.stderr.write("faux docker : dossier de travail non absolu : %s\n" % cwd)
        return 125
    os.chdir(cwd)
    return conteneur(a[i + 1:])


def cyrius(a):
    d = apres(a, "-o")
    with open(os.path.join(d, "x.tsv"), "w", encoding="utf-8") as fh:
        fh.write("Sample\tGenotype\tFilter\n")
        fh.write("x\t%s\t%s\n" % (os.environ.get("FAUX_D6", "*1/*2"),
                                   os.environ.get("FAUX_D6_FILTRE", "PASS")))
    if PANNE == "d6_illisible":
        ecrire(os.path.join(d, "x.tsv"), b"Sample\tGenotype\tFilter\nx\t\xff\xfe\tPASS\n")
    ecrire(os.path.join(d, "x.json"), b"{}")
    return 0


def main():
    role = os.environ.get("FAUX_ROLE") or os.path.basename(sys.argv[0])
    a = sys.argv[1:]
    journal(role, a)
    if role == "samtools":
        return samtools(a)
    if role == "docker":
        return docker(a)
    if role == "star_caller.py":
        return cyrius(a)
    sys.stderr.write("doublure inconnue : %s\n" % role)
    return 64


if __name__ == "__main__":
    sys.exit(main())
