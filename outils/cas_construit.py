# -*- coding: utf-8 -*-
"""Fabrique un fichier de variants portant un allele donne, pour exercer un
chemin que la cohorte de validation ne contient pas.

Pourquoi c'est necessaire. Sur deux cent quarante-trois genomes publics, les
trois positions de MT-RNR1 retenues par le panel sont de reference partout :
aucun porteur. Or MT-RNR1 est le seul gene de l'etage de typage complementaire
a porter un phenotype et une recommandation. L'appel negatif est donc verifie
deux cent quarante-deux fois, et l'appel positif pas une seule — ni le typage,
ni l'appel externe vers l'interpreteur, ni la restitution de la consigne au
compte rendu.

Ce que cela vaut, et ce que cela ne vaut pas. Un variant injecte dans le
fichier de variants valide la chaine logicielle : le typeur nomme-t-il
l'allele, l'interpreteur rend-il le phenotype, la consigne arrive-t-elle au
document. Il ne valide pas la detection : l'alignement ne porte pas le variant,
et aucune conclusion sur la sensibilite analytique ne peut en etre tiree. Un
materiau de reference porteur reste necessaire pour cela.

Usage :
  python3 outils/cas_construit.py --variants E.vcf.gz --sortie E_mtrnr1.vcf.gz \\
      --allele m.1555A>G
"""
import argparse
import gzip
import os
import subprocess
import sys

# Les trois positions de classe 1 du panel pour MT-RNR1, telles que le
# complement RNPGx les nomme.
ALLELES = {
    "m.1555A>G": ("chrM", 1555, "A", "G", "rs267606617"),
    "m.1494C>T": ("chrM", 1494, "C", "T", "rs267606619"),
    "m.1095T>C": ("chrM", 1095, "T", "C", "rs267606618"),
}


def ouvrir(chemin):
    if chemin.endswith(".gz"):
        return gzip.open(chemin, "rt", encoding="utf-8")
    return open(chemin, encoding="utf-8")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--variants", required=True)
    p.add_argument("--sortie", required=True)
    p.add_argument("--allele", required=True, choices=sorted(ALLELES))
    p.add_argument("--genotype", default="1/1",
                   help="l'ADN mitochondrial est haploide ; 1/1 pour un "
                        "porteur homoplasmique")
    p.add_argument("--gq", type=int, default=99)
    a = p.parse_args()

    contig, pos, ref, alt, rsid = ALLELES[a.allele]
    lignes, entete, vu = [], [], False
    with ouvrir(a.variants) as fh:
        for l in fh:
            if l.startswith("#"):
                entete.append(l)
                continue
            c = l.split("\t", 2)
            if c[0] == contig and int(c[1]) == pos:
                vu = True
                continue      # on remplace l'enregistrement existant
            lignes.append(l)
    if not entete:
        sys.exit("fichier de variants sans entete")

    neuf = "\t".join([contig, str(pos), rsid, ref, alt, ".", "PASS", ".",
                      "GT:GQ:DP", "%s:%d:200" % (a.genotype, a.gq)]) + "\n"
    lignes.append(neuf)

    # bcftools sort remet l'enregistrement a sa place et indexe : un VCF non
    # trie est refuse par les etages suivants.
    tmp = a.sortie + ".brut.vcf"
    with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
        fh.writelines(entete)
        fh.writelines(lignes)
    r = subprocess.run(["bcftools", "sort", "-Oz", "-o", a.sortie, tmp],
                       capture_output=True)
    os.unlink(tmp)
    if r.returncode != 0:
        sys.exit(r.stderr.decode("utf-8", "replace")[-300:])
    subprocess.run(["bcftools", "index", "-f", "-t", a.sortie], check=False)
    print("%s : %s injecte en %s (%s%s)"
          % (os.path.basename(a.sortie), a.allele, a.genotype, rsid,
             ", enregistrement existant remplace" if vu else ""))


if __name__ == "__main__":
    main()
