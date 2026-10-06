# -*- coding: utf-8 -*-
"""Tests de l'etage des appels externes.

C'est l'etage le plus discret et l'un des plus dangereux : une faute d'ecriture
n'y leve aucune erreur, l'interpreteur rend simplement « indetermine » et le
compte rendu perd une recommandation sans que rien ne le signale.
"""
import io
import os
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))

import appels_externes as a  # noqa: E402

ABSENT = os.path.join(tempfile.gettempdir(), "absent_xyz_appels.tsv")


def fichier(contenu):
    f = tempfile.mktemp(suffix=".tsv")
    with io.open(f, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(contenu)
    return f


class Cyp2d6(unittest.TestCase):
    def lis(self, genotype):
        f = fichier("Genotype\tFilter\n%s\tPASS\n" % genotype)
        try:
            return a.cyp2d6(f)
        finally:
            os.unlink(f)

    def test_diplotype_simple(self):
        self.assertEqual(self.lis("*1/*4"), ("*1/*4", None))

    def test_tandem_entoure_d_espaces_autour_du_plus(self):
        # « *36+*10 » fait rendre « indetermine » a l'interpreteur, sans erreur ;
        # « *36 + *10 » est interprete. La faute ne se voit qu'au compte rendu.
        g, err = self.lis("*36+*10/*4")
        self.assertIsNone(err)
        self.assertEqual(g, "*36 + *10/*4")

    def test_espaces_deja_presents_ne_sont_pas_doubles(self):
        self.assertEqual(self.lis("*36 + *10/*4")[0], "*36 + *10/*4")

    def test_refuse_un_diplotype_ambigu(self):
        # Cyrius separe par un point-virgule les lectures entre lesquelles il
        # n'a pas tranche ; en choisir une serait affirmer ce qu'il n'a pas dit,
        # et l'interpreteur s'arrete sur une telle ligne, perdant le genome.
        for brut in ("*1/*2;*1/*41", "*1/*2,*1/*41", "*1/*2/*3"):
            g, err = self.lis(brut)
            self.assertIsNone(g, brut)
            self.assertIn("plusieurs diplotypes", err)

    def test_aucun_diplotype(self):
        for brut in ("", "None"):
            g, err = self.lis(brut)
            self.assertIsNone(g)
            self.assertIn("aucun diplotype", err)

    def test_fichier_absent_ou_vide(self):
        g, err = a.cyp2d6(ABSENT)
        self.assertIsNone(g)
        self.assertIn("absent", err)
        f = fichier("Genotype\tFilter\n")
        try:
            g, err = a.cyp2d6(f)
        finally:
            os.unlink(f)
        self.assertIsNone(g)
        self.assertIn("vide", err)


class Hla(unittest.TestCase):
    def lis(self, a1, a2, b1, b2):
        f = fichier("A1\tA2\tB1\tB2\n%s\t%s\t%s\t%s\n" % (a1, a2, b1, b2))
        try:
            return a.hla(f)
        finally:
            os.unlink(f)

    def test_deux_genes(self):
        out, err = self.lis("A*31:01", "A*02:01", "B*57:01", "B*07:02")
        self.assertIsNone(err)
        self.assertEqual(out["HLA-A"], "*31:01/*02:01")
        self.assertEqual(out["HLA-B"], "*57:01/*07:02")

    def test_tronque_a_deux_champs(self):
        # L'interpreteur attend deux champs ; un troisieme le ferait echouer,
        # et le troisieme champ n'a aucune consequence de prescription.
        out, err = self.lis("A*31:01:02", "A*02:01:01:01",
                            "B*57:01:01", "B*07:02:01")
        self.assertIsNone(err)
        self.assertEqual(out["HLA-A"], "*31:01/*02:01")
        self.assertEqual(out["HLA-B"], "*57:01/*07:02")

    def test_un_allele_illisible_fait_tomber_son_gene_seul(self):
        # Rendre une ligne a demi remplie vaudrait mieux que rien pour un
        # humain, mais l'interpreteur ne lit pas une paire incomplete.
        out, err = self.lis("A*31:01", "", "B*57:01", "B*07:02")
        self.assertNotIn("HLA-A", out)
        self.assertIn("HLA-B", out)
        self.assertIn("HLA-A", err)

    def test_fichier_absent(self):
        out, err = a.hla(ABSENT)
        self.assertEqual(out, {})
        self.assertIn("absent", err)


class MtRnr1(unittest.TestCase):
    ENTETE = "gene\tclasse_rnpgx\tstatut\tallele_externe\tmotif\n"

    def lis(self, lignes):
        f = fichier(self.ENTETE + lignes)
        try:
            return a.mtrnr1(f)
        finally:
            os.unlink(f)

    def test_allele_rendu(self):
        # Gene haploide : un allele seul, jamais une paire.
        g, err = self.lis("MT-RNR1\t1\trendu\tm.1555A>G\t\n")
        self.assertIsNone(err)
        self.assertEqual(g, "m.1555A>G")
        self.assertNotIn("/", g)

    def test_gene_non_conclusif_n_est_pas_repris(self):
        # Reprendre un gene non conclusif ferait rendre « risque normal de
        # surdite sous aminoside » sur un gene que l'etage a refuse de typer.
        g, err = self.lis("MT-RNR1\t1\tnon conclusif\t\tcontig chrM absent\n")
        self.assertIsNone(g)
        self.assertIn("chrM", err)

    def test_statut_rendu_mais_allele_vide(self):
        g, err = self.lis("MT-RNR1\t1\trendu\t\t\n")
        self.assertIsNone(g)
        self.assertTrue(err)

    def test_gene_absent_du_typage(self):
        g, err = self.lis("BCHE\t1\trendu\t\t\n")
        self.assertIsNone(g)
        self.assertIn("absent", err)

    def test_fichier_absent(self):
        # L'etage 4b est optionnel : son fichier peut legitimement manquer, et
        # l'absence doit etre un motif, pas une exception.
        g, err = a.mtrnr1(ABSENT)
        self.assertIsNone(g)
        self.assertIn("absent", err)


if __name__ == "__main__":
    unittest.main(verbosity=2)
