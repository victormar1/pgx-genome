# -*- coding: utf-8 -*-
"""Tests de la regle qui choisit l'allele CYP2D6 de la verite par assemblage.

deCYPher peut donner plusieurs alleles dont toutes les variations definissantes
sont presentes dans une copie. En retenir un au hasard, ou le premier, change le
phenotype de la verite : c'est la verite elle-meme qui devient fausse, et le
module est juge contre elle. Deux regles plus simples echouent chacune sur un
cas reel ; elles sont testees ici.
"""
import io
import os
import shutil
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bench", "orthogonale"))

import verite_cyp2d6 as v6  # noqa: E402

# Definitions de deCYPher, telles qu'elles sont dans ses VCF de reference.
DEFINITIONS = {
    "*2": ["7870 C T", "9200 G C"],
    "*4": ["6866 G A"],
    "*10": ["5119 C T", "9200 G C"],
    "*49": ["5119 C T", "6631 T A", "9200 G C"],
    "*56": ["8221 C T", "9200 G C"],
    "*84": ["7594 C A", "7870 C T", "9200 G C"],
}


class Regle(unittest.TestCase):
    """La regle seule, sans passer par une sortie de deCYPher."""

    def setUp(self):
        self.trav = tempfile.mkdtemp()
        d = os.path.join(self.trav, "outils", "deCYPher", "ref", "CYP2D6", "RefSeqGene")
        os.makedirs(d)
        for etoile, lignes in DEFINITIONS.items():
            with io.open(os.path.join(d, "CYP2D6.%s.vcf" % etoile[1:]), "w", newline="\n") as f:
                f.write("##fileformat=VCFv4.2\n")
                for l in lignes:
                    pos, ref, alt = l.split()
                    f.write("NG_008376.4\t%s\t.\t%s\t%s\t.\t.\t.\n" % (pos, ref, alt))

    def tearDown(self):
        shutil.rmtree(self.trav, ignore_errors=True)

    def retenir(self, candidats, rendu):
        return v6._retenir(candidats, rendu, self.trav)

    def test_allele_le_plus_specifique_gagne(self):
        """*49 porte toutes les variations de *10 et une de plus. Retenir *10
        rendrait une copie a 0,25 d'activite au lieu de 0,5."""
        self.assertEqual(self.retenir(["*10", "*49"], "*49"), "*49")

    def test_le_plus_specifique_gagne_meme_sans_activite_chiffree(self):
        """*84 contient *2 et n'a pas de valeur CPIC. Choisir par la seule
        activite le remplacait par *2."""
        self.assertEqual(self.retenir(["*2", "*84"], "*84"), "*84")

    def test_a_specificite_egale_le_moins_actif_gagne(self):
        """*56 et *10 ne se contiennent pas. *56 est sans fonction : c'est la
        variation la plus deletere qui decide, comme dans PharmVar."""
        self.assertEqual(self.retenir(["*10", "*56"], "*10"), "*56")

    def test_definition_plus_courte_ne_perd_pas(self):
        """*4 n'a qu'une variation definissante chez deCYPher, *10 en a deux.
        Une regle qui compte les variations rendait *10, un allele a fonction
        diminuee, la ou la copie est sans fonction."""
        self.assertEqual(self.retenir(["*10", "*4"], "*4"), "*4")

    def test_appel_de_decypher_garde_a_egalite(self):
        self.assertEqual(self.retenir(["*4", "*56"], "*56"), "*56")

    def test_candidat_inconnu_du_dossier_de_definitions(self):
        """Un allele sans fichier de definition ne doit pas etre ecarte en
        silence : il reste candidat."""
        self.assertIn(self.retenir(["*10", "*999"], "*999"), ("*999", "*10"))


class Ecriture(unittest.TestCase):
    def test_duplication(self):
        self.assertEqual(v6.haplotype(["*4", "*4"]), "*4x2")

    def test_copie_unique(self):
        self.assertEqual(v6.haplotype(["*4"]), "*4")

    def test_tandem_hybride_en_premier(self):
        self.assertEqual(v6.haplotype(["*10", "*36"]), "*36+*10")

    def test_haplotype_vide(self):
        self.assertIsNone(v6.haplotype(None))


class SansPhase(unittest.TestCase):
    """Deux ecritures d'un meme jeu de copies ne sont pas un desaccord."""

    def test_tandem_contre_duplication(self):
        self.assertEqual(v6.multiensemble("*10/*36x2"), v6.multiensemble("*36/*36+*10"))

    def test_phases_differentes(self):
        self.assertEqual(v6.multiensemble("*2/*36+*10"), v6.multiensemble("*36+*2/*10"))

    def test_copies_differentes(self):
        self.assertNotEqual(v6.multiensemble("*1/*4"), v6.multiensemble("*1/*10"))

    def test_sous_allele_ramene_a_son_allele(self):
        self.assertEqual(v6.multiensemble("*4.021/*10"), ["*10", "*4"])

    def test_ecriture_illisible(self):
        self.assertIsNone(v6.multiensemble("indetermine"))


class Activite(unittest.TestCase):
    def test_table_chargee(self):
        self.assertEqual(v6.ACTIVITE["*1"], 1.0)
        self.assertEqual(v6.ACTIVITE["*4"], 0.0)

    def test_allele_sans_valeur_absent_de_la_table(self):
        """« n/a » ne doit pas devenir 0 : un allele de fonction inconnue
        deviendrait le plus severe de tous."""
        self.assertNotIn("*84", v6.ACTIVITE)


if __name__ == "__main__":
    unittest.main()
