# -*- coding: utf-8 -*-
"""Tests du denominateur de la mesure du banc orthogonal.

Un couple gene x genome ne doit compter que si la verite tranche. Juger cela
sur la vacuite de la colonne d'assemblage est faux : l'outil de verite ecrit
« Unknown » quand il ne conclut pas, ce qui remplit la colonne sans rien dire.
Le couple restait alors au denominateur, comptait comme un echec, et affichait
son gene a 0,00 % dans le tableau publie. Le tri se fait sur le verdict.
"""
import os
import sys
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bench", "orthogonale"))

import mesurer  # noqa: E402


def ligne(**kw):
    d = {"perimetre": "oui", "echantillon": "E1", "gene": "G1",
         "assemblage": "*1/*1", "module": "*1/*1",
         "module_vs_assemblage": "concordant", "note": ""}
    d.update(kw)
    return d


class TestRetenus(unittest.TestCase):

    def test_verite_conclusive_retenue(self):
        gardes, ecartes = mesurer.retenus([ligne()])
        self.assertEqual(len(gardes), 1)
        self.assertFalse(ecartes)

    def test_verite_unknown_ecartee_malgre_colonne_remplie(self):
        """Le cas qui faisait afficher un gene a 0,00 %."""
        l = ligne(assemblage="Unknown", module_vs_assemblage="pas de verite",
                  note="verite non conclusive")
        gardes, ecartes = mesurer.retenus([l])
        self.assertEqual(gardes, [])
        self.assertEqual(sum(ecartes.values()), 1)

    def test_un_unknown_ne_fait_pas_chuter_le_taux(self):
        """Avec le defaut, 1 concordant + 1 Unknown donnaient 1/2."""
        lignes = [ligne(), ligne(echantillon="E2", assemblage="Unknown",
                                 module_vs_assemblage="pas de verite")]
        gardes, _ = mesurer.retenus(lignes)
        self.assertEqual(len(gardes), 1)

    def test_colonne_vide_toujours_ecartee(self):
        gardes, ecartes = mesurer.retenus([ligne(assemblage="")])
        self.assertEqual(gardes, [])
        self.assertEqual(sum(ecartes.values()), 1)

    def test_hors_perimetre_ni_garde_ni_compte(self):
        gardes, ecartes = mesurer.retenus([ligne(perimetre="non")])
        self.assertEqual(gardes, [])
        self.assertFalse(ecartes)

    def test_module_sans_compte_rendu_ecarte(self):
        gardes, ecartes = mesurer.retenus([ligne(module="sans compte rendu")])
        self.assertEqual(gardes, [])
        self.assertEqual(ecartes["genome sans compte rendu"], 1)

    def test_locus_douteux_ecarte(self):
        gardes, ecartes = mesurer.retenus([ligne()], douteux={("E1", "G1")})
        self.assertEqual(gardes, [])
        self.assertEqual(ecartes["locus a chaine douteuse"], 1)


if __name__ == "__main__":
    unittest.main()
