# -*- coding: utf-8 -*-
"""Tests du rejeu des seuils.

L'outil rejoue la regle de l'etage 2b a d'autres seuils de GQ et de profondeur,
pour mesurer ce que les seuils par defaut coutent. S'il ne reproduisait pas
exactement la regle, la mesure de sensibilite porterait sur autre chose que le
module.
"""
import os
import sys
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))

import sensibilite_seuils as s  # noqa: E402


def pos(genes, gq, prof, statut, p="chr1:100"):
    return (genes, gq, prof, statut, p)


class Rejeu(unittest.TestCase):
    def test_une_position_mesuree_qui_tient_les_deux_seuils(self):
        att, ret = s.rejoue([pos(["TPMT"], "99", "30", "mesuree")], 20, 10)
        self.assertEqual((att["TPMT"], ret["TPMT"]), (1, 1))

    def test_profondeur_sous_le_seuil(self):
        att, ret = s.rejoue([pos(["TPMT"], "99", "4", "mesuree")], 20, 10)
        self.assertEqual((att["TPMT"], ret["TPMT"]), (1, 0))

    def test_qualite_sous_le_seuil(self):
        att, ret = s.rejoue([pos(["TPMT"], "7", "30", "mesuree")], 20, 10)
        self.assertEqual(ret["TPMT"], 0)

    def test_une_qualite_absente_ecarte_la_position(self):
        # Un seuil qu'on ne peut pas appliquer n'est pas un seuil satisfait :
        # c'est la regle de l'etage 2b, et le rejeu doit la suivre.
        for gq in ("", ".", None, "abc"):
            att, ret = s.rejoue([pos(["TPMT"], gq, "30", "mesuree")], 20, 10)
            self.assertEqual(ret["TPMT"], 0, repr(gq))

    def test_une_profondeur_illisible_ecarte_la_position(self):
        att, ret = s.rejoue([pos(["TPMT"], "99", ".", "mesuree")], 20, 10)
        self.assertEqual((att["TPMT"], ret["TPMT"]), (1, 0))

    def test_une_reference_lue_ne_depend_pas_du_gq(self):
        # Une position de reference n'a pas de genotype appele : exiger un GQ
        # ferait perdre les positions conformes a la reference, c'est-a-dire la
        # grande majorite du perimetre.
        att, ret = s.rejoue([pos(["TPMT"], "", "30", "reference, lue")], 20, 10)
        self.assertEqual(ret["TPMT"], 1)

    def test_une_reference_lue_depend_de_la_profondeur(self):
        att, ret = s.rejoue([pos(["TPMT"], "", "4", "reference, lue")], 20, 10)
        self.assertEqual(ret["TPMT"], 0)

    def test_qualite_insuffisante_se_rejoue_sur_le_gq(self):
        # Le statut du fichier est le resultat du filtrage aux seuils par
        # defaut : juger sur le statut et non sur le GQ rendrait le rejeu
        # insensible aux seuils, donc sans objet.
        ligne = [pos(["TPMT"], "15", "30", "qualite insuffisante")]
        self.assertEqual(s.rejoue(ligne, 20, 10)[1]["TPMT"], 0)
        self.assertEqual(s.rejoue(ligne, 10, 10)[1]["TPMT"], 1)

    def test_un_statut_hors_liste_n_est_jamais_retenu(self):
        for statut in ("profondeur insuffisante", "non sequencee",
                       "couverture seule", ""):
            att, ret = s.rejoue([pos(["TPMT"], "99", "30", statut)], 20, 10)
            self.assertEqual(ret["TPMT"], 0, statut)

    def test_une_position_partagee_compte_pour_chaque_gene(self):
        att, ret = s.rejoue(
            [pos(["CYP2C9", "CYP2C19"], "99", "30", "mesuree")], 20, 10)
        self.assertEqual((att["CYP2C9"], att["CYP2C19"]), (1, 1))
        self.assertEqual((ret["CYP2C9"], ret["CYP2C19"]), (1, 1))

    def test_une_position_sans_gene_ne_compte_pour_personne(self):
        att, ret = s.rejoue([pos([], "99", "30", "mesuree")], 20, 10)
        self.assertEqual(sum(att.values()), 0)
        self.assertEqual(sum(ret.values()), 0)

    def test_les_attendues_ne_dependent_pas_des_seuils(self):
        # Le denominateur doit rester le meme d'un seuil a l'autre, sinon les
        # taux comparent des choses differentes.
        lignes = [pos(["TPMT"], "99", "30", "mesuree"),
                  pos(["TPMT"], "5", "3", "mesuree")]
        a1, _ = s.rejoue(lignes, 20, 10)
        a2, _ = s.rejoue(lignes, 50, 30)
        self.assertEqual(a1["TPMT"], a2["TPMT"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
