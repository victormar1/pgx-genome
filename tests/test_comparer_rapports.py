# -*- coding: utf-8 -*-
"""Tests de l'outil de comparaison de deux executions.

C'est lui qui dira si un changement du module a modifie un resultat rendu :
s'il manque un ecart, la revalidation conclura a tort que rien n'a bouge. Les
cas qui comptent sont donc ceux ou il doit voir une difference, et celui ou il
doit la qualifier plutot que la compter comme un desaccord.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "outils"))

import comparer_rapports as c  # noqa: E402


def rapport(diplotypes, medicaments=()):
    """Un rapport d'interpreteur reduit a ce que l'outil en lit."""
    genes = {}
    for g, (lab, phe) in diplotypes.items():
        genes[g] = {"recommendationDiplotypes": [
            {"label": lab, "phenotypes": [phe] if phe else []}]}
    annotations = {}
    for nom, texte in medicaments:
        annotations.setdefault(nom, {"guidelines": [{"annotations": []}]})
        annotations[nom]["guidelines"][0]["annotations"].append(
            {"classification": "Strong", "drugRecommendation": texte})
    return {"genes": genes,
            "drugs": {"CPIC Guideline Annotation": annotations}}


class Dossiers(unittest.TestCase):
    def setUp(self):
        self.racine = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.racine, ignore_errors=True)

    def ecris(self, cote, ech, rep):
        d = os.path.join(self.racine, cote, ech, "sortie")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, ech + ".report.json"), "w",
                  encoding="utf-8", newline="\n") as fh:
            json.dump(rep, fh)

    def compare(self, perimetre=("CYP2C19",)):
        return c.compare(os.path.join(self.racine, "avant"),
                         os.path.join(self.racine, "apres"), list(perimetre))

    def test_rien_ne_change(self):
        r = rapport({"CYP2C19": ("*1/*2", "Intermediate")},
                    [("clopidogrel", "autre antiagregant")])
        self.ecris("avant", "E1", r)
        self.ecris("apres", "E1", r)
        communs, genes, ecarts, cond = self.compare()
        self.assertEqual((len(communs), genes, ecarts, cond), (1, 1, [], []))

    def test_un_diplotype_change(self):
        self.ecris("avant", "E1", rapport({"CYP2C19": ("*1/*1", "Normal")}))
        self.ecris("apres", "E1", rapport({"CYP2C19": ("*1/*2", "Intermediate")}))
        _, _, ecarts, _ = self.compare()
        self.assertEqual(len(ecarts), 1)
        self.assertFalse(ecarts[0]["ambiguite_tranchee"])

    def test_une_ambiguite_tranchee_est_qualifiee(self):
        # Un appel unique pris dans l'ambiguite d'avant, meme phenotype : les
        # deux cotes disent la meme chose. L'ecart est compte, mais nomme.
        self.ecris("avant", "E1", rapport({"CYP2C19": ("*1/*2;*1/*3", "Intermediate")}))
        self.ecris("apres", "E1", rapport({"CYP2C19": ("*1/*2", "Intermediate")}))
        _, _, ecarts, _ = self.compare()
        self.assertEqual(len(ecarts), 1)
        self.assertTrue(ecarts[0]["ambiguite_tranchee"])

    def test_un_appel_hors_de_l_ambiguite_n_est_pas_qualifie(self):
        self.ecris("avant", "E1", rapport({"CYP2C19": ("*1/*2;*1/*3", "Intermediate")}))
        self.ecris("apres", "E1", rapport({"CYP2C19": ("*1/*17", "Rapid")}))
        _, _, ecarts, _ = self.compare()
        self.assertFalse(ecarts[0]["ambiguite_tranchee"])

    def test_un_phenotype_qui_change_n_est_jamais_une_ambiguite_tranchee(self):
        # Meme ensemble, mais la classe de fonction bouge : c'est un ecart
        # a consequence, et le qualifier de tranche le ferait passer pour benin.
        self.ecris("avant", "E1", rapport({"CYP2C19": ("*1/*2;*1/*3", "Intermediate")}))
        self.ecris("apres", "E1", rapport({"CYP2C19": ("*1/*2", "Poor")}))
        _, _, ecarts, _ = self.compare()
        self.assertFalse(ecarts[0]["ambiguite_tranchee"])

    def test_une_recommandation_apparait(self):
        self.ecris("avant", "E1", rapport({"CYP2C19": ("*1/*1", "Normal")}))
        self.ecris("apres", "E1", rapport({"CYP2C19": ("*1/*1", "Normal")},
                                          [("rosuvastatine", "dose souhaitee")]))
        _, _, ecarts, cond = self.compare()
        self.assertEqual(ecarts, [])
        self.assertEqual(cond[0]["apparus"], ["rosuvastatine"])

    def test_un_texte_de_recommandation_qui_change(self):
        self.ecris("avant", "E1", rapport({"CYP2C19": ("*1/*1", "Normal")},
                                          [("clopidogrel", "dose standard")]))
        self.ecris("apres", "E1", rapport({"CYP2C19": ("*1/*1", "Normal")},
                                          [("clopidogrel", "autre antiagregant")]))
        _, _, _, cond = self.compare()
        self.assertEqual(cond[0]["texte_change"], ["clopidogrel"])

    def test_un_echantillon_absent_d_un_cote_est_ignore(self):
        # Comparer un echantillon a rien produirait un faux ecart.
        self.ecris("avant", "E1", rapport({"CYP2C19": ("*1/*1", "Normal")}))
        self.ecris("apres", "E2", rapport({"CYP2C19": ("*1/*2", "Intermediate")}))
        communs, genes, ecarts, cond = self.compare()
        self.assertEqual((communs, genes, ecarts, cond), ([], 0, [], []))

    def test_un_gene_absent_des_deux_cotes_ne_compte_pas(self):
        r = rapport({"CYP2C19": ("*1/*1", "Normal")})
        self.ecris("avant", "E1", r)
        self.ecris("apres", "E1", r)
        _, genes, _, _ = self.compare(perimetre=("CYP2C19", "POR"))
        self.assertEqual(genes, 1)


class Inclus(unittest.TestCase):
    def test_appel_unique_dans_l_ambiguite(self):
        self.assertTrue(c.inclus("*1/*2;*1/*3", "*1/*2"))

    def test_appel_unique_hors_de_l_ambiguite(self):
        self.assertFalse(c.inclus("*1/*2;*1/*3", "*1/*4"))

    def test_deux_appels_uniques(self):
        self.assertFalse(c.inclus("*1/*2", "*1/*3"))

    def test_ambiguite_qui_s_elargit(self):
        self.assertFalse(c.inclus("*1/*2", "*1/*2;*1/*3"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
