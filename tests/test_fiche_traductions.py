# -*- coding: utf-8 -*-
"""Tests de la fiche de relecture des traductions.

La fiche sert a faire valider les traductions par un pharmacologue : si elle
annonce comme non traduit un texte que le compte rendu rend en francais, la
relecture part sur une fausse piste. C'est arrive deux fois a l'ecriture de
l'outil — une fois en comptant des medicaments que le compte rendu ne restitue
pas, une fois en comparant des chaines brutes la ou la table est indexee sur un
texte normalise. Les deux cas sont tenus ici.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "outils"))

import fiche_traductions as f  # noqa: E402

SOURCE = "Avoid abacavir"
FR = "Ne pas prescrire l'abacavir"
TABLE = {"version": "1.0", "etabli_le": "2026-01-01", "statut": "proposition",
         "medicaments_restitues": {"abacavir": "abacavir"},
         "traductions": {SOURCE: {"fr": FR, "type": "contre-indication"}}}


def rapport(medicaments):
    """medicaments : [(nom, texte, classification, genes)]"""
    d = {}
    for nom, texte, classe, genes in medicaments:
        d.setdefault(nom, {"guidelines": [{"annotations": []}]})
        d[nom]["guidelines"][0]["annotations"].append(
            {"classification": classe, "drugRecommendation": texte,
             "lookupKey": [{g: "x" for g in genes}]})
    return {"drugs": {"CPIC Guideline Annotation": d}}


class Usages(unittest.TestCase):
    def setUp(self):
        self.racine = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.racine, ignore_errors=True)

    def ecris(self, ech, rep):
        d = os.path.join(self.racine, ech, "sortie")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, ech + ".report.json"), "w",
                  encoding="utf-8", newline="\n") as fh:
            json.dump(rep, fh)

    def test_compte_les_echantillons_et_les_genes(self):
        self.ecris("E1", rapport([("abacavir", SOURCE, "Strong", ["HLA-B"])]))
        self.ecris("E2", rapport([("abacavir", SOURCE, "Strong", ["HLA-B"])]))
        u = f.usages(self.racine, {"abacavir"})
        self.assertEqual(len(u[SOURCE]["echantillons"]), 2)
        self.assertEqual(u[SOURCE]["genes"], {"HLA-B"})

    def test_un_medicament_non_restitue_ne_compte_pas(self):
        # Il existe dans le rapport de l'interpreteur et n'atteint jamais le
        # document : le compter ferait annoncer un texte non traduit qui n'est
        # pas rendu.
        self.ecris("E1", rapport([("lovastatine", "Other text", "Strong", ["SLCO1B1"])]))
        self.assertEqual(dict(f.usages(self.racine, {"abacavir"})), {})

    def test_une_classification_faible_ne_compte_pas(self):
        self.ecris("E1", rapport([("abacavir", SOURCE, "Optional", ["HLA-B"])]))
        self.assertEqual(dict(f.usages(self.racine, {"abacavir"})), {})

    def test_le_texte_est_normalise_comme_dans_le_module(self):
        # Entites HTML et espaces multiples : la table est indexee sur le texte
        # normalise, et c'est ainsi que le module la consulte.
        self.ecris("E1", rapport([("abacavir", "Limit dose to &lt;20mg", "Strong", ["X"])]))
        u = f.usages(self.racine, {"abacavir"})
        self.assertIn("Limit dose to <20mg", u)

    def test_sans_restriction_tout_compte(self):
        self.ecris("E1", rapport([("lovastatine", "Other", "Strong", ["X"])]))
        self.assertEqual(len(f.usages(self.racine, None)), 1)


class Fiche(unittest.TestCase):
    def test_classe_par_usage_decroissant(self):
        table = {"traductions": {"A": {"fr": "a", "type": ""},
                                 "B": {"fr": "b", "type": ""}}}
        u = {"B": {"echantillons": {1, 2}, "medicaments": set(), "genes": set(),
                   "classes": set()},
             "A": {"echantillons": {1}, "medicaments": set(), "genes": set(),
                   "classes": set()}}
        lignes, _, _ = f.fiche(table, u)
        self.assertEqual([x["source"] for x in lignes], ["B", "A"])

    def test_un_texte_rendu_et_non_traduit_est_signale(self):
        table = {"traductions": {"A": {"fr": "a", "type": ""}}}
        u = {"Z": {"echantillons": {1}, "medicaments": set(), "genes": set(),
                   "classes": set()}}
        _, absents, _ = f.fiche(table, u)
        self.assertEqual(absents, ["Z"])

    def test_un_texte_normalise_n_est_pas_signale_comme_absent(self):
        # Le defaut qui annoncait quatorze textes non traduits : la cle de la
        # table portait des entites HTML, le rapport aussi, et la comparaison
        # brute ne les rapprochait pas.
        table = {"traductions": {"Limit dose to &lt;20mg": {"fr": "a", "type": ""}}}
        u = {"Limit dose to <20mg": {"echantillons": {1}, "medicaments": set(),
                                     "genes": set(), "classes": set()}}
        _, absents, _ = f.fiche(table, u)
        self.assertEqual(absents, [])

    def test_une_traduction_inemployee_est_signalee(self):
        table = {"traductions": {"A": {"fr": "a", "type": ""}}}
        _, _, inutilisees = f.fiche(table, {})
        self.assertEqual(inutilisees, ["A"])

    def test_la_fiche_s_ecrit_et_porte_la_version(self):
        d = tempfile.mkdtemp()
        try:
            p = os.path.join(d, "fiche.md")
            lignes, absents, inutilisees = f.fiche(
                TABLE, {SOURCE: {"echantillons": {1}, "medicaments": {"abacavir"},
                                 "genes": {"HLA-B"}, "classes": {"Strong"}}})
            f.ecris(p, TABLE, lignes, absents, inutilisees, True)
            with open(p, encoding="utf-8") as fh:
                t = fh.read()
            self.assertIn("1.0", t)
            self.assertIn("proposition", t)
            self.assertIn(FR, t)
            self.assertIn("Visa", t)
        finally:
            shutil.rmtree(d, ignore_errors=True)


class Tronque(unittest.TestCase):
    def test_reduit_les_espaces(self):
        self.assertEqual(f.tronque("a   b\n c"), "a b c")

    def test_coupe_et_marque_la_coupe(self):
        t = f.tronque("x" * 50, n=10)
        self.assertEqual(len(t), 10)
        self.assertTrue(t.endswith("…"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
