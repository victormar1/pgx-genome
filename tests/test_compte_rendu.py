# -*- coding: utf-8 -*-
"""Le compte rendu, execute de bout en bout.

Le regroupement des medicaments est decide et teste dans cr_logique. Ce test-ci
verifie que le document le suit : une statine seule y garde son nom, sans quoi
la consigne sort anonyme et la conclusion dit qu'il n'y a aucune recommandation.
"""
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(RACINE, "bin", "compte_rendu.py")


@unittest.skipUnless(importlib.util.find_spec("reportlab"), "reportlab absent")
class Regroupement(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="pgx_cr_")
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)

    def rendre(self, medicaments):
        consigne = lambda m: {
            "classification": "Strong", "lookupKey": [{"ABCG2": "Poor Function"}],
            "drugRecommendation": "Prescribe <=20mg of %s." % m,
            "dosingInformation": True}
        rapport = {
            "pharmcatVersion": "3.4.0",
            "genes": {"ABCG2": {"recommendationDiplotypes": [
                {"label": "rs2231142 variant (T)/rs2231142 variant (T)",
                 "phenotypes": ["Poor Function"]}]}},
            "drugs": {"CPIC Guideline Annotation": {
                m: {"citations": [{"pmid": "35152405", "year": 2022}],
                    "guidelines": [{"annotations": [consigne(m)]}]}
                for m in medicaments}}}
        perimetre = {
            "genes": {"ABCG2": {"statut": "complet", "positions_attendues": 1,
                                "positions_perdues": 0}},
            "perimetre_clinique": ["ABCG2"], "couverture_mediane_retenue": 31}
        chemins = {}
        for nom, contenu in (("rapport", rapport), ("perimetre", perimetre)):
            chemins[nom] = os.path.join(self.d, nom + ".json")
            with open(chemins[nom], "w", encoding="utf-8") as fh:
                json.dump(contenu, fh)
        r = subprocess.run(
            [sys.executable, "-B", SCRIPT, "--rapport", chemins["rapport"],
             "--perimetre", chemins["perimetre"], "--echantillon", "ECH",
             "--sortie", os.path.join(self.d, "cr.pdf")],
            capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return r.stdout

    def test_une_statine_seule_garde_sa_ligne_et_son_nom(self):
        self.assertIn("1 medicaments restitues", self.rendre(["rosuvastatin"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
