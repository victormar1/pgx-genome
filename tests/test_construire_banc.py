# -*- coding: utf-8 -*-
"""La liste des prelevements du banc se deduit du lot, et le statut de chaque
ligne de sa trace. Un refus compte comme un refus, pas comme un echec, et un
echec d'etage n'est pas un prelevement analyse."""
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "outils"))

import construire_banc as b  # noqa: E402


class Statut(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def trace(self, nom, **champs):
        p = os.path.join(self.d, nom, "sortie")
        os.makedirs(p)
        with io.open(os.path.join(p, "provenance.json"), "w", encoding="utf-8") as fh:
            json.dump(champs, fh)
        return os.path.join(p, "provenance.json")

    def test_reussite_complete_est_analyse(self):
        p = self.trace("E1", reussite_complete=True, etages={})
        self.assertEqual(b.statut(p), "analysé")

    def test_recevabilite_en_echec_est_un_refus_avec_son_motif(self):
        p = self.trace("E2", reussite_complete=False,
                       etages={"recevabilite": {"etat": "ECHEC", "detail": "alignement tronque ; "}})
        self.assertEqual(b.statut(p), "refusé : alignement tronque")

    def test_echec_d_etage_nomme_l_etage(self):
        p = self.trace("E3", reussite_complete=False, etages_en_echec=["cyp2d6"],
                       etages={"recevabilite": {"etat": "OK"}, "cyp2d6": {"etat": "ECHEC"}})
        self.assertEqual(b.statut(p), "échec : cyp2d6")

    def test_une_ligne_par_dossier_du_lot_avec_ses_verites(self):
        self.trace("E1", reussite_complete=True, etages={})
        self.trace("E2", reussite_complete=True, etages={})
        os.makedirs(os.path.join(self.d, "sans_trace"))
        rows = b.lignes(self.d, getrm={"E1"}, hla={"E2"}, pop={"E1": "GBR"})
        self.assertEqual([r[0] for r in rows], ["E1", "E2"])
        self.assertEqual(rows[0][2:], ["GBR", "analysé", "oui", "non"])
        self.assertEqual(rows[1][2:], ["", "analysé", "non", "oui"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
