# -*- coding: utf-8 -*-
"""Tests de la trace d'execution.

La provenance porte le champ qu'une plateforme surveille, `reussite_complete`.
Il ne vaut vrai que si les neuf etages requis ont abouti : le typage
complementaire, optionnel, n'y entre pas, et un etat tolere ne doit pas se
confondre avec un echec.

Le calcul de ce champ est verifie en executant le module, et non en relisant
son source : un test qui cherche un mot dans le code survit a presque toutes
les fautes de logique.
"""
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))

import provenance as p  # noqa: E402

REQUIS = ["recevabilite", "filtre", "tranche", "qc", "cyp2d6", "hla",
          "appels", "pharmcat", "rendu"]


class Empreinte(unittest.TestCase):
    def fichier(self, octets):
        f = tempfile.mktemp()
        with open(f, "wb") as fh:
            fh.write(octets)
        return f

    def test_empreinte_complete(self):
        f = self.fichier(b"abc")
        try:
            self.assertEqual(p.sha256(f), hashlib.sha256(b"abc").hexdigest())
        finally:
            os.unlink(f)

    def test_empreinte_partielle_se_declare_comme_telle(self):
        # Un CRAM de cinquante gigaoctets ne se hache pas entierement : la
        # limite doit apparaitre dans la valeur, sinon deux empreintes
        # incomparables se ressemblent.
        f = self.fichier(b"x" * (1 << 21))
        try:
            v = p.sha256(f, limite=1 << 20)
        finally:
            os.unlink(f)
        self.assertIn("+partiel:", v)

    def test_fichier_absent_rend_une_chaine_vide(self):
        # Pas d'exception : la provenance s'ecrit toujours, meme sur un genome
        # refuse a l'entree.
        self.assertEqual(p.sha256(None), "")
        self.assertEqual(p.sha256(os.path.join(tempfile.gettempdir(),
                                               "absent_xyz_prov")), "")


class VersionsExterieures(unittest.TestCase):
    def test_un_outil_absent_ne_leve_pas(self):
        # Les versions sont relevees par sous-processus ; l'absence de l'outil
        # ne doit pas interrompre l'ecriture de la trace.
        self.assertEqual(p.version_pypgx(""), {})
        d = p.version_pypgx("commande_qui_n_existe_pas_xyz")
        self.assertEqual(d["commande"], "commande_qui_n_existe_pas_xyz")
        self.assertEqual(d["version"], "")

    def test_une_image_introuvable_rend_son_etiquette(self):
        nom = "depot/image_inexistante_xyz:1.0"
        self.assertEqual(p.version_image(nom), nom)

    def test_un_depot_git_absent_ne_leve_pas(self):
        self.assertIsInstance(p.version_git(tempfile.gettempdir()), (str, dict))


class ReussiteComplete(unittest.TestCase):
    """Le champ que surveille une plateforme, mesure en executant le module."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.d, "travail"))
        self.entrees = {}
        for nom in ("cram", "vcf"):
            c = os.path.join(self.d, nom + ".bin")
            with open(c, "wb") as fh:
                fh.write(b"x")
            self.entrees[nom] = c

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def trace(self, etats):
        chemin = os.path.join(self.d, "travail", "etats.tsv")
        with io.open(chemin, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("etage\tetat\tsecondes\tdetail\n")
            for etage, etat in etats:
                fh.write("%s\t%s\t1\tessai\n" % (etage, etat))
        r = subprocess.run(
            [sys.executable, os.path.join(RACINE, "bin", "provenance.py"),
             "--sortie", self.d, "--echantillon", "ECH",
             "--cram", self.entrees["cram"], "--vcf", self.entrees["vcf"],
             "--ressources", os.path.join(RACINE, "ressources"),
             "--gq", "20", "--profondeur", "10"],
            capture_output=True)
        p_json = os.path.join(self.d, "sortie", "provenance.json")
        if not os.path.exists(p_json):
            self.fail("provenance.json non ecrit : %s"
                      % r.stderr.decode("utf-8", "replace")[-400:])
        with io.open(p_json, encoding="utf-8") as fh:
            return json.load(fh)

    def test_les_neuf_etages_ok_donnent_une_reussite_complete(self):
        d = self.trace([(e, "OK") for e in REQUIS])
        self.assertTrue(d["reussite_complete"])

    def test_un_etage_en_echec_la_fait_tomber(self):
        for perdu in REQUIS:
            etats = [(e, "ECHEC" if e == perdu else "OK") for e in REQUIS]
            d = self.trace(etats)
            self.assertFalse(d["reussite_complete"], perdu)

    def test_un_etage_manquant_la_fait_tomber(self):
        # Un etage qui n'a pas inscrit son etat n'a pas abouti : l'absence ne
        # doit pas se lire comme un succes.
        for absent in REQUIS:
            etats = [(e, "OK") for e in REQUIS if e != absent]
            d = self.trace(etats)
            self.assertFalse(d["reussite_complete"], absent)

    def test_sans_resultat_n_est_pas_un_echec(self):
        # Cyrius qui ne tranche pas entre deux diplotypes n'est pas une faute
        # du module : le genome reste exploitable, le gene est signale.
        etats = [(e, "SANS_RESULTAT" if e == "cyp2d6" else "OK") for e in REQUIS]
        d = self.trace(etats)
        self.assertTrue(d["reussite_complete"])

    def test_l_etage_optionnel_ignore_ne_compte_pas(self):
        # Sans PyPGx installe, l'etage sort IGNORE ; le genome doit rester en
        # reussite complete, sinon le module serait declare en echec pour une
        # dependance qu'il annonce facultative.
        d = self.trace([(e, "OK") for e in REQUIS] + [("pypgx", "IGNORE")])
        self.assertTrue(d["reussite_complete"])
        self.assertIn("pypgx", d.get("etages_ignores") or [])

    def test_l_etage_optionnel_en_echec_fait_tomber_le_genome(self):
        # Non arme, il ne pese pas ; arme et en echec, il compte. On arme un
        # etage pour qu'il aboutisse, et un typage complementaire manquant
        # retire du compte rendu un gene de classe 1.
        d = self.trace([(e, "OK") for e in REQUIS] + [("pypgx", "ECHEC")])
        self.assertFalse(d["reussite_complete"])

    def test_la_trace_porte_la_version_du_module(self):
        # Le module tracait la version de chacune de ses dependances et pas la
        # sienne. C'est la premiere chose qu'un auditeur demande d'un compte
        # rendu : quelle version du logiciel l'a produit.
        d = self.trace([(e, "OK") for e in REQUIS])
        self.assertIn("module", d)
        self.assertTrue(d["module"].get("version"),
                        "la version du module est vide")

    def test_la_trace_porte_les_seuils_et_les_empreintes(self):
        d = self.trace([(e, "OK") for e in REQUIS])
        texte = json.dumps(d)
        self.assertIn("20", texte)          # GQ
        self.assertIn("10", texte)          # profondeur
        self.assertTrue(any("sha256" in k or "empreinte" in k
                            for k in json.dumps(d).split('"')), d.keys())


if __name__ == "__main__":
    unittest.main(verbosity=2)
