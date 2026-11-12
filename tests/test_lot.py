# -*- coding: utf-8 -*-
"""Le lot : ce que le manifeste laisse passer, et ce que dit le tableau de bord.

Un identifiant devient un nom de dossier, et le tableau de bord est ce qu'un
exploitant lit pour savoir si un lot de cent genomes a abouti : un echantillon
qui en disparait est un echec que personne ne voit.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))

import lot  # noqa: E402


class Manifeste(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="pgx_lot_")
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.cram = self.fichier("a.cram")
        self.vcf = self.fichier("a.vcf.gz")

    def fichier(self, nom):
        p = os.path.join(self.d, nom)
        with open(p, "w", encoding="utf-8") as fh:
            fh.write("x")
        return p

    def verifier(self, *lignes):
        m = os.path.join(self.d, "manifeste.tsv")
        with open(m, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(lignes) + "\n")
        return lot.verifier(m)

    def ligne(self, ech, cram=None, vcf=None):
        return "\t".join([ech, cram or self.cram, vcf or self.vcf])

    def test_une_ligne_recevable(self):
        ok, refus = self.verifier("# commentaire", "", self.ligne("ECH-1.a_b"))
        self.assertEqual(ok, [("ECH-1.a_b", self.cram, self.vcf)])
        self.assertEqual(refus, [])

    def test_moins_de_trois_colonnes(self):
        ok, refus = self.verifier("ECH\t" + self.cram)
        self.assertEqual(ok, [])
        self.assertEqual(refus, ["ligne 1 : moins de trois colonnes"])

    def test_un_identifiant_qui_sortirait_du_dossier_est_refuse(self):
        # L'identifiant sert de nom de dossier sous la sortie du lot.
        for ech in ("../ailleurs", "a/b", "avec espace", "", "x" * 65):
            ok, refus = self.verifier(self.ligne(ech))
            self.assertEqual(ok, [], ech)
            self.assertIn("identifiant inutilisable", refus[0], ech)

    def test_un_identifiant_en_double_ne_passe_qu_une_fois(self):
        # Deux lignes sous le meme nom ecriraient dans le meme dossier.
        ok, refus = self.verifier(self.ligne("ECH"), self.ligne("ECH"))
        self.assertEqual(len(ok), 1)
        self.assertEqual(refus, ["ligne 2 : identifiant en double : ECH"])

    def test_un_fichier_introuvable_est_refuse_et_nomme(self):
        ok, refus = self.verifier(self.ligne("A", cram="/absent/a.cram"),
                                  self.ligne("B", vcf="/absent/b.vcf.gz"))
        self.assertEqual(ok, [])
        self.assertEqual(refus, ["ligne 1 : alignement introuvable : /absent/a.cram",
                                 "ligne 2 : fichier de variants introuvable : /absent/b.vcf.gz"])

    def test_un_identifiant_refuse_ne_reserve_pas_son_nom(self):
        # La premiere ligne est refusee pour son fichier : la seconde, valide,
        # ne doit pas etre prise pour un doublon.
        ok, refus = self.verifier(self.ligne("ECH", cram="/absent/a.cram"),
                                  self.ligne("ECH"))
        self.assertEqual([x[0] for x in ok], ["ECH"])
        self.assertEqual(len(refus), 1)


class TableauDeBord(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="pgx_lot_")
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)

    def lot(self, **traces):
        """Un echantillon par argument ; None = aucune sortie."""
        with open(os.path.join(self.d, "manifeste_verifie.tsv"), "w",
                  encoding="utf-8", newline="\n") as fh:
            for ech in traces:
                fh.write("%s\tx.cram\tx.vcf.gz\n" % ech)
        for ech, etages in traces.items():
            if etages is None:
                continue
            s = os.path.join(self.d, ech, "sortie")
            os.makedirs(s)
            echec = any(v not in ("OK", "SANS_RESULTAT", "IGNORE") for v in etages.values())
            with open(os.path.join(s, "provenance.json"), "w", encoding="utf-8") as fh:
                json.dump({"reussite_complete": not echec, "date": "2026-01-01",
                           "etages": {e: {"etat": v, "secondes": 3}
                                      for e, v in etages.items()}}, fh)
        return lot.tableau(self.d)

    @staticmethod
    def tous(etat="OK", **autres):
        d = {e: etat for e in lot.ATTENDUS}
        d["pypgx"] = "IGNORE"
        d.update(autres)
        return d

    def test_un_echantillon_sans_sortie_garde_sa_ligne(self):
        lignes = self.lot(A=self.tous(), B=None)
        self.assertEqual([l["echantillon"] for l in lignes], ["A", "B"])
        self.assertFalse(lignes[1]["reussite_complete"])
        self.assertEqual(lignes[1]["hla"], "AUCUNE SORTIE")
        self.assertIn("  etage hla           en echec chez 1 : B", lot.resume(lignes))

    def test_un_etage_absent_de_la_trace_est_dit_absent(self):
        etages = self.tous()
        del etages["rendu"]
        lignes = self.lot(A=etages)
        self.assertEqual(lignes[0]["rendu"], "ABSENT")
        self.assertIn("  etage rendu         en echec chez 1 : A", lot.resume(lignes))

    def test_sans_resultat_n_est_pas_compte_en_echec(self):
        r = lot.resume(self.lot(A=self.tous(cyp2d6="SANS_RESULTAT")))
        self.assertEqual(r[0], "echantillons du manifeste : 1, reussite complete : 1")
        self.assertIn("  etage cyp2d6        sans resultat chez 1 : A", r)
        self.assertFalse([l for l in r if "en echec" in l])

    def test_l_etage_optionnel_est_compte_a_part(self):
        r = lot.resume(self.lot(A=self.tous(), B=self.tous(pypgx="OK")))
        self.assertEqual(r[-1], "  etage pypgx         (optionnel) : IGNORE chez 1, OK chez 1")

    def test_le_fichier_porte_une_colonne_par_etage_et_sa_duree(self):
        lignes = self.lot(A=self.tous())
        with open(lot.ecrire(self.d, lignes), encoding="utf-8") as fh:
            entete, ligne = [l.rstrip("\n").split("\t") for l in fh]
        self.assertEqual(entete[:3], ["echantillon", "reussite_complete", "date"])
        for e in lot.ATTENDUS + lot.OPTIONNELS:
            self.assertIn(e, entete)
            self.assertIn(e + "_s", entete)
        self.assertEqual(ligne[entete.index("qc_s")], "3")

    def test_la_liste_des_etages_est_celle_de_la_trace(self):
        import provenance
        self.assertIs(lot.ATTENDUS, provenance.ATTENDUS)


if __name__ == "__main__":
    unittest.main(verbosity=2)
