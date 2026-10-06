# -*- coding: utf-8 -*-
"""Tests des outils du depot : le validateur de messages de commit, et le
generateur de cas construits.

Les garde-fous doivent etre tenus par des tests comme le reste : le validateur
a deja laisse passer un message qu'il devait refuser, et un cas construit dont
les positions ne correspondraient pas a la ressource exercerait un allele qui
n'existe pas dans le perimetre.
"""
import io
import os
import re
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "outils"))

import valider_message as V  # noqa: E402
import cas_construit as C  # noqa: E402


class ValidateurDeMessage(unittest.TestCase):
    BON = "doc: ajoute une note\n\nUn corps bref et conforme.\n"

    def test_un_message_conforme_passe(self):
        self.assertEqual(V.verifie(self.BON), [])

    def test_sujet_seul(self):
        self.assertEqual(V.verifie("doc: ajoute une note\n"), [])

    def test_sujet_sans_domaine(self):
        self.assertTrue(V.verifie("Ajoute une note\n\nCorps.\n"))

    def test_domaine_inconnu(self):
        f = V.verifie("divers: ajoute une note\n\nCorps.\n")
        self.assertTrue(any("inconnu" in x for x in f))

    def test_majuscule_apres_le_deux_points(self):
        f = V.verifie("doc: Ajoute une note\n\nCorps.\n")
        self.assertTrue(any("majuscule" in x for x in f))

    def test_point_final(self):
        f = V.verifie("doc: ajoute une note.\n\nCorps.\n")
        self.assertTrue(any("point" in x for x in f))

    def test_sujet_trop_long(self):
        f = V.verifie("doc: " + "a" * 80 + "\n\nCorps.\n")
        self.assertTrue(any("caracteres" in x for x in f))

    def test_la_longueur_se_compte_en_caracteres_et_non_en_octets(self):
        # Un accent pese deux octets en UTF-8. Un controle qui compte des
        # octets declarait trop longs des sujets qui ne le sont pas.
        sujet = "doc: " + "é" * 60       # 65 caracteres, 125 octets
        self.assertEqual(len(sujet), 65)
        self.assertEqual(V.verifie(sujet + "\n"), [])

    def test_corps_trop_large(self):
        f = V.verifie("doc: ajoute une note\n\n" + "mot " * 25 + "fin\n")
        self.assertTrue(any("corps" in x for x in f))

    def test_une_ligne_de_fin_ne_se_replie_pas(self):
        long_trailer = ("Reviewed-by: Quelqu'un avec un nom tres long "
                        "<adresse.tres.longue@exemple.invalide>")
        self.assertGreater(len(long_trailer), 72)
        self.assertEqual(V.verifie("doc: ajoute une note\n\nCorps.\n\n"
                                   + long_trailer + "\n"), [])

    def test_une_url_ne_se_replie_pas(self):
        url = "  https://exemple.invalide/" + "a" * 60
        self.assertEqual(V.verifie("doc: ajoute une note\n\n" + url + "\n"), [])

    def test_casse_d_une_ligne_de_fin(self):
        # Le projet Git ne met en majuscule que l'initiale du jeton, et la
        # forme fautive est celle qu'insere l'outillage.
        f = V.verifie("doc: ajoute une note\n\nCorps.\n\n"
                      "Co-Authored-By: Machin <x@y.z>\n")
        self.assertTrue(any("initiale" in x for x in f))
        self.assertEqual(V.verifie("doc: ajoute une note\n\nCorps.\n\n"
                                   "Co-authored-by: Machin <x@y.z>\n"), [])

    def test_pas_de_ligne_vide_apres_le_sujet(self):
        f = V.verifie("doc: ajoute une note\nCorps colle.\n")
        self.assertTrue(any("ligne vide" in x for x in f))

    def test_les_commentaires_de_git_sont_ignores(self):
        msg = (self.BON + "# Please enter the commit message for your changes."
               + "\n# " + "a" * 200 + "\n")
        self.assertEqual(V.verifie(msg), [])

    def test_motifs_interdits(self):
        cas = [
            "Voir C:\\Users\\quelqu_un\\travail.",
            "Voir /mnt/d/travail/ECH.",
            "Mesure sur HG00188.",
            "Mesure sur NA12878.",
            "Finalement corrige.",
            "Traductions a valider par un pharmacologue.",
            "TODO : revoir ce point.",
        ]
        for corps in cas:
            f = V.verifie("doc: ajoute une note\n\n" + corps + "\n")
            self.assertTrue(f, corps)

    def test_un_message_vide_est_refuse(self):
        self.assertTrue(V.verifie("\n\n"))
        self.assertTrue(V.verifie("# que des commentaires\n"))


class CasConstruit(unittest.TestCase):
    def setUp(self):
        with io.open(os.path.join(RACINE, "ressources",
                                  "rnpgx_complement.vcf"),
                     encoding="utf-8") as fh:
            self.complement = fh.read()

    def test_les_alleles_proposes_existent_dans_la_ressource(self):
        # Un cas construit sur une position absente du perimetre exercerait un
        # allele que le module ne mesure pas.
        for nom, (contig, pos, ref, alt, rsid) in C.ALLELES.items():
            motif = r"^%s\t%d\t%s\t%s\t" % (re.escape(contig), pos,
                                            re.escape(rsid), re.escape(ref))
            self.assertRegex(self.complement, re.compile(motif, re.M),
                             "%s introuvable dans le complement" % nom)

    def test_l_allele_alternatif_correspond(self):
        for nom, (contig, pos, ref, alt, rsid) in C.ALLELES.items():
            for ligne in self.complement.splitlines():
                c = ligne.split("\t")
                if len(c) > 4 and c[0] == contig and c[1] == str(pos):
                    self.assertIn(alt, c[4].split(","),
                                  "%s : allele alternatif absent" % nom)
                    break
            else:
                self.fail("%s : position absente" % nom)

    def test_les_trois_positions_sont_de_classe_un(self):
        # Ce sont les positions a consequence : celles pour lesquelles un cas
        # construit vaut la peine.
        for nom, (contig, pos, _, _, _) in C.ALLELES.items():
            for ligne in self.complement.splitlines():
                c = ligne.split("\t")
                if len(c) > 7 and c[0] == contig and c[1] == str(pos):
                    self.assertIn("PXCLASSE=1", c[7], nom)
                    break


if __name__ == "__main__":
    unittest.main(verbosity=2)
