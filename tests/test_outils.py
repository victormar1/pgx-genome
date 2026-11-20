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
import verifier_depot as D  # noqa: E402
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
        # octets declarerait trop longs des sujets qui ne le sont pas.
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

    def test_trailer_d_assistance_sous_sa_seule_forme(self):
        # Un nom de modele change a chaque version : une forme unique garde
        # l'historique homogene.
        debut = "doc: ajoute une note\n\nCorps.\n\n"
        self.assertEqual(V.verifie(debut + V.CO_AUTEUR_IA + "\n"), [])
        for variante in (
                "Co-authored-by: Claude Opus 5 <noreply@anthropic.com>",
                "Co-authored-by: Claude Fable 5.1 <noreply@anthropic.com>",
                "Co-authored-by: claude <noreply@anthropic.com>",
                "Co-authored-by: Claude <claude@anthropic.com>"):
            f = V.verifie(debut + variante + "\n")
            self.assertTrue(any("assistance" in x for x in f), variante)

    def test_un_co_auteur_humain_reste_libre(self):
        self.assertEqual(V.verifie(
            "doc: ajoute une note\n\nCorps.\n\n"
            "Co-authored-by: Marie Curie <marie@exemple.invalide>\n"), [])

    def test_le_domaine_ci_est_reconnu(self):
        self.assertEqual(V.verifie("ci: ajoute une note\n"), [])

    def test_un_co_auteur_humain_prenomme_claude_reste_libre(self):
        # Le trailer d'assistance se reconnait a son adresse, pas au prenom.
        self.assertEqual(V.verifie(
            "doc: ajoute une note\n\nCorps.\n\n"
            "Co-authored-by: Claude Dupont <claude@exemple.invalide>\n"), [])

    def test_un_chemin_windows_en_minuscule_ou_a_barre_oblique(self):
        for chemin in ("c:" + chr(92) + "Users", "D:" + "/" + "Victor"):
            f = V.verifie("doc: ajoute une note\n\nVoir " + chemin + " ici.\n")
            self.assertTrue(any("Windows" in x for x in f), chemin)

    def test_pas_de_ligne_vide_apres_le_sujet(self):
        f = V.verifie("doc: ajoute une note\nCorps colle.\n")
        self.assertTrue(any("ligne vide" in x for x in f))

    def test_les_commentaires_de_git_sont_ignores(self):
        msg = (self.BON + "# Please enter the commit message for your changes."
               + "\n# " + "a" * 200 + "\n")
        self.assertEqual(V.verifie(msg), [])

    def test_motifs_interdits(self):
        # Les chaines sont assemblees par morceaux a dessein : ecrites en clair,
        # elles feraient echouer le controle d'integrite du depot, qui cherche
        # ces motifs dans les fichiers suivis et ne sait pas distinguer un
        # exemple d'une fuite. Exempter ce fichier creerait un trou.
        cas = [
            "Voir " + "C:" + chr(92) + "Users" + chr(92) + "quelqu_un.",
            "Voir " + "/mnt" + "/d/" + "travail/ECH.",
            "Mesure sur " + "HG" + "00188.",
            "Mesure sur " + "NA" + "12878.",
            "Finalement corrige.",
            "Seuil a valider par le laboratoire.",
            "TODO : revoir ce point.",
        ]
        for corps in cas:
            f = V.verifie("doc: ajoute une note\n\n" + corps + "\n")
            self.assertTrue(f, corps)

    def test_un_message_vide_est_refuse(self):
        self.assertTrue(V.verifie("\n\n"))
        self.assertTrue(V.verifie("# que des commentaires\n"))


class MotifsDuDepot(unittest.TestCase):
    """Ce que le controle du depot cherche dans les fichiers suivis."""

    def trouve(self, texte):
        return [quoi for motif, quoi in D.INTERDITS if re.search(motif, texte)]

    def test_un_chemin_windows_sous_toutes_ses_ecritures(self):
        # Le motif exigeait deux antislashs, l'ecriture d'une chaine echappee :
        # un chemin ecrit normalement dans un document passait.
        barre = chr(92)
        for sep in (barre, barre * 2, "/"):
            for lecteur in ("C:", "c:", "D:"):
                chemin = lecteur + sep + "Users" + sep + "quelqu_un"
                self.assertIn("chemin Windows d'un poste", self.trouve("voir " + chemin), chemin)

    def test_un_texte_ordinaire_ne_declenche_rien(self):
        for texte in ("GT:GQ:DP", "https://example.org/Users/guide", "ratio 1:2/Users",
                      "les Users du service", "https:" + "//Users.example.org"):
            self.assertEqual(self.trouve(texte), [], texte)


class DocumentationAJour(unittest.TestCase):
    """Le code ne part pas sans sa documentation : ce que le controle du depot
    lit pour le dire."""

    SCRIPT = ('while [ $# -gt 0 ]; do\n  case "$1" in\n'
              '    --cram) CRAM="$2"; shift 2;;\n'
              '    --si-incomplet) SI=1; shift;;\n'
              '    -f|--fils) FILS="$2"; shift 2;;\n'
              '    -h|--help) exit 0;;\n'
              '    *) exit 2;;\n  esac\ndone\n'
              'outil --pas-une-option) fin\n')

    def test_les_options_d_un_script_se_lisent_dans_son_aiguillage(self):
        # Une option qui porte aussi une forme courte se lit sous sa forme longue.
        self.assertEqual(D.options_acceptees(self.SCRIPT),
                         ["--cram", "--fils", "--si-incomplet"])

    def test_l_aide_n_a_pas_a_etre_documentee(self):
        self.assertNotIn("--help", D.options_acceptees(self.SCRIPT))

    def test_le_journal_se_lit_a_sa_premiere_entree(self):
        journal = "# Journal\n\nLes versions suivent `M.m.c`.\n\n## 0.6.1\n\nTexte.\n\n## 0.6.0\n"
        self.assertEqual(D.version_du_journal(journal), "0.6.1")
        self.assertEqual(D.version_du_journal("# Journal\n\n## 0.6.1.1\n"), "0.6.1.1")

    def test_un_journal_sans_entree_ne_donne_aucune_version(self):
        self.assertEqual(D.version_du_journal("# Journal\n\n## Preambule\n"), "")

    def test_du_code_sans_document_est_une_faute(self):
        for code in ("bin/pgx_genome.sh", "ressources/perimetre_rnpgx.json",
                     "env/requirements.txt", "exemples/pgx_genome.nf"):
            self.assertTrue(D.sans_documentation([code, "tests/test_x.py"]), code)

    def test_du_code_avec_un_document_passe(self):
        for doc in ("doc/MODE_EMPLOI.md", "README.md", "CHANGELOG.md"):
            self.assertFalse(D.sans_documentation(["bin/pgx_genome.sh", doc]), doc)

    def test_sans_code_modifie_rien_n_est_exige(self):
        self.assertFalse(D.sans_documentation(["tests/test_x.py", "outils/x.py"]))
        self.assertFalse(D.sans_documentation([]))

    def test_un_fichier_au_nom_proche_n_est_pas_un_document(self):
        self.assertTrue(D.sans_documentation(["bin/a.py", "bench/README.md", "docs.txt"]))

    def test_le_depot_tient_sa_propre_regle(self):
        self.assertEqual(D.version_du_journal(D.lire("CHANGELOG.md")),
                         D.lire("VERSION").strip())
        mode = D.lire("doc", "MODE_EMPLOI.md")
        for script in ("pgx_genome.sh", "pgx_lot.sh"):
            for option in D.options_acceptees(D.lire("bin", script)):
                self.assertIn(option, mode, script)


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
