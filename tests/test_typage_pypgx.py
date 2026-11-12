# -*- coding: utf-8 -*-
"""Tests de l'etage 4b. Chaque test encode un garde-fou ou une regression deja
constatee, et porte en commentaire ce qu'il empeche de revenir.

Aucune dependance externe : les fonctions testees ne font pas appel au typeur,
qui est invoque en sous-processus ailleurs dans le module.
"""
import gzip
import io
import json
import os
import sys
import tempfile
import unittest
import zipfile
from unittest import mock

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))

import typage_pypgx as t  # noqa: E402


def ligne(h1, h2, genotype="", donnees=""):
    return {"Haplotype1": h1, "Haplotype2": h2,
            "Genotype": genotype, "VariantData": donnees}


BCHE = {
    "positions_panel": ["chr3:165773492", "chr3:165830741"],
    "alleles_positions": {"rs1803274": ["chr3:165773492"],
                          "rs1799807": ["chr3:165830741"]},
    "allele_reference": "Reference",
}
POR = {
    "positions_panel": ["chr7:75985688"],
    "alleles_positions": {"*28": ["chr7:75985688"],
                          "*10": ["chr7:75985688", "chr7:75990000"],
                          "*5": ["chr7:75900000"]},
    "allele_reference": "*1",
}


class Haplotypes(unittest.TestCase):
    def test_lit_tous_les_alleles_de_chaque_haplotype(self):
        # La colonne Genotype de PyPGx n'en garde qu'un par haplotype : sur BCHE
        # elle taisait l'allele atypique rs1799807 chez quatre porteurs.
        self.assertEqual(t.haplotypes(ligne("rs1803274;rs1799807;", "Reference;")),
                         [["rs1803274", "rs1799807"], ["Reference"]])

    def test_normalise_les_espaces_autour_des_noms(self):
        # Un nom entoure d'espaces ne correspondrait a aucune cle du catalogue
        # et sortirait sous une forme differente d'un patient a l'autre.
        self.assertEqual(t.haplotypes(ligne(" *28 ; ", "*1;")),
                         [["*28"], ["*1"]])

    def test_haplotype_vide(self):
        self.assertEqual(t.haplotypes(ligne("", "")), [[], []])


class GenotypeComplet(unittest.TestCase):
    def test_garde_les_deux_alleles_d_un_haplotype(self):
        self.assertEqual(
            t.genotype_complet(ligne("rs1803274;rs1799807;", "Reference;")),
            "Reference/rs1799807+rs1803274")

    def test_ordre_stable_quel_que_soit_le_rang_des_colonnes(self):
        # PyPGx n'attribue pas Haplotype1 au meme chromosome d'un genome a
        # l'autre : sans tri, le meme porteur sortirait « *1/*28 » ici et
        # « *28/*1 » la, et rien ne serait comparable entre patients.
        a = t.genotype_complet(ligne("*1;", "*28;"))
        b = t.genotype_complet(ligne("*28;", "*1;"))
        self.assertEqual(a, b)

    def test_rend_vide_si_aucun_haplotype(self):
        self.assertEqual(t.genotype_complet(ligne("", "")), "")


class GenotypePanel(unittest.TestCase):
    def test_rend_toujours_trois_valeurs(self):
        # L'appelant deballe trois valeurs : un chemin qui n'en rendrait que
        # deux leverait une ValueError sur tout genome valide.
        for cas in (ligne("rs1803274;", "Reference;"),
                    ligne("inconnu;", "Reference;"),
                    ligne("", "")):
            self.assertEqual(len(t.genotype_panel(cas, BCHE)), 3)
        self.assertEqual(len(t.genotype_panel(ligne("*1;", "*1;"), {})), 3)

    def test_retient_les_alleles_du_panel(self):
        g, ecartes, souci = t.genotype_panel(
            ligne("rs1803274;rs1799807;", "Reference;"), BCHE)
        self.assertEqual(g, "Reference/rs1799807+rs1803274")
        self.assertEqual(ecartes, 0)
        self.assertIsNone(souci)

    def test_ecarte_les_alleles_hors_panel_et_les_compte(self):
        g, ecartes, souci = t.genotype_panel(ligne("*5;", "*28;"), POR)
        self.assertEqual(g, "*1/*28")
        self.assertEqual(ecartes, 1)
        self.assertIsNone(souci)

    def test_nomme_la_reference_comme_le_gene_la_nomme(self):
        # Rendre « Reference/Reference » la ou le typeur dit « *1/*1 » ferait
        # croire a un autre resultat.
        g, _, _ = t.genotype_panel(ligne("*1;", "*1;"), POR)
        self.assertEqual(g, "*1/*1")

    def test_refuse_un_catalogue_incomplet(self):
        # Sans refus, le rendu retomberait sur le genotype complet de PyPGx
        # comme s'il etait restreint au perimetre clinique.
        for manque in ("positions_panel", "alleles_positions", "allele_reference"):
            info = dict(BCHE)
            del info[manque]
            g, _, souci = t.genotype_panel(ligne("rs1803274;", "Reference;"), info)
            self.assertIsNone(g, manque)
            self.assertTrue(souci)

    def test_refuse_un_allele_inconnu_du_catalogue(self):
        # Le compter hors panel reviendrait a decider sans savoir s'il touche
        # une position du panel ; c'est le cas d'une montee de version du typeur.
        g, _, souci = t.genotype_panel(ligne("rs9999999;", "Reference;"), BCHE)
        self.assertIsNone(g)
        self.assertIn("rs9999999", souci)


class AlleleUnique(unittest.TestCase):
    def test_homozygote(self):
        self.assertEqual(t.allele_unique(ligne("Reference;", "Reference;")),
                         ("Reference", False, None))

    def test_heteroplasmie_garde_l_allele_non_reference(self):
        a, hetero, souci = t.allele_unique(ligne("m.1555A>G;", "Reference;"))
        self.assertEqual(a, "m.1555A>G")
        self.assertTrue(hetero)
        self.assertIsNone(souci)

    def test_refuse_deux_alleles_non_reference_differents(self):
        # Garder le premier des deux reproduirait la perte silencieuse sur le
        # seul gene de l'etage qui porte un phenotype et une recommandation.
        a, _, souci = t.allele_unique(ligne("m.1555A>G;", "m.1494C>T;"))
        self.assertIsNone(a)
        self.assertTrue(souci)

    def test_refuse_un_haplotype_a_plusieurs_alleles(self):
        a, _, souci = t.allele_unique(ligne("m.1555A>G;m.1494C>T;", "Reference;"))
        self.assertIsNone(a)
        self.assertTrue(souci)


def enregistrement(pos, gt="0/1", gq=99, ad="15,15", filtre="PASS", contig="chr3"):
    """Une ligne de VCF a dix colonnes, comme les lit le controle qualite."""
    cles, vals = ["GT"], [gt]
    for nom, v in (("GQ", gq), ("AD", ad)):
        if v is not None:
            cles.append(nom)
            vals.append(str(v))
    return [contig, str(pos), ".", "C", "T", ".", filtre, ".", ":".join(cles), ":".join(vals)]


class Recevabilite(unittest.TestCase):
    INFO = {"contig": "chr3", "positions": ["chr3:100:C:T", "chr3:200:T:C"]}
    PROF = {("chr3", 100): 30, ("chr3", 200): 25}

    def recevable(self, *enregistrements, **kw):
        lignes = {}
        for c in enregistrements or (enregistrement(100),):
            lignes.setdefault((c[0], int(c[1])), []).append(c)
        a = {"info": self.INFO, "prof": self.PROF, "vus": {"chr3"},
             "prof_min": 10, "gq_min": 20, "ignorer_filtre": False}
        a.update(kw)
        return t.recevabilite(a["info"], a["prof"], a["vus"], lignes,
                              a["prof_min"], a["gq_min"], a["ignorer_filtre"])

    def test_cas_nominal(self):
        ok, lues, total, motif = self.recevable()
        self.assertTrue(ok)
        self.assertEqual((lues, total, motif), (2, 2, ""))

    def test_refuse_un_contig_absent_du_fichier_de_variants(self):
        # Sans ce refus, un genome dont le chrM n'a jamais ete appele
        # ressortirait « risque normal de surdite sous aminoside ».
        ok, _, _, motif = self.recevable(vus={"chr1"})
        self.assertFalse(ok)
        self.assertIn("chr3", motif)

    def test_refuse_une_liste_de_contigs_vide(self):
        ok, _, _, motif = self.recevable(vus=set())
        self.assertFalse(ok)
        self.assertIn("aucun contig", motif)

    def test_refuse_une_profondeur_non_mesuree(self):
        ok, _, _, motif = self.recevable(prof={})
        self.assertFalse(ok)
        self.assertIn("profondeur", motif)

    def test_refuse_une_position_sous_le_seuil(self):
        ok, lues, total, motif = self.recevable(
            prof={("chr3", 100): 30, ("chr3", 200): 4})
        self.assertFalse(ok)
        self.assertEqual((lues, total), (1, 2))
        self.assertIn("1 position", motif)

    def test_refuse_un_genotype_sans_gq(self):
        ok, _, _, motif = self.recevable(enregistrement(100, gq=None))
        self.assertFalse(ok)
        self.assertIn("chr3:100 (qualite non mesurable)", motif)

    def test_refuse_un_gq_sous_le_seuil_et_passe_au_seuil(self):
        self.assertFalse(self.recevable(enregistrement(100, gq=19))[0])
        self.assertTrue(self.recevable(enregistrement(100, gq=20))[0])

    def test_refuse_un_filtre_autre_que_pass(self):
        ok, _, _, motif = self.recevable(enregistrement(100, filtre="LowQual"))
        self.assertFalse(ok)
        self.assertIn("site filtre", motif)

    def test_ignorer_filtre_laisse_passer_le_filtre_mais_pas_le_gq(self):
        ok, _, _, _ = self.recevable(enregistrement(100, filtre="LowQual"),
                                     ignorer_filtre=True)
        self.assertTrue(ok)
        ok, _, _, motif = self.recevable(enregistrement(100, filtre="LowQual", gq=3),
                                         ignorer_filtre=True)
        self.assertFalse(ok)
        self.assertIn("qualite insuffisante", motif)

    # Les trois regles ci-dessous sont celles du controle qualite que le typage
    # n'appliquait pas : il nommait un allele sur un genotype que le reste du
    # module avait ecarte.

    def test_refuse_un_appel_sur_trop_peu_de_lectures(self):
        ok, _, _, motif = self.recevable(enregistrement(100, ad="5,4"))
        self.assertFalse(ok)
        self.assertIn("profondeur d appel insuffisante", motif)
        self.assertTrue(self.recevable(enregistrement(100, ad="5,5"))[0])

    def test_refuse_un_heterozygote_desequilibre(self):
        ok, _, _, motif = self.recevable(enregistrement(100, ad="76,24"))
        self.assertFalse(ok)
        self.assertIn("desequilibre allelique", motif)
        self.assertTrue(self.recevable(enregistrement(100, ad="75,25"))[0])

    def test_refuse_un_genotype_incomplet(self):
        ok, _, _, motif = self.recevable(enregistrement(100, gt="0/."))
        self.assertFalse(ok)
        self.assertIn("non appelee", motif)

    def test_une_heteroplasmie_n_est_pas_un_desequilibre(self):
        # Sur l'ADN mitochondrial, un allele a 10 % est le resultat a rendre,
        # et c'est le seul gene de l'etage qui porte une recommandation.
        info = dict(self.INFO, haploide=True)
        self.assertTrue(self.recevable(enregistrement(100, ad="90,10"), info=info)[0])

    def test_un_second_enregistrement_douteux_suffit(self):
        ok, _, _, motif = self.recevable(enregistrement(100), enregistrement(100, gq=3))
        self.assertFalse(ok)
        self.assertIn("qualite insuffisante", motif)

    def test_la_seconde_position_est_jugee_aussi(self):
        ok, _, _, motif = self.recevable(enregistrement(100), enregistrement(200, gq=3))
        self.assertFalse(ok)
        self.assertIn("chr3:200", motif)


class LectureDesEntrees(unittest.TestCase):
    def test_option_non_fournie_rend_none(self):
        self.assertIsNone(t.contigs(None))

    def test_un_fichier_demande_mais_introuvable_arrete_l_etage(self):
        # Sinon le garde-fou devient contournable par simple suppression du
        # fichier de controle, et un genotype douteux passe en silence.
        with self.assertRaises(SystemExit):
            t.exige(os.path.join(tempfile.gettempdir(), "absent_xyz.tsv"), "contigs")


ENTETE_VCF = ("##fileformat=VCFv4.2\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n")


def ecrire_vcf(chemin, enregistrements, echantillons="ECH"):
    with gzip.open(chemin, "wt", encoding="utf-8", newline="\n") as fh:
        fh.write(ENTETE_VCF % echantillons)
        for c in enregistrements:
            fh.write("\t".join(c) + "\n")


class TypeurEnPanne(unittest.TestCase):
    """Un typeur qui ne tourne pas est une panne de l'installation : l'etage doit
    echouer. Rendu en « non conclusif », il passerait en reussite complete, et
    MT-RNR1 perdrait ses recommandations sans un mot."""

    def execute(self, enregistrements, echantillons="ECH"):
        d = tempfile.mkdtemp()
        p = lambda n: os.path.join(d, n)
        cat = {"genes": {"POR": {"classe_rnpgx": "2", "contig": "chr7",
                                 "positions": ["chr7:75985688"],
                                 "positions_panel": ["chr7:75985688"],
                                 "alleles_positions": {"*28": ["chr7:75985688"]},
                                 "allele_reference": "*1", "phenotype": False,
                                 "appel_externe_pharmcat": False, "haploide": False}}}
        with open(p("cat.json"), "w", encoding="utf-8") as fh:
            json.dump(cat, fh)
        with open(p("prof.txt"), "w", encoding="utf-8") as fh:
            fh.write("chr7\t75985688\t40\n")
        with open(p("contigs.txt"), "w", encoding="utf-8") as fh:
            fh.write("chr7\n")
        ecrire_vcf(p("v.vcf.gz"), enregistrements, echantillons)
        with open(p("panne.py"), "w", encoding="utf-8") as fh:
            fh.write("import sys\nsys.stderr.write('FATAL: image introuvable\\n')\nsys.exit(255)\n")
        argv = ["typage_pypgx.py", "--catalogue", p("cat.json"), "--variants", p("v.vcf.gz"),
                "--profondeur", p("prof.txt"), "--contigs-vcf", p("contigs.txt"),
                "--travail", d, "--sortie", p("sortie.tsv"),
                "--pypgx", "%s %s" % (sys.executable, p("panne.py"))]
        sortie = io.StringIO()
        with mock.patch.object(sys, "argv", argv), mock.patch("sys.stdout", sortie), \
                mock.patch("sys.stderr", io.StringIO()):
            code = t.main()
        if not os.path.exists(p("sortie.tsv")):
            return code, ""
        with open(p("sortie.tsv"), encoding="utf-8") as fh:
            return code, fh.read()

    def test_un_typeur_qui_ne_tourne_pas_fait_echouer_l_etage(self):
        code, lignes = self.execute([enregistrement(75985688, contig="chr7")])
        self.assertEqual(code, 1)
        self.assertIn("\techec\t", lignes)
        self.assertNotIn("non conclusif", lignes)

    def test_un_genotype_douteux_n_atteint_pas_le_typeur(self):
        # Le fichier de variants est juge par l'etage lui-meme : le gene sort non
        # conclusif avec son motif, et le typeur, ici en panne, n'est pas lance.
        code, lignes = self.execute([enregistrement(75985688, contig="chr7", ad="76,24")])
        self.assertEqual(code, 0)
        self.assertIn("\tnon conclusif\t", lignes)
        self.assertIn("chr7:75985688 (desequilibre allelique)", lignes)

    def test_deux_echantillons_arretent_l_etage(self):
        # Le jugement lirait le premier et le typeur les deux.
        code, lignes = self.execute(
            [enregistrement(75985688, contig="chr7") + ["0/1:99:15,15"]],
            echantillons="ECH\tAUTRE")
        self.assertEqual((code, lignes), (1, ""))


class ResultatDuTypeur(unittest.TestCase):
    def archive(self, entrees):
        d = tempfile.mkdtemp()
        with zipfile.ZipFile(os.path.join(d, "results.zip"), "w") as z:
            for nom, contenu in entrees:
                z.writestr(nom, contenu)
        return d

    TABLE = "\tGenotype\tHaplotype1\tHaplotype2\nECH\t*1/*28\t%s\t%s\n"

    def test_lit_la_table(self):
        d = self.archive([("x/data.tsv", self.TABLE % ("*1;", "*28;"))])
        r, err = t.lis_resultat(d)
        self.assertIsNone(err)
        self.assertEqual(r["Genotype"], "*1/*28")

    def test_refuse_un_haplotype_vide(self):
        # Un haplotype vide se lisait comme l'allele de reference : un typeur
        # qui n'a rien rendu sortait « *1/*1 », statut rendu.
        for h1, h2 in (("", ""), ("*28;", ""), ("", "*28;"), (";", "*28;")):
            r, err = t.lis_resultat(self.archive([("x/data.tsv", self.TABLE % (h1, h2))]))
            self.assertIsNone(r, (h1, h2))
            self.assertIn("haplotype", err)

    def test_refuse_une_table_sans_colonnes_d_haplotype(self):
        r, err = t.lis_resultat(self.archive([("x/data.tsv", "\tGenotype\nECH\t*1/*28\n")]))
        self.assertIsNone(r)
        self.assertIn("haplotype", err)

    def test_refuse_plusieurs_tables(self):
        # Prendre la premiere reviendrait a tirer au sort le resultat.
        d = self.archive([("a/data.tsv", "\tGenotype\nE\t*1\n"),
                          ("b/data.tsv", "\tGenotype\nE\t*2\n")])
        r, err = t.lis_resultat(d)
        self.assertIsNone(r)
        self.assertIn("data.tsv", err)

    def test_refuse_plusieurs_echantillons(self):
        d = self.archive([("a/data.tsv", "\tGenotype\nE1\t*1\nE2\t*2\n")])
        r, err = t.lis_resultat(d)
        self.assertIsNone(r)
        self.assertIn("echantillons", err)

    def test_refuse_une_archive_absente(self):
        r, err = t.lis_resultat(tempfile.mkdtemp())
        self.assertIsNone(r)
        self.assertIn("results.zip", err)


class Fraction(unittest.TestCase):
    def test_lit_la_fraction_de_l_allele(self):
        l = ligne("m.1555A>G;", "Reference;",
                  donnees="m.1555A>G:M-1555-A-G:0.62;Reference:default;")
        self.assertAlmostEqual(t.fraction(l, "m.1555A>G"), 0.62)

    def test_allele_absent_des_donnees(self):
        self.assertIsNone(t.fraction(ligne("a;", "b;"), "a"))


class CatalogueDuDepot(unittest.TestCase):
    """Le catalogue livre doit porter ce que le module lui demande."""

    def setUp(self):
        with open(os.path.join(RACINE, "ressources", "pypgx_genes.json"),
                  encoding="utf-8") as fh:
            self.genes = json.load(fh)["genes"]

    def test_les_quatre_genes_sont_la(self):
        self.assertEqual(sorted(self.genes),
                         ["BCHE", "MT-RNR1", "MTHFR", "POR"])

    def test_chaque_gene_porte_les_champs_requis(self):
        for nom, info in self.genes.items():
            for champ in ("classe_rnpgx", "contig", "positions",
                          "positions_panel", "alleles_positions",
                          "allele_reference", "haploide",
                          "appel_externe_pharmcat"):
                self.assertIn(champ, info, "%s sans %s" % (nom, champ))

    def test_chaque_position_du_panel_est_portee_par_un_allele(self):
        # Une position du panel qu'aucun allele ne definit ne serait jamais
        # rendue : le gene sortirait toujours de reference.
        for nom, info in self.genes.items():
            portees = set()
            for p in info["alleles_positions"].values():
                portees.update(p)
            for p in info["positions_panel"]:
                self.assertIn(p, portees, "%s : %s sans allele" % (nom, p))

    def test_seul_mt_rnr1_est_haploide_et_porte_un_phenotype(self):
        for nom, info in self.genes.items():
            attendu = (nom == "MT-RNR1")
            self.assertEqual(bool(info["haploide"]), attendu, nom)
            self.assertEqual(bool(info["appel_externe_pharmcat"]), attendu, nom)


if __name__ == "__main__":
    unittest.main(verbosity=2)
