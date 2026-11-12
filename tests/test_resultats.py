# -*- coding: utf-8 -*-
"""Tests des sorties structurees : le genotype d'un cote, l'annotation de l'autre.

Le genotype est ce qui a ete mesure ; l'annotation est ce qu'un referentiel en
dit aujourd'hui. Les deux fichiers se rejoignent par l'identifiant HGNC, et
doivent dire la meme chose que le compte rendu : ils lisent la meme logique.
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))
sys.path.insert(0, os.path.join(RACINE, "outils"))

import resultats as R  # noqa: E402

HGNC = {"CYP2C19": "HGNC:2621", "CYP2D6": "HGNC:2625", "HLA-B": "HGNC:4932",
        "TPMT": "HGNC:12014", "POR": "HGNC:9208", "BCHE": "HGNC:983",
        "MT-RNR1": "HGNC:7470", "ABCG2": "HGNC:74"}


def gene(labels=(), phenotypes=(), source=None, score=None, variants=()):
    d = {"recommendationDiplotypes":
         [{"label": l, "phenotypes": list(phenotypes), "activityScore": score}
          for l in labels],
         "variants": list(variants)}
    if source:
        d["callSource"] = source
    return d


def qc(statut="complet", attendues=3, perdues=0):
    return {"statut": statut, "positions_attendues": attendues,
            "positions_retenues": attendues - perdues, "positions_perdues": perdues}


def position(gene, pos, rsid, ref, alt, gt, retenue="oui", statut="mesuree",
             chrom="chr10"):
    return {"chrom": chrom, "pos": str(pos), "rsid": rsid, "ref": ref, "alt": alt,
            "genes": gene, "genotype": gt, "retenue": retenue, "statut": statut}


def par_gene(lignes):
    return {l["gene"]: l for l in lignes}


class Identifiant(unittest.TestCase):
    def test_le_symbole_est_traduit(self):
        self.assertEqual(R.hgnc(HGNC, "CYP2D6"), "HGNC:2625")

    def test_un_gene_sans_identifiant_arrete_tout(self):
        # Une cellule vide casserait la jointure sans le dire.
        with self.assertRaises(ValueError):
            R.hgnc(HGNC, "CYP9Z9")


class Genotypes(unittest.TestCase):
    def rendre(self, genes, meta_genes, clinique, typage=(), hors=()):
        return par_gene(R.genotypes(
            {"genes": genes},
            {"genes": meta_genes, "perimetre_clinique": clinique,
             "perimetre_clinique_hors_interpreteur": list(hors)},
            list(typage), HGNC))

    def test_un_gene_complet_porte_son_identifiant_et_son_diplotype(self):
        g = self.rendre({"CYP2C19": gene(["*1/*17"], ["Rapid Metabolizer"])},
                        {"CYP2C19": qc(attendues=35)}, ["CYP2C19"])["CYP2C19"]
        self.assertEqual(
            (g["hgnc_id"], g["genotype"], g["statut"], g["source"],
             g["positions_attendues"], g["positions_retenues"]),
            ("HGNC:2621", "*1/*17", "rendu", "PharmCAT", 35, 35))

    def test_le_genotype_ne_porte_aucune_annotation(self):
        g = self.rendre({"CYP2C19": gene(["*1/*17"], ["Rapid Metabolizer"])},
                        {"CYP2C19": qc()}, ["CYP2C19"])["CYP2C19"]
        self.assertNotIn("Rapid", " ".join(str(v) for v in g.values()))
        self.assertNotIn("phenotype", g)

    def test_un_gene_partiel_est_rendu_sous_reserve_avec_son_denominateur(self):
        g = self.rendre({"TPMT": gene(["*1/*1"], ["Normal Metabolizer"])},
                        {"TPMT": qc("partiel", 45, 2)}, ["TPMT"])["TPMT"]
        self.assertEqual((g["statut"], g["positions_retenues"]), ("reserve", 43))
        self.assertEqual(g["motif"], "2 positions non lues sur 45")

    def test_une_seule_position_perdue_s_ecrit_au_singulier(self):
        g = self.rendre({"TPMT": gene(["*1/*1"])}, {"TPMT": qc("partiel", 45, 1)},
                        ["TPMT"])["TPMT"]
        self.assertEqual(g["motif"], "1 position non lue sur 45")

    def test_un_diplotype_ambigu_garde_toutes_ses_lectures(self):
        g = self.rendre({"TPMT": gene(["*1/*3A", "*3B/*3C"], ["Intermediate Metabolizer"])},
                        {"TPMT": qc()}, ["TPMT"])["TPMT"]
        self.assertEqual((g["statut"], g["genotype"]), ("ambigu", "*1/*3A | *3B/*3C"))

    def test_un_gene_non_rendu_n_a_pas_de_genotype(self):
        g = self.rendre({"TPMT": gene([])}, {"TPMT": qc("absent", 45, 45)},
                        ["TPMT"])["TPMT"]
        self.assertEqual((g["statut"], g["genotype"]), ("non_rendu", ""))
        self.assertTrue(g["motif"])

    def test_un_appel_externe_nomme_son_outil_et_non_des_positions(self):
        r = self.rendre(
            {"CYP2D6": gene(["*1/*4"], ["Intermediate Metabolizer"], "OUTSIDE"),
             "HLA-B": gene(["*07:02/*58:01"], ["*58:01 positive"], "OUTSIDE")},
            {"CYP2D6": qc(attendues=157)}, ["CYP2D6", "HLA-B"])
        self.assertEqual((r["CYP2D6"]["source"], r["HLA-B"]["source"]),
                         ("Cyrius", "OptiType"))
        self.assertEqual(r["CYP2D6"]["positions_attendues"], "")

    def test_un_gene_hors_interpreteur_est_rendu_par_le_typage_complementaire(self):
        typage = [{"gene": "POR", "statut": "rendu", "genotype": "*1/*28",
                   "positions_lues": "46", "positions_attendues": "46", "motif": ""}]
        g = self.rendre({}, {}, ["POR"], typage, hors=["POR"])["POR"]
        self.assertEqual((g["genotype"], g["statut"], g["source"], g["positions_retenues"]),
                         ("*1/*28", "rendu", "PyPGx", "46"))

    def test_mt_rnr1_passe_par_l_interpreteur_garde_sa_source_et_ses_positions(self):
        typage = [{"gene": "MT-RNR1", "statut": "rendu", "genotype": "Reference/Reference",
                   "positions_lues": "23", "positions_attendues": "23", "motif": ""}]
        g = self.rendre({"MT-RNR1": gene(["Reference"], ["normal risk"], "OUTSIDE")},
                        {}, ["MT-RNR1"], typage)["MT-RNR1"]
        self.assertEqual((g["genotype"], g["source"], g["positions_attendues"]),
                         ("Reference", "PyPGx", "23"))

    def test_un_gene_type_hors_perimetre_clinique_est_la(self):
        typage = [{"gene": "BCHE", "statut": "rendu", "genotype": "Reference/rs1803274",
                   "positions_lues": "2", "positions_attendues": "2", "motif": ""}]
        g = self.rendre({}, {}, ["CYP2C19"], typage)["BCHE"]
        self.assertEqual((g["hgnc_id"], g["statut"], g["genotype"]),
                         ("HGNC:983", "rendu", "Reference/rs1803274"))

    def test_un_typage_non_conclusif_n_est_pas_un_genotype(self):
        typage = [{"gene": "BCHE", "statut": "non conclusif", "genotype": "",
                   "positions_lues": "1", "positions_attendues": "2",
                   "motif": "1 position non lue"}]
        g = self.rendre({}, {}, [], typage)["BCHE"]
        self.assertEqual((g["statut"], g["genotype"], g["motif"]),
                         ("non_rendu", "", "1 position non lue"))

    def test_un_gene_hors_perimetre_clinique_est_ignore(self):
        self.assertEqual(self.rendre({"CYP2C19": gene(["*1/*1"])},
                                     {"CYP2C19": qc()}, ["TPMT"]).get("CYP2C19"), None)

    def test_les_lignes_sont_triees(self):
        lignes = R.genotypes(
            {"genes": {"TPMT": gene(["*1/*1"]), "CYP2C19": gene(["*1/*1"])}},
            {"genes": {"TPMT": qc(), "CYP2C19": qc()},
             "perimetre_clinique": ["TPMT", "CYP2C19"]}, [], HGNC)
        self.assertEqual([l["gene"] for l in lignes], ["CYP2C19", "TPMT"])


class Variants(unittest.TestCase):
    GENO = [{"gene": "CYP2C19", "source": "PharmCAT", "statut": "rendu"},
            {"gene": "CYP2D6", "source": "Cyrius", "statut": "rendu"},
            {"gene": "POR", "source": "PyPGx", "statut": "rendu"}]

    def rendre(self, positions, rapport=None):
        return R.variants(positions, rapport or {"genes": {}}, self.GENO, HGNC)

    def test_un_variant_observe_sort_avec_son_rs_et_ses_alleles(self):
        v = self.rendre([position("CYP2C19", 94761900, "rs12248560", "C", "T", "0/1")])
        self.assertEqual(len(v), 1)
        self.assertEqual(
            (v[0]["hgnc_id"], v[0]["gene"], v[0]["chrom"], v[0]["pos"], v[0]["rsid"],
             v[0]["ref"], v[0]["alt"], v[0]["genotype"], v[0]["retenue"]),
            ("HGNC:2621", "CYP2C19", "chr10", 94761900, "rs12248560", "C", "T",
             "C/T", "oui"))

    def test_une_position_de_reference_n_est_pas_un_variant(self):
        self.assertEqual(self.rendre([
            position("CYP2C19", 1, "rs1", "C", "T", "0/0"),
            position("CYP2C19", 2, "rs2", "C", "T", "", statut="reference, lue")]), [])

    def test_un_site_multiallelique_nomme_l_allele_porte(self):
        v = self.rendre([position("CYP2C19", 5, "rs5", "G", "A,C", "0/2")])
        self.assertEqual(v[0]["genotype"], "G/C")

    def test_le_phasage_est_conserve(self):
        v = self.rendre([position("CYP2C19", 5, "rs5", "G", "A", "1|0")])
        self.assertEqual(v[0]["genotype"], "A|G")

    def test_un_genotype_haploide(self):
        v = self.rendre([position("POR", 5, "rs5", "A", "G", "1", chrom="chr7")])
        self.assertEqual(v[0]["genotype"], "G")

    def test_une_position_sans_rs_garde_ses_coordonnees(self):
        v = self.rendre([position("CYP2C19", 7, ".", "C", "T", "1/1")])
        self.assertEqual((v[0]["rsid"], v[0]["pos"], v[0]["genotype"]), ("", 7, "T/T"))

    def test_un_variant_ecarte_par_le_controle_qualite_reste_visible(self):
        # Vu puis ecarte n'est pas « reference » : le lecteur doit le savoir.
        v = self.rendre([position("CYP2C19", 7, "rs7", "C", "T", "0/1", "non",
                                  "qualite insuffisante")])
        self.assertEqual((v[0]["retenue"], v[0]["statut"]),
                         ("non", "qualite insuffisante"))

    def test_un_gene_appele_par_un_outil_dedie_n_a_pas_de_variants(self):
        # Le diplotype de CYP2D6 vient de Cyrius : les petits variants du
        # fichier de variants ne l'ont pas fonde, et peuvent le contredire.
        self.assertEqual(self.rendre([position("CYP2D6", 9, "rs9", "C", "T", "0/1",
                                               chrom="chr22")]), [])

    def test_un_gene_absent_des_genotypes_n_a_pas_de_variants(self):
        self.assertEqual(self.rendre([position("TPMT", 9, "rs9", "C", "T", "0/1")]), [])

    def test_les_alleles_definis_viennent_de_l_interpreteur(self):
        rapport = {"genes": {"CYP2C19": gene(variants=[
            {"chromosome": "chr10", "position": 94761900, "alleles": ["*4", "*17"]}])}}
        v = self.rendre([position("CYP2C19", 94761900, "rs12248560", "C", "T", "0/1")],
                        rapport)
        self.assertEqual(v[0]["alleles_definis"], "*4;*17")

    def test_les_variants_sont_tries_par_gene_puis_par_position(self):
        v = self.rendre([position("POR", 30, "rs3", "C", "T", "0/1", chrom="chr7"),
                         position("CYP2C19", 200, "rs2", "C", "T", "0/1"),
                         position("CYP2C19", 30, "rs1", "C", "T", "0/1")])
        self.assertEqual([x["rsid"] for x in v], ["rs1", "rs2", "rs3"])


TR = {"version": "9.9",
      "medicaments_restitues": {"clopidogrel": "clopidogrel", "codeine": "codéine"},
      "traductions": {"Avoid.": {"fr": "Ne pas prescrire.", "type": "contre-indication"}}}


def rapport_annote(genes, consignes):
    return {"pharmcatVersion": "3.4.0", "dataVersion": "2026-01-01", "genes": genes,
            "drugs": {"CPIC Guideline Annotation": {
                nom: {"citations": [{"pmid": "123", "year": 2020}],
                      "guidelines": [{"annotations": [
                          {"classification": "Strong", "lookupKey": [cle],
                           "drugRecommendation": texte}]}]}
                for nom, (cle, texte) in consignes.items()}}}


class Annotations(unittest.TestCase):
    def rendre(self, genes, consignes, clinique):
        meta = {"genes": {g: qc() for g in genes}, "perimetre_clinique": clinique}
        return R.annotations(rapport_annote(genes, consignes), meta, TR, HGNC)

    def test_une_consigne_sort_avec_son_gene_son_phenotype_et_sa_source(self):
        a = self.rendre(
            {"CYP2C19": gene(["*2/*2"], ["Poor Metabolizer"], score="0.0")},
            {"clopidogrel": ({"CYP2C19": "Poor Metabolizer"}, "Avoid.")}, ["CYP2C19"])
        self.assertEqual(len(a), 1)
        self.assertEqual(
            (a[0]["hgnc_id"], a[0]["gene"], a[0]["phenotype"], a[0]["phenotype_fr"],
             a[0]["score_activite"], a[0]["medicament"], a[0]["medicament_fr"],
             a[0]["type"], a[0]["consigne"], a[0]["source"], a[0]["pmid"],
             a[0]["traduite"]),
            ("HGNC:2621", "CYP2C19", "Poor Metabolizer", "métaboliseur lent", "0.0",
             "clopidogrel", "clopidogrel", "contre-indication", "Ne pas prescrire.",
             "CPIC, forte", "123", "oui"))

    def test_l_annotation_dit_de_quel_referentiel_elle_vient(self):
        a = self.rendre({"CYP2C19": gene(["*1/*1"], ["Normal Metabolizer"])}, {},
                        ["CYP2C19"])
        self.assertEqual(a[0]["referentiel"],
                         "PharmCAT 3.4.0 ; CPIC 2026-01-01 ; traductions 9.9")

    def test_l_annotation_ne_repete_pas_le_genotype(self):
        a = self.rendre({"CYP2C19": gene(["*2/*2"], ["Poor Metabolizer"])}, {},
                        ["CYP2C19"])
        self.assertNotIn("genotype", a[0])
        self.assertNotIn("*2/*2", " ".join(str(v) for v in a[0].values()))

    def test_un_gene_sans_consigne_garde_son_phenotype(self):
        a = self.rendre({"CYP2C19": gene(["*1/*1"], ["Normal Metabolizer"])}, {},
                        ["CYP2C19"])
        self.assertEqual((len(a), a[0]["phenotype"], a[0]["medicament"]),
                         (1, "Normal Metabolizer", ""))
        self.assertEqual(a[0]["score_activite"], "")

    def test_une_consigne_ne_s_attache_qu_au_gene_qui_la_declenche(self):
        a = self.rendre(
            {"CYP2C19": gene(["*2/*2"], ["Poor Metabolizer"]),
             "CYP2D6": gene(["*1/*1"], ["Normal Metabolizer"], "OUTSIDE")},
            {"clopidogrel": ({"CYP2C19": "Poor Metabolizer"}, "Avoid.")},
            ["CYP2C19", "CYP2D6"])
        self.assertEqual([(x["gene"], x["medicament"]) for x in a],
                         [("CYP2C19", "clopidogrel"), ("CYP2D6", "")])

    def test_un_texte_non_traduit_est_signale(self):
        a = self.rendre(
            {"CYP2D6": gene(["*4/*4"], ["Poor Metabolizer"], "OUTSIDE")},
            {"codeine": ({"CYP2D6": "Poor Metabolizer"}, "Texte inconnu.")}, ["CYP2D6"])
        self.assertEqual((a[0]["consigne"], a[0]["traduite"]), ("Texte inconnu.", "non"))

    def test_un_gene_ambigu_n_a_aucune_consigne(self):
        a = self.rendre(
            {"CYP2C19": gene(["*2/*2", "*1/*2"], ["Poor Metabolizer"])},
            {"clopidogrel": ({"CYP2C19": "Poor Metabolizer"}, "Avoid.")}, ["CYP2C19"])
        self.assertEqual([x["medicament"] for x in a], [""])

    def test_un_gene_non_rendu_n_a_pas_d_annotation(self):
        self.assertEqual(self.rendre({"CYP2C19": gene([])}, {}, ["CYP2C19"]), [])

    def test_une_consigne_a_deux_genes_donne_une_ligne_par_gene(self):
        a = self.rendre(
            {"CYP2C19": gene(["*2/*2"], ["Poor Metabolizer"]),
             "CYP2D6": gene(["*4/*4"], ["Poor Metabolizer"], "OUTSIDE")},
            {"clopidogrel": ({"CYP2C19": "Poor Metabolizer",
                              "CYP2D6": "Poor Metabolizer"}, "Avoid.")},
            ["CYP2C19", "CYP2D6"])
        self.assertEqual([(x["gene"], x["medicament"]) for x in a],
                         [("CYP2C19", "clopidogrel"), ("CYP2D6", "clopidogrel")])


class TableDuDepot(unittest.TestCase):
    def test_tout_gene_des_ressources_a_son_identifiant(self):
        import construire_hgnc as C
        res = os.path.join(RACINE, "ressources")
        with open(os.path.join(res, "genes_hgnc.json"), encoding="utf-8") as fh:
            table = json.load(fh)["genes"]
        self.assertEqual(sorted(table), C.genes_des_ressources(res))
        for sym, ident in table.items():
            self.assertRegex(ident, r"^HGNC:\d+$", sym)
        self.assertEqual(len(set(table.values())), len(table))

    def test_une_etiquette_qui_n_est_pas_un_gene_n_y_entre_pas(self):
        import construire_hgnc as C
        self.assertNotIn("CYP2C_cluster",
                         C.genes_des_ressources(os.path.join(RACINE, "ressources")))


class Execution(unittest.TestCase):
    """Le script, sur le genome de qualification du depot."""

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="pgx_res_")
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.travail = os.path.join(self.d, "travail")
        os.makedirs(self.travail)
        q = os.path.join(RACINE, "bench", "qualification", "HG00699")
        self.rapport = os.path.join(q, "sortie", "HG00699.report.json")
        shutil.copy(os.path.join(q, "travail", "pypgx.tsv"), self.travail)
        with open(self.rapport, encoding="utf-8") as fh:
            genes = json.load(fh)["genes"]
        clinique = ["CYP2C19", "CYP2D6", "HLA-B", "POR", "TPMT"]
        with open(os.path.join(self.travail, "perimetre.json"), "w",
                  encoding="utf-8") as fh:
            json.dump({"genes": {g: qc() for g in clinique if g in genes},
                       "perimetre_clinique": clinique,
                       "perimetre_clinique_hors_interpreteur": ["POR"]}, fh)
        self.positions([position("CYP2C19", 94761900, "rs12248560", "C", "T", "0/1")])

    def positions(self, lignes):
        with open(os.path.join(self.travail, "qc_positions.tsv"), "w",
                  encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(position("", 0, "", "", "", "")),
                               delimiter="\t")
            w.writeheader()
            w.writerows(lignes)

    def lancer(self, hgnc=None):
        return subprocess.run(
            [sys.executable, "-B", os.path.join(RACINE, "bin", "resultats.py"),
             "--rapport", self.rapport,
             "--perimetre", os.path.join(self.travail, "perimetre.json"),
             "--echantillon", "ECH", "--sortie", self.d]
            + (["--hgnc", hgnc] if hgnc else []),
            capture_output=True, text=True, encoding="utf-8")

    def lire(self, nom):
        with open(os.path.join(self.d, "ECH.%s.tsv" % nom), encoding="utf-8",
                  newline="") as fh:
            brut = fh.read()
        self.assertNotIn("\r", brut)
        return list(csv.DictReader(brut.splitlines(), delimiter="\t"))

    def test_les_trois_fichiers_se_rejoignent_par_l_identifiant(self):
        r = self.lancer()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        g, v, a = self.lire("genotypes"), self.lire("variants"), self.lire("annotations")
        self.assertEqual([x["gene"] for x in g],
                         ["BCHE", "CYP2C19", "CYP2D6", "HLA-B", "MT-RNR1", "MTHFR",
                          "POR", "TPMT"])
        ids = {x["hgnc_id"] for x in g}
        self.assertTrue({x["hgnc_id"] for x in v} <= ids)
        self.assertTrue({x["hgnc_id"] for x in a} <= ids)
        self.assertEqual(v[0]["rsid"], "rs12248560")
        self.assertEqual(list(g[0]), R.COLONNES["genotypes"])
        self.assertEqual(list(v[0]), R.COLONNES["variants"])
        self.assertEqual(list(a[0]), R.COLONNES["annotations"])

    def test_deux_executions_rendent_les_memes_octets(self):
        def octets():
            out = {}
            for n in sorted(os.listdir(self.d)):
                if n.endswith(".tsv"):
                    with open(os.path.join(self.d, n), "rb") as fh:
                        out[n] = fh.read()
            return out
        self.lancer()
        avant = octets()
        self.lancer()
        self.assertEqual(len(avant), 3)
        self.assertEqual(octets(), avant)

    def test_un_gene_sans_identifiant_fait_echouer_le_script(self):
        table = os.path.join(self.d, "hgnc.json")
        with open(table, "w", encoding="utf-8") as fh:
            json.dump({"genes": {"CYP2C19": "HGNC:2621"}}, fh)
        r = self.lancer(table)
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertIn("sorties structurees refusees", r.stderr)
        self.assertIn("HGNC", r.stderr)
        self.assertNotIn("Traceback", r.stderr)
        self.assertEqual([n for n in os.listdir(self.d) if n.endswith(".tsv")], [])

    def test_sans_mesure_des_positions_le_script_echoue(self):
        os.remove(os.path.join(self.travail, "qc_positions.tsv"))
        self.assertNotEqual(self.lancer().returncode, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
