# -*- coding: utf-8 -*-
"""Tests de la logique du compte rendu.

C'est la partie qui decide ce que le clinicien lit. Elle a deja perdu quatre
conduites a tenir sur neuf, dont la codeine et le tramadol, et elle n'etait pas
verifiable : le module de rendu lit ses arguments a l'import.
"""
import os
import sys
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))

import cr_logique as L  # noqa: E402


def gene(labels=(), phenotypes=(), source=None):
    d = {"recommendationDiplotypes":
         [{"label": l, "phenotypes": list(phenotypes)} for l in labels]}
    if source:
        d["callSource"] = source
    return d


def annotation(classification, cle, texte, **kw):
    a = {"classification": classification,
         "lookupKey": [cle],
         "drugRecommendation": texte}
    a.update(kw)
    return a


def rapport(annotations_par_med):
    return {"drugs": {"CPIC Guideline Annotation": {
        nom: {"citations": [{"pmid": "123", "year": 2020}],
              "guidelines": [{"annotations": ann}]}
        for nom, ann in annotations_par_med.items()}}}


TR = {
    "traductions": {
        "Avoid codeine use.": {"fr": "Ne pas prescrire.", "type": "contre-indication"},
        "Use at standard dose.": {"fr": "Dose standard.", "type": "standard"},
        "Reduce starting dose.": {"fr": "Réduire la dose.", "type": "adaptation"},
        "Monitor.": {"fr": "Surveiller.", "type": "vigilance"},
    },
    "regles_hla_rnpgx": {
        "HLA-B*15:11": {"fr": "Éviter la carbamazépine.", "type": "éviter",
                        "source": "RNPGx", "medicaments": ["carbamazepine"]},
    },
}
MEDS = {"codeine": "codéine", "tramadol": "tramadol",
        "carbamazepine": "carbamazépine", "simvastatin": "simvastatine",
        "atorvastatin": "atorvastatine", "gentamicin": "gentamicine",
        "amikacin": "amikacine"}


class Libelle(unittest.TestCase):
    def test_traduit_les_termes_de_l_interpreteur(self):
        g = {"X": gene(labels=["Reference/Unknown"])}
        self.assertEqual(L.libelle(g, "X"), "référence/non déterminé")

    def test_ambiguite_annoncee_comme_telle(self):
        # Deux diplotypes possibles : le compte rendu doit le dire, et aucune
        # conduite ne doit en decouler.
        g = {"X": gene(labels=["*1/*2", "*1/*3"])}
        self.assertTrue(L.libelle(g, "X").startswith("ambigu : "))
        self.assertEqual(L.ambigus(g, ["X"]), {"X"})

    def test_abcg2_recompose_le_rsid_par_allele(self):
        g = {"ABCG2": gene(labels=["rs2231142 reference (G)/rs2231142 variant (T)"])}
        self.assertEqual(L.libelle(g, "ABCG2"), "rs2231142 G/T")

    def test_gene_sans_diplotype(self):
        self.assertEqual(L.libelle({}, "X"), "—")


class Phenotype(unittest.TestCase):
    def test_traduit_un_metabolisme(self):
        g = {"CYP2D6": gene(labels=["*1/*4"], phenotypes=["Poor Metabolizer"])}
        self.assertEqual(L.phenotype(g, "CYP2D6"), "métaboliseur lent")

    def test_cyp3a5_porte_la_mention_d_expression(self):
        # Un metaboliseur normal de CYP3A5 est un expresseur : sans la mention,
        # le lecteur conclut l'inverse de ce qui compte pour le tacrolimus.
        g = {"CYP3A5": gene(labels=["*1/*1"], phenotypes=["Normal Metabolizer"])}
        self.assertIn("expresseur", L.phenotype(g, "CYP3A5"))

    def test_mt_rnr1_rend_un_risque_et_non_un_metabolisme(self):
        g = {"MT-RNR1": gene(
            labels=["Reference"],
            phenotypes=["increased risk of aminoglycoside-induced hearing loss"])}
        self.assertEqual(L.phenotype(g, "MT-RNR1"),
                         "risque augmenté de surdité sous aminoside")

    def test_hla_present(self):
        g = {"HLA-B": gene(labels=["*57:01/*07:02"],
                           phenotypes=["B*57:01 positive"])}
        self.assertEqual(L.phenotype(g, "HLA-B"), "B*57:01 présent")

    def test_hla_absent_au_pluriel(self):
        g = {"HLA-B": gene(labels=["*07:02/*08:01"],
                           phenotypes=["B*57:01 negative", "B*58:01 negative"])}
        self.assertIn("absents", L.phenotype(g, "HLA-B"))

    def test_phenotype_inconnu_rendu_tel_quel(self):
        # Un terme que la table ne connait pas ne doit pas disparaitre.
        g = {"X": gene(labels=["*1/*1"], phenotypes=["Terme Inedit"])}
        self.assertEqual(L.phenotype(g, "X"), "Terme Inedit")


class GeneNormal(unittest.TestCase):
    def test_phenotype_normal(self):
        g = {"X": gene(labels=["*1/*1"], phenotypes=["Normal Metabolizer"])}
        self.assertTrue(L.gene_normal(g, "X"))

    def test_phenotype_anormal(self):
        g = {"X": gene(labels=["*1/*4"], phenotypes=["Poor Metabolizer"])}
        self.assertFalse(L.gene_normal(g, "X"))

    def test_hla_negatif_compte_comme_normal(self):
        g = {"HLA-B": gene(labels=["*07:02/*08:01"],
                           phenotypes=["B*57:01 negative"])}
        self.assertTrue(L.gene_normal(g, "HLA-B"))

    def test_sans_phenotype_n_est_pas_normal(self):
        # Sans phenotype, rien ne permet de dire que le gene ne declenche pas
        # de consigne ; l'affirmer retirerait la consigne.
        self.assertFalse(L.gene_normal({"X": gene(labels=["*1/*1"])}, "X"))


class Perimetre(unittest.TestCase):
    def test_gene_complet_est_rendu(self):
        rap = {"genes": {"TPMT": gene(labels=["*1/*1"])}}
        meta = {"perimetre_clinique": ["TPMT"],
                "genes": {"TPMT": {"statut": "complet"}}}
        rendus, absents, reserves = L.perimetre(rap, meta)
        self.assertEqual(rendus, ["TPMT"])
        self.assertEqual(absents, [])

    def test_gene_partiel_est_rendu_avec_reserve(self):
        rap = {"genes": {"TPMT": gene(labels=["*1/*1"])}}
        meta = {"perimetre_clinique": ["TPMT"],
                "genes": {"TPMT": {"statut": "partiel",
                                   "positions_perdues": 4,
                                   "positions_attendues": 7}}}
        rendus, _, reserves = L.perimetre(rap, meta)
        self.assertEqual(rendus, ["TPMT"])
        # La reserve porte son denominateur : quatre positions perdues sur
        # sept et quatre sur quarante-six ne se lisent pas de la meme facon,
        # et la titration a montre qu'un faux diplotype sort dans le premier
        # cas et pas dans le second.
        self.assertEqual(reserves["TPMT"], (4, 7))

    def test_aucune_position_lue_n_est_pas_rendu(self):
        # Le coeur de la garantie : un diplotype calcule sur des positions non
        # sequencees ne doit pas sortir.
        rap = {"genes": {"TPMT": gene(labels=["*1/*1"])}}
        meta = {"perimetre_clinique": ["TPMT"],
                "genes": {"TPMT": {"statut": "absent"}}}
        rendus, absents, _ = L.perimetre(rap, meta)
        self.assertEqual(rendus, [])
        self.assertEqual(absents[0][0], "TPMT")

    def test_un_appel_externe_echappe_au_jugement_sur_les_positions(self):
        # CYP2D6, HLA et MT-RNR1 ne viennent pas du fichier de variants : les
        # juger sur ses positions les refuserait a tort.
        rap = {"genes": {"MT-RNR1": gene(labels=["Reference"], source="OUTSIDE")}}
        meta = {"perimetre_clinique": ["MT-RNR1"],
                "genes": {"MT-RNR1": {"statut": "absent"}}}
        rendus, _, _ = L.perimetre(rap, meta)
        self.assertEqual(rendus, ["MT-RNR1"])

    def test_diplotype_inconnu_n_est_pas_rendu(self):
        rap = {"genes": {"TPMT": gene(labels=["Unknown/Unknown"])}}
        meta = {"perimetre_clinique": ["TPMT"],
                "genes": {"TPMT": {"statut": "complet"}}}
        rendus, absents, _ = L.perimetre(rap, meta)
        self.assertEqual(rendus, [])

    def test_gene_hors_tables_de_l_interpreteur(self):
        rap = {"genes": {}}
        meta = {"perimetre_clinique": ["POR"],
                "perimetre_clinique_hors_interpreteur": ["POR"]}
        _, absents, _ = L.perimetre(rap, meta)
        self.assertEqual(absents[0][0], "POR")
        self.assertIn("absent des tables", absents[0][1])

    def test_un_gene_hors_perimetre_clinique_est_ignore(self):
        rap = {"genes": {"CYP4F2": gene(labels=["*1/*1"])}}
        meta = {"perimetre_clinique": ["TPMT"], "genes": {}}
        rendus, absents, _ = L.perimetre(rap, meta)
        self.assertNotIn("CYP4F2", rendus)
        self.assertNotIn("CYP4F2", [a[0] for a in absents])


class Recommandations(unittest.TestCase):
    GENES = {"CYP2D6": gene(labels=["*1/*4"], phenotypes=["Poor Metabolizer"])}

    def selection(self, rep, genes=None, mesures=("CYP2D6",), ambigus=()):
        return L.recommandations(rep, genes or self.GENES, set(mesures),
                                 set(ambigus), TR, MEDS)

    def test_une_contre_indication_forte_est_restituee(self):
        rep = rapport({"codeine": [annotation(
            "Strong", {"CYP2D6": "Poor Metabolizer"}, "Avoid codeine use.")]})
        b = self.selection(rep)
        self.assertEqual(len(b), 1)
        self.assertEqual(b[0][1], "codeine")
        self.assertEqual(b[0][3], "contre-indication")

    def test_une_recommandation_faible_est_ecartee(self):
        # Seules les fortes et les moderees modifient la prise en charge.
        rep = rapport({"codeine": [annotation(
            "Optional", {"CYP2D6": "Poor Metabolizer"}, "Avoid codeine use.")]})
        self.assertEqual(self.selection(rep), [])

    def test_une_consigne_de_dose_standard_n_est_pas_restituee(self):
        # « Dose standard » n'est pas une conduite a tenir : la restituer
        # noierait les quatre lignes qui comptent.
        rep = rapport({"codeine": [annotation(
            "Strong", {"CYP2D6": "Normal Metabolizer"}, "Use at standard dose.")]})
        self.assertEqual(self.selection(rep), [])

    def test_un_gene_hors_perimetre_n_en_declenche_aucune(self):
        rep = rapport({"codeine": [annotation(
            "Strong", {"CYP2C19": "Poor Metabolizer"}, "Avoid codeine use.")]})
        self.assertEqual(self.selection(rep), [])

    def test_un_gene_ambigu_n_en_declenche_aucune(self):
        rep = rapport({"codeine": [annotation(
            "Strong", {"CYP2D6": "Poor Metabolizer"}, "Avoid codeine use.")]})
        self.assertEqual(self.selection(rep, ambigus=("CYP2D6",)), [])

    def test_un_texte_inconnu_de_la_ressource_est_rendu_tel_quel(self):
        # Jamais tu : un texte non traduit sort en anglais et le compte rendu
        # le signale, plutot que de perdre la conduite.
        rep = rapport({"codeine": [annotation(
            "Strong", {"CYP2D6": "Poor Metabolizer"}, "Texte inedit de CPIC.",
            dosingInformation=True)]})
        b = self.selection(rep)
        self.assertEqual(len(b), 1)
        self.assertEqual(b[0][2], "Texte inedit de CPIC.")
        self.assertFalse(b[0][7])          # marque comme non traduit

    def test_la_regle_hla_hors_interpreteur_s_applique(self):
        # HLA-B*15:11 est de classe 1 au RNPGx et l'interpreteur ne l'evalue
        # pas : sans cette regle, un porteur ne recevrait aucune consigne.
        genes = {"HLA-B": gene(labels=["*15:11/*07:02"],
                               phenotypes=["B*15:02 negative"])}
        b = L.recommandations(rapport({}), genes, {"HLA-B"}, set(), TR, MEDS)
        self.assertEqual(len(b), 1)
        self.assertEqual(b[0][1], "carbamazepine")
        self.assertEqual(b[0][4], "RNPGx")

    def test_le_gene_cite_est_celui_qui_declenche(self):
        # Avec deux genes dans la cle dont un normal, c'est l'anormal qu'il
        # faut nommer ; citer le gene normal ferait chercher au mauvais endroit.
        genes = {"CYP2D6": gene(labels=["*1/*4"], phenotypes=["Poor Metabolizer"]),
                 "CYP2C19": gene(labels=["*1/*1"], phenotypes=["Normal Metabolizer"])}
        rep = rapport({"codeine": [annotation(
            "Strong", {"CYP2D6": "Poor Metabolizer", "CYP2C19": "Normal Metabolizer"},
            "Avoid codeine use.")]})
        b = L.recommandations(rep, genes, {"CYP2D6", "CYP2C19"}, set(), TR, MEDS)
        self.assertEqual(b[0][5], ("CYP2D6",))


class ParMedicament(unittest.TestCase):
    def brut(self, sev, med, fr, genes=("X",), source="CPIC, forte"):
        return (sev, med, fr, "contre-indication", source, tuple(genes),
                "123", True)

    def test_la_consigne_la_plus_severe_l_emporte(self):
        d = L.par_medicament([self.brut(2, "codeine", "Adapter"),
                              self.brut(0, "codeine", "Ne pas prescrire")])
        self.assertEqual(d["codeine"]["sev"], 0)
        self.assertEqual(d["codeine"]["fr"], ["Ne pas prescrire"])

    def test_a_severite_egale_les_textes_distincts_s_additionnent(self):
        # Deux genes peuvent concourir au meme medicament pour des raisons
        # differentes ; n'en garder qu'une retirerait une conduite du document.
        d = L.par_medicament([self.brut(0, "codeine", "Raison A", genes=("X",)),
                              self.brut(0, "codeine", "Raison B", genes=("Y",))])
        self.assertEqual(d["codeine"]["fr"], ["Raison A", "Raison B"])
        self.assertEqual(d["codeine"]["genes"], {"X", "Y"})

    def test_un_texte_identique_n_est_pas_duplique(self):
        d = L.par_medicament([self.brut(0, "codeine", "Idem", genes=("X",)),
                              self.brut(0, "codeine", "Idem", genes=("Y",))])
        self.assertEqual(d["codeine"]["fr"], ["Idem"])
        self.assertEqual(d["codeine"]["genes"], {"X", "Y"})

    def test_un_seul_texte_non_traduit_marque_le_medicament(self):
        a = self.brut(0, "codeine", "A")
        b = list(self.brut(0, "codeine", "B"))
        b[7] = False
        d = L.par_medicament([a, tuple(b)])
        self.assertFalse(d["codeine"]["ok"])

    def test_chaque_medicament_est_traite_separement(self):
        d = L.par_medicament([self.brut(0, "codeine", "A"),
                              self.brut(3, "tramadol", "B")])
        self.assertEqual(sorted(d), ["codeine", "tramadol"])


class TablesDuDepot(unittest.TestCase):
    def test_les_phenotypes_normaux_sont_dans_la_table_de_traduction(self):
        for p in L.NORMAUX:
            self.assertIn(p, L.PHENO_FR, p)

    def test_les_niveaux_de_restitution_sont_ordonnes_du_plus_grave(self):
        self.assertEqual(L.RESTITUES["contre-indication"], 0)
        self.assertLess(L.RESTITUES["contre-indication"], L.RESTITUES["éviter"])
        self.assertLess(L.RESTITUES["éviter"], L.RESTITUES["adaptation"])
        self.assertLess(L.RESTITUES["adaptation"], L.RESTITUES["vigilance"])

    def test_les_onze_aminosides_de_la_recommandation_sont_la(self):
        self.assertEqual(len(L.AMINOSIDES), 11)


class Couverture(unittest.TestCase):
    """La couverture mesuree est une condition d'emploi : elle figure sur le
    document, avec le seuil applique, ou elle n'y figure pas du tout."""

    def test_couverture_et_domaine(self):
        t = L.couverture({"couverture_mediane_retenue": 36}, 18)
        self.assertIn("36", t)
        self.assertIn("18", t)

    def test_sans_seuil_le_domaine_n_est_pas_annonce(self):
        # --couverture-min 0 : la porte est levee, et le compte rendu ne doit pas
        # pretendre qu'un domaine a ete applique.
        t = L.couverture({"couverture_mediane_retenue": 36}, 0)
        self.assertIn("36", t)
        self.assertNotIn("domaine", t)

    def test_sans_mesure_rien_n_est_affirme(self):
        self.assertEqual(L.couverture({}, 18), "")
        self.assertEqual(L.couverture({"couverture_mediane_retenue": None}, 18), "")

    def test_la_valeur_est_celle_sur_laquelle_la_porte_a_statue(self):
        # Le rendu ne refait pas le choix entre mediane clinique et globale : il
        # lit celle que l'etage a retenue, sans quoi les deux pourraient differer.
        meta = {"couverture_mediane_retenue": 18, "couverture_mediane": 40,
                "couverture_mediane_clinique": 18}
        self.assertIn("18", L.couverture(meta, 18))
        self.assertNotIn("40", L.couverture(meta, 18))

if __name__ == "__main__":
    unittest.main(verbosity=2)
