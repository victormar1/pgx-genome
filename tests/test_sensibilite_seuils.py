# -*- coding: utf-8 -*-
"""Tests du rejeu des seuils.

L'outil rejoue la regle de l'etage 2b a d'autres seuils de GQ et de profondeur,
pour mesurer ce que les seuils par defaut coutent. S'il ne reproduisait pas
exactement la regle, la mesure de sensibilite porterait sur autre chose que le
module : il appelle donc la fonction de l'etage, et le premier test verifie
qu'aux seuils par defaut il rend le perimetre que l'etage a rendu.
"""
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))

import qc_perimetre as q  # noqa: E402
import sensibilite_seuils as s  # noqa: E402


def pos(genes, gt="0/1", gq="99", lu="30", ad="15,15", filtre="PASS", dp=""):
    return s.position(genes, gt, gq, dp, ad, filtre, lu, "100")


def retenues(ligne, gq_min=20, prof_min=10, **kw):
    return s.rejoue([ligne], gq_min, prof_min, **kw)[1]["TPMT"]


# Un enregistrement par statut que l'etage sait rendre : (position, gene,
# profondeur de l'alignement, colonnes FILTER, FORMAT et echantillon).
CAS = [
    (100, "TPMT", 30, None),                                    # reference, lue
    (200, "TPMT", 4, None),                                     # non lue
    (300, "TPMT", 30, ("PASS", "GT:GQ:AD", "0/1:99:15,15")),    # mesuree
    (400, "NUDT15", 30, ("PASS", "GT:GQ:AD", "0/1:15:15,15")),  # qualite insuffisante
    (500, "NUDT15", 30, ("PASS", "GT:AD", "0/1:15,15")),        # qualite non mesurable
    (600, "DPYD", 6, ("PASS", "GT:GQ:AD", "0/1:99:3,3")),       # profondeur insuffisante
    (700, "DPYD", 30, ("LowQual", "GT:GQ:AD", "0/1:99:15,15")),  # site filtre
    (800, "UGT1A1", 30, ("PASS", "GT:GQ:AD", "0/1:99:8,7")),    # mesuree, quinze lectures
    (900, "UGT1A1", 30, ("PASS", "GT:GQ:AD", "0/1:99:5,4")),    # profondeur d appel
    (1000, "CYP2C9", 30, ("PASS", "GT:GQ:AD", "0/1:99:80,20")),  # desequilibre
    (1100, "CYP2C9", 30, ("PASS", "GT:GQ:AD", "./.:99:0,0")),   # non appelee
    (1200, "CYP2D6", 30, ("PASS", "GT:GQ:AD", "0/1:99:80,20")),  # exempte d'equilibre
    (1300, "VKORC1", 30, ("PASS", "GT:GQ:DP", "1/1:35:22")),    # mesuree, sans AD
]


class Fidelite(unittest.TestCase):
    """Le rejeu doit rendre ce que l'etage a rendu, statut par statut."""

    def setUp(self):
        self.lot = tempfile.mkdtemp(prefix="pgx_sens_")
        self.addCleanup(shutil.rmtree, self.lot, ignore_errors=True)
        self.travail = os.path.join(self.lot, "ECH", "travail")
        os.makedirs(self.travail)
        p = lambda n: os.path.join(self.travail, n)
        with open(p("positions.vcf"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for n, gene, _, _ in CAS:
                fh.write("chr1\t%d\trs%d\tC\tT\t.\tPASS\tPX=%s\n" % (n, n, gene))
        with open(p("complement.vcf"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("chr1\t5000\trs5000\tC\tT\t.\tPASS\tPX=TPMT;PXCLASSE=1;PXTYPE=snv\n")
        with open(p("ech.vcf"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("##fileformat=VCFv4.2\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tECH\n")
            for n, _, _, v in CAS:
                if v:
                    fh.write("chr1\t%d\t.\tC\tT\t.\t%s\t.\t%s\t%s\n" % ((n,) + v))
        with open(p("profondeur.txt"), "w", encoding="utf-8", newline="\n") as fh:
            for n, _, lu, _ in CAS:
                fh.write("chr1\t%d\t%d\n" % (n, lu))
            fh.write("chr1\t5000\t2\n")
        argv = ["qc_perimetre.py", "--vcf", p("ech.vcf"), "--profondeur", p("profondeur.txt"),
                "--positions", p("positions.vcf"), "--positions-complement", p("complement.vcf"),
                "--sortie", self.travail]
        with mock.patch.object(sys, "argv", argv), \
                contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(q.main(), 0)
        with open(p("perimetre.json"), encoding="utf-8") as fh:
            self.perimetre = json.load(fh)["genes"]

    def test_les_statuts_du_jeu_sont_bien_tous_la(self):
        # Sans quoi le test suivant ne prouverait rien sur les statuts absents.
        vus = set()
        for g in self.perimetre.values():
            vus |= set(g["detail"])
        self.assertEqual(vus, {
            "reference, lue", "non lue", "mesuree", "qualite insuffisante",
            "qualite non mesurable", "profondeur insuffisante", "site filtre",
            "profondeur d appel insuffisante", "desequilibre allelique", "non appelee"})

    def test_aux_seuils_par_defaut_le_rejeu_rend_le_perimetre_de_l_etage(self):
        att, ret = s.rejoue(s.charge(self.lot)["ECH"], 20, 10)
        self.assertEqual(set(att), set(self.perimetre))
        for gene, p in self.perimetre.items():
            self.assertEqual((att[gene], ret[gene]),
                             (p["positions_attendues"], p["positions_retenues"]), gene)

    def test_une_position_complementaire_ne_pese_pas_sur_le_gene(self):
        # Elle est a deux lectures : comptee, elle rendrait TPMT partiel de plus.
        att, _ = s.rejoue(s.charge(self.lot)["ECH"], 20, 10)
        self.assertEqual(att["TPMT"], 3)

    def test_baisser_un_seuil_rend_des_positions(self):
        # Le rejeu lisait le statut deja etabli aux seuils par defaut : une
        # position ecartee pour sa profondeur le restait a tout seuil.
        lignes = s.charge(self.lot)["ECH"]
        self.assertEqual(s.rejoue(lignes, 20, 10)[1]["DPYD"], 0)
        self.assertEqual(s.rejoue(lignes, 20, 5)[1]["DPYD"], 1)
        self.assertEqual(s.rejoue(lignes, 10, 10)[1]["NUDT15"], 1)

    def test_monter_la_profondeur_juge_aussi_celle_de_l_appel(self):
        # A quinze lectures alleliques, l'appel ne tient pas un seuil de vingt,
        # meme sous un alignement a trente.
        lignes = s.charge(self.lot)["ECH"]
        self.assertEqual(s.rejoue(lignes, 20, 10)[1]["UGT1A1"], 1)
        self.assertEqual(s.rejoue(lignes, 20, 20)[1]["UGT1A1"], 0)
        # Sans lectures alleliques, c'est la profondeur que declare l'appel.
        self.assertEqual(s.rejoue(lignes, 20, 22)[1]["VKORC1"], 1)
        self.assertEqual(s.rejoue(lignes, 20, 23)[1]["VKORC1"], 0)

    def test_le_filtre_se_leve_sur_demande(self):
        lignes = s.charge(self.lot)["ECH"]
        self.assertEqual(s.rejoue(lignes, 20, 10)[1]["DPYD"], 0)
        self.assertEqual(s.rejoue(lignes, 20, 10, ignorer_filtre=True)[1]["DPYD"], 1)

    def test_un_fichier_d_une_version_anterieure_est_refuse(self):
        # Sans les lectures alleliques ni le filtre, le rejeu rendrait un
        # perimetre plus large que celui de l'etage, sans le dire.
        f = os.path.join(self.travail, "qc_positions.tsv")
        with open(f, encoding="utf-8") as fh:
            lignes = [l.rstrip("\n").split("\t") for l in fh]
        i = lignes[0].index("lectures_alleliques")
        with open(f, "w", encoding="utf-8", newline="\n") as fh:
            for l in lignes:
                fh.write("\t".join(l[:i] + l[i + 1:]) + "\n")
        with self.assertRaises(ValueError):
            s.charge(self.lot)


class Rejeu(unittest.TestCase):
    def test_une_position_mesuree_qui_tient_les_deux_seuils(self):
        att, ret = s.rejoue([pos(["TPMT"])], 20, 10)
        self.assertEqual((att["TPMT"], ret["TPMT"]), (1, 1))

    def test_les_bornes_des_deux_seuils(self):
        self.assertEqual(retenues(pos(["TPMT"], lu="10")), 1)
        self.assertEqual(retenues(pos(["TPMT"], lu="9")), 0)
        self.assertEqual(retenues(pos(["TPMT"], gq="20")), 1)
        self.assertEqual(retenues(pos(["TPMT"], gq="19")), 0)

    def test_une_qualite_absente_ecarte_la_position(self):
        # Un seuil qu'on ne peut pas appliquer n'est pas un seuil satisfait :
        # c'est la regle de l'etage 2b, et le rejeu doit la suivre.
        for gq in ("", ".", None, "abc"):
            self.assertEqual(retenues(pos(["TPMT"], gq=gq)), 0, repr(gq))

    def test_une_profondeur_illisible_ecarte_la_position(self):
        for lu in ("", ".", None):
            self.assertEqual(retenues(pos(["TPMT"], lu=lu)), 0, repr(lu))

    def test_une_reference_lue_ne_depend_pas_du_gq(self):
        # Une position de reference n'a pas de genotype appele : exiger un GQ
        # ferait perdre les positions conformes a la reference, c'est-a-dire la
        # grande majorite du perimetre.
        self.assertEqual(retenues(pos(["TPMT"], gt="", gq=""), gq_min=40), 1)

    def test_une_reference_lue_depend_de_la_profondeur(self):
        self.assertEqual(retenues(pos(["TPMT"], gt="", gq="", lu="4")), 0)

    def test_une_position_partagee_compte_pour_chaque_gene(self):
        att, ret = s.rejoue([pos(["CYP2C9", "CYP2C19"])], 20, 10)
        self.assertEqual((att["CYP2C9"], att["CYP2C19"]), (1, 1))
        self.assertEqual((ret["CYP2C9"], ret["CYP2C19"]), (1, 1))

    def test_une_position_sans_gene_ne_compte_pour_personne(self):
        att, ret = s.rejoue([pos([])], 20, 10)
        self.assertEqual(sum(att.values()), 0)
        self.assertEqual(sum(ret.values()), 0)

    def test_les_attendues_ne_dependent_pas_des_seuils(self):
        # Le denominateur doit rester le meme d'un seuil a l'autre, sinon les
        # taux comparent des choses differentes.
        lignes = [pos(["TPMT"]), pos(["TPMT"], gq="5", lu="3")]
        a1, _ = s.rejoue(lignes, 20, 10)
        a2, _ = s.rejoue(lignes, 50, 30)
        self.assertEqual(a1["TPMT"], a2["TPMT"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
