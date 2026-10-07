# -*- coding: utf-8 -*-
"""Tests du controle qualite, etage 2b. C'est l'etage qui distingue « position
lue, identique a la reference » de « position non sequencee » : ses regles sont
celles qui empechent un faux genotype de reference de sortir.

Chaque test porte en commentaire le defaut qu'il empeche de revenir ; plusieurs
viennent d'un audit position par position qui avait releve cent quarante et un
ecarts, dont de faux CYP4F2*17 dans une region ou trois genes se ressemblent.
"""
import gzip
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))

import qc_perimetre as q  # noqa: E402


def champs(gt, gq=None, dp=None, ad=None, ref="C", alt="T"):
    """Une ligne de VCF a dix colonnes, FORMAT et echantillon compris."""
    cles, vals = ["GT"], [gt]
    for nom, v in (("GQ", gq), ("DP", dp), ("AD", ad)):
        if v is not None:
            cles.append(nom)
            vals.append(v)
    return ["chr1", "100", "rs1", ref, alt, ".", "PASS", ".",
            ":".join(cles), ":".join(str(x) for x in vals)]


class Entier(unittest.TestCase):
    def test_valeurs_lisibles(self):
        self.assertEqual(q.entier("12"), 12)

    def test_valeurs_illisibles_rendent_none(self):
        # Un point, une chaine vide ou None sont les trois formes sous
        # lesquelles un VCF annonce une valeur absente.
        for x in (".", "", None, "abc", "1.5"):
            self.assertIsNone(q.entier(x), repr(x))


class Genotype(unittest.TestCase):
    def test_lit_gt_gq_dp(self):
        self.assertEqual(q.genotype(champs("0/1", gq=99, dp=30)),
                         ("0/1", 99, 30))

    def test_gq_absent_rend_none_et_non_zero(self):
        # Rendre zero ferait echouer le seuil pour la mauvaise raison, et
        # rendre une valeur par defaut le ferait passer a tort.
        gt, gq, dp = q.genotype(champs("0/1", dp=30))
        self.assertEqual(gt, "0/1")
        self.assertIsNone(gq)

    def test_ligne_sans_echantillon(self):
        self.assertEqual(q.genotype(["chr1", "100", "rs1", "C", "T"]),
                         (None, None, None))


class LecturesAlleliques(unittest.TestCase):
    def test_lit_ad(self):
        self.assertEqual(q.lectures_alleliques(champs("0/1", ad="12,18")),
                         [12, 18])

    def test_ad_absent_ou_illisible(self):
        self.assertIsNone(q.lectures_alleliques(champs("0/1")))
        self.assertIsNone(q.lectures_alleliques(champs("0/1", ad="12,x")))
        self.assertIsNone(q.lectures_alleliques(None))


class PartMinimale(unittest.TestCase):
    def test_heterozygote_equilibre(self):
        self.assertAlmostEqual(q.part_minimale("0/1", [15, 15]), 0.5)

    def test_heterozygote_desequilibre(self):
        # 0,10 a 0,23 est la plage ou se trouvaient les faux CYP4F2*17 de la
        # region paralogue : c'est cette mesure qui les ecarte.
        self.assertAlmostEqual(q.part_minimale("0/1", [27, 3]), 0.1)

    def test_homozygote_n_a_pas_de_part_minimale(self):
        # La regle ne s'applique pas ; elle ne doit pas se deviner.
        self.assertIsNone(q.part_minimale("1/1", [0, 30]))
        self.assertIsNone(q.part_minimale("0/0", [30, 0]))

    def test_phase_traitee_comme_non_phase(self):
        self.assertAlmostEqual(q.part_minimale("0|1", [15, 15]), 0.5)

    def test_sans_lectures_ou_sans_genotype(self):
        self.assertIsNone(q.part_minimale("0/1", None))
        self.assertIsNone(q.part_minimale(None, [15, 15]))
        self.assertIsNone(q.part_minimale("0/1", [0, 0]))

    def test_allele_hors_du_tableau_de_lectures(self):
        # Un AD plus court que les alleles appeles ne doit pas lever
        # d'exception ni fabriquer une part.
        self.assertIsNone(q.part_minimale("1/2", [10, 10]))

    def test_multiallelique(self):
        self.assertAlmostEqual(q.part_minimale("1/2", [0, 8, 12]), 0.4)


class Appele(unittest.TestCase):
    def test_genotypes_appeles(self):
        for gt in ("0/0", "0/1", "1/1", "1|2"):
            self.assertTrue(q.appele(gt), gt)

    def test_genotypes_non_appeles(self):
        # Un demi-appel « ./1 » n'est pas un appel : le prendre pour tel
        # fabriquerait un homozygote ou un heterozygote imaginaire.
        for gt in (".", "./.", ".|.", "./1", "1/.", "", None):
            self.assertFalse(q.appele(gt), repr(gt))


class LirePositions(unittest.TestCase):
    ENTETE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

    def ecris(self, lignes):
        f = tempfile.mktemp(suffix=".vcf")
        with open(f, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(self.ENTETE)
            fh.writelines(lignes)
        return f

    def test_lit_les_annotations_du_complement(self):
        f = self.ecris(["chr3\t165773492\trs1803274\tC\tT\t.\tPASS\t"
                        "PX=BCHE;PXCLASSE=1;PXTYPE=snv;PXINTERP=aucun;"
                        "PXNOTE=p.Ala539Thr,_variante_K\n"])
        try:
            d = q.lire_positions(f, complement=True)
        finally:
            os.unlink(f)
        e = d[("chr3", 165773492)]
        self.assertEqual(e["genes"], ["BCHE"])
        self.assertEqual(e["classe"], 1)
        self.assertEqual(e["type"], "snv")
        self.assertEqual(e["interpreteur"], "aucun")
        self.assertIn("variante K", e["note"])      # les tirets bas redeviennent des espaces
        self.assertTrue(e["complement"])

    def test_defauts_quand_les_annotations_manquent(self):
        f = self.ecris(["chr1\t100\trs1\tC\tT\t.\tPASS\tPX=GENE\n"])
        try:
            d = q.lire_positions(f)
        finally:
            os.unlink(f)
        e = d[("chr1", 100)]
        self.assertEqual(e["type"], "snv")
        self.assertEqual(e["interpreteur"], "pharmcat")
        self.assertIsNone(e["classe"])
        self.assertFalse(e["complement"])

    def test_fin_deduite_de_la_longueur_de_la_reference(self):
        # Un indel couvre plusieurs bases : sa fin ne se confond pas avec sa
        # position de depart, sans quoi la couverture serait jugee trop court.
        f = self.ecris(["chr1\t100\trs1\tCTTT\tC\t.\tPASS\tPX=GENE\n",
                        "chr18\t657646\trs2\tC\t<VNTR>\t.\tPASS\t"
                        "PX=TYMS;PXTYPE=vntr;END=657712\n"])
        try:
            d = q.lire_positions(f, complement=True)
        finally:
            os.unlink(f)
        self.assertEqual(d[("chr1", 100)]["fin"], 103)
        self.assertEqual(d[("chr18", 657646)]["fin"], 657712)   # END fait foi


class LireVcf(unittest.TestCase):
    CORPS = ("##fileformat=VCFv4.2\n"
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tECH\n"
             "chr1\t100\trs1\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
             "chr1\t100\trs1\tC\tG\t.\tPASS\t.\tGT\t0/1\n"
             "chr1\t200\trs2\tA\tG\t.\tPASS\t.\tGT\t1/1\n")

    def test_groupe_les_lignes_par_position(self):
        # Deux enregistrements a la meme position arrivent d'un VCF eclate par
        # allele : les perdre reviendrait a oublier un allele.
        f = tempfile.mktemp(suffix=".vcf")
        with open(f, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(self.CORPS)
        try:
            entete, lignes, n = q.lire_vcf(f)
        finally:
            os.unlink(f)
        self.assertEqual(n, 1)
        self.assertEqual(len(lignes[("chr1", 100)]), 2)
        self.assertEqual(len(lignes[("chr1", 200)]), 1)
        self.assertTrue(any(l.startswith("#CHROM") for l in entete))

    def test_lit_un_fichier_comprime(self):
        f = tempfile.mktemp(suffix=".vcf.gz")
        with gzip.open(f, "wt", encoding="utf-8") as fh:
            fh.write(self.CORPS)
        try:
            _, lignes, n = q.lire_vcf(f)
        finally:
            os.unlink(f)
        self.assertEqual(n, 1)
        self.assertEqual(len(lignes), 2)

    def test_compte_plusieurs_echantillons(self):
        # Le module refuse un VCF multi-echantillon ; encore faut-il les compter.
        f = tempfile.mktemp(suffix=".vcf")
        with open(f, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
                     "\tA\tB\n")
        try:
            _, _, n = q.lire_vcf(f)
        finally:
            os.unlink(f)
        self.assertEqual(n, 2)


class LireProfondeur(unittest.TestCase):
    def test_lit_samtools_depth(self):
        f = tempfile.mktemp(suffix=".txt")
        with open(f, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("chr1\t100\t30\nchr1\t101\t0\n")
        try:
            p = q.lire_profondeur(f)
        finally:
            os.unlink(f)
        self.assertEqual(p[("chr1", 100)], 30)
        self.assertEqual(p[("chr1", 101)], 0)

    def test_fichier_absent_rend_un_dictionnaire_vide(self):
        self.assertEqual(q.lire_profondeur(None), {})
        self.assertEqual(q.lire_profondeur("/absent/xyz.txt"), {})


class Seuils(unittest.TestCase):
    def test_les_seuils_du_module_sont_ceux_documentes(self):
        # Changer ces valeurs change les resultats et doit etre remesure ; le
        # test est la pour que la modification soit explicite.
        self.assertAlmostEqual(q.EQUILIBRE_MIN, 0.25)
        self.assertEqual(q.SANS_EQUILIBRE, {"CYP2D6"})


class Mediane(unittest.TestCase):
    def test_nombre_impair_de_valeurs(self):
        self.assertEqual(q.mediane([30, 10, 20]), 20)

    def test_nombre_pair_de_valeurs(self):
        # La moyenne des deux valeurs centrales, et un entier quand elle en est
        # un : une couverture de 24,0 fois se lit mal sur un compte rendu.
        self.assertEqual(q.mediane([10, 20, 30, 40]), 25)
        self.assertNotIsInstance(q.mediane([10, 20, 30, 40]), float)
        self.assertEqual(q.mediane([10, 11]), 10.5)

    def test_sans_valeur(self):
        # Rendre zero ferait refuser le genome pour une mesure qui n'existe pas.
        self.assertIsNone(q.mediane([]))

    def test_une_seule_valeur(self):
        self.assertEqual(q.mediane([17]), 17)


class HorsDomaine(unittest.TestCase):
    """Le verdict qui refuse un genome. Sous le seuil mesure, un allele variant
    peut etre lu comme reference et l'interpreteur l'affirme sans reserve."""

    def test_sous_le_seuil(self):
        self.assertTrue(q.hors_domaine({"couverture_mediane_clinique": 11}, 17))

    def test_au_seuil_exactement(self):
        # Le domaine inclut sa borne : elle a ete mesuree, et elle passe.
        self.assertFalse(q.hors_domaine({"couverture_mediane_clinique": 17}, 17))

    def test_au_dessus_du_seuil(self):
        self.assertFalse(q.hors_domaine({"couverture_mediane_clinique": 34}, 17))

    def test_seuil_nul_leve_la_porte(self):
        self.assertFalse(q.hors_domaine({"couverture_mediane_clinique": 2}, 0))

    def test_couverture_inconnue_n_est_pas_un_refus(self):
        # L'etage a echoue : c'est un echec d'etage, pas un genome hors domaine.
        self.assertFalse(q.hors_domaine({}, 17))
        self.assertFalse(q.hors_domaine({"couverture_mediane_clinique": None,
                                         "couverture_mediane": None}, 17))

    def test_le_perimetre_clinique_decide_quand_il_existe(self):
        # Les genes rendus sans verite ne doivent pas porter la decision : ici
        # la mediane globale passerait, la clinique non.
        meta = {"couverture_mediane": 40, "couverture_mediane_clinique": 11}
        self.assertTrue(q.hors_domaine(meta, 17))
        self.assertEqual(q.couverture_retenue(meta), 11)

    def test_sans_perimetre_clinique_la_mediane_globale_decide(self):
        meta = {"couverture_mediane": 11, "couverture_mediane_clinique": None}
        self.assertTrue(q.hors_domaine(meta, 17))
        self.assertEqual(q.couverture_retenue(meta), 11)


class CouvertureMesuree(unittest.TestCase):
    """La mediane sort de l'etage complet : ce que mesurent les tests ci-dessus
    ne vaut que si l'etage leur donne les bonnes positions."""

    POSITIONS = [("chr10", 100, "CYP2C19"), ("chr10", 200, "CYP2C19"),
                 ("chr7", 300, "POR"), ("chr7", 400, "POR")]

    def execute(self, profondeurs, couv_min=17):
        d = tempfile.mkdtemp()
        pos = os.path.join(d, "positions.vcf")
        with open(pos, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for i, (c, p, g) in enumerate(self.POSITIONS):
                fh.write("%s\t%d\trs%d\tC\tT\t.\tPASS\tPX=%s\n" % (c, p, i, g))
        vcf = os.path.join(d, "ech.vcf")
        with open(vcf, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("##fileformat=VCFv4.2\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tECH\n")
        prof = os.path.join(d, "profondeur.txt")
        with open(prof, "w", encoding="utf-8", newline="\n") as fh:
            for (c, p, _), v in zip(self.POSITIONS, profondeurs):
                if v is not None:
                    fh.write("%s\t%d\t%d\n" % (c, p, v))
        clin = os.path.join(d, "clinique.json")
        with open(clin, "w", encoding="utf-8", newline="\n") as fh:
            json.dump({"genes": {"CYP2C19": {}}}, fh)
        argv = sys.argv
        sys.argv = ["qc_perimetre.py", "--vcf", vcf, "--profondeur", prof,
                    "--positions", pos, "--sortie", d, "--echantillon", "ECH",
                    "--perimetre-clinique", clin,
                    "--couverture-min", str(couv_min)]
        try:
            self.assertEqual(q.main(), 0)
            with io.open(os.path.join(d, "perimetre.json"), encoding="utf-8") as fh:
                return json.load(fh)
        finally:
            sys.argv = argv
            shutil.rmtree(d, ignore_errors=True)

    def test_la_mediane_clinique_ignore_les_autres_genes(self):
        m = self.execute([30, 40, 4, 6])
        self.assertEqual(m["couverture_mediane_clinique"], 35)
        self.assertEqual(m["couverture_mediane"], 18)
        # c'est la valeur clinique qui est retenue, et c'est elle que le compte
        # rendu affichera
        self.assertEqual(m["couverture_mediane_retenue"], 35)
        self.assertFalse(m["couverture_hors_domaine"])

    def test_le_perimetre_clinique_fait_refuser_seul(self):
        # La mediane globale passerait le seuil ; celle du perimetre non.
        m = self.execute([10, 12, 90, 90])
        self.assertEqual(m["couverture_mediane_clinique"], 11)
        self.assertTrue(m["couverture_hors_domaine"])

    def test_une_position_absente_compte_pour_zero(self):
        # samtools depth -a rend toute position du BED : une absence est une
        # position hors de l'alignement, et la passer sous silence ferait
        # paraitre bien couvert un alignement tronque.
        m = self.execute([None, 40, 90, 90])
        self.assertEqual(m["couverture_mediane_clinique"], 20)

    def test_le_seuil_est_trace(self):
        m = self.execute([30, 40, 30, 40], couv_min=17)
        self.assertEqual(m["couverture_mediane_min"], 17)

if __name__ == "__main__":
    unittest.main(verbosity=2)
