# -*- coding: utf-8 -*-
"""Tests du controle qualite, etage 2b. C'est l'etage qui distingue « position
lue, identique a la reference » de « position non sequencee » : ses regles sont
celles qui empechent un faux genotype de reference de sortir.

Chaque test porte en commentaire le defaut qu'il empeche, par exemple un faux
CYP4F2*17 dans une region ou trois genes se ressemblent.
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

    def test_une_profondeur_illisible_est_une_erreur(self):
        # Le typage complementaire lit la profondeur par cette meme fonction :
        # une ligne abimee ne doit pas devenir une position « non mesuree ».
        f = tempfile.mktemp(suffix=".txt")
        with open(f, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("chr1\t100\t30\nchr1\t101\tNaN\n")
        try:
            with self.assertRaises(ValueError):
                q.lire_profondeur(f)
        finally:
            os.unlink(f)


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

    def execute(self, profondeurs, couv_min=17, melange_max=0.0, variants=()):
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
            for i, (gt, ad) in enumerate(variants):
                fh.write("chr1\t%d\trs%d\tC\tT\t.\tPASS\t.\tGT:GQ:AD\t%s:99:%s\n"
                         % (900000 + i, i, gt, ad))
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
                    "--couverture-min", str(couv_min),
                    "--melange-max", str(melange_max)]
        try:
            self.assertEqual(q.main(), 0)
            with io.open(os.path.join(d, "perimetre.json"), encoding="utf-8") as fh:
                return json.load(fh)
        finally:
            sys.argv = argv
            shutil.rmtree(d, ignore_errors=True)

    def test_un_perimetre_clinique_introuvable_arrete_l_etage(self):
        # Demande et absent, il laissait le perimetre clinique vide : la porte
        # de couverture jugeait alors tous les genes, et le compte rendu en
        # rendait vingt-trois au lieu de treize.
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        p = lambda n: os.path.join(d, n)
        with open(p("positions.vcf"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     "chr10\t100\trs1\tC\tT\t.\tPASS\tPX=CYP2C19\n")
        with open(p("ech.vcf"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tECH\n")
        with open(p("profondeur.txt"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("chr10\t100\t30\n")
        base = ["qc_perimetre.py", "--vcf", p("ech.vcf"), "--profondeur", p("profondeur.txt"),
                "--positions", p("positions.vcf"), "--sortie", d]
        for argv, code in ((base, 0), (base + ["--perimetre-clinique", p("absent.json")], 1)):
            ancien, sys.argv = sys.argv, argv
            sortie, erreur = io.StringIO(), io.StringIO()
            a_out, a_err, sys.stdout, sys.stderr = sys.stdout, sys.stderr, sortie, erreur
            try:
                self.assertEqual(q.main(), code)
            finally:
                sys.argv, sys.stdout, sys.stderr = ancien, a_out, a_err
        self.assertIn("perimetre clinique introuvable", erreur.getvalue())

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

    def test_le_verdict_de_melange_est_inscrit(self):
        # Trente heterozygotes dont vingt-cinq desequilibres : la part depasse
        # tout ce qu'un echantillon pur porte, et l'etage doit le dire.
        v = [("0/1", "20,4")] * 25 + [("0/1", "20,20")] * 5
        m = self.execute([30, 40, 30, 40], melange_max=35.0, variants=v)
        self.assertEqual(m["melange_heterozygotes"], 30)
        self.assertEqual(m["melange_desequilibrees"], 25)
        self.assertAlmostEqual(m["melange_part"], 83.3, places=1)
        self.assertTrue(m["melange_suspect"])

    def test_un_echantillon_equilibre_n_est_pas_suspect(self):
        v = [("0/1", "20,20")] * 30
        m = self.execute([30, 40, 30, 40], melange_max=35.0, variants=v)
        self.assertEqual(m["melange_desequilibrees"], 0)
        self.assertFalse(m["melange_suspect"])

    def test_le_seuil_est_trace(self):
        m = self.execute([30, 40, 30, 40], couv_min=17)
        self.assertEqual(m["couverture_mediane_min"], 17)

class IndicateurMelange(unittest.TestCase):
    """Un melange de deux individus fabrique de faux heterozygotes a faible
    part. L'indicateur compte ceux-la, et rien d'autre."""

    def lignes(self, *champs):
        return {("chr1", i): [c] for i, c in enumerate(champs)}

    def test_heterozygote_equilibre_compte_sans_etre_suspect(self):
        het, desq = q.indicateur_melange(
            self.lignes(champs("0/1", ad="20,20"), champs("0/1", ad="18,22")))
        self.assertEqual((het, desq), (2, 0))

    def test_heterozygote_desequilibre(self):
        # 4 sur 24, soit un sixieme : la forme d'un faux heterozygote.
        het, desq = q.indicateur_melange(self.lignes(champs("0/1", ad="20,4")))
        self.assertEqual((het, desq), (1, 1))

    def test_homozygote_hors_du_compte(self):
        # La part minimale ne veut rien dire pour un homozygote : l'y compter
        # diluerait l'indicateur par ce que porte le genome, pas le melange.
        het, desq = q.indicateur_melange(
            self.lignes(champs("0/0", ad="30,0"), champs("1/1", ad="0,30")))
        self.assertEqual((het, desq), (0, 0))

    def test_profondeur_insuffisante_hors_du_compte(self):
        # Sur quatre lectures, une part ne mesure rien.
        het, desq = q.indicateur_melange(self.lignes(champs("0/1", ad="3,1")))
        self.assertEqual((het, desq), (0, 0))

    def test_genotype_incomplet_hors_du_compte(self):
        # part_minimale exige deux alleles appeles : « ./. » et « 0/. » sont
        # ecartes par elle, et non par un filtre supplementaire.
        for gt in ("./.", "0/.", "./1", "."):
            het, desq = q.indicateur_melange(self.lignes(champs(gt, ad="20,4")))
            self.assertEqual((het, desq), (0, 0), gt)

    def test_lectures_illisibles_hors_du_compte(self):
        het, desq = q.indicateur_melange(self.lignes(champs("0/1")))
        self.assertEqual((het, desq), (0, 0))

    def test_la_borne_est_incluse_dans_l_equilibre(self):
        # Exactement trois dixiemes : la position n'est pas desequilibree.
        het, desq = q.indicateur_melange(self.lignes(champs("0/1", ad="21,9")))
        self.assertEqual((het, desq), (1, 0))


class SuspectMelange(unittest.TestCase):
    def test_au_dela_du_maximum(self):
        self.assertTrue(q.suspect_melange(
            {"melange_part": 41.5, "melange_heterozygotes": 53}, 35))

    def test_au_maximum_exactement(self):
        # Le seuil est une borne de ce qui passe : 35 n'est pas un refus.
        self.assertFalse(q.suspect_melange(
            {"melange_part": 35.0, "melange_heterozygotes": 53}, 35))

    def test_en_dessous(self):
        self.assertFalse(q.suspect_melange(
            {"melange_part": 25.0, "melange_heterozygotes": 44}, 35))

    def test_trop_peu_d_heterozygotes(self):
        # Deux positions sur trois font soixante-six pour cent, et ne disent
        # rien : une proportion ne s'etablit pas sur trois mesures.
        self.assertFalse(q.suspect_melange(
            {"melange_part": 66.7, "melange_heterozygotes": 3}, 35))

    def test_part_non_mesuree(self):
        self.assertFalse(q.suspect_melange({}, 35))
        self.assertFalse(q.suspect_melange(
            {"melange_part": None, "melange_heterozygotes": 50}, 35))

    def test_seuil_nul_leve_la_porte(self):
        self.assertFalse(q.suspect_melange(
            {"melange_part": 90.0, "melange_heterozygotes": 50}, 0))

    def test_les_constantes_sont_figees(self):
        # Les changer change un refus : la modification doit etre explicite.
        self.assertAlmostEqual(q.MELANGE_PART, 0.30)
        self.assertEqual(q.MELANGE_LECTURES, 10)
        self.assertEqual(q.MELANGE_HETEROZYGOTES_MIN, 20)

class Choix(unittest.TestCase):
    def test_sans_enregistrement(self):
        self.assertIsNone(q.choisir([]))

    def test_le_genotype_appele_l_emporte_sur_le_premier(self):
        vide, plein = champs("./."), champs("0/1", gq=99)
        self.assertIs(q.choisir([vide, plein]), plein)

    def test_sinon_le_premier(self):
        a, b = champs("./."), champs("./.")
        self.assertIs(q.choisir([a, b]), a)


def filtre(c, valeur):
    c = list(c)
    c[6] = valeur
    return c


class Classement(unittest.TestCase):
    """La decision par position : c'est elle qui laisse passer ou masque un
    genotype avant l'interpreteur. Chaque borne est tenue des deux cotes."""

    SNV = {"genes": ["CYP2C19"], "type": "snv"}
    D6 = {"genes": ["CYP2D6"], "type": "snv"}
    VNTR = {"genes": ["TYMS"], "type": "vntr"}
    BON = dict(gq=99, ad="15,15")

    CAS = [
        # (nom, info, champs, profondeur lue, ignorer le filtre, attendu)
        ("reference lue", SNV, None, 30, False, ("reference, lue", False)),
        ("reference au seuil", SNV, None, 10, False, ("reference, lue", False)),
        ("reference sous le seuil", SNV, None, 9, False, ("non lue", False)),
        ("reference sans profondeur", SNV, None, None, False,
         ("profondeur inconnue", False)),
        ("non appelee", SNV, champs("./.", gq=99), 30, False, ("non appelee", False)),
        ("a moitie appelee", SNV, champs("0/.", gq=99), 30, False,
         ("non appelee", False)),
        ("sans qualite", SNV, champs("0/1", ad="15,15"), 30, False,
         ("qualite non mesurable", True)),
        ("qualite sous le seuil", SNV, champs("0/1", gq=19, ad="15,15"), 30, False,
         ("qualite insuffisante", True)),
        ("qualite au seuil", SNV, champs("0/1", gq=20, ad="15,15"), 30, False,
         ("mesuree", False)),
        ("appelee sans profondeur", SNV, champs("0/1", gq=99), None, False,
         ("profondeur inconnue", True)),
        ("profondeur sous le seuil", SNV, champs("0/1", **BON), 9, False,
         ("profondeur insuffisante", True)),
        ("profondeur au seuil", SNV, champs("0/1", **BON), 10, False,
         ("mesuree", False)),
        ("site filtre", SNV, filtre(champs("0/1", **BON), "LowQual"), 30, False,
         ("site filtre", True)),
        ("site filtre, filtre ignore", SNV, filtre(champs("0/1", **BON), "LowQual"),
         30, True, ("mesuree", False)),
        ("filtre absent", SNV, filtre(champs("0/1", **BON), "."), 30, False,
         ("mesuree", False)),
        # Le filtre juge le variant du site, pas l'individu : s'il ne le porte
        # pas, il ne dit rien de lui. Un fichier de cohorte donne une ligne a
        # chacun la ou un seul varie.
        ("site filtre, individu a la reference", SNV,
         filtre(champs("0/0", gq=99, ad="30,0"), "LowQual"), 30, False,
         ("mesuree", False)),
        ("site filtre, reference phasee", SNV,
         filtre(champs("0|0", gq=99, ad="30,0"), "LowQual"), 30, False,
         ("mesuree", False)),
        ("site filtre, un seul allele alternatif suffit", SNV,
         filtre(champs("0/2", gq=99, ad="15,0,15"), "LowQual"), 30, False,
         ("site filtre", True)),
        ("site filtre, sa propre qualite reste jugee", SNV,
         filtre(champs("0/0", gq=5, ad="30,0"), "LowQual"), 30, False,
         ("qualite insuffisante", True)),
        ("site filtre, sa propre profondeur reste jugee", SNV,
         filtre(champs("0/0", gq=99, ad="6,0"), "LowQual"), 30, False,
         ("profondeur d appel insuffisante", True)),
        ("la qualite se juge avant le filtre", SNV,
         filtre(champs("0/1", gq=5, ad="15,15"), "LowQual"), 30, False,
         ("qualite insuffisante", True)),
        ("appel sur neuf lectures", SNV, champs("0/1", gq=99, ad="5,4"), 30, False,
         ("profondeur d appel insuffisante", True)),
        ("appel sur dix lectures", SNV, champs("0/1", gq=99, ad="5,5"), 30, False,
         ("mesuree", False)),
        ("sans AD, DP de l'appel a neuf", SNV, champs("0/1", gq=99, dp=9), 30, False,
         ("profondeur d appel insuffisante", True)),
        ("allele mineur a 24 %", SNV, champs("0/1", gq=99, ad="76,24"), 30, False,
         ("desequilibre allelique", True)),
        ("allele mineur a 25 %", SNV, champs("0/1", gq=99, ad="75,25"), 30, False,
         ("mesuree", False)),
        ("CYP2D6 exempte d'equilibre", D6, champs("0/1", gq=99, ad="76,24"), 30, False,
         ("mesuree", False)),
        ("homozygote sans equilibre a juger", SNV, champs("1/1", gq=99, ad="0,30"),
         30, False, ("mesuree", False)),
        ("VNTR couvert", VNTR, None, 10, False, ("couverture seule", False)),
        ("VNTR sous le seuil", VNTR, None, 9, False, ("non lue", False)),
        ("VNTR sans profondeur", VNTR, None, None, False,
         ("profondeur inconnue", False)),
        ("VNTR : le genotype n'est pas juge", VNTR, champs("0/1", gq=3, ad="2,1"),
         30, False, ("couverture seule", False)),
        # La profondeur que declare le fichier de variants ne remplace pas
        # celle de l'alignement : un DP confortable sur une position que
        # l'alignement n'a pas mesuree ferait passer une position non verifiee.
        ("le DP du fichier ne vaut pas mesure", SNV,
         champs("0/1", gq=99, dp=40, ad="20,20"), None, False,
         ("profondeur inconnue", True)),
        ("VNTR sous le seuil, genotype present", VNTR, champs("0/1", gq=99, ad="15,15"),
         9, False, ("non lue", False)),
    ]

    def test_chaque_statut_et_chaque_borne(self):
        for nom, info, c, lu, ignorer, attendu in self.CAS:
            with self.subTest(nom):
                self.assertEqual(q.classer(info, c, lu, 20, 10, ignorer), attendu[0])
                statut, porteur, masquer = q.juger(info, [c] if c else [], lu, 20, 10, ignorer)
                self.assertEqual((statut, masquer), attendu)
                self.assertIs(porteur, c)

    def test_un_gene_haploide_n_a_pas_d_equilibre_a_tenir(self):
        # Sur l'ADN mitochondrial, un allele a 10 % est une heteroplasmie.
        c = champs("0/1", gq=99, ad="90,10")
        self.assertEqual(q.classer(self.SNV, c, 100, 20, 10), "desequilibre allelique")
        self.assertEqual(q.classer(self.SNV, c, 100, 20, 10, equilibre=False), "mesuree")

    def test_les_statuts_retenus(self):
        self.assertEqual(q.RETENUS, ("mesuree", "reference, lue", "couverture seule"))


class PlusieursEnregistrements(unittest.TestCase):
    """L'interpreteur lit tous les enregistrements d'une position : un seul
    douteux suffit a fausser son appel, et c'est donc le pire qui decide."""

    SNV = {"genes": ["CYP2C19"], "type": "snv"}
    BON = champs("0/1", gq=99, ad="15,15")

    def juge(self, *groupe, **options):
        return q.juger(self.SNV, list(groupe), 30, 20, 10, **options)

    def test_un_enregistrement_sans_genotype_ne_cache_pas_celui_qui_en_porte(self):
        # Huit colonnes : un site sans echantillon, qui se lit comme une
        # reference. Le statut de la position est celui du genotype appele.
        self.assertEqual(self.juge(self.BON[:8], self.BON), ("mesuree", self.BON, False))

    def test_deux_bons_enregistrements(self):
        autre = champs("0/1", gq=80, ad="12,14", alt="G")
        self.assertEqual(self.juge(self.BON, autre), ("mesuree", self.BON, False))

    def test_un_douteux_ecarte_la_position_quel_que_soit_son_rang(self):
        douteux = champs("0/1", gq=5, ad="15,15", alt="G")
        for groupe in ((self.BON, douteux), (douteux, self.BON)):
            statut, porteur, masquer = self.juge(*groupe)
            self.assertEqual(statut, "qualite insuffisante")
            self.assertIs(porteur, douteux)
            self.assertTrue(masquer)

    def test_un_non_appele_ecarte_la_position_et_masque_l_autre(self):
        # Le genotype appele ne doit pas atteindre l'interpreteur seul : il lui
        # ferait rendre un appel sur la moitie de ce que porte la position.
        vide = champs("./.", alt="G")
        statut, porteur, masquer = self.juge(self.BON, vide)
        self.assertEqual(statut, "non appelee")
        self.assertIs(porteur, vide)
        self.assertTrue(masquer)

    def test_rien_a_masquer_quand_rien_n_est_appele(self):
        self.assertEqual(self.juge(champs("./."), champs("./."))[2], False)

    def test_le_premier_ecarte_donne_le_statut(self):
        a = champs("0/1", gq=5, ad="15,15")
        b = champs("0/1", gq=99, ad="3,3", alt="G")
        self.assertEqual(self.juge(a, b)[0], "qualite insuffisante")
        self.assertEqual(self.juge(b, a)[0], "profondeur d appel insuffisante")

    def test_le_filtre_ignore_vaut_pour_chaque_enregistrement(self):
        filtre_ = filtre(champs("0/1", gq=99, ad="15,15", alt="G"), "LowQual")
        self.assertEqual(self.juge(self.BON, filtre_)[0], "site filtre")
        self.assertEqual(self.juge(self.BON, filtre_, ignorer_filtre=True)[0], "mesuree")


class Verdict(unittest.TestCase):
    """La ligne que l'orchestration relit pour refuser un genome."""

    META = {"couverture_mediane": 40, "couverture_mediane_clinique": 12,
            "couverture_hors_domaine": True, "melange_part": 41.5,
            "melange_suspect": True, "positions_retenues": 7, "positions_visees": 9,
            "positions_sans_GQ": 2,
            "genes": {"A": {"statut": "complet"}, "B": {"statut": "partiel"},
                      "C": {"statut": "complet"}}}

    def test_les_neuf_champs_dans_l_ordre(self):
        self.assertEqual(q.verdict(self.META),
                         ["12", "1", "41.5", "1", "7", "9", "2", "3", "2"])

    def test_une_mesure_absente_ne_laisse_aucun_champ_vide(self):
        # Un champ vide decalerait les suivants : un melange suspect serait lu
        # comme une couverture.
        meta = dict(self.META, couverture_mediane=None,
                    couverture_mediane_clinique=None, melange_part=None,
                    couverture_hors_domaine=False, melange_suspect=False)
        v = q.verdict(meta)
        self.assertEqual(v[:4], ["-", "0", "-", "0"])
        self.assertEqual(len(v), 9)
        self.assertNotIn("", v)


class VcfQualifie(unittest.TestCase):
    """Ce que l'interpreteur recoit : une position ecartee ne lui laisse aucun
    genotype appele, sur aucun de ses enregistrements."""

    ENTETE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

    def execute(self, enregistrements, profondeurs, complement=()):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        p = lambda n: os.path.join(d, n)
        with open(p("complement.vcf"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(self.ENTETE + "\n")
            for pos in complement:
                fh.write("chr1\t%d\trs%d\tC\tT\t.\tPASS\tPX=TPMT;PXCLASSE=1;PXTYPE=snv\n"
                         % (pos, pos))
        with open(p("positions.vcf"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(self.ENTETE + "\n")
            for pos in sorted(profondeurs):
                fh.write("chr1\t%d\trs%d\tC\tT\t.\tPASS\tPX=TPMT\n" % (pos, pos))
        with open(p("ech.vcf"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("##fileformat=VCFv4.2\n" + self.ENTETE + "\tFORMAT\tECH\n")
            for pos, alt, echantillon in enregistrements:
                fh.write("chr1\t%d\t.\tC\t%s\t.\tPASS\t.\tGT:GQ:DP:AD\t%s\n"
                         % (pos, alt, echantillon))
        with open(p("profondeur.txt"), "w", encoding="utf-8", newline="\n") as fh:
            for pos, lu in profondeurs.items():
                if lu is not None:
                    fh.write("chr1\t%d\t%d\n" % (pos, lu))
        argv = ["qc_perimetre.py", "--vcf", p("ech.vcf"), "--profondeur", p("profondeur.txt"),
                "--positions", p("positions.vcf"), "--sortie", d]
        if complement:
            argv += ["--positions-complement", p("complement.vcf")]
        sortie = io.StringIO()
        ancien, sys.argv, sys.stdout = (sys.argv, sys.stdout), argv, sortie
        try:
            self.assertEqual(q.main(), 0)
        finally:
            sys.argv, sys.stdout = ancien
        with open(p("qualifie.vcf"), encoding="utf-8") as fh:
            gts = [(int(c[1]), c[4], c[9]) for c in
                   (l.rstrip("\n").split("\t") for l in fh if not l.startswith("#"))]
        with open(p("qc_positions.tsv"), encoding="utf-8") as fh:
            lignes = [l.rstrip("\n").split("\t") for l in fh]
        tsv = {int(r[lignes[0].index("pos")]): dict(zip(lignes[0], r)) for r in lignes[1:]}
        return gts, tsv

    def test_un_enregistrement_douteux_masque_toute_la_position(self):
        gts, tsv = self.execute(
            [(100, "T", "0/1:99:30:15,15"), (100, "G", "0/1:5:30:15,15"),
             (200, "T", "0/1:99:30:15,15"), (200, "G", "0/1:80:30:14,16")],
            {100: 30, 200: 30})
        self.assertEqual(gts, [(100, "T", "./."), (100, "G", "./."),
                               (200, "T", "0/1"), (200, "G", "0/1")])
        self.assertEqual(tsv[100]["statut"], "qualite insuffisante")
        self.assertEqual(tsv[100]["retenue"], "non")
        self.assertEqual(tsv[100]["enregistrements"], "2")
        # la ligne porte l'enregistrement qui a decide, de quoi rejouer la decision
        self.assertEqual((tsv[100]["GQ"], tsv[100]["filtre"], tsv[100]["lectures_alleliques"]),
                         ("5", "PASS", "15,15"))
        self.assertEqual(tsv[200]["retenue"], "oui")

    def test_une_position_complementaire_ecartee_n_est_pas_masquee(self):
        # L'interpreteur ne s'en sert pas : la masquer ne changerait pas son
        # appel, mais modifierait son entree. Elle est ecartee dans le suivi,
        # et laissee intacte dans le fichier.
        gts, tsv = self.execute([(100, "T", "0/1:99:30:15,15"), (300, "T", "0/1:5:30:15,15")],
                                {100: 30}, complement=(300,))
        self.assertEqual(gts, [(100, "T", "0/1"), (300, "T", "0/1")])
        self.assertEqual((tsv[300]["statut"], tsv[300]["retenue"], tsv[300]["source"]),
                         ("qualite insuffisante", "non", "RNPGx"))

    def test_la_profondeur_du_fichier_ne_sauve_pas_une_position_non_mesuree(self):
        # La position 100 manque a la mesure de l'alignement ; le fichier de
        # variants y declare quarante lectures.
        gts, tsv = self.execute([(100, "T", "0/1:99:40:20,20")], {100: None, 200: 30})
        self.assertEqual(gts, [(100, "T", "./.")])
        self.assertEqual(tsv[100]["statut"], "profondeur inconnue")
        self.assertEqual(tsv[100]["profondeur_VCF"], "40")


if __name__ == "__main__":
    unittest.main(verbosity=2)
