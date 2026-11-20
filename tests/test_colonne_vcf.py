# -*- coding: utf-8 -*-
"""La colonne d'un individu, extraite d'un fichier a plusieurs echantillons.

Une fusion de fichiers individuels ajoute a chaque individu ce qu'elle ne sait
pas de lui : une ligne vide a chaque site ou un autre varie, et, sur ses
propres lignes, les alleles des autres sans lectures. Lus tels quels, les
premieres font ecarter des positions qu'il a lues, les seconds font echouer le
typeur. Un appel conjoint, lui, a mesure l'individu : ses lignes ne changent pas.
"""
import io
import os
import sys
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RACINE, "bin"))
import colonne_vcf as C  # noqa: E402

NOMBRES = {"GT": "1", "AD": "R", "AF": "A", "DP": "1", "GQ": "1", "PL": "G",
           "F1R2": "R", "SB": "4", "PS": "1"}


def ligne(fmt, valeur, alt="G", info="AC=2;AN=8"):
    return "\t".join(["chr1", "100", ".", "A", alt, "50", "PASS", info, fmt, valeur]) + "\n"


def champs(fmt, valeur, alt="G"):
    return ligne(fmt, valeur, alt).rstrip("\n").split("\t")


ENTETE = ("##fileformat=VCFv4.2\n"
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
          '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">\n'
          '##FORMAT=<ID=AF,Number=A,Type=Float,Description="f">\n'
          '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="p">\n'
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tE\n")


class SansValeur(unittest.TestCase):
    def test_une_ligne_ajoutee_par_une_fusion_ne_porte_rien(self):
        self.assertTrue(C.sans_valeur("GT:AD:DP:GQ", "./.:.:.:."))
        self.assertTrue(C.sans_valeur("GT:AD:DP:GQ", ".|.:.:.:."))
        self.assertTrue(C.sans_valeur("GT:AD:DP:GQ", ".:.:.:."))

    def test_des_champs_omis_en_fin_de_colonne_sont_des_champs_vides(self):
        # La norme permet de ne pas ecrire les derniers champs manquants.
        self.assertTrue(C.sans_valeur("GT:AD:DP:GQ", "./."))

    def test_une_reference_affirmee_sans_aucune_valeur_ne_porte_rien(self):
        # Ce que laisse une fusion qui ramene les absents a la reference.
        self.assertTrue(C.sans_valeur("GT:AD:DP:GQ", "0/0:.:.:."))
        self.assertTrue(C.sans_valeur("GT", "0/0"))

    def test_un_non_appele_qui_porte_des_lectures_est_une_mesure(self):
        # Un appel conjoint : l'appelant a lu l'individu et n'a pas conclu.
        self.assertFalse(C.sans_valeur("GT:AD:DP:GQ:PL", "./.:12,0:12:.:."))
        self.assertFalse(C.sans_valeur("GT:AD:DP:GQ", "./.:.:0:."))

    def test_une_reference_mesuree_est_gardee(self):
        # Le controle qualite la juge, comme dans un fichier individuel.
        self.assertFalse(C.sans_valeur("GT:AD:DP:GQ", "0/0:29,0:29:81"))
        self.assertFalse(C.sans_valeur("GT:AD:DP:GQ", "0/0:.:.:12"))

    def test_un_variant_sans_valeur_est_garde_pour_etre_ecarte_plus_loin(self):
        # Un allele alternatif affirme sans qualite n'est pas une absence : le
        # controle qualite doit le voir, et ecarter la position.
        self.assertFalse(C.sans_valeur("GT:AD:DP:GQ", "0/1:.:.:."))
        self.assertFalse(C.sans_valeur("GT:AD:DP:GQ", "./1:.:.:."))
        self.assertFalse(C.sans_valeur("GT:AD:DP:GQ", "1/.:.:.:."))

    def test_sans_genotype_en_tete_rien_n_est_presume(self):
        self.assertFalse(C.sans_valeur("AD:DP", ".:."))


class AllelesSansLectures(unittest.TestCase):
    def alleger(self, fmt, valeur, alt, nombres=NOMBRES):
        c = champs(fmt, valeur, alt)
        change = C.retirer_alleles_sans_lectures(c, nombres)
        return change, c

    def test_l_allele_d_un_autre_sort_de_tous_les_champs_par_allele(self):
        # Une ligne reelle de fusion : l'individu porte le premier allele, le
        # second vient d'un autre fichier.
        change, c = self.alleger(
            "GT:AD:AF:DP:F1R2:GQ:PL:SB",
            "1/1:0,5,.:1,.:5:0,0,.:36:100,39,0,.,.,.:0,0,0,5", "TCC,TCCC")
        self.assertTrue(change)
        self.assertEqual(c[4], "TCC")
        self.assertEqual(c[9], "1/1:0,5:1:5:0,0:36:100,39,0:0,0,0,5")
        self.assertEqual(c[7], ".")

    def test_les_rangs_du_genotype_suivent_les_alleles_gardes(self):
        # Il porte le troisieme allele : retirer les deux premiers le renumerote.
        change, c = self.alleger("GT:AD:AF:PL", "0|3:9,.,.,11:.,.,0.55:"
                                 "40,.,.,.,.,.,0,.,.,60", "C,G,T")
        self.assertTrue(change)
        self.assertEqual(c[4], "T")
        self.assertEqual(c[9], "0|1:9,11:0.55:40,0,60")

    def test_un_allele_non_porte_mais_compte_reste(self):
        # Un appel conjoint a mesure cet allele chez l'individu : zero lecture
        # est une mesure, la ligne ne change pas.
        c = champs("GT:AD:DP:GQ", "0/1:10,12,0:22:99", "C,G")
        avant = list(c)
        self.assertFalse(C.retirer_alleles_sans_lectures(c, NOMBRES))
        self.assertEqual(c, avant)

    def test_seuls_les_alleles_sans_lectures_sortent(self):
        change, c = self.alleger("GT:AD", "0/1:10,12,0,.", "C,G,T")
        self.assertTrue(change)
        self.assertEqual((c[4], c[9]), ("C,G", "0/1:10,12,0"))

    def test_les_lectures_se_lisent_par_allele_quoi_que_dise_l_en_tete(self):
        # D'anciens appelants declarent AD de taille libre.
        change, c = self.alleger("GT:AD", "1/1:0,7,.", "C,G", dict(NOMBRES, AD="."))
        self.assertTrue(change)
        self.assertEqual((c[4], c[9]), ("C", "1/1:0,7"))

    def test_un_champ_de_taille_fixe_ou_inconnue_ne_change_pas(self):
        change, c = self.alleger("GT:AD:SB:XX", "1/1:0,7,.:1,2,3,4:5,6,7", "C,G")
        self.assertTrue(change)
        self.assertEqual(c[9], "1/1:0,7:1,2,3,4:5,6,7")

    def test_un_genotype_haploide_se_reduit_par_allele(self):
        change, c = self.alleger("GT:AD:PL", "1:0,30,.:90,0,.", "C,G")
        self.assertTrue(change)
        self.assertEqual((c[4], c[9]), ("C", "1:0,30:90,0"))

    def test_une_reference_dont_tous_les_alleles_sont_d_autres_perd_son_alternatif(self):
        change, c = self.alleger("GT:AD:DP", "0/0:30,.:30", "C")
        self.assertTrue(change)
        self.assertEqual((c[4], c[9]), (".", "0/0:30:30"))

    def test_un_allele_porte_reste_meme_sans_lectures(self):
        # Il le porte : le retirer changerait son genotype. Les lectures
        # absentes sont l'affaire du controle qualite.
        change, c = self.alleger("GT:AD", "0/2:9,.,.", "C,G")
        self.assertTrue(change)
        self.assertEqual((c[4], c[9]), ("G", "0/1:9,."))

    def test_trois_alleles_gardes_gardent_l_ordre_des_genotypes(self):
        # Rangs de la norme : (0,0) (0,1) (1,1) (0,2) (1,2) (2,2), puis ceux
        # du quatrieme allele, qui sortent.
        change, c = self.alleger("GT:AD:PL", "1/2:0,10,12,.:90,91,92,93,94,95,.,.,.,.", "C,G,T")
        self.assertTrue(change)
        self.assertEqual((c[4], c[9]), ("C,G", "1/2:0,10,12:90,91,92,93,94,95"))

    def test_un_champ_par_allele_d_une_autre_taille_n_est_pas_touche(self):
        # Declare par allele mais ecrit autrement : rien ne dit quelle valeur
        # va avec quel allele.
        change, c = self.alleger("GT:AD:F1R2:AF", "1/1:0,7,.:1,2,3,4:1,.,.", "C,G")
        self.assertTrue(change)
        self.assertEqual(c[9], "1/1:0,7:1,2,3,4:1,.,.")

    def test_sans_lectures_par_allele_rien_ne_se_decide(self):
        for fmt, valeur in (("GT:DP", "0/1:20"), ("GT:AD", "0/1"), ("AD:GT", "9,9,.:0/1")):
            c = champs(fmt, valeur, "C,G")
            self.assertFalse(C.retirer_alleles_sans_lectures(c, NOMBRES), fmt)

    def test_des_lectures_qui_ne_suivent_pas_les_alleles_ne_sont_pas_touchees(self):
        c = champs("GT:AD", "0/1:10,.", "C,G")
        avant = list(c)
        self.assertFalse(C.retirer_alleles_sans_lectures(c, NOMBRES))
        self.assertEqual(c, avant)


class Filtrer(unittest.TestCase):
    def passer(self, texte):
        sortie = io.StringIO()
        n = C.filtrer(io.StringIO(texte), sortie)
        return n, sortie.getvalue()

    def test_seul_ce_que_la_fusion_a_ajoute_change(self):
        vide = ligne("GT:AD:AF:PL", "./.:.:.:.", alt="C")
        conjoint = ligne("GT:AD:PL", "0/1:10,12,0:99,0,99,120,130,200", alt="C,G")
        texte = (ENTETE
                 + ligne("GT:AD:AF:PL", "0/1:10,12:0.55:99,0,99")
                 + vide
                 + ligne("GT:AD:PL", "0/0:30,0:0,90,900")
                 + conjoint
                 + ligne("GT:AD:AF:PL", "1/1:0,7,.:1,.:70,21,0,.,.,.", alt="C,G")
                 + ligne("GT:AD:PL", "./.:12,0:.", alt="C"))
        (retirees, allegees), rendu = self.passer(texte)
        self.assertEqual((retirees, allegees), (1, 1))
        attendu = texte.replace(vide, "").replace(
            ligne("GT:AD:AF:PL", "1/1:0,7,.:1,.:70,21,0,.,.,.", alt="C,G"),
            ligne("GT:AD:AF:PL", "1/1:0,7:1:70,21,0", alt="C", info="."))
        self.assertEqual(rendu, attendu)
        self.assertIn(conjoint, rendu)

    def test_une_ligne_trop_courte_passe_telle_quelle(self):
        for courte in ("chr1\t100\t.\tA\tG\t50\tPASS\t.\n",
                       "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:AD\n"):
            n, rendu = self.passer(ENTETE + courte)
            self.assertEqual(n, (0, 0))
            self.assertEqual(rendu, ENTETE + courte)

    def test_plus_d_une_colonne_est_un_refus(self):
        # Le filtre juge une colonne : sur deux, il jugerait la premiere.
        texte = ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
                 + ligne("GT", "0/1")[:-1] + "\t0/0\n")
        with self.assertRaises(ValueError):
            self.passer(texte)

    def test_un_flux_sans_en_tete_de_colonnes_est_un_refus(self):
        with self.assertRaises(ValueError):
            self.passer(ligne("GT", "0/1"))
        with self.assertRaises(ValueError):
            self.passer("##fileformat=VCFv4.2\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
