# -*- coding: utf-8 -*-
"""La chaine executee de bout en bout, sur des doublures des outils externes.

Les autres tests du decoupage lisent le texte du script. Ceux-ci l'executent :
un refus qui n'arrete rien, une reprise qui retrouve un fichier tronque ou un
code de retour qui ment ne se voient pas a la lecture. Les doublures sont dans
tests/faux/outils.py ; elles ne calculent rien, et seule l'orchestration se
teste ici.

Linux seulement : le module ne s'execute pas ailleurs.
"""
import json
import os
import shutil
import stat
import subprocess
import tempfile
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(RACINE, "bin", "pgx_genome.sh")
LOT = os.path.join(RACINE, "bin", "pgx_lot.sh")
DOUBLURES = os.path.join(RACINE, "tests", "faux", "outils.py")
VCF = os.path.join(RACINE, "bench", "qualification", "HG00699.vcf.gz")


def executable():
    if os.name == "nt":
        return "le module ne s'execute que sous Linux"
    for outil in ("bash", "python3"):
        if not shutil.which(outil):
            return outil + " absent"
    r = subprocess.run(["python3", "-c", "import reportlab"], capture_output=True)
    if r.returncode != 0:
        return "reportlab absent de python3"
    return ""


TABLES = ["HG00699.annotations.tsv", "HG00699.genotypes.tsv"]


@unittest.skipIf(executable(), executable())
class Execution(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="pgx_essai_")
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        faux = os.path.join(self.d, "faux")
        os.makedirs(faux)
        for nom in ("samtools", "docker"):
            p = os.path.join(faux, nom)
            with open(p, "w", encoding="utf-8", newline="\n") as fh:
                fh.write('#!/bin/sh\nFAUX_ROLE=%s exec python3 "%s" "$@"\n'
                         % (nom, DOUBLURES))
            os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR)
        self.cyrius = os.path.join(self.d, "cyrius")
        os.makedirs(self.cyrius)
        with open(os.path.join(self.cyrius, "star_caller.py"), "w",
                  encoding="utf-8", newline="\n") as fh:
            fh.write("import os, runpy\nos.environ['FAUX_ROLE'] = 'star_caller.py'\n"
                     "runpy.run_path(%r, run_name='__main__')\n" % DOUBLURES)
        self.cram = os.path.join(self.d, "HG00699.cram")
        with open(self.cram, "wb") as fh:
            fh.write(b"alignement de substitution")
        self.sortie = os.path.join(self.d, "resultat")
        self.appels = os.path.join(self.d, "appels.jsonl")
        self.env = {k: v for k, v in os.environ.items()
                    if not k.startswith(("PGX_", "FAUX_", "PYPGX"))}
        self.env["PATH"] = faux + os.pathsep + self.env["PATH"]
        self.env["PGX_CYRIUS"] = self.cyrius
        self.env["FAUX_JOURNAL"] = self.appels

    def lancer(self, *options, **reglages):
        cwd = reglages.pop("cwd", None)
        sortie = reglages.pop("sortie", self.sortie)
        cram = reglages.pop("cram", self.cram)
        env = dict(self.env, **reglages)
        return subprocess.run(
            ["bash", SCRIPT, "--cram", cram, "--vcf", VCF, "--sortie", sortie,
             "--echantillon", "HG00699"] + list(options),
            env=env, cwd=cwd, capture_output=True, text=True)

    def trace(self, sortie=None):
        p = os.path.join(sortie or self.sortie, "sortie", "provenance.json")
        with open(p, encoding="utf-8") as fh:
            return json.load(fh)

    def etat(self, etage, sortie=None):
        return self.trace(sortie)["etages"][etage]["etat"]

    def appele(self, mot):
        if not os.path.exists(self.appels):
            return 0
        with open(self.appels, encoding="utf-8") as fh:
            return sum(1 for l in fh if mot in l)

    def livrables(self):
        d = os.path.join(self.sortie, "sortie")
        return sorted(os.listdir(d)) if os.path.isdir(d) else []

    # ------------------------------------------------------------ chaine entiere
    def test_une_execution_complete_rend_ses_livrables(self):
        r = self.lancer()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        t = self.trace()
        self.assertTrue(t["reussite_complete"])
        self.assertEqual(t["etages_en_echec"], [])
        self.assertEqual(self.livrables(), ["CR_HG00699.pdf"] + TABLES + [
            "HG00699.report.json", "HG00699.variants.tsv", "provenance.json"])
        with open(os.path.join(self.sortie, "travail", "appels_externes.tsv"),
                  encoding="utf-8") as fh:
            genes = [l.split("\t")[0] for l in fh if l.strip()]
        self.assertEqual(genes, ["CYP2D6", "HLA-A", "HLA-B"])

    def test_l_en_tete_de_l_alignement_n_est_lu_qu_une_fois(self):
        # Sur un alignement distant, chaque lecture d'en-tete est un aller-retour.
        self.assertEqual(self.lancer().returncode, 0)
        with open(self.appels, encoding="utf-8") as fh:
            lectures = [l for l in fh if '"-H"' in l and "HG00699.cram" in l]
        self.assertEqual(len(lectures), 1, lectures)

    def test_un_dossier_de_sortie_relatif_est_accepte(self):
        # Les conteneurs s'executent dans leur propre dossier de travail : un
        # chemin relatif au dossier d'appel n'y designe plus rien.
        r = self.lancer(sortie="resultat", cram="HG00699.cram", cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(self.trace()["reussite_complete"])

    # ------------------------------------------------------------------- refus
    def test_une_couverture_hors_domaine_ne_rend_aucun_resultat(self):
        r = self.lancer(FAUX_PROFONDEUR="12")
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertEqual(self.livrables(), ["provenance.json"])
        t = self.trace()
        self.assertFalse(t["reussite_complete"])
        self.assertEqual(t["etages"]["qc"]["etat"], "OK")
        for e in ("cyp2d6", "hla", "appels", "pharmcat", "rendu"):
            self.assertEqual(t["etages"][e]["etat"], "ECHEC", e)
        self.assertEqual(self.appele("star_caller"), 0)
        self.assertEqual(self.appele("pharmcat.jar"), 0)

    def test_un_etage_lance_apres_un_refus_ne_tourne_pas(self):
        self.lancer(FAUX_PROFONDEUR="12")
        r = self.lancer("--etage", "pharmcat", FAUX_PROFONDEUR="12")
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertEqual(self.appele("pharmcat.jar"), 0)

    def test_des_identites_divergentes_refusent_le_genome(self):
        r = self.lancer("--echantillon", "AUTRE", FAUX_SM="INCONNU")
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        t = self.trace()
        self.assertEqual(t["etages"]["recevabilite"]["etat"], "ECHEC")
        self.assertEqual(self.livrables(), ["provenance.json"])
        self.assertEqual(self.appele('"view", "-b"'), 0)

    def refus_a_l_entree(self, detail, **reglages):
        """Refuse a la recevabilite : le motif est dans l'etat, rien n'est
        produit, et l'alignement n'est pas meme parcouru."""
        r = self.lancer(**reglages)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        e = self.trace()["etages"]["recevabilite"]
        self.assertEqual(e["etat"], "ECHEC")
        self.assertIn(detail, e["detail"])
        self.assertEqual(self.livrables(), ["provenance.json"])
        self.assertEqual(self.appele('"view", "-b"'), 0)

    def test_un_vcf_illisible_n_est_pas_un_vcf_vide(self):
        # bcftools n'a rien lu : dire « zero echantillon » enverrait chercher
        # la faute dans les donnees du patient.
        self.refus_a_l_entree("lecture du VCF impossible", FAUX_PANNE="vcf_illisible")

    def test_un_vcf_a_deux_echantillons_est_refuse(self):
        self.refus_a_l_entree("2 echantillons", FAUX_VCF="deux_echantillons")

    def test_un_vcf_sans_ad_est_refuse_a_l_entree(self):
        # Sans lectures par allele, ni l'equilibre allelique ni le melange ne
        # se mesurent : le controle qualite les sauterait sans le dire.
        self.refus_a_l_entree("VCF sans AD", FAUX_VCF="sans_ad")

    def test_un_alignement_tronque_est_refuse(self):
        self.refus_a_l_entree("alignement tronque", FAUX_PANNE="quickcheck")

    def test_des_identites_divergentes_acceptees_sont_dites_comme_telles(self):
        # Le champ SM d'un alignement peut etre un identifiant interne : la
        # demande nomme le VCF, la divergence est acceptee, jamais presentee
        # comme une concordance.
        r = self.lancer(FAUX_SM="INTERNE")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        e = self.trace()["etages"]["recevabilite"]
        self.assertEqual(e["etat"], "OK")
        self.assertIn("identites divergentes acceptees", e["detail"])
        self.assertIn("alignement=INTERNE", e["detail"])
        self.assertNotIn("identites concordantes", r.stdout)

    def test_la_trace_d_un_refus_porte_les_seuils_appliques(self):
        # Un genome refuse a l'entree n'a pas ete juge sans seuil : la trace
        # doit dire ceux de l'execution, pas zero.
        self.lancer("--echantillon", "AUTRE", FAUX_SM="INCONNU")
        seuils = self.trace()["seuils"]
        self.assertEqual(seuils["couverture_mediane_min"], 18)
        self.assertEqual(seuils["melange_part_max"], 35.0)

    # ------------------------------------------------------- typages sans appel
    def test_un_cyp2d6_ambigu_n_est_pas_transmis(self):
        r = self.lancer(FAUX_D6="*1/*46;*43/*45")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        t = self.trace()
        self.assertEqual(t["etages"]["cyp2d6"]["etat"], "SANS_RESULTAT")
        self.assertTrue(t["reussite_complete"])
        with open(os.path.join(self.sortie, "travail", "appels_externes.tsv"),
                  encoding="utf-8") as fh:
            self.assertNotIn("CYP2D6", fh.read())

    def test_un_cyp2d6_que_cyrius_n_a_pas_valide_n_est_pas_transmis(self):
        # L'etage et celui des appels externes jugent par la meme fonction : un
        # diplotype marque par l'outil est un etage sans resultat, et non un
        # etage reussi dont la ligne manquerait ensuite.
        for d6, filtre in (("*1/*2x3", "LowQ_high_CN"), ("*1/*2/*3", "PASS"), ("*1", "PASS")):
            r = self.lancer(FAUX_D6=d6, FAUX_D6_FILTRE=filtre)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            t = self.trace()
            self.assertEqual(t["etages"]["cyp2d6"]["etat"], "SANS_RESULTAT", d6)
            self.assertEqual(t["etages"]["appels"]["etat"], "OK", d6)
            self.assertTrue(t["reussite_complete"])
            with open(os.path.join(self.sortie, "travail", "appels_externes.tsv"),
                      encoding="utf-8") as fh:
                self.assertNotIn("CYP2D6", fh.read())
        self.assertIn("LowQ_high_CN", self.trace_de("*1/*2x3", "LowQ_high_CN"))

    def trace_de(self, d6, filtre):
        self.lancer(FAUX_D6=d6, FAUX_D6_FILTRE=filtre)
        return self.trace()["etages"]["cyp2d6"]["detail"]

    def test_un_fichier_de_cyrius_illisible_est_un_echec_et_non_une_abstention(self):
        # Le jugement plante sur ce fichier : le prendre pour « l'outil n'a rien
        # conclu » rendrait le gene non analyse, module en reussite complete.
        r = self.lancer(FAUX_PANNE="d6_illisible")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        t = self.trace()
        self.assertEqual(t["etages"]["cyp2d6"]["etat"], "ECHEC")
        self.assertIn("cyp2d6", t["etages_en_echec"])
        self.assertFalse(t["reussite_complete"])

    def test_un_tandem_valide_est_transmis(self):
        r = self.lancer(FAUX_D6="*36+*10/*4")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.etat("cyp2d6"), "OK")
        with open(os.path.join(self.sortie, "travail", "appels_externes.tsv"),
                  encoding="utf-8") as fh:
            self.assertIn("CYP2D6\t*36 + *10/*4\n", fh.read())

    def test_un_fichier_d_identite_sans_patient_ne_rend_pas_de_compte_rendu(self):
        ident = os.path.join(self.d, "identite.json")
        with open(ident, "w", encoding="utf-8") as fh:
            fh.write("{}")
        r = self.lancer("--identite", ident)
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertEqual(self.trace()["etages_en_echec"], ["rendu"])
        self.assertNotIn("CR_HG00699.pdf", self.livrables())
        with open(ident, "w", encoding="utf-8") as fh:
            json.dump({"patient": "ESSAI Marie", "naissance": "01/02/1960"}, fh)
        r = self.lancer("--identite", ident)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("CR_HG00699.pdf", self.livrables())

    def test_sans_les_tables_le_compte_rendu_ne_sort_pas_seul(self):
        # Le document annonce des resultats « disponibles au format
        # structure » : il ne sort pas sans eux, et ceux d'avant ne restent pas.
        self.assertEqual(self.lancer().returncode, 0)
        os.remove(os.path.join(self.sortie, "travail", "qc_positions.tsv"))
        r = self.lancer("--etage", "rendu")
        self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.livrables(), ["HG00699.report.json", "provenance.json"])
        with open(os.path.join(self.sortie, "travail", "etats.tsv"),
                  encoding="utf-8") as fh:
            self.assertIn("rendu	ECHEC", fh.read().split("rendu	OK")[-1])

    def test_un_rendu_qui_echoue_ne_laisse_pas_ceux_d_avant(self):
        self.assertEqual(self.lancer().returncode, 0)
        ident = os.path.join(self.d, "identite.json")
        with open(ident, "w", encoding="utf-8") as fh:
            fh.write("{}")
        r = self.lancer("--etage", "rendu", "--identite", ident)
        self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.livrables(), ["HG00699.report.json", "provenance.json"])

    def test_aucun_appel_externe_n_est_pas_un_echec_des_appels(self):
        # Cyrius ne conclut pas et le typage HLA echoue : il n'y a aucune ligne
        # a transmettre, et c'est le typage HLA qui a failli, pas cet etage.
        r = self.lancer(FAUX_D6="*1/*46;*43/*45", FAUX_PAIRES="10")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        t = self.trace()
        self.assertEqual(t["etages"]["hla"]["etat"], "ECHEC")
        self.assertEqual(t["etages"]["appels"]["etat"], "OK")
        self.assertEqual(t["etages_en_echec"], ["hla"])
        self.assertNotIn("integer expression", r.stderr)

    def test_une_extraction_hla_interrompue_n_est_pas_typee(self):
        # Le flux s'arrete en route : les fichiers de lectures existent, assez
        # gros pour passer le seuil, et ne portent qu'une partie du complexe.
        r = self.lancer(FAUX_PANNE="collate")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertEqual(self.etat("hla"), "ECHEC")
        self.assertEqual(self.appele('"OptiTypePipeline.py", "-i"'), 0)

    # ------------------------------------------------------------------ reprise
    def test_la_reprise_ne_type_pas_des_lectures_tronquees(self):
        self.lancer(FAUX_PANNE="collate")
        r = self.lancer("--reprise")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        fq = os.path.join(self.sortie, "travail", "mhc_1.fq")
        with open(fq, encoding="ascii") as fh:
            self.assertEqual(sum(1 for _ in fh) // 4, 500)

    def test_la_reprise_survit_a_un_tri_interrompu(self):
        r = self.lancer(FAUX_PANNE="sort_tue")
        self.assertEqual(r.returncode, 143, r.stdout + r.stderr)
        r = self.lancer("--reprise")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(self.sortie, "travail", "tranche.bam"), "rb") as fh:
            self.assertEqual(fh.read(), b"alignement trie")

    def test_la_reprise_n_emploie_pas_un_index_interrompu(self):
        r = self.lancer(FAUX_PANNE="index_tue")
        self.assertEqual(r.returncode, 143, r.stdout + r.stderr)
        r = self.lancer("--reprise")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(self.sortie, "travail", "tranche.bam.bai"), "rb") as fh:
            self.assertEqual(fh.read(), b"index")

    def test_la_reprise_reutilise_la_tranche_et_le_typage_hla(self):
        self.lancer()
        os.remove(self.appels)
        r = self.lancer("--reprise")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.appele('"view", "-b"'), 0)
        self.assertEqual(self.appele('"OptiTypePipeline.py", "-i"'), 0)
        self.assertIn("repris", self.trace()["etages"]["hla"]["detail"])

    def test_sans_reprise_rien_n_est_reutilise(self):
        self.lancer()
        os.remove(self.appels)
        r = self.lancer()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.appele('"view", "-b"'), 1)
        self.assertEqual(self.appele('"OptiTypePipeline.py", "-i"'), 1)

    # ------------------------------------------------------------------- etages
    def test_les_etages_lances_un_a_un_rendent_la_meme_trace(self):
        for e in ("recevabilite", "filtre", "tranche", "qc", "cyp2d6", "hla",
                  "pypgx", "appels", "pharmcat", "rendu"):
            r = self.lancer("--etage", e)
            self.assertEqual(r.returncode, 0, e + "\n" + r.stdout + r.stderr)
        r = self.lancer("--etage", "provenance")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        t = self.trace()
        self.assertTrue(t["reussite_complete"])
        self.assertEqual(t["etages"]["pypgx"]["etat"], "IGNORE")
        self.assertEqual(t["seuils"]["couverture_mediane_min"], 18)
        self.assertNotEqual(t["couverture_mediane_perimetre"], "")

    def test_une_trace_perimee_ne_decide_pas_du_code_de_retour(self):
        # La trace ne s'ecrit pas : le code de retour ne doit pas se lire dans
        # celle, complete, de l'execution anterieure.
        self.assertEqual(self.lancer().returncode, 0)
        with open(os.path.join(self.sortie, "travail", "etats.tsv"), "a",
                  encoding="utf-8", newline="\n") as fh:
            fh.write("rendu\tOK\tillisible\tligne abimee\n")
        r = self.lancer("--etage", "provenance")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertEqual(self.livrables(), ["CR_HG00699.pdf"] + TABLES + [
            "HG00699.report.json", "HG00699.variants.tsv"])

    def test_un_etage_seul_refuse_d_autres_entrees(self):
        self.lancer()
        r = self.lancer("--etage", "pharmcat", "--gq", "30")
        self.assertEqual(r.returncode, 2)
        self.assertIn("recevabilite", r.stderr)

    def test_une_seconde_execution_sur_le_meme_dossier_est_refusee(self):
        os.makedirs(os.path.join(self.sortie, ".verrou"))
        with open(os.path.join(self.sortie, ".verrou", "pid"), "w") as fh:
            fh.write(str(os.getpid()))
        r = self.lancer()
        self.assertEqual(r.returncode, 3, r.stdout + r.stderr)
        self.assertEqual(self.appele("samtools"), 0)

    def test_des_intermediaires_sans_signature_ne_sont_jamais_repris(self):
        # Une tranche dont rien ne dit d'ou elle vient : aucune option ne doit
        # la faire passer pour celle de ces entrees.
        t = os.path.join(self.sortie, "travail")
        os.makedirs(t)
        for nom in ("tranche.bam", "tranche.bam.bai"):
            with open(os.path.join(t, nom), "wb") as fh:
                fh.write(b"origine inconnue")
        for options in ((), ("--reprise",)):
            r = self.lancer(*options)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            with open(os.path.join(t, "tranche.bam"), "rb") as fh:
                self.assertEqual(fh.read(), b"alignement trie", options)
            shutil.rmtree(t)
            os.makedirs(t)
            for nom in ("tranche.bam", "tranche.bam.bai"):
                with open(os.path.join(t, nom), "wb") as fh:
                    fh.write(b"origine inconnue")

    # ---------------------------------------------------------------------- lot
    def lot(self, cram):
        manifeste = os.path.join(self.d, "manifeste.tsv")
        with open(manifeste, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("HG00699\t%s\t%s\n" % (cram, VCF))
        return subprocess.run(
            ["bash", LOT, "--manifeste", manifeste, "--sortie",
             os.path.join(self.d, "lot"), "--parallele", "1"],
            env=self.env, capture_output=True, text=True)

    def test_un_lot_relance_saute_l_echantillon_deja_complet(self):
        r = self.lot(self.cram)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("reussite complete : 1", r.stdout)
        os.remove(self.appels)
        r = self.lot(self.cram)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.appele("samtools"), 0)
        self.assertIn("reussite complete : 1", r.stdout)
        with open(os.path.join(self.d, "lot", "lot.log"), encoding="utf-8") as fh:
            self.assertIn("deja complet", fh.read())

    def test_un_lot_relance_sur_un_autre_alignement_le_retraite(self):
        # Le manifeste est corrige : le meme identifiant designe un autre
        # fichier. Le resultat du premier ne doit pas rester sous ce nom.
        self.lot(self.cram)
        autre = os.path.join(self.d, "corrige.cram")
        with open(autre, "wb") as fh:
            fh.write(b"un autre alignement")
        os.remove(self.appels)
        r = self.lot(autre)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.appele('"view", "-b"'), 1)
        t = self.trace(os.path.join(self.d, "lot", "HG00699"))
        self.assertEqual(t["entrees"]["alignement"], autre)


if __name__ == "__main__":
    unittest.main(verbosity=2)
