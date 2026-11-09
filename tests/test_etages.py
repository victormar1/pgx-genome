# -*- coding: utf-8 -*-
"""Le decoupage en etages est un contrat entre trois endroits : la liste des
noms que le module accepte, les blocs qui les executent, et la trace qui les
attend. Un etage ajoute ou renomme dans l'un et pas dans les autres serait
saute sans un mot — accepte en ligne de commande, jamais execute.

Ces tests lisent le script sans l'executer : ils tiennent sur tout poste, avec
ou sans moteur de conteneurs.
"""
import io
import os
import re
import shutil
import subprocess
import unittest

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(RACINE, "bin", "pgx_genome.sh")
TRACE = os.path.join(RACINE, "bin", "provenance.py")


def lire(chemin):
    with io.open(chemin, encoding="utf-8") as fh:
        return fh.read()


class Contrat(unittest.TestCase):
    def setUp(self):
        self.s = lire(SCRIPT)

    def connus(self):
        m = re.search(r'^CONNUS="([^"]+)"', self.s, re.M)
        self.assertIsNotNone(m, "liste CONNUS introuvable")
        return m.group(1).split()

    def blocs(self):
        """Les blocs d'etage, dans l'ordre du script, sans repetition."""
        vus = []
        for nom in re.findall(r"^if a_faire (\w+); then$", self.s, re.M):
            if nom not in vus:
                vus.append(nom)
        return vus

    def test_chaque_nom_accepte_a_son_bloc(self):
        # Un nom accepte sans bloc serait un etage annonce et jamais execute.
        manquants = set(self.connus()) - set(self.blocs())
        self.assertEqual(manquants, set())

    def test_chaque_bloc_porte_un_nom_accepte(self):
        # Un bloc sans nom accepte ne pourrait etre lance qu'en execution
        # complete : on ne saurait plus le relancer seul.
        orphelins = set(self.blocs()) - set(self.connus())
        self.assertEqual(orphelins, set())

    def test_les_blocs_suivent_l_ordre_de_la_chaine(self):
        # L'ordre des blocs est l'ordre de l'execution complete ; la liste
        # publiee dans l'aide doit le reproduire, sans quoi l'exploitant
        # enchainerait les etages a contresens.
        self.assertEqual(self.blocs(), self.connus())

    def test_la_trace_n_attend_que_des_etages_connus(self):
        m = re.search(r"ATTENDUS = \[([^\]]+)\]", lire(TRACE), re.S)
        self.assertIsNotNone(m)
        attendus = re.findall(r'"(\w+)"', m.group(1))
        for e in attendus:
            self.assertIn(e, self.connus(), e)

    def test_l_aide_nomme_chaque_etage(self):
        aide = "\n".join(l for l in self.s.split("\n")[:60] if l.startswith("#"))
        for e in self.connus():
            self.assertIn(e, aide, e)

    def test_la_plage_d_aide_couvre_tout_l_en_tete(self):
        # L'aide s'imprime par plage de lignes ; chaque option ajoutee allonge
        # l'en-tete, et la plage a deja du etre corrigee a la main deux fois.
        # Elle s'arrete sur la ligne vide qui ferme l'en-tete : plus court,
        # l'aide perd ses dernieres options ; plus long, elle imprime du code.
        m = re.search(r"sed -n '2,(\d+)p'", self.s)
        self.assertIsNotNone(m)
        lignes = self.s.split("\n")
        derniere = 1          # numero, compte a partir de un, du dernier commentaire
        while derniere < len(lignes) and lignes[derniere].startswith("#"):
            derniere += 1
        borne = int(m.group(1))
        self.assertEqual(borne, derniere + 1,
                         "la plage d'aide s'arrete a %d, l'en-tete se ferme a %d"
                         % (borne, derniere + 1))
        self.assertEqual(lignes[borne - 1].strip(), "",
                         "la ligne qui ferme l'en-tete n'est pas vide")


class Garanties(unittest.TestCase):
    """Ce que le decoupage doit preserver, lu dans le script."""

    def setUp(self):
        self.s = lire(SCRIPT)

    def test_seule_la_recevabilite_ouvre_le_dossier(self):
        # La purge du dossier de travail et la remise a zero de la trace ne
        # doivent dependre que de la recevabilite : un etage lance seul qui
        # purgerait effacerait le travail des autres.
        i = self.s.index("rm -rf \"$T\"; mkdir -p \"$T\"")
        avant = self.s[:i]
        self.assertGreater(avant.rfind("if a_faire recevabilite; then"),
                           avant.rfind("\nfi\n"),
                           "la purge n'est pas gardee par la recevabilite")

    def test_un_etage_seul_verifie_la_signature(self):
        self.assertIn('"$(cat "$SIG_FICHIER")" != "$SIGNATURE"', self.s)

    def test_un_refus_est_inscrit_pour_les_etages_suivants(self):
        # Les deux chemins de refus doivent laisser la marque que les etages
        # suivants lisent.
        self.assertEqual(self.s.count('> "$ARRET"'), 2)

    def test_le_code_de_retour_ne_lit_que_ses_propres_etages(self):
        # Trois etages tournent en parallele et ecrivent dans la meme table ;
        # l'echec de l'un ne doit pas devenir celui des autres.
        self.assertIn('index(e, "," $1 ",")', self.s)

    def test_l_ajout_a_la_table_passe_par_un_verrou(self):
        self.assertRegex(self.s, r"etat\(\) \{ inscrire")

    def test_une_image_designee_par_un_chemin_doit_exister(self):
        # Un .sif absent deviendrait un nom de depot invalide, et la
        # recevabilite conclurait a un fichier de variants sans echantillon.
        self.assertIn('[ -f "$i" ] || { echo "pgx_genome : image introuvable', self.s)

    def test_une_lecture_du_vcf_en_echec_n_est_pas_un_vcf_vide(self):
        i = self.s.index("NOMS_VCF=$(BCF query -l")
        bloc = self.s[i:i + 900]
        self.assertIn("RC_VCF=$?", bloc)
        self.assertLess(bloc.index('if [ "$RC_VCF" != 0 ]'), bloc.index('"$NB_VCF" != "1"'))

    def test_le_rendu_recoit_l_identite_et_le_contexte(self):
        # Sans eux, tout compte rendu sortirait en prototype, sans preindication :
        # une plateforme ne pourrait pas produire un document de patient.
        i = self.s.index("compte_rendu.py")
        bloc = self.s[i:self.s.index("\n\n", i)]
        self.assertIn('--identite "$IDENTITE"', bloc)
        self.assertIn('--contexte "$CONTEXTE"', bloc)

    def test_un_vcf_sans_ad_est_refuse_a_l_entree(self):
        # Sans lectures par allele, l'equilibre allelique de chaque heterozygote
        # et le critere de melange ne se mesureraient pas : le controle qualite
        # les sauterait sans le dire.
        i = self.s.index("NOMS_VCF=$(BCF query -l")
        bloc = self.s[i:self.s.index("\nfi\n", i)]
        j = bloc.find("grep -c '^##FORMAT=<ID=AD,'")
        self.assertGreater(j, 0, "la recevabilite ne verifie plus AD")
        self.assertIn("RECEVABLE=0", bloc[j:bloc.index("\nel", j)])

    def test_aucun_grep_q_derriere_un_tuyau(self):
        # Le script tourne sous pipefail : grep -q ferme le tuyau des la premiere
        # correspondance, la commande en amont meurt de SIGPIPE, et le tuyau
        # entier passe pour un echec. Un VCF conforme serait ainsi refuse.
        self.assertIn("set -uo pipefail", self.s)
        self.assertNotRegex(self.s, r"\|\s*grep -q")

    def test_des_identites_divergentes_ne_sont_pas_dites_concordantes(self):
        # Le champ SM d'un alignement peut etre un identifiant interne : la
        # divergence est acceptee, mais elle est dite et tracee comme telle.
        i = self.s.index("NOMS_VCF=$(BCF query -l")
        bloc = self.s[i:self.s.index("\nfi\n", i)]
        j = bloc.index('[ "$SM_CRAM" != "$NOMS_VCF" ]; then')
        suite = bloc[j:bloc.index("\nelse\n", j)]
        self.assertIn("identites divergentes acceptees", suite)
        self.assertNotIn("identites concordantes", suite)

    def test_un_signal_arrete_l_execution(self):
        # Sans exit dans le piege, bash reprend apres le nettoyage : l'execution
        # continue sans verrou, et une seconde peut demarrer sur le meme dossier.
        if not shutil.which("bash"):
            self.skipTest("bash absent")
        pieges = re.findall(r"^trap .*$", self.s, re.M)
        self.assertTrue(pieges)
        script = ("nettoyer() { :; }\n" + "\n".join(pieges)
                  + "\nkill -TERM $$\necho apres-le-signal\n")
        r = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
        self.assertNotIn("apres-le-signal", r.stdout)
        self.assertEqual(r.returncode, 143)

    def test_la_signature_de_reprise_porte_la_version_du_module(self):
        # Sans elle, --reprise reutiliserait apres une mise a jour les
        # intermediaires produits par la version precedente.
        m = re.search(r'^SIGNATURE="(.*)"$', self.s, re.M)
        self.assertIsNotNone(m)
        self.assertIn('version=$(cat "$RACINE/VERSION"', m.group(1))

    def test_la_recevabilite_purge_les_livrables_precedents(self):
        # Un etage en echec laisserait sinon le compte rendu ou la trace d'une
        # autre execution dans le dossier de sortie.
        i = self.s.index("if a_faire recevabilite; then")
        bloc = self.s[i:self.s.index("\nelse\n", i)]
        self.assertIn('rm -rf "$SORTIE/sortie"; mkdir -p "$SORTIE/sortie"', bloc)

    def test_une_trace_perimee_ne_decide_pas_du_code_de_retour(self):
        # Si la trace ne s'ecrit pas, le code de retour ne doit pas se lire dans
        # celle d'une execution anterieure.
        i = self.s.index("REUSSI=$(python3")
        avant = self.s[self.s.rindex("if a_faire provenance; then", 0, i):i]
        j = avant.index('rm -f "$SORTIE/sortie/provenance.json"')
        self.assertLess(j, avant.index("bin/provenance.py"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
