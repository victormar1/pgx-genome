# -*- coding: utf-8 -*-
"""Etage 4b : typage des genes du core panel que l'interpreteur ne sait pas appeler.

PharmCAT connait vingt-trois genes. Quatre genes de classe 1 et 2 du core panel
RNPGx lui echappent, et PyPGx sait les typer : MT-RNR1, BCHE, MTHFR, POR. Cet
etage les lui confie, puis rend le seul resultat que l'on peut justifier.

Le point dur n'est pas le typage, c'est la recevabilite. Sur un fichier de
variants qui ne couvre pas la region, PyPGx rend « Reference/Reference » sans
lever la moindre erreur : un genome dont le chrM n'a jamais ete appele
ressortirait « risque normal de surdite sous aminoside ». Un appel n'est donc
rendu que si ses positions definissantes ont ete lues sur l'alignement et si le
contig figure bien dans le fichier de variants. Sinon le gene est non conclusif,
et le compte rendu le dit.

MT-RNR1 est le seul a porter un phenotype. PharmCAT embarque sa table de
recommandations sans appeleur : le resultat lui revient en appel externe, par le
meme chemin que HLA et CYP2D6. Les trois autres sortent avec un allele nomme et
sans interpretation, faute de source clinique pour la leur donner.
"""
import argparse, csv, json, os, shutil, subprocess, sys, tempfile, zipfile


def profondeurs(chemin):
    """samtools depth -a : contig, position, profondeur."""
    d = {}
    if not chemin or not os.path.exists(chemin):
        return d
    with open(chemin, encoding="utf-8") as fh:
        for ligne in fh:
            c = ligne.rstrip("\n").split("\t")
            if len(c) >= 3:
                try:
                    d[(c[0], int(c[1]))] = int(c[2])
                except ValueError:
                    pass
    return d


def exige(chemin, quoi):
    """Un fichier de controle demande mais introuvable n'est pas une absence
    d'information : c'est une entree cassee. S'en passer en silence rendrait le
    garde-fou contournable par simple suppression du fichier — l'etage doit
    echouer, pas devenir indulgent."""
    if chemin and not os.path.exists(chemin):
        sys.stderr.write("%s demande mais introuvable : %s\n" % (quoi, chemin))
        sys.exit(3)


def contigs(chemin):
    if not chemin:
        return None   # option non fournie
    with open(chemin, encoding="utf-8") as fh:
        return {l.strip() for l in fh if l.strip()}


def qualites(chemin):
    """bcftools query sur les positions definissantes : contig, position, FILTER,
    GQ, DP. Les positions absentes du fichier n'y figurent pas, et c'est la
    profondeur d'alignement qui tranche leur sort."""
    d = {}
    if not chemin:
        return None   # option non fournie
    with open(chemin, encoding="utf-8") as fh:
        for ligne in fh:
            c = ligne.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            nb = lambda x: None if x in ("", ".", None) else (
                int(float(x)) if str(x).replace(".", "", 1).isdigit() else None)
            try:
                pos = int(c[1])
            except ValueError:
                continue   # ligne illisible : la profondeur tranchera
            d[(c[0], pos)] = {"filtre": c[2],
                                    "gq": nb(c[3]) if len(c) > 3 else None,
                                    "dp": nb(c[4]) if len(c) > 4 else None}
    return d


def recevabilite(gene, info, prof, vus, qual, prof_min, gq_min, ignorer_filtre):
    """Le gene peut-il etre type ? Renvoie (recevable, lues, attendues, motif).

    Deux refus distincts, et il faut les deux. La profondeur d'alignement traite
    la position absente du fichier de variants : sans elle, une position jamais
    lue passerait pour une reference. La qualite traite la position presente
    mais douteuse : le typeur s'en servirait sans le dire, alors que le controle
    qualite l'a ecartee du reste du module. Un seul critere pour tout le module,
    ou le compte rendu se contredit d'une section a l'autre."""
    if vus is not None and not vus:
        return False, 0, len(info["positions"]), \
            "aucun contig releve dans le fichier de variants"
    if vus is not None and info["contig"] not in vus:
        return False, 0, len(info["positions"]), \
            "contig %s absent du fichier de variants" % info["contig"]
    lues, inconnues, douteuses = 0, 0, []
    for p in info["positions"]:
        contig, pos = p.split(":")[0], int(p.split(":")[1])
        d = prof.get((contig, pos))
        if d is None:
            inconnues += 1
        elif d >= prof_min:
            lues += 1
        q = (qual or {}).get((contig, pos))
        if q is None:
            continue   # position absente du fichier : la profondeur a deja tranche
        if not ignorer_filtre and q["filtre"] not in ("PASS", ".", ""):
            douteuses.append("%s:%d (%s)" % (contig, pos, q["filtre"]))
        elif q["gq"] is None:
            douteuses.append("%s:%d (sans GQ)" % (contig, pos))
        elif q["gq"] < gq_min:
            douteuses.append("%s:%d (GQ %d)" % (contig, pos, q["gq"]))
    total = len(info["positions"])
    if inconnues == total:
        return False, lues, total, "profondeur non mesuree sur les positions du gene"
    if lues < total:
        return False, lues, total, \
            "%d position(s) sur %d sous %dx ou non lue(s)" % (total - lues, total, prof_min)
    if douteuses:
        return False, lues, total, "%d genotype(s) ecarte(s) sur leur qualite : %s%s" % (
            len(douteuses), ", ".join(douteuses[:3]), " ..." if len(douteuses) > 3 else "")
    return True, lues, total, ""


def lance(pypgx_cmd, gene, variants, dossier):
    """run-ngs-pipeline sur le seul fichier de variants : aucun de ces quatre
    genes ne porte de variation de structure, ni profondeur ni statistiques de
    controle ne sont donc requises."""
    if os.path.isdir(dossier):
        shutil.rmtree(dossier)
    cmd = pypgx_cmd + ["run-ngs-pipeline", gene, dossier, "--variants", variants,
                       "--assembly", "GRCh38",
                       "--do-not-plot-copy-number", "--do-not-plot-allele-fraction"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        derniere = [l for l in (r.stdout + r.stderr).splitlines() if l.strip()]
        return None, (derniere[-1][:160] if derniere else "code %d" % r.returncode)
    return dossier, None


def lis_resultat(dossier):
    """La table de resultats est dans results.zip, sous un nom qui depend du gene."""
    z = os.path.join(dossier, "results.zip")
    if not os.path.exists(z):
        return None, "results.zip absent"
    with zipfile.ZipFile(z) as arc:
        noms = [n for n in arc.namelist() if n.endswith("data.tsv")]
        if not noms:
            return None, "data.tsv absent de results.zip"
        # Prendre le premier d'une archive qui en contient plusieurs reviendrait
        # a tirer au sort le resultat ; de meme pour une table a plusieurs
        # lignes, qui signalerait un appel multi-echantillon.
        if len(noms) > 1:
            return None, "%d data.tsv dans results.zip" % len(noms)
        texte = arc.read(noms[0]).decode("utf-8")
    lignes = list(csv.DictReader(texte.splitlines(), delimiter="\t"))
    if not lignes:
        return None, "table de resultats vide"
    if len(lignes) > 1:
        return None, "%d echantillons dans la table de resultats" % len(lignes)
    return lignes[0], None


def haplotypes(ligne):
    """Les colonnes Haplotype1 et Haplotype2 portent TOUS les alleles nommes sur
    chaque chromosome. La colonne Genotype, elle, n'en retient qu'un par
    haplotype, et le choix suit l'ordre de la table interne, pas la gravite.

    Sur BCHE, un porteur de rs1803274 et rs1799807 sur le meme chromosome
    ressortait « Reference/rs1803274 » : l'allele atypique — celui de l'apnee
    prolongee a la succinylcholine — disparaissait du genotype alors qu'il
    figurait encore, une colonne plus loin, dans l'haplotype. Quatre porteurs
    sur deux cent quarante-trois etaient ainsi tus, tous dans le meme sens.
    On lit donc les haplotypes, jamais le genotype tout fait."""
    return [[x.strip() for x in (ligne.get(c) or "").split(";") if x.strip()]
            for c in ("Haplotype1", "Haplotype2")]


def genotype_complet(ligne):
    """Un haplotype qui porte plusieurs alleles les garde tous, joints par « + »,
    la notation dont PharmCAT se sert deja pour les duplications.

    L'ordre est fixe par tri, et non par le rang des colonnes : PyPGx n'attribue
    pas Haplotype1 au meme chromosome d'un genome a l'autre, et le meme porteur
    sortait « *1/*28 » ici et « *28/*1 » la. Deux genomes identiques doivent
    rendre la meme chaine, sinon rien n'est comparable — ni entre patients, ni
    d'une version du pipeline a la suivante."""
    h = haplotypes(ligne)
    if not any(h):
        return ""
    return "/".join(sorted("+".join(sorted(x)) if x else "Reference" for x in h))


def genotype_panel(ligne, info):
    """Le genotype restreint aux alleles que le panel RNPGx nomme.

    Le typeur en connait bien plus que le panel n'en retient : neuf alleles sur
    MTHFR contre deux, quarante-sept sur POR contre ceux qui portent *28. Tout
    rendre donnerait « rs1476413+rs17367504+rs1801131+rs2274976+rs3737967+
    rs4846051/... » sur trois genomes sur quatre — exact, illisible, et sans un
    seul element actionnable de plus.

    La restriction n'est pas l'ecrasement que l'on corrige : elle suit le
    perimetre clinique deja retenu pour tout le produit, elle est la meme pour
    tous les genomes, et le nombre d'alleles mis de cote est rendu a cote. Sur
    BCHE elle ne retire rien : les deux positions sont au panel, de classe 1."""
    pan = set(info.get("positions_panel") or [])
    carte = info.get("alleles_positions") or {}
    ref = info.get("allele_reference")
    # Un catalogue incomplet ne doit pas faire retomber le rendu sur le genotype
    # complet comme s'il etait restreint : c'est une ressource a regenerer, et
    # l'etage doit le dire.
    if not pan or not carte or not ref:
        return None, 0, "catalogue sans positions du panel, carte d'alleles ou allele de reference"
    retenus, ecartes = [], 0
    for h in haplotypes(ligne):
        garde = []
        for a in h:
            if a in ("Reference", ref):
                continue
            if a not in carte:
                # Le typeur nomme un allele que le catalogue ne connait pas : le
                # compter hors panel reviendrait a decider sans savoir.
                return None, 0, "allele %s absent du catalogue (version de PyPGx ?)" % a
            if pan & set(carte[a]):
                garde.append(a)
            else:
                ecartes += 1
        retenus.append("+".join(sorted(garde)) if garde else ref)
    return "/".join(sorted(retenus)), ecartes, None


def allele_unique(ligne):
    """Pour un gene haploide, l'interpreteur attend un allele et non une paire.
    Deux haplotypes differents sur l'ADN mitochondrial signalent une
    heteroplasmie : on retient alors l'allele non reference, parce qu'un porteur
    heteroplasmique reste expose, et on le signale.

    Un haplotype porteur de plusieurs alleles nommes n'a pas de traduction en
    appel externe : PharmCAT attend un nom, pas une combinaison. Plutot que d'en
    choisir un — c'est exactement la faute que l'on corrige ici — on le declare
    ambigu et le gene sort non conclusif. Il en va de meme pour deux alleles non
    reference differents : garder le premier des deux reproduirait la perte
    silencieuse sur le gene ou elle coute le plus cher."""
    h1, h2 = haplotypes(ligne)
    if len(h1) > 1 or len(h2) > 1:
        return None, False, "+".join(h1) + "/" + "+".join(h2)
    h1 = h1[0] if h1 else ""
    h2 = h2[0] if h2 else ""
    if h1 == h2:
        return h1, False, None
    non_ref = [h for h in (h1, h2) if h and h != "Reference"]
    if len(non_ref) > 1:
        return None, True, h1 + "/" + h2
    return (non_ref[0] if non_ref else h1), True, None


def fraction(ligne, allele):
    """PyPGx joint a chaque allele la fraction allelique observee ; sur l'ADN
    mitochondrial c'est le taux d'heteroplasmie."""
    for bloc in (ligne.get("VariantData") or "").split(";"):
        p = bloc.split(":")
        if len(p) >= 3 and p[0] == allele:
            try:
                return float(p[-1])
            except ValueError:
                return None
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--catalogue", required=True)
    p.add_argument("--variants", required=True)
    p.add_argument("--profondeur")
    p.add_argument("--contigs-vcf")
    p.add_argument("--qualites", help="TSV contig/position/FILTER/GQ/DP sur les positions definissantes")
    p.add_argument("--profondeur-min", type=int, default=10)
    p.add_argument("--gq", type=int, default=20)
    p.add_argument("--ignorer-filtre", action="store_true")
    p.add_argument("--travail", required=True)
    p.add_argument("--sortie", required=True)
    p.add_argument("--json")
    p.add_argument("--pypgx", default="pypgx")
    a = p.parse_args()

    cat = json.load(open(a.catalogue, encoding="utf-8"))
    for chemin, quoi in ((a.profondeur, "profondeur"),
                         (a.contigs_vcf, "liste des contigs"),
                         (a.qualites, "qualites")):
        exige(chemin, quoi)
    prof = profondeurs(a.profondeur)
    vus = contigs(a.contigs_vcf)
    qual = qualites(a.qualites)
    pypgx_cmd = a.pypgx.split()
    travail = os.path.join(a.travail, "pypgx")
    os.makedirs(travail, exist_ok=True)

    resultats, rendus = [], 0
    for gene in sorted(cat["genes"]):
        info = cat["genes"][gene]
        ligne = {"gene": gene, "classe_rnpgx": info["classe_rnpgx"],
                 "genotype": "", "genotype_complet": "",
                 "genotype_pypgx": "", "alleles_hors_panel": 0,
                 "phenotype": "",
                 "allele_externe": "",
                 "heteroplasmie": False, "fraction": None}
        ok, lues, total, motif = recevabilite(gene, info, prof, vus, qual,
                                              a.profondeur_min, a.gq, a.ignorer_filtre)
        ligne.update({"positions_lues": lues, "positions_attendues": total})
        if not ok:
            ligne.update({"statut": "non conclusif", "motif": motif})
            resultats.append(ligne)
            continue
        dossier, err = lance(pypgx_cmd, gene, a.variants, os.path.join(travail, gene))
        if err:
            ligne.update({"statut": "non conclusif", "motif": "PyPGx : " + err})
            resultats.append(ligne)
            continue
        brut, err = lis_resultat(dossier)
        if err:
            ligne.update({"statut": "non conclusif", "motif": "PyPGx : " + err})
            resultats.append(ligne)
            continue
        complet = genotype_complet(brut) or (brut.get("Genotype") or "").strip()
        panel, ecartes, souci = genotype_panel(brut, info)
        ligne["genotype_complet"] = complet
        if panel is None:
            ligne.update({"statut": "non conclusif", "motif": souci})
            resultats.append(ligne)
            continue
        # Le compte rendu lit « genotype » : c'est la version du panel qui y va,
        # et la version complete reste a cote pour qui veut tout voir.
        ligne["genotype"] = panel
        ligne["alleles_hors_panel"] = ecartes
        # Le genotype tel que PyPGx le rend est conserve a cote : c'est la trace
        # de l'ecart, et le typeur reste verifiable sur sa propre sortie.
        ligne["genotype_pypgx"] = (brut.get("Genotype") or "").strip()
        phen = (brut.get("Phenotype") or "").strip()
        # « Indeterminate » n'est pas un phenotype : c'est PyPGx qui dit qu'il
        # n'a pas de table pour ce gene. On ne le restitue pas comme un resultat.
        ligne["phenotype"] = phen if info["phenotype"] and phen != "Indeterminate" else ""
        if info["appel_externe_pharmcat"] and info["haploide"]:
            allele, hetero, ambigu = allele_unique(brut)
            if ambigu:
                ligne.update({"statut": "non conclusif",
                              "motif": "plusieurs alleles nommes, sans traduction"
                                       " en appel externe : " + ambigu})
                resultats.append(ligne)
                continue
            ligne["allele_externe"] = allele
            ligne["heteroplasmie"] = hetero
            ligne["fraction"] = fraction(brut, allele)
        ligne.update({"statut": "rendu", "motif": ""})
        rendus += 1
        resultats.append(ligne)

    champs = ["gene", "classe_rnpgx", "statut", "genotype", "genotype_complet",
              "genotype_pypgx", "alleles_hors_panel", "phenotype", "allele_externe",
              "heteroplasmie", "fraction", "positions_lues", "positions_attendues", "motif"]
    with open(a.sortie, "w", newline="\n", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=champs, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in resultats:
            w.writerow({k: ("" if r.get(k) is None else r.get(k, "")) for k in champs})
    if a.json:
        with open(a.json, "w", newline="\n", encoding="utf-8") as fh:
            json.dump({"source": cat.get("source"), "profondeur_min": a.profondeur_min,
                       "genes": resultats}, fh, ensure_ascii=False, indent=1)
            fh.write("\n")

    for r in resultats:
        detail = r["genotype"] or r["motif"]
        if r["phenotype"]:
            detail += " -> " + r["phenotype"]
        print("  %-9s %-14s %s" % (r["gene"], r["statut"], detail))
    print("typage complementaire : %d rendu(s) sur %d" % (rendus, len(resultats)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
