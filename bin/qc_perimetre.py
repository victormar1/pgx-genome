# -*- coding: utf-8 -*-
"""Etage 2b : controle qualite et perimetre reellement mesure.

Un fichier de variants ne distingue pas une position identique a la reference
d une position non lue. L alignement, lui, le dit. Cet etage confronte les deux
et produit le perimetre : la liste des genes dont toutes les positions
diagnostiques ont ete lues a une profondeur et une qualite suffisantes. Ce
perimetre commande ensuite ce que le compte rendu a le droit d affirmer.

Regles, toutes destinees a ne jamais transformer une absence de mesure en
resultat normal :
  - une qualite de genotype absente ou illisible ecarte la position, elle ne la
    laisse pas passer ;
  - un genotype partiellement appele est traite comme non appele ;
  - plusieurs enregistrements a la meme position sont tous juges, et un seul
    ecarte suffit a ecarter la position ;
  - la profondeur est celle de l alignement, jamais celle que declare le
    fichier de variants ;
  - l etage echoue si le fichier porte plus d un echantillon ;
  - un site rejete par le filtrage de l appelant est ecarte si l individu y
    porte un allele alternatif, sauf demande explicite par --ignorer-filtre :
    le filtre juge le variant du site, et ne dit rien de qui ne le porte pas ;
  - un genotype calcule sur moins de lectures que le seuil de profondeur est
    ecarte, meme si l alignement en porte davantage : c est l appel qu on juge ;
  - un heterozygote dont un allele porte moins du quart des lectures est ecarte,
    sauf sur les genes dont le nombre de copies varie.

Les trois dernieres regles viennent de l audit des 141 ecarts du lot de 105
genomes, et chacune a ete mesuree avant d etre ecrite.

Sorties :
  qc_positions.tsv   une ligne par position, ce qui a ete lu et ce qui a ete retenu
  perimetre.json     par gene : attendu, retenu, perdu, statut
  qualifie.vcf       le meme fichier, genotypes insuffisants ramenes a ./., GT seul
"""
import argparse, collections, csv, gzip, json, os, sys


def lire_positions(chemin, complement=False):
    """Positions diagnostiques attendues, depuis un VCF de positions.

    Le fichier complementaire du core panel RNPGx porte en plus PXCLASSE (1 ou 2),
    PXTYPE (snv, indel, vntr), PXINTERP et PXNOTE. Un VNTR n'est pas genotypable en
    lectures courtes : sa couverture est mesuree, son genotype n'est pas juge.
    """
    attendu = {}
    with open(chemin, encoding="utf-8") as fh:
        for l in fh:
            if l.startswith("#"):
                continue
            c = l.rstrip("\n").split("\t")
            info = {}
            for x in c[7].split(";"):
                k, _, v = x.partition("=")
                info[k] = v
            attendu[(c[0], int(c[1]))] = {
                "rsid": c[2], "ref": c[3], "alt": c[4],
                "genes": info.get("PX", "").split(",") if info.get("PX") else [],
                "complement": complement,
                "classe": entier(info.get("PXCLASSE")),
                "type": info.get("PXTYPE", "snv"),
                "interpreteur": info.get("PXINTERP", "pharmcat"),
                "note": info.get("PXNOTE", "").replace("_", " "),
                "fin": entier(info.get("END")) or int(c[1]) + len(c[3]) - 1,
            }
    return attendu


def lire_profondeur(chemin):
    """samtools depth -a : contig, position, profondeur. Une profondeur
    illisible leve une erreur : la sauter ferait d une position mesuree une
    position inconnue, et l etage doit echouer plutot que le taire."""
    prof = {}
    if not chemin or not os.path.exists(chemin):
        return prof
    with open(chemin, encoding="utf-8") as fh:
        for l in fh:
            c = l.split()
            if len(c) >= 3:
                prof[(c[0], int(c[1]))] = int(c[2])
    return prof


def mediane(valeurs):
    """La couverture mediane des positions du panel porte le domaine de validite
    du module : sous le seuil mesure, un allele variant peut etre lu comme
    reference, et l interpreteur l affirme alors sans aucune reserve."""
    v = sorted(valeurs)
    if not v:
        return None
    m = len(v) // 2
    if len(v) % 2:
        return v[m]
    moy = (v[m - 1] + v[m]) / 2.0
    return int(moy) if moy == int(moy) else moy


# Un melange de deux individus fabrique de faux heterozygotes a faible part.
# Sous cette part, un heterozygote appele est dit desequilibre ; en dessous de
# cette profondeur, la part n'a pas de sens ; et sous ce nombre
# d'heterozygotes, une proportion ne s'etablit pas.
MELANGE_PART = 0.30
MELANGE_LECTURES = 10
MELANGE_HETEROZYGOTES_MIN = 20


def indicateur_melange(lignes):
    """Part des heterozygotes appeles dont l'allele mineur est trop peu soutenu.

    La mesure porte sur tout le fichier de variants filtre et non sur les seules
    positions du panel : une vingtaine d'heterozygotes ne suffit pas a etablir
    une proportion, et les regions en portent une cinquantaine.
    """
    het = desq = 0
    for groupe in lignes.values():
        for champs in groupe:
            gt, _, _ = genotype(champs)
            ad = lectures_alleliques(champs)
            if not ad or sum(ad) < MELANGE_LECTURES:
                continue
            # part_minimale exige deux alleles appeles distincts : elle
            # ecarte d'elle-meme un homozygote et un genotype incomplet.
            part = part_minimale(gt, ad)
            if part is None:
                continue
            het += 1
            if part < MELANGE_PART:
                desq += 1
    return het, desq


def suspect_melange(meta, maximum):
    """Au-dela du maximum, l'echantillon n'en est pas un seul.

    Une proportion etablie sur trop peu d'heterozygotes ne vaut rien : elle ne
    declenche alors aucun soupcon, et l'absence de mesure est dite comme telle
    plutot que prise pour un resultat.
    """
    p = meta.get("melange_part")
    n = meta.get("melange_heterozygotes") or 0
    # Un maximum nul leve la porte, et il faut l'ecrire : contrairement a un
    # seuil plancher, qu'une comparaison a zero neutralise d'elle-meme, « au-dela
    # de zero pour cent » serait vrai de presque tout echantillon.
    if not maximum or p is None or n < MELANGE_HETEROZYGOTES_MIN:
        return False
    return p > maximum


def couverture_retenue(meta):
    """Le perimetre clinique porte la decision quand il est connu : c'est sur lui
    que le domaine de validite a ete mesure, et non sur les genes rendus sans
    verite."""
    c = meta.get("couverture_mediane_clinique")
    return meta.get("couverture_mediane") if c is None else c


def hors_domaine(meta, minimum):
    """Sous le seuil, un allele variant peut etre lu comme reference : le module
    sort de son domaine et le genome doit etre refuse. Une couverture inconnue
    n'est pas un refus, c'est un echec de l'etage, deja signale par ailleurs."""
    c = couverture_retenue(meta)
    return c is not None and c < minimum


def verdict(meta):
    """Ce que l orchestration relit de l etage, en neuf champs jamais vides :
    couverture retenue, hors domaine, part de melange, melange suspect,
    positions retenues, positions visees, genes complets, genes, positions
    sans qualite. Un champ vide decalerait les suivants a la lecture."""
    genes = meta.get("genes", {})
    champs = [couverture_retenue(meta), int(bool(meta.get("couverture_hors_domaine"))),
              meta.get("melange_part"), int(bool(meta.get("melange_suspect"))),
              meta["positions_retenues"], meta["positions_visees"],
              sum(1 for v in genes.values() if v["statut"] == "complet"), len(genes),
              meta.get("positions_sans_GQ", 0)]
    return ["-" if c is None else str(c) for c in champs]


def ouvrir(chemin):
    return gzip.open(chemin, "rt", encoding="utf-8") if chemin.endswith(".gz") \
        else open(chemin, encoding="utf-8")


def lire_vcf(chemin):
    """entete, lignes groupees par position, nombre d echantillons."""
    lignes = collections.defaultdict(list)
    entete, n_ech = [], 0
    with ouvrir(chemin) as fh:
        for l in fh:
            if l.startswith("##"):
                entete.append(l)
                continue
            if l.startswith("#CHROM"):
                entete.append(l)
                n_ech = max(0, len(l.rstrip("\n").split("\t")) - 9)
                continue
            c = l.rstrip("\n").split("\t")
            lignes[(c[0], int(c[1]))].append(c)
    return entete, lignes, n_ech


# Part minimale des lectures pour chacun des alleles d un heterozygote. Mesure
# sur 105 genomes : a 0,25, aucun heterozygote du perimetre clinique n est perdu
# hors CYP2D6, et les faux CYP4F2*17 de la region paralogue CYP4F sont a 0,10-0,23.
EQUILIBRE_MIN = 0.25
# Genes dont le nombre de copies varie : une duplication deplace legitimement la
# fraction allelique, jusqu a 0,15 sur un CYP2D6 a onze copies. Leur diplotype
# vient d un outil dedie, pas de ces positions.
SANS_EQUILIBRE = {"CYP2D6"}


def entier(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def genotype(champs):
    """GT, GQ, DP. GQ et DP valent None des qu ils sont absents ou illisibles."""
    if len(champs) < 10:
        return None, None, None
    d = dict(zip(champs[8].split(":"), champs[9].split(":")))
    return d.get("GT"), entier(d.get("GQ")), entier(d.get("DP"))


def lectures_alleliques(champs):
    """Lectures par allele, champ AD. None si absent ou illisible."""
    if not champs or len(champs) < 10:
        return None
    d = dict(zip(champs[8].split(":"), champs[9].split(":")))
    brut = d.get("AD")
    if not brut:
        return None
    try:
        return [int(x) for x in brut.split(",")]
    except ValueError:
        return None


def part_minimale(gt, ad):
    """Part de lectures de l allele le moins soutenu, parmi les alleles appeles.

    None pour un homozygote, ou quand les lectures ne sont pas lisibles : la
    regle ne s applique alors pas, elle ne se devine pas.
    """
    if not ad or not gt:
        return None
    idx = sorted({int(x) for x in gt.replace("|", "/").split("/") if x.isdigit()})
    if len(idx) < 2 or any(i >= len(ad) for i in idx):
        return None
    total = sum(ad)
    if total <= 0:
        return None
    return min(ad[i] for i in idx) / float(total)


def appele(gt):
    """Un genotype n est appele que si tous ses alleles le sont."""
    if not gt or gt in (".", "./.", ".|."):
        return False
    return "." not in gt.replace("|", "/").split("/")


def porte_un_alternatif(gt):
    """Vrai si un genotype appele porte un allele autre que la reference."""
    return any(a != "0" for a in gt.replace("|", "/").split("/"))


# Les statuts qui valent mesure : la position compte comme retenue.
RETENUS = ("mesuree", "reference, lue", "couverture seule")


def choisir(groupe):
    """L enregistrement le plus informatif d une position : celui qui porte un
    genotype appele, sinon le premier."""
    for c in groupe:
        if appele(genotype(c)[0]):
            return c
    return groupe[0] if groupe else None


def classer(info, champs, lu, gq_min, prof_min, ignorer_filtre=False, equilibre=True):
    """Le statut d un enregistrement, ou d une position sans enregistrement.

    `lu` est la profondeur de l alignement, None si elle n a pas ete mesuree :
    celle que declare le fichier de variants ne la remplace pas. `equilibre`
    se leve pour un gene haploide, ou une fraction allelique est une
    heteroplasmie et non un defaut.
    """
    gt, gq, dp = genotype(champs) if champs else (None, None, None)
    if info.get("type") == "vntr" or gt is None:
        # Un VNTR ne se genotype pas en lectures courtes : on mesure sa
        # couverture et on s arrete la. Le pretendre appele serait un resultat
        # invente.
        if lu is None:
            return "profondeur inconnue"
        if lu < prof_min:
            return "non lue"
        return "couverture seule" if info.get("type") == "vntr" else "reference, lue"
    if not appele(gt):
        return "non appelee"
    if gq is None:
        return "qualite non mesurable"
    if gq < gq_min:
        return "qualite insuffisante"
    if lu is None:
        return "profondeur inconnue"
    if lu < prof_min:
        return "profondeur insuffisante"
    # Le filtre de l appelant est ecrit une fois par site : il juge le variant,
    # pas l individu. Il ne compte donc que pour qui porte un allele
    # alternatif. Un fichier de cohorte donne une ligne a chacun la ou un seul
    # varie, et le filtre du variant d un autre ecarterait une reference bien
    # lue. La qualite et les lectures de l individu restent jugees.
    filtre = (champs[6] if len(champs) > 6 else "") or ""
    if not ignorer_filtre and filtre not in (".", "PASS", "") and porte_un_alternatif(gt):
        return "site filtre"
    # la profondeur de l appel est celle sur laquelle le genotype a ete
    # calcule ; l alignement peut en porter davantage que l appelant n a lu
    ad = lectures_alleliques(champs)
    dp_appel = sum(ad) if ad else dp
    if dp_appel is not None and dp_appel < prof_min:
        return "profondeur d appel insuffisante"
    if equilibre and not any(g in SANS_EQUILIBRE for g in info["genes"]):
        part = part_minimale(gt, ad)
        if part is not None and part < EQUILIBRE_MIN:
            return "desequilibre allelique"
    return "mesuree"


def juger(info, groupe, lu, gq_min, prof_min, ignorer_filtre=False, equilibre=True):
    """Le statut d une position, l enregistrement qui le porte, et s il faut
    masquer ses genotypes.

    Tous les enregistrements sont juges et le premier ecarte donne son statut a
    la position : l interpreteur les lit tous, un seul douteux suffit a fausser
    son appel. Une position ecartee ne lui laisse aucun genotype appele, y
    compris ceux de ses enregistrements qui, seuls, auraient ete retenus.
    """
    for c in groupe:
        statut = classer(info, c, lu, gq_min, prof_min, ignorer_filtre, equilibre)
        if statut not in RETENUS:
            return statut, c, (info.get("type") != "vntr"
                               and any(appele(genotype(x)[0]) for x in groupe))
    c = choisir(groupe)
    return classer(info, c, lu, gq_min, prof_min, ignorer_filtre, equilibre), c, False


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--profondeur", default=None)
    p.add_argument("--positions", required=True)
    p.add_argument("--positions-complement", default=None, dest="complement",
                   help="VCF des positions RNPGx absentes des definitions PharmCAT")
    p.add_argument("--sortie", required=True)
    p.add_argument("--gq", type=int, default=20)
    p.add_argument("--profondeur-min", type=int, default=10, dest="prof_min")
    p.add_argument("--couverture-min", type=int, default=0, dest="couv_min",
                   help="couverture mediane en deca de laquelle le genome sort du domaine")
    p.add_argument("--melange-max", type=float, default=0.0, dest="melange_max",
                   help="part d'heterozygotes desequilibres au-dela de laquelle "
                        "l'echantillon est tenu pour melange")
    p.add_argument("--echantillon", default="")
    p.add_argument("--ignorer-filtre", action="store_true", dest="ignorer_filtre",
                   help="accepte les sites que l appelant a rejetes")
    p.add_argument("--perimetre-clinique", default=None, dest="clinique",
                   help="JSON des genes du perimetre clinique a distinguer")
    a = p.parse_args()

    os.makedirs(a.sortie, exist_ok=True)
    attendu = lire_positions(a.positions)
    if not attendu:
        print("fichier de positions vide ou illisible :", a.positions, file=sys.stderr)
        return 1
    if a.complement:
        if not os.path.exists(a.complement):
            print("fichier complementaire introuvable :", a.complement, file=sys.stderr)
            return 1
        comp = lire_positions(a.complement, complement=True)
        if not comp:
            print("fichier complementaire vide :", a.complement, file=sys.stderr)
            return 1
        # une position deja suivie par PharmCAT reste sous son propre regime
        doublons = sorted(set(comp) & set(attendu))
        attendu.update({k: v for k, v in comp.items() if k not in attendu})
        if doublons:
            print("positions complementaires deja suivies, ignorees : %d" % len(doublons))
    clinique = set()
    if a.clinique:
        if not os.path.exists(a.clinique):
            print("perimetre clinique introuvable :", a.clinique, file=sys.stderr)
            return 1
        with open(a.clinique, encoding="utf-8") as fh:
            clinique = set(json.load(fh).get("genes", {}))
    prof = lire_profondeur(a.profondeur)
    if not prof:
        print("profondeur absente : le perimetre ne peut pas etre etabli", file=sys.stderr)
        return 1
    entete, lignes, n_ech = lire_vcf(a.vcf)
    if n_ech != 1:
        print("le fichier de variants porte %d echantillons, un seul est accepte" % n_ech,
              file=sys.stderr)
        return 1

    lignes_qc, complement_qc = [], []
    couv_toutes, couv_clinique = [], []
    par_gene = collections.defaultdict(collections.Counter)
    a_masquer = set()
    sans_gq = 0

    for cle in sorted(attendu, key=lambda k: (k[0], k[1])):
        info = attendu[cle]
        groupe = lignes.get(cle, [])
        lu = prof.get(cle)
        statut, champs, masquer = juger(info, groupe, lu, a.gq, a.prof_min, a.ignorer_filtre)
        gt, gq, dp = genotype(champs) if champs else (None, None, None)
        if masquer:
            a_masquer.add(cle)
        sans_gq += statut == "qualite non mesurable"
        retenue = statut in RETENUS
        lignes_qc.append({
            "chrom": cle[0], "pos": cle[1], "rsid": info["rsid"],
            "ref": info["ref"], "alt": info["alt"], "genes": ",".join(info["genes"]),
            "enregistrements": len(groupe),
            "genotype": gt or "", "GQ": "" if gq is None else gq,
            "profondeur_alignement": "" if lu is None else lu,
            "profondeur_VCF": "" if dp is None else dp,
            "statut": statut, "retenue": "oui" if retenue else "non",
            "source": "RNPGx" if info.get("complement") else "PharmCAT",
            "classe_rnpgx": info.get("classe") or "",
            # De quoi rejouer la decision a d autres seuils sans relire le
            # fichier de variants : c est l enregistrement qui a decide.
            "filtre": champs[6] if champs and len(champs) > 6 else "",
            "lectures_alleliques": ",".join(map(str, lectures_alleliques(champs) or [])),
        })
        if info.get("complement"):
            genotype_rendu = ("non génotypé (VNTR)" if info.get("type") == "vntr"
                              else (gt if retenue and gt else ("référence" if statut == "reference, lue" else "")))
            complement_qc.append({
                "gene": info["genes"][0] if info["genes"] else "", "rsid": info["rsid"],
                "chrom": cle[0], "pos": cle[1], "classe": info.get("classe"),
                "type": info.get("type"), "interpreteur": info.get("interpreteur"),
                "genotype": genotype_rendu, "statut": statut,
                "profondeur": "" if lu is None else lu, "note": info.get("note", ""),
            })
        # Le statut d'un gene ne compte que les positions dont l'interpreteur se sert.
        # Une position complementaire ratee ne doit pas faire passer le gene en
        # « partiel » : elle est suivie a part, dans complement_rnpgx.
        if not info.get("complement"):
            # La couverture est celle de l alignement, et non celle du VCF : elle
            # existe pour une position sur laquelle l appelant n a rien ecrit,
            # qui est justement le cas ou le module se trompe.
            couv_toutes.append(lu or 0)
            if clinique & set(info["genes"]):
                couv_clinique.append(lu or 0)
            for g in info["genes"]:
                par_gene[g]["attendu"] += 1
                par_gene[g]["retenue" if retenue else "perdue"] += 1
                par_gene[g][statut] += 1

    # Une position complementaire n'est pas utilisee par l'interpreteur : la masquer
    # ne changerait rien a son appel, mais modifierait son entree. On la laisse donc
    # intacte, et le VCF qualifie reste identique a ce qu'il serait sans complement.
    a_masquer -= {c for c in a_masquer if attendu[c].get("complement")}

    with open(os.path.join(a.sortie, "qc_positions.tsv"), "w", newline="\n", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(lignes_qc[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(lignes_qc)

    perimetre = {}
    for g, c in sorted(par_gene.items()):
        att, ret = c["attendu"], c["retenue"]
        perimetre[g] = {
            "positions_attendues": att,
            "positions_retenues": ret,
            "positions_perdues": att - ret,
            "statut": "complet" if ret == att else ("partiel" if ret else "absent"),
            "clinique": (g in clinique) if clinique else None,
            "detail": {k: v for k, v in c.items() if k not in ("attendu", "retenue", "perdue")},
        }
    meta = {
        "echantillon": a.echantillon,
        "seuil_GQ": a.gq,
        "seuil_profondeur": a.prof_min,
        "positions_visees": sum(1 for x in lignes_qc if x["source"] == "PharmCAT"),
        "positions_retenues": sum(1 for x in lignes_qc if x["source"] == "PharmCAT" and x["retenue"] == "oui"),
        "positions_visees_toutes": len(attendu),
        "positions_sans_GQ": sans_gq,
        "positions_a_plusieurs_enregistrements": sum(1 for x in lignes_qc if x["enregistrements"] > 1),
        "profondeur_lue_depuis_alignement": True,
        "couverture_mediane": mediane(couv_toutes),
        "couverture_mediane_clinique": mediane(couv_clinique),
        "perimetre_clinique": sorted(clinique) if clinique else [],
        "perimetre_clinique_hors_interpreteur": sorted(clinique - set(perimetre)) if clinique else [],
        "genes": perimetre,
        "complement_rnpgx": sorted(complement_qc, key=lambda x: (x["gene"], x["pos"])),
        "complement_rnpgx_resume": {
            "positions": len(complement_qc),
            "retenues": sum(1 for x in complement_qc if x["statut"] in RETENUS),
            "classe_1": sum(1 for x in complement_qc if x["classe"] == 1),
            "classe_2": sum(1 for x in complement_qc if x["classe"] == 2),
        } if complement_qc else None,
    }
    het, desq = indicateur_melange(lignes)
    meta["melange_heterozygotes"] = het
    meta["melange_desequilibrees"] = desq
    meta["melange_part"] = round(100.0 * desq / het, 1) if het else None
    meta["melange_part_max"] = a.melange_max
    meta["melange_suspect"] = suspect_melange(meta, a.melange_max)
    meta["couverture_mediane_retenue"] = couverture_retenue(meta)
    meta["couverture_mediane_min"] = a.couv_min
    meta["couverture_hors_domaine"] = hors_domaine(meta, a.couv_min)
    with open(os.path.join(a.sortie, "perimetre.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(a.sortie, "verdict.tsv"), "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\t".join(verdict(meta)) + "\n")

    # VCF qualifie : tous les enregistrements sont reecrits, GT seul.
    # Certains champs de format font echouer le pretraitement sans message clair,
    # et l interpreteur ne lit que le genotype.
    lues = sum(len(v) for v in lignes.values())
    ecrites = 0
    sortie_vcf = os.path.join(a.sortie, "qualifie.vcf")
    with open(sortie_vcf, "w", newline="\n", encoding="utf-8") as fh:
        for l in entete:
            fh.write(l)
        for cle in sorted(lignes, key=lambda k: (k[0], k[1])):
            for c in lignes[cle]:
                c = list(c)
                if len(c) >= 10:
                    gt = dict(zip(c[8].split(":"), c[9].split(":"))).get("GT", "./.")
                    if cle in a_masquer:
                        gt = "./." if "/" in gt or "|" in gt else "."
                    c[8], c[9] = "GT", gt
                fh.write("\t".join(c) + "\n")
                ecrites += 1
    if ecrites != lues:
        print("perte d enregistrements : %d lus, %d ecrits" % (lues, ecrites), file=sys.stderr)
        return 1

    complets = sum(1 for v in perimetre.values() if v["statut"] == "complet")
    print("positions visees %d, retenues %d, sans qualite %d"
          % (meta["positions_visees"], meta["positions_retenues"], sans_gq))
    print("genes complets %d sur %d" % (complets, len(perimetre)))
    for g, v in sorted(perimetre.items()):
        if v["statut"] != "complet":
            print("  %-10s %s : %d positions perdues sur %d"
                  % (g, v["statut"], v["positions_perdues"], v["positions_attendues"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
