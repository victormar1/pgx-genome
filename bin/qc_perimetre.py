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
  - plusieurs enregistrements a la meme position sont tous conserves ;
  - l etage echoue si le fichier porte plus d un echantillon ;
  - un site rejete par le filtrage de l appelant est ecarte, sauf demande
    explicite par --ignorer-filtre ;
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
    prof = {}
    if not chemin or not os.path.exists(chemin):
        return prof
    with open(chemin, encoding="utf-8") as fh:
        for l in fh:
            c = l.split()
            if len(c) >= 3:
                prof[(c[0], int(c[1]))] = int(c[2])
    return prof


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
    if a.clinique and os.path.exists(a.clinique):
        clinique = set(json.load(open(a.clinique, encoding="utf-8")).get("genes", {}))
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
    par_gene = collections.defaultdict(collections.Counter)
    a_masquer = set()
    sans_gq = 0

    for cle in sorted(attendu, key=lambda k: (k[0], k[1])):
        info = attendu[cle]
        groupe = lignes.get(cle, [])
        lu = prof.get(cle)
        # on retient l enregistrement le plus informatif : celui qui porte un
        # genotype appele, sinon le premier
        champs = None
        for c in groupe:
            gt, _, _ = genotype(c)
            if appele(gt):
                champs = c
                break
        if champs is None and groupe:
            champs = groupe[0]
        gt, gq, dp = genotype(champs) if champs else (None, None, None)
        profondeur = lu if lu is not None else dp
        filtre = (champs[6] if champs and len(champs) > 6 else "") or ""
        ad = lectures_alleliques(champs)
        # la profondeur de l appel est celle sur laquelle le genotype a ete
        # calcule ; l alignement peut en porter davantage que l appelant n a lu
        dp_appel = sum(ad) if ad else dp
        gene_sans_equilibre = any(g in SANS_EQUILIBRE for g in info["genes"])
        part = None if gene_sans_equilibre else (
            part_minimale(gt, ad) if gt and appele(gt) else None)

        if info.get("type") == "vntr":
            # Un VNTR ne se genotype pas en lectures courtes : on mesure sa couverture
            # et on s'arrete la. Le pretendre appele serait un resultat invente.
            statut = ("couverture seule" if profondeur is not None and profondeur >= a.prof_min
                      else ("non lue" if profondeur is not None else "profondeur inconnue"))
        elif gt is None:
            if profondeur is not None and profondeur >= a.prof_min:
                statut = "reference, lue"
            elif profondeur is not None:
                statut = "non lue"
            else:
                statut = "profondeur inconnue"
        elif not appele(gt):
            statut = "non appelee"
        elif gq is None:
            statut = "qualite non mesurable"
            sans_gq += 1
            a_masquer.add(cle)
        elif gq < a.gq:
            statut = "qualite insuffisante"
            a_masquer.add(cle)
        elif profondeur is None:
            statut = "profondeur inconnue"
            a_masquer.add(cle)
        elif profondeur < a.prof_min:
            statut = "profondeur insuffisante"
            a_masquer.add(cle)
        elif not a.ignorer_filtre and filtre not in (".", "PASS", ""):
            statut = "site filtre"
            a_masquer.add(cle)
        elif dp_appel is not None and dp_appel < a.prof_min:
            statut = "profondeur d appel insuffisante"
            a_masquer.add(cle)
        elif part is not None and part < EQUILIBRE_MIN:
            statut = "desequilibre allelique"
            a_masquer.add(cle)
        else:
            statut = "mesuree"

        retenue = statut in ("mesuree", "reference, lue", "couverture seule")
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
        })
        if info.get("complement"):
            genotype_rendu = ("non génotypé (VNTR)" if info.get("type") == "vntr"
                              else (gt if retenue and gt else ("référence" if statut == "reference, lue" else "")))
            complement_qc.append({
                "gene": info["genes"][0] if info["genes"] else "", "rsid": info["rsid"],
                "chrom": cle[0], "pos": cle[1], "classe": info.get("classe"),
                "type": info.get("type"), "interpreteur": info.get("interpreteur"),
                "genotype": genotype_rendu, "statut": statut,
                "profondeur": "" if profondeur is None else profondeur, "note": info.get("note", ""),
            })
        # Le statut d'un gene ne compte que les positions dont l'interpreteur se sert.
        # Une position complementaire ratee ne doit pas faire passer le gene en
        # « partiel » : elle est suivie a part, dans complement_rnpgx.
        if not info.get("complement"):
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
        "perimetre_clinique": sorted(clinique) if clinique else [],
        "perimetre_clinique_hors_interpreteur": sorted(clinique - set(perimetre)) if clinique else [],
        "genes": perimetre,
        "complement_rnpgx": sorted(complement_qc, key=lambda x: (x["gene"], x["pos"])),
        "complement_rnpgx_resume": {
            "positions": len(complement_qc),
            "retenues": sum(1 for x in complement_qc if x["statut"] in ("mesuree", "reference, lue", "couverture seule")),
            "classe_1": sum(1 for x in complement_qc if x["classe"] == 1),
            "classe_2": sum(1 for x in complement_qc if x["classe"] == 2),
        } if complement_qc else None,
    }
    with open(os.path.join(a.sortie, "perimetre.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=1)

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
