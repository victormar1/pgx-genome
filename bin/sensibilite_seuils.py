# -*- coding: utf-8 -*-
"""Sensibilite du perimetre aux deux seuils de qualite.

Les seuils par defaut, GQ 20 et profondeur 10, sont des choix. Une plateforme
qui met le module en production demande a voir ce qu'ils changent : rejouer le
controle qualite a d'autres seuils, sur les memes donnees, et compter combien de
genes basculent.

Tout se calcule sur les qc_positions.tsv deja produits : aucune relecture
d'alignement, aucune nouvelle execution. La decision est celle de l'etage 2b,
par sa propre fonction : le rejeu aux seuils par defaut rend le perimetre que
l'etage a rendu. Seule limite : une position a plusieurs enregistrements est
rejouee sur celui qui avait decide.
"""
import argparse, collections, csv, glob, json, os, sys

from qc_perimetre import RETENUS, classer, entier, genotype

SEP = chr(9)


def position(genes, gt, gq, dp, ad, filtre, lu, pos=""):
    """Une ligne de qc_positions.tsv, remise sous la forme que juge l'etage 2b :
    (genes, description de la position, enregistrement ou None, profondeur lue,
    position). Une valeur absente ou illisible le reste : c'est a la regle de
    dire ce qu'elle en fait."""
    champs = None
    if gt:
        champs = ["", str(pos), "", "", "", "", filtre or "", "", "GT:GQ:DP:AD",
                  ":".join(str(x) if x not in (None, "") else "." for x in (gt, gq, dp, ad))]
    return genes, {"genes": genes}, champs, entier(lu), pos


def charge(lot):
    """{echantillon: [position]}, sur les positions dont l'interpreteur se sert.

    Les positions complementaires sont suivies a part par l'etage et ne pesent
    pas sur le statut d'un gene : elles ne pesent pas davantage ici. Ce sont
    aussi les seules a porter un VNTR, dont le genotype ne se juge pas."""
    out = {}
    for f in sorted(glob.glob(os.path.join(lot, "*", "travail", "qc_positions.tsv"))):
        ech = os.path.basename(os.path.dirname(os.path.dirname(f)))
        lignes = []
        with open(f, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter=SEP):
                if "lectures_alleliques" not in r:
                    raise ValueError("%s vient d'une version anterieure du module et ne porte "
                                     "pas de quoi rejouer la decision : relancer l'etage qc" % f)
                if r.get("source") == "RNPGx":
                    continue
                genes = [g.strip() for g in (r.get("genes") or "").split(",") if g.strip()]
                lignes.append(position(genes, r.get("genotype"), r.get("GQ"),
                                       r.get("profondeur_VCF"), r.get("lectures_alleliques"),
                                       r.get("filtre"), r.get("profondeur_alignement"),
                                       r.get("pos")))
        if lignes:
            out[ech] = lignes
    return out


def rejoue(lignes, gq_min, prof_min, ignorer_filtre=False):
    """Rend, par gene, (attendues, retenues) aux seuils donnes."""
    att = collections.Counter()
    ret = collections.Counter()
    for genes, info, champs, lu, _ in lignes:
        retenue = classer(info, champs, lu, gq_min, prof_min, ignorer_filtre) in RETENUS
        for g in genes:
            att[g] += 1
            ret[g] += retenue
    return att, ret


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lot", required=True)
    p.add_argument("--sortie", required=True)
    p.add_argument("--ignorer-filtre", action="store_true", dest="ignorer_filtre",
                   help="a passer si le lot a ete produit avec cette option")
    a = p.parse_args()
    os.makedirs(a.sortie, exist_ok=True)

    try:
        data = charge(a.lot)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    print("genomes lus : %d" % len(data))
    if not data:
        return 1

    GQ = [0, 10, 20, 30, 40]
    PROF = [5, 10, 15, 20]
    table = []
    for gq in GQ:
        for pr in PROF:
            complets = partiels = absents = 0
            perdues = 0
            for ech, lignes in data.items():
                att, ret = rejoue(lignes, gq, pr, a.ignorer_filtre)
                for g in att:
                    if ret[g] == att[g]:
                        complets += 1
                    elif ret[g] == 0:
                        absents += 1
                    else:
                        partiels += 1
                    perdues += att[g] - ret[g]
            tot = complets + partiels + absents
            table.append({"gq": gq, "profondeur": pr, "genes": tot,
                          "complets": complets, "partiels": partiels, "absents": absents,
                          "part_complets": round(100 * complets / tot, 2) if tot else None,
                          "positions_perdues": perdues})

    # quels genes basculent le plus quand on serre le seuil de qualite
    bascule = collections.Counter()
    for ech, lignes in data.items():
        a20, r20 = rejoue(lignes, 20, 10, a.ignorer_filtre)
        a30, r30 = rejoue(lignes, 30, 10, a.ignorer_filtre)
        for g in a20:
            if r20[g] == a20[g] and r30[g] != a30[g]:
                bascule[g] += 1

    # positions qui tombent en premier quand on serre : elles nomment les zones fragiles
    fragiles = collections.Counter()
    for ech, lignes in data.items():
        for genes, info, champs, lu, pos in lignes:
            gq = genotype(champs)[1] if champs else None
            if gq is None or not 20 <= gq < 40:
                continue
            if classer(info, champs, lu, 20, 10, a.ignorer_filtre) == "mesuree":
                for g in genes:
                    fragiles["%s %s" % (g, pos)] += 1

    res = {"genomes": len(data), "grille": table,
           "genes_qui_basculent_de_20_a_30": bascule.most_common(15),
           "positions_entre_20_et_40": fragiles.most_common(15)}
    with open(os.path.join(a.sortie, "sensibilite_seuils.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)

    print()
    print("%-5s %-6s %8s %9s %8s %10s" % ("GQ", "prof", "complets", "partiels", "absents", "part %"))
    for t in table:
        marque = "  <- defaut" if (t["gq"], t["profondeur"]) == (20, 10) else ""
        print("%-5d %-6d %8d %9d %8d %9.2f%s"
              % (t["gq"], t["profondeur"], t["complets"], t["partiels"], t["absents"],
                 t["part_complets"], marque))
    print()
    print("genes qui perdent leur statut complet en passant de GQ 20 a GQ 30 :")
    for g, n in bascule.most_common(10):
        print("  %-10s %d genomes" % (g, n))
    print()
    print("positions mesurees dont la qualite est entre 20 et 40 (les premieres a tomber) :")
    for k, n in fragiles.most_common(10):
        print("  %-22s %d genomes" % (k, n))
    print()
    print("ecrit :", os.path.join(a.sortie, "sensibilite_seuils.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
