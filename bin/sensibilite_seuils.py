# -*- coding: utf-8 -*-
"""Sensibilite du perimetre aux deux seuils de qualite.

Les seuils par defaut, GQ 20 et profondeur 10, sont des choix. Une plateforme
qui met le module en production demande a voir ce qu'ils changent : rejouer le
controle qualite a d'autres seuils, sur les memes donnees, et compter combien de
genes basculent.

Tout se calcule sur les qc_positions.tsv deja produits : aucune relecture
d'alignement, aucune nouvelle execution.
"""
import argparse, collections, csv, glob, json, os, statistics, sys

SEP = chr(9)


def charge(lot):
    """[(echantillon, [positions])] ; une position = (gene, GQ, profondeur, statut)."""
    out = {}
    for f in sorted(glob.glob(os.path.join(lot, "*", "travail", "qc_positions.tsv"))):
        ech = os.path.basename(os.path.dirname(os.path.dirname(f)))
        lignes = []
        for r in csv.DictReader(open(f, encoding="utf-8"), delimiter=SEP):
            genes = [g.strip() for g in (r.get("genes") or "").split(",") if g.strip()]
            gq = r.get("GQ")
            pa = r.get("profondeur_alignement")
            lignes.append((genes, gq, pa, r.get("statut"), r.get("pos")))
        if lignes:
            out[ech] = lignes
    return out


def rejoue(lignes, gq_min, prof_min):
    """Rend, par gene, (attendues, retenues). Reproduit la regle de qc_perimetre.

    Une position est retenue si la profondeur tient le seuil, et, lorsqu'un
    genotype est appele, si sa qualite le tient aussi. Une qualite absente ou
    illisible ecarte la position : un seuil qu'on ne peut pas appliquer n'est pas
    un seuil satisfait.
    """
    att = collections.Counter()
    ret = collections.Counter()
    for genes, gq, pa, statut, _ in lignes:
        for g in genes:
            att[g] += 1
        try:
            prof = int(pa)
        except (TypeError, ValueError):
            continue
        if prof < prof_min:
            continue
        # Le statut du fichier est deja le resultat du filtrage aux seuils par
        # defaut : "qualite insuffisante" designe une position appelee qu'ils
        # ont ecartee. Le seuil se rejoue donc sur le GQ, pas sur le statut.
        if statut in ("mesuree", "qualite insuffisante"):
            try:
                q = int(gq)
            except (TypeError, ValueError):
                continue
            if q < gq_min:
                continue
        elif statut not in ("reference, lue",):
            continue
        for g in genes:
            ret[g] += 1
    return att, ret


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lot", required=True)
    p.add_argument("--sortie", required=True)
    a = p.parse_args()
    os.makedirs(a.sortie, exist_ok=True)

    data = charge(a.lot)
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
                att, ret = rejoue(lignes, gq, pr)
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
    ref = {(20, 10)}
    bascule = collections.Counter()
    for ech, lignes in data.items():
        a20, r20 = rejoue(lignes, 20, 10)
        a30, r30 = rejoue(lignes, 30, 10)
        for g in a20:
            if r20[g] == a20[g] and r30[g] != a30[g]:
                bascule[g] += 1

    # positions qui tombent en premier quand on serre : elles nomment les zones fragiles
    fragiles = collections.Counter()
    for ech, lignes in data.items():
        for genes, gq, pa, statut, pos in lignes:
            if statut not in ("mesuree", "qualite insuffisante"):
                continue
            try:
                q = int(gq)
            except (TypeError, ValueError):
                continue
            if 20 <= q < 40:
                for g in genes:
                    fragiles["%s %s" % (g, pos)] += 1

    res = {"genomes": len(data), "grille": table,
           "genes_qui_basculent_de_20_a_30": bascule.most_common(15),
           "positions_entre_20_et_40": fragiles.most_common(15)}
    json.dump(res, open(os.path.join(a.sortie, "sensibilite_seuils.json"), "w",
                        encoding="utf-8"), ensure_ascii=False, indent=1)

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
