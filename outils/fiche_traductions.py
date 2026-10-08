# -*- coding: utf-8 -*-
"""Produit la fiche de relecture des traductions des recommandations CPIC.

Les traductions sont declarees « proposition, a valider » : la validation
demande un pharmacologue, pas un outil. Ce que l'outil peut faire, c'est la
rendre possible — mettre chaque texte d'origine en face de sa traduction, dire
quel medicament et quel gene elle concerne, et combien de comptes rendus s'en
servent, pour que la relecture commence par ce qui est effectivement rendu.

Il releve aussi deux incoherences que personne ne voit a l'oeil :

  - un texte present dans les rapports et absent de la table : le compte rendu
    le laisse en anglais, ce qui est le repli voulu, mais il faut le savoir ;
  - une traduction que plus aucun rapport n'emploie : elle a survecu a un
    changement de referentiel et sera relue pour rien.

Usage :
    fiche_traductions.py --sortie FICHE.md [--lot DOSSIER]

Sans --lot, la fiche porte la table seule, sans comptage d'usage.
"""
import argparse
import collections
import json
import os
import sys

# La cle de la table est un texte normalise : espaces reduits et entites HTML
# decodees. Comparer des chaines brutes faisait annoncer comme non traduits des
# textes que le compte rendu rend bel et bien en francais. On emploie donc la
# fonction du module, et non une seconde copie de la regle.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
from cr_logique import norm  # noqa: E402


def lire(chemin):
    with open(chemin, encoding="utf-8") as fh:
        return json.load(fh)


def rapports(dossier):
    """Les rapports d'interpreteur d'un dossier de lot, un par echantillon."""
    if not os.path.isdir(dossier):
        return
    for ech in sorted(os.listdir(dossier)):
        d = os.path.join(dossier, ech, "sortie")
        if not os.path.isdir(d):
            continue
        for nom in sorted(os.listdir(d)):
            if nom.endswith(".report.json"):
                try:
                    yield ech, lire(os.path.join(d, nom))
                except ValueError:
                    pass
                break


def usages(dossier, restitues=None):
    """Pour chaque texte d'origine : echantillons, medicaments, genes, classes.

    Le comptage ne porte que sur les medicaments que le compte rendu restitue.
    Les autres existent dans le rapport de l'interpreteur sans jamais atteindre
    le document : les compter ferait annoncer des dizaines de textes « rendus et
    non traduits » qui ne sont pas rendus du tout.
    """
    u = collections.defaultdict(
        lambda: {"echantillons": set(), "medicaments": set(), "genes": set(),
                 "classes": set()})
    for ech, rep in rapports(dossier):
        for nom, m in (rep.get("drugs", {})
                       .get("CPIC Guideline Annotation", {})).items():
            if restitues is not None and nom.lower() not in restitues:
                continue
            for gl in m.get("guidelines") or []:
                for a in gl.get("annotations") or []:
                    if a.get("classification") not in ("Strong", "Moderate"):
                        continue
                    src = norm(a.get("drugRecommendation") or "")
                    if not src:
                        continue
                    e = u[src]
                    e["echantillons"].add(ech)
                    e["medicaments"].add(nom)
                    e["classes"].add(a.get("classification"))
                    for cle in (a.get("lookupKey") or []):
                        for g in cle:
                            e["genes"].add(g)
    return u


def tronque(t, n=300):
    t = " ".join((t or "").split())
    return t if len(t) <= n else t[:n - 1] + "…"


def fiche(table, u):
    """La fiche, par ordre d'usage decroissant : on relit d'abord ce qui sert."""
    tr = table.get("traductions", {})
    lignes = []
    for src, cible in tr.items():
        fr = cible.get("fr", "") if isinstance(cible, dict) else cible
        typ = cible.get("type", "") if isinstance(cible, dict) else ""
        e = u.get(norm(src))
        lignes.append({
            "source": src, "fr": fr, "type": typ,
            "echantillons": len(e["echantillons"]) if e else 0,
            "medicaments": sorted(e["medicaments"]) if e else [],
            "genes": sorted(e["genes"]) if e else [],
            "classes": sorted(e["classes"]) if e else [],
        })
    lignes.sort(key=lambda x: (-x["echantillons"], x["source"]))
    absents = sorted(set(u) - {norm(s) for s in tr})
    inutilisees = [x["source"] for x in lignes if x["echantillons"] == 0]
    return lignes, absents, inutilisees


def ecris(chemin, table, lignes, absents, inutilisees, lot):
    o = ["# Fiche de relecture des traductions CPIC", "",
         "Version de la table : **%s**, etablie le %s, statut **%s**."
         % (table.get("version", "?"), table.get("etabli_le", "?"),
            table.get("statut", "?")), "",
         "Chaque ligne porte le texte d'origine du referentiel et la traduction "
         "proposee. La colonne *visa* est a remplir par le relecteur ; une "
         "traduction sans visa reste une proposition.", ""]
    if lot:
        o += ["Les lignes sont classees par nombre de comptes rendus qui les "
              "emploient : la relecture commence par ce qui est effectivement "
              "rendu.", ""]
    o += ["| Comptes rendus | Medicaments | Genes | Classe | Texte d'origine | "
          "Traduction proposee | Type | Visa |",
          "|---|---|---|---|---|---|---|---|"]
    for x in lignes:
        o.append("| %d | %s | %s | %s | %s | %s | %s |  |"
                 % (x["echantillons"], ", ".join(x["medicaments"]) or "—",
                    ", ".join(x["genes"]) or "—",
                    ", ".join(x["classes"]) or "—",
                    tronque(x["source"]), tronque(x["fr"]), x["type"] or "—"))
    if absents:
        o += ["", "## Textes rendus et non traduits", "",
              "Le compte rendu les laisse en anglais, ce qui est le repli voulu.",
              ""]
        o += ["- %s" % tronque(s) for s in absents]
    if inutilisees:
        o += ["", "## Traductions qu'aucun compte rendu du lot n'emploie", "",
              "Une entree inemployee ne signale pas une erreur : elle porte une "
              "recommandation qu'aucun genome du banc ne declenche. Deux causes "
              "attendues, et une a surveiller.", "",
              "- l'allele qui la declenche est absent de la cohorte — c'est le "
              "cas de l'avertissement sur les aminosides, faute de porteur de "
              "MT-RNR1 ;",
              "- le medicament est hors de la liste restituee, par exemple une "
              "statine ecartee du perimetre ;",
              "- ou bien le texte d'origine a change dans le referentiel, et "
              "l'entree ne correspond plus a rien. C'est la seule des trois qui "
              "demande une correction.", ""]
        o += ["- %s" % tronque(s) for s in inutilisees]
    with open(chemin, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(o) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--sortie", required=True)
    p.add_argument("--lot", default="", help="dossier d'un lot, pour compter l'usage")
    p.add_argument("--table", default="",
                   help="table des traductions ; par defaut celle du depot")
    a = p.parse_args()

    chemin = a.table or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "ressources", "traductions_cpic_fr.json")
    table = lire(chemin)
    restitues = {x.lower() for x in (table.get("medicaments_restitues") or {})}
    u = usages(a.lot, restitues or None) if a.lot else {}
    lignes, absents, inutilisees = fiche(table, u)
    ecris(a.sortie, table, lignes, absents, inutilisees, bool(a.lot))
    print("%d traductions, %d rendues au moins une fois, %d textes non traduits, "
          "%d traductions inemployees"
          % (len(lignes), sum(1 for x in lignes if x["echantillons"]),
             len(absents), len(inutilisees)))
    print("fiche ecrite : %s" % a.sortie)
    return 0


if __name__ == "__main__":
    sys.exit(main())
