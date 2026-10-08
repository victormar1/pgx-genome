# -*- coding: utf-8 -*-
"""Compare ce que deux executions du module ont rendu.

Le protocole de revalidation demande de rejouer des genomes apres un changement
et de verifier que rien n'a bouge. Sans outil, cette verification se fait a
l'oeil sur des fichiers de plusieurs milliers de lignes, et un ecart s'y perd.

Trois niveaux, du nom vers la conduite a tenir, parce qu'ils ne se valent pas :

  - le diplotype : un nom peut changer sans consequence, par exemple quand une
    ambiguite est tranchee en l'un de ses propres membres ;
  - le phenotype : c'est lui qui porte la classe de fonction ;
  - l'ensemble des recommandations fortes et moderees du document : c'est ce
    que le patient recevrait.

Usage :
    comparer_rapports.py --avant DOSSIER --apres DOSSIER [--perimetre FICHIER]

Chaque dossier contient un sous-dossier par echantillon, comme les sorties du
module. Le code de retour vaut 0 si rien n'a change, 1 sinon : l'outil peut
servir de porte.
"""
import argparse
import json
import os
import sys


def lire(chemin):
    with open(chemin, encoding="utf-8") as fh:
        return json.load(fh)


def rapport(dossier, ech):
    """Le rapport de l'interpreteur, quel que soit le nom qu'il porte."""
    d = os.path.join(dossier, ech, "sortie")
    if not os.path.isdir(d):
        return None
    for nom in sorted(os.listdir(d)):
        if nom.endswith(".report.json"):
            return lire(os.path.join(d, nom))
    return None


def echantillons(dossier):
    if not os.path.isdir(dossier):
        return []
    return sorted(x for x in os.listdir(dossier)
                  if os.path.isdir(os.path.join(dossier, x, "sortie")))


def etat(rep, gene):
    """Diplotype et phenotype rendus pour un gene, sous forme comparable."""
    g = (rep.get("genes") or {}).get(gene) or {}
    x = g.get("recommendationDiplotypes") or g.get("sourceDiplotypes") or []
    return (";".join(sorted({y.get("label", "") for y in x if y.get("label")})),
            ";".join(sorted({p for y in x for p in (y.get("phenotypes") or [])})))


def conduites(rep):
    """Les recommandations fortes et moderees, par medicament.

    Elles ne sont pas attribuees a un gene : une recommandation qui porte sur
    deux genes change de cle quand l'un devient indetermine, et l'attribuer
    ferait compter comme apparue une recommandation qui etait deja la.
    """
    out = {}
    for nom, m in (rep.get("drugs", {}).get("CPIC Guideline Annotation", {})).items():
        for gl in m.get("guidelines") or []:
            for a in gl.get("annotations") or []:
                if a.get("classification") in ("Strong", "Moderate"):
                    out.setdefault(nom, set()).add(
                        (a.get("drugRecommendation") or "").strip())
    return out


def inclus(avant, apres):
    """Un appel unique pris dans l'ambiguite d'avant n'est pas un desaccord :
    les deux cotes disent la meme chose, l'un ayant tranche."""
    a = {x for x in avant.split(";") if x}
    b = {x for x in apres.split(";") if x}
    return len(a) > 1 and len(b) == 1 and b <= a


def compare(avant, apres, perimetre):
    genes, ecarts, conduites_ecarts = 0, [], []
    communs = [e for e in echantillons(apres) if rapport(avant, e)]
    for e in communs:
        a, b = rapport(avant, e), rapport(apres, e)
        for g in perimetre:
            da, pa = etat(a, g)
            db, pb = etat(b, g)
            if not da and not db:
                continue
            genes += 1
            if (da, pa) != (db, pb):
                ecarts.append({"echantillon": e, "gene": g, "avant": da,
                               "apres": db, "phenotype_avant": pa,
                               "phenotype_apres": pb,
                               "ambiguite_tranchee": inclus(da, db) and pa == pb})
        ca, cb = conduites(a), conduites(b)
        app = sorted(set(cb) - set(ca))
        dis = sorted(set(ca) - set(cb))
        chg = sorted(m for m in set(ca) & set(cb) if ca[m] != cb[m])
        if app or dis or chg:
            conduites_ecarts.append({"echantillon": e, "apparus": app,
                                     "disparus": dis, "texte_change": chg})
    return communs, genes, ecarts, conduites_ecarts


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--avant", required=True)
    p.add_argument("--apres", required=True)
    p.add_argument("--perimetre", default="",
                   help="JSON du perimetre clinique ; par defaut, celui du depot")
    p.add_argument("--json", default="", help="ecrit le detail dans ce fichier")
    a = p.parse_args()

    chemin = a.perimetre or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "ressources", "perimetre_rnpgx.json")
    perimetre = sorted(lire(chemin).get("genes", {}))

    communs, genes, ecarts, cond = compare(a.avant, a.apres, perimetre)
    if not communs:
        print("aucun echantillon commun aux deux dossiers", file=sys.stderr)
        return 2

    tranchees = sum(1 for x in ecarts if x["ambiguite_tranchee"])
    print("echantillons compares : %d, couples gene x echantillon : %d"
          % (len(communs), genes))
    print("ecarts de diplotype ou de phenotype : %d, dont %d ambiguites "
          "tranchees sans changement de phenotype" % (len(ecarts), tranchees))
    print("echantillons dont les recommandations changent : %d" % len(cond))
    for x in ecarts:
        print("  %s %s : %s -> %s%s"
              % (x["echantillon"], x["gene"], x["avant"] or "-",
                 x["apres"] or "-",
                 "  (ambiguite tranchee)" if x["ambiguite_tranchee"] else ""))
    for x in cond:
        print("  %s : apparus %s ; disparus %s ; texte change %s"
              % (x["echantillon"], ", ".join(x["apparus"]) or "aucun",
                 ", ".join(x["disparus"]) or "aucun",
                 ", ".join(x["texte_change"]) or "aucun"))
    if a.json:
        with open(a.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump({"echantillons": communs, "couples": genes,
                       "ecarts": ecarts, "recommandations": cond}, fh,
                      ensure_ascii=False, indent=1)
    return 1 if (ecarts or cond) else 0


if __name__ == "__main__":
    sys.exit(main())
