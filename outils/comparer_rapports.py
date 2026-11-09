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
    """Les recommandations fortes et moderees, par medicament, avec leur force :
    une consigne qui passe de forte a moderee change le compte rendu, meme a
    texte egal.

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
                        (a.get("classification"),
                         (a.get("drugRecommendation") or "").strip()))
    return out


def inclus(avant, apres):
    """Un appel unique pris dans l'ambiguite d'avant n'est pas un desaccord :
    les deux cotes disent la meme chose, l'un ayant tranche."""
    a = {x for x in avant.split(";") if x}
    b = {x for x in apres.split(";") if x}
    return len(a) > 1 and len(b) == 1 and b <= a


def rendus_sans_reference(avant, apres):
    """Les genomes que la reference n'a pas rendus et que l'execution comparee
    rend : un refus ou un echec qui se leve change le comportement du module.
    Un genome absent de la reference reste ignore, comme dans compare."""
    connus = set(echantillons(avant))
    return [e for e in echantillons(apres)
            if e in connus and not rapport(avant, e) and rapport(apres, e)]


def compare(avant, apres, perimetre):
    genes, ecarts, conduites_ecarts, sans_rapport, communs = 0, [], [], [], []
    for e in echantillons(apres):
        a = rapport(avant, e)
        if not a:
            continue
        b = rapport(apres, e)
        if not b:
            # Le genome a ete repasse et n'a rien rendu : refus, echec ou
            # interruption. C'est le plus grave des ecarts, pas un absent.
            sans_rapport.append(e)
            continue
        communs.append(e)
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
    return communs, genes, ecarts, conduites_ecarts, sans_rapport


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

    communs, genes, ecarts, cond, sans = compare(a.avant, a.apres, perimetre)
    leves = rendus_sans_reference(a.avant, a.apres)
    if not communs and not sans and not leves:
        print("aucun echantillon commun aux deux dossiers", file=sys.stderr)
        return 2

    tranchees = sum(1 for x in ecarts if x["ambiguite_tranchee"])
    print("echantillons compares : %d, couples gene x echantillon : %d"
          % (len(communs), genes))
    print("ecarts de diplotype ou de phenotype : %d, dont %d ambiguites "
          "tranchees sans changement de phenotype" % (len(ecarts), tranchees))
    print("echantillons dont les recommandations changent : %d" % len(cond))
    print("echantillons repasses sans rapport : %d%s"
          % (len(sans), " (" + ", ".join(sans) + ")" if sans else ""))
    print("echantillons rendus sans rapport de reference : %d%s"
          % (len(leves), " (" + ", ".join(leves) + ")" if leves else ""))
    for x in ecarts:
        print("  %s %s : %s -> %s%s"
              % (x["echantillon"], x["gene"], x["avant"] or "-",
                 x["apres"] or "-",
                 "  (ambiguite tranchee)" if x["ambiguite_tranchee"] else ""))
    for x in cond:
        print("  %s : apparus %s ; disparus %s ; texte ou force change %s"
              % (x["echantillon"], ", ".join(x["apparus"]) or "aucun",
                 ", ".join(x["disparus"]) or "aucun",
                 ", ".join(x["texte_change"]) or "aucun"))
    if a.json:
        with open(a.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump({"echantillons": communs, "couples": genes,
                       "ecarts": ecarts, "recommandations": cond,
                       "sans_rapport": sans,
                       "rendus_sans_reference": leves}, fh,
                      ensure_ascii=False, indent=1)
    return 1 if (ecarts or cond or sans or leves) else 0


if __name__ == "__main__":
    sys.exit(main())
