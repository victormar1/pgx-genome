# -*- coding: utf-8 -*-
"""Sorties structurees : le genotype d'un cote, l'annotation de l'autre.

Trois tables, a plat, qui se rejoignent par l'identifiant HGNC du gene :

  <ECH>.genotypes.tsv    ce qui a ete mesure, un gene par ligne
  <ECH>.variants.tsv     les variants observes qui fondent ces genotypes,
                         avec leur rsID quand il existe
  <ECH>.annotations.tsv  ce que le referentiel en dit : phenotype et consignes

Le genotype ne change pas ; l'annotation suit la version du referentiel et se
refait sans reprendre la mesure. Les separer permet de conserver l'un et de
remplacer l'autre.

Rien n'est decide ici : le perimetre, les reserves, les ambiguites et les
consignes viennent de `cr_logique`, que lit aussi le compte rendu. Les deux
sorties ne peuvent donc pas se contredire.
"""
import argparse
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cr_logique as L  # noqa: E402

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLONNES = {
    "genotypes": ["hgnc_id", "gene", "genotype", "statut", "motif", "source",
                  "positions_attendues", "positions_retenues"],
    "variants": ["hgnc_id", "gene", "chrom", "pos", "rsid", "ref", "alt",
                 "genotype", "retenue", "statut", "alleles_definis"],
    "annotations": ["hgnc_id", "gene", "phenotype", "phenotype_fr",
                    "score_activite", "medicament", "medicament_fr", "type",
                    "consigne", "source", "pmid", "traduite", "referentiel"],
}
# L'outil qui a rendu un appel externe. Tout autre gene remis a l'interpreteur
# de cette facon vient du typage complementaire.
OUTILS = {"CYP2D6": "Cyrius", "HLA-A": "OptiType", "HLA-B": "OptiType"}
# Les sources dont le genotype se lit dans le fichier de variants : pour les
# autres, les positions de ce fichier n'ont pas fonde l'appel.
SUR_VARIANTS = {"PharmCAT", "PyPGx"}


def hgnc(table, symbole):
    """Un gene sans identifiant arrete tout : une cellule vide casserait la
    jointure entre les tables sans que rien ne le signale."""
    if symbole not in table:
        raise ValueError("gene sans identifiant HGNC : %s" % symbole)
    return table[symbole]


def _pl(n):
    return "s" if n > 1 else ""


# ------------------------------------------------------------------ genotype
def genotypes(rapport, meta, typage, table):
    """Un gene par ligne : son genotype, d'ou il vient, sur quoi il repose.

    statut : rendu, reserve (des positions n'ont pas ete lues), ambigu
    (plusieurs lectures possibles, toutes donnees), non_rendu.
    """
    genes, qc = rapport.get("genes") or {}, meta.get("genes") or {}
    mesures, absents, reserves = L.perimetre(rapport, meta)
    les_ambigus = L.ambigus(genes, mesures)
    types = {r["gene"]: r for r in typage}
    lignes = {}

    def ligne(sym, **champs):
        d = dict.fromkeys(COLONNES["genotypes"], "")
        d.update(hgnc_id=hgnc(table, sym), gene=sym, **champs)
        lignes[sym] = d

    def du_typage(sym):
        r = types[sym]
        compte = {"positions_attendues": r.get("positions_attendues", ""),
                  "positions_retenues": r.get("positions_lues", "")}
        if r.get("statut") == "rendu":
            ligne(sym, genotype=r.get("genotype") or "", statut="rendu",
                  source="PyPGx", **compte)
        else:
            ligne(sym, statut="non_rendu", source="PyPGx",
                  motif=r.get("motif") or "échec du typeur", **compte)

    for sym in mesures:
        externe = genes[sym].get("callSource") == "OUTSIDE"
        champs = {"genotype": " | ".join(L.diplotypes(genes, sym)[0]),
                  "statut": "rendu",
                  "source": OUTILS.get(sym, "PyPGx") if externe else "PharmCAT"}
        if sym in les_ambigus:
            champs["statut"] = "ambigu"
        elif sym in reserves:
            n, tot = reserves[sym]
            champs.update(statut="reserve", motif="%d position%s non lue%s sur %d"
                          % (n, _pl(n), _pl(n), tot))
        if not externe:
            champs.update(positions_attendues=qc[sym]["positions_attendues"],
                          positions_retenues=qc[sym]["positions_retenues"])
        elif sym in types:
            champs.update(positions_attendues=types[sym].get("positions_attendues", ""),
                          positions_retenues=types[sym].get("positions_lues", ""))
        ligne(sym, **champs)
    for sym, motif in absents:
        if sym in types:
            du_typage(sym)
        else:
            ligne(sym, statut="non_rendu", motif=motif)
    for sym in types:
        if sym not in lignes:
            du_typage(sym)
    return [lignes[s] for s in sorted(lignes)]


# ------------------------------------------------------------------ variants
def _alleles(gt, ref, alt):
    """Le genotype en alleles : sur un site a plusieurs alleles, l'indice seul
    ne dit pas lequel est porte."""
    noms = [ref] + alt.split(",")
    return "".join(noms[int(x)] if x.isdigit() and int(x) < len(noms) else x
                   for x in re.split(r"([/|])", gt))


def variants(positions, rapport, lignes_genotypes, table):
    """Les variants observes des genes dont le genotype se lit dans le fichier
    de variants. Une position ecartee par le controle qualite reste listee,
    marquee : vue puis ecartee, elle n'est pas « reference »."""
    retenus = {l["gene"] for l in lignes_genotypes if l["source"] in SUR_VARIANTS}
    definis = {(v.get("chromosome"), str(v.get("position"))): v.get("alleles") or []
               for e in (rapport.get("genes") or {}).values()
               for v in e.get("variants") or []}
    out = []
    for p in positions:
        gt = p.get("genotype") or ""
        if not any(x.isdigit() and int(x) > 0 for x in re.split(r"[/|]", gt)):
            continue
        for sym in (p.get("genes") or "").split(","):
            if sym not in retenus:
                continue
            rsid = p.get("rsid") or ""
            out.append({
                "hgnc_id": hgnc(table, sym), "gene": sym, "chrom": p["chrom"],
                "pos": int(p["pos"]), "rsid": "" if rsid == "." else rsid,
                "ref": p["ref"], "alt": p["alt"],
                "genotype": _alleles(gt, p["ref"], p["alt"]),
                "retenue": p.get("retenue", ""), "statut": p.get("statut", ""),
                "alleles_definis": ";".join(definis.get((p["chrom"], p["pos"]), []))})
    return sorted(out, key=lambda v: (v["gene"], v["chrom"], v["pos"]))


# ---------------------------------------------------------------- annotation
def annotations(rapport, meta, tr, table):
    """Le phenotype de chaque gene rendu, et une ligne par consigne restituee
    et par gene qui la declenche. Un gene sans consigne garde une ligne : son
    phenotype est une annotation a lui seul."""
    genes = rapport.get("genes") or {}
    meds = tr["medicaments_restitues"]
    mesures, _, _ = L.perimetre(rapport, meta)
    par_med = L.par_medicament(L.recommandations(
        rapport, genes, set(mesures), L.ambigus(genes, mesures), tr, meds))
    referentiel = "PharmCAT %s ; CPIC %s ; traductions %s" % (
        rapport.get("pharmcatVersion", "?"), rapport.get("dataVersion", "?"),
        tr.get("version", "?"))
    out = []
    for sym in mesures:
        e = genes[sym]
        dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
        base = dict.fromkeys(COLONNES["annotations"], "")
        base.update(
            hgnc_id=hgnc(table, sym), gene=sym, referentiel=referentiel,
            phenotype=" | ".join(L.diplotypes(genes, sym)[1]),
            phenotype_fr=L.phenotype(genes, sym),
            score_activite=" | ".join(sorted(
                {str(d["activityScore"]) for d in dips
                 if d.get("activityScore") is not None})))
        consignes = sorted(n for n, d in par_med.items() if sym in d["genes"])
        for nom in consignes:
            d = par_med[nom]
            out.append(dict(
                base, medicament=nom, medicament_fr=meds[nom], type=d["type"],
                consigne=" ".join(d["fr"]), source=" ; ".join(d["source"]),
                pmid=";".join(sorted(d["pmid"])),
                traduite="oui" if d["ok"] else "non"))
        if not consignes:
            out.append(base)
    return out


# -------------------------------------------------------------------- script
def _table(chemin):
    with open(chemin, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    p = argparse.ArgumentParser(description="Genotype et annotation, en tables.")
    p.add_argument("--rapport", required=True, help="le JSON produit par PharmCAT")
    p.add_argument("--perimetre", required=True,
                   help="perimetre.json ; qc_positions.tsv et pypgx.tsv sont lus a cote")
    p.add_argument("--sortie", required=True, help="dossier des livrables")
    p.add_argument("--echantillon", required=True)
    p.add_argument("--traductions",
                   default=os.path.join(RACINE, "ressources", "traductions_cpic_fr.json"))
    p.add_argument("--hgnc", default=os.path.join(RACINE, "ressources", "genes_hgnc.json"))
    a = p.parse_args()
    travail = os.path.dirname(os.path.abspath(a.perimetre))
    def charge(chemin):
        with open(chemin, encoding="utf-8") as fh:
            return json.load(fh)
    rapport, meta, tr = charge(a.rapport), charge(a.perimetre), charge(a.traductions)
    table = charge(a.hgnc)["genes"]
    typage = os.path.join(travail, "pypgx.tsv")
    typage = _table(typage) if os.path.exists(typage) else []
    # Tout est calcule avant la premiere ecriture : un echec ne laisse pas une
    # table sans les deux autres.
    try:
        g = genotypes(rapport, meta, typage, table)
        tables = {"genotypes": g,
                  "variants": variants(_table(os.path.join(travail, "qc_positions.tsv")),
                                       rapport, g, table),
                  "annotations": annotations(rapport, meta, tr, table)}
    except ValueError as e:
        sys.exit("sorties structurees refusees : %s" % e)
    os.makedirs(a.sortie, exist_ok=True)
    for nom, lignes in tables.items():
        with open(os.path.join(a.sortie, "%s.%s.tsv" % (a.echantillon, nom)), "w",
                  encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLONNES[nom], delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(lignes)
    print("%s : %d genes, %d variants, %d annotations" % (
        a.echantillon, len(g), len(tables["variants"]), len(tables["annotations"])))


if __name__ == "__main__":
    main()
