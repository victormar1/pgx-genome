# -*- coding: utf-8 -*-
"""Verite CYP2D6 d'un haplotype d'assemblage, a partir de la sortie de deCYPher.

Quand plusieurs alleles etoiles voient toutes leurs variations definissantes
presentes dans une copie :
1. on ecarte ceux dont la definition est contenue dans celle d'un autre
   candidat : l'allele le plus specifique l'emporte, comme dans PharmVar
   (*49 contient *10, *84 contient *2) ;
2. entre les candidats restants, que rien n'ordonne, la verite retient le plus
   severe selon la table d'activite CPIC de PharmCAT 3.4.0 : la variation la
   plus deletere decide de la fonction (*56 contre *10, *4 contre *10). A
   egalite, l'appel de deCYPher est garde.

deCYPher retient le premier allele de son alignement : il a rendu *10 pour une
copie portant aussi R344X, donc *56 (un genome du pont). L'etape 1 est un
amendement du 9 novembre 2026, pose avant tout resultat du module sur les 211 :
la seule etape 2 remplacait un allele specifique sans activite chiffree par un
allele plus general (*84 en *2) ou un *49 par *10. Un decompte des variations
definissantes ne convient pas non plus : chez deCYPher, la definition de *4 est
plus courte que celle de *10.
"""
import json
import os
import re

import chemins

HYBRIDES = {"*13", "*36", "*61", "*63", "*68", "*83"}
ACTIVITE = {}
with open(os.path.join(chemins.VERITES, "activite_cyp2d6.json"), encoding="utf-8") as _f:
    for _a, _v in json.load(_f).items():
        try:
            ACTIVITE[_a] = float(_v)
        except (TypeError, ValueError):
            pass      # « n/a » : allele sans valeur d'activite CPIC


def _refseq(trav):
    return os.path.join(trav, "outils", "deCYPher", "ref", "CYP2D6", "RefSeqGene")


def _definition(etoile, trav):
    f = os.path.join(_refseq(trav), f"CYP2D6.{etoile[1:]}.vcf")
    if not os.path.exists(f):
        return None
    with open(f, encoding="utf-8") as fh:
        champs = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]
    return frozenset((c[1], c[3], c[4]) for c in champs if len(c) > 4)


def _maximaux(candidats, trav):
    defs = {c: _definition(c, trav) for c in candidats}
    return [c for c in candidats if defs[c] is None or not any(
        defs[o] is not None and defs[c] < defs[o] for o in candidats if o != c)]


def _plus_severe(candidats, rendu):
    notes = [c for c in candidats if c in ACTIVITE]
    if not notes:
        return rendu if rendu in candidats else None
    mini = min(ACTIVITE[c] for c in notes)
    egaux = [c for c in notes if ACTIVITE[c] == mini]
    return rendu if rendu in egaux else sorted(egaux)[0]


def _retenir(candidats, rendu, trav):
    maxi = _maximaux(candidats, trav)
    return maxi[0] if len(maxi) == 1 else _plus_severe(maxi, rendu)


def _etoiles_completes(txt):
    """Alleles dont deCYPher dit, dans son compte rendu, qu'ils correspondent a
    100 % a la copie."""
    etoiles, dans_liste = [], False
    with open(txt, encoding="utf-8") as fh:
        for ligne in fh:
            if "100% matching" in ligne:
                dans_liste = True
                continue
            m = re.match(r"\s*-\s*CYP2D6\.([0-9]+)\.vcf", ligne)
            if dans_liste and m:
                etoiles.append("*" + m.group(1))
    return etoiles


def copies(p, trav=None):
    """Copies de CYP2D6 de l'haplotype p, apres les deux regles, et les
    corrections appliquees (liste de (rendu deCYPher, retenu))."""
    trav = chemins.travail(trav)
    f = os.path.join(trav, "decypher", p, p + ".json")
    if not os.path.exists(f):
        return None, []
    with open(f, encoding="utf-8") as fh:
        d = json.load(fh).get("CYP2D6") or {}
    if "note" in d or not d:
        return None, []
    rendus = [v for k, v in sorted(d.items(), key=lambda kv: int(kv[0][6:])) if k.startswith("allele")]
    if rendus == ["*5"]:
        return ["*5"], []           # gene entierement delete : rien a departager
    res = os.path.join(trav, "decypher", p, "CYP2D6", "result")
    txts = sorted(t for t in os.listdir(res) if t.endswith(".txt")) if os.path.isdir(res) else []
    retenus, corrections = [], []
    for i, rendu in enumerate(rendus):
        if i < len(txts) and rendu not in ("*13",):
            candidats = _etoiles_completes(os.path.join(res, txts[i]))
            if candidats:
                meilleur = _retenir(candidats, rendu, trav)
                if meilleur is None:
                    return None, corrections
                if meilleur != rendu:
                    corrections.append((rendu, meilleur))
                    rendu = meilleur
        retenus.append(rendu)
    if any(not str(c).startswith("*") for c in retenus):
        return None, corrections
    return retenus, corrections


def haplotype(cps):
    """Ecriture d'un haplotype : *Xxn pour une duplication, les hybrides
    d'abord dans un tandem."""
    if not cps:
        return None
    if len(set(cps)) == 1:
        return cps[0] + (f"x{len(cps)}" if len(cps) > 1 else "")
    return "+".join(sorted(cps, key=lambda c: (c not in HYBRIDES, c)))


def multiensemble(diplotype):
    """Copies d'un diplotype ecrit '*10/*36x2' ou '*2/*36+*10', sans la phase."""
    sortie = []
    for hap in re.split(r"\s*/\s*", diplotype.strip()):
        for c in re.split(r"\s*\+\s*", hap):
            m = re.match(r"(\*[0-9A-Za-z.]+?)(?:x([0-9]+))?$", c.strip())
            if not m:
                return None
            sortie += [m.group(1).split(".")[0]] * int(m.group(2) or 1)
    return sorted(sortie)
