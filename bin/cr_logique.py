# -*- coding: utf-8 -*-
"""Logique du compte rendu, sans mise en page.

Ce module ne contient que ce qui decide : quels genes sont rendus, comment un
diplotype se lit en francais, quelles recommandations sont restituees et
laquelle l'emporte quand plusieurs portent sur le meme medicament. La mise en
page reste dans `compte_rendu.py`.

La separation existe pour une seule raison : rien de tout cela n'etait
verifiable. Le module de rendu lit ses arguments et ses ressources a l'import,
donc aucune de ses fonctions ne pouvait etre appelee par un test — et c'est
pourtant ici qu'on avait deja perdu quatre conduites a tenir sur neuf, dont la
codeine et le tramadol.
"""
import html
import re

# ------------------------------------------------------------------ libelles
PHENO_FR = {
    "Normal Metabolizer": "métaboliseur normal",
    "Intermediate Metabolizer": "métaboliseur intermédiaire",
    "Poor Metabolizer": "métaboliseur lent",
    "Rapid Metabolizer": "métaboliseur rapide",
    "Ultrarapid Metabolizer": "métaboliseur ultrarapide",
    "Likely Intermediate Metabolizer": "métaboliseur intermédiaire probable",
    "Likely Poor Metabolizer": "métaboliseur lent probable",
    "Normal Function": "fonction normale",
    "Decreased Function": "fonction diminuée",
    "Possible Decreased Function": "fonction possiblement diminuée",
    "Poor Function": "fonction faible",
    "Increased Function": "fonction augmentée",
    "Indeterminate": "non interprétable",
    "No Result": "non appelé", "n/a": "—",
    # MT-RNR1 : l'interpreteur rend un niveau de risque, pas un metabolisme.
    "normal risk of aminoglycoside-induced hearing loss":
        "risque normal de surdité sous aminoside",
    "uncertain risk of aminoglycoside-induced hearing loss":
        "risque incertain de surdité sous aminoside",
    "increased risk of aminoglycoside-induced hearing loss":
        "risque augmenté de surdité sous aminoside",
}
CYP3A5_EXPR = {"Normal Metabolizer": "expresseur",
               "Intermediate Metabolizer": "expresseur",
               "Poor Metabolizer": "non-expresseur"}
FORCE_FR = {"Strong": "forte", "Moderate": "modérée"}
RESTITUES = {"contre-indication": 0, "éviter": 1, "adaptation": 2, "vigilance": 3}
STATINES = {"simvastatin", "atorvastatin", "rosuvastatin", "pravastatin",
            "fluvastatin"}
# Les onze aminosides de la recommandation MT-RNR1 portent la meme consigne :
# les nommer un par un remplirait la page sans rien ajouter.
AMINOSIDES = {"amikacin", "dibekacin", "gentamicin", "kanamycin", "neomycin",
              "netilmicin", "paromomycin", "plazomicin", "ribostamycin",
              "streptomycin", "tobramycin"}
NORMAUX = {"Normal Metabolizer", "Normal Function",
           "normal risk of aminoglycoside-induced hearing loss"}


def e_(t):
    return html.escape(str(t), quote=False)


def norm(t):
    return re.sub(r"\s+", " ", html.unescape(str(t))).strip()


# ------------------------------------------------------------------ perimetre
def perimetre(rapport, meta):
    """Croise ce que le controle qualite a mesure et ce que l'interpreteur a rendu.

      complet  toutes les positions diagnostiques lues -> rendu sans reserve
      partiel  au moins une position non lue           -> rendu avec reserve
      absent   aucune position lue, ou aucun diplotype -> non rendu

    Les genes venus d'un outil dedie (Cyrius, OptiType, typage complementaire)
    ont le perimetre de cet outil : leur appel ne vient pas du fichier de
    variants, et le juger sur les positions de celui-ci le refuserait a tort.
    """
    qc = meta.get("genes", {})
    clinique = set(meta.get("perimetre_clinique") or [])
    hors_interpreteur = list(meta.get("perimetre_clinique_hors_interpreteur") or [])
    rendus, reserves, absents = [], {}, []
    for sym, e in rapport.get("genes", {}).items():
        if clinique and sym not in clinique:
            continue
        dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
        lab = ";".join(sorted({x.get("label", "") for x in dips if x.get("label")}))
        appele = bool(lab) and "Unknown" not in lab
        v = qc.get(sym)
        statut = (v or {}).get("statut")
        if not appele:
            absents.append((sym, "aucun diplotype rendu" if v
                            else "non appelable depuis le fichier de variants"))
        elif e.get("callSource") == "OUTSIDE" or v is None or statut == "complet":
            rendus.append(sym)
        elif statut == "partiel":
            rendus.append(sym)
            # La reserve porte son denominateur : perdre quatre positions sur
            # sept n'est pas perdre trois positions sur quarante-six, et c'est
            # la proportion qui dit si le diplotype est utilisable.
            reserves[sym] = (v.get("positions_perdues"),
                             v.get("positions_attendues"))
        else:
            absents.append((sym, "aucune position lue"))
    for sym in sorted(clinique - set(rapport.get("genes", {}))):
        absents.append((sym, "gène absent des tables de l'interpréteur"
                        if sym in hors_interpreteur
                        else "non rendu par l'interpréteur"))
    return sorted(rendus), sorted(absents), reserves


# ------------------------------------------------------------- lecture d'un gene
def diplotypes(genes, sym):
    e = genes.get(sym) or {}
    dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
    labs = sorted({x.get("label", "") for x in dips if x.get("label")})
    phen = sorted({p for x in dips for p in (x.get("phenotypes") or [])})
    return labs, phen


def libelle(genes, sym):
    labs, _ = diplotypes(genes, sym)
    out = []
    for lab in labs:
        # ABCG2 : l'interpreteur nomme le rsID par allele, pas le diplotype.
        m = re.findall(r"(rs\d+) (?:reference|variant) \((\w+)\)", lab)
        if m and len(m) == 2:
            lab = "%s %s/%s" % (m[0][0], m[0][1], m[1][1])
        lab = lab.replace("Reference", "référence").replace("Unknown", "non déterminé")
        out.append(lab)
    if len(out) > 1:
        return "ambigu : " + " ou ".join(out)
    return out[0] if out else "—"


def phenotype(genes, sym):
    _, phen = diplotypes(genes, sym)
    if sym.startswith("HLA-"):
        pos = [p.replace(" positive", "") for p in phen if p.endswith("positive")]
        neg = [p.replace(" negative", "") for p in phen if p.endswith("negative")]
        if pos:
            return ", ".join("%s présent" % a for a in pos)
        if neg:
            return ", ".join(neg) + (" absents" if len(neg) > 1 else " absent")
        return "—"
    fr = []
    for p in phen:
        t = PHENO_FR.get(p, p)
        if sym == "CYP3A5" and p in CYP3A5_EXPR:
            t += " (%s)" % CYP3A5_EXPR[p]
        fr.append(t)
    return ", ".join(fr) or "—"


def ambigus(genes, mesures):
    return {s for s in mesures if len(diplotypes(genes, s)[0]) > 1}


def gene_normal(genes, sym):
    """Phenotype normal, ou aucun allele HLA a risque : ce gene ne declenche pas
    la consigne. La cle de l'interpreteur porte parfois un score d'activite,
    d'ou ce jugement sur le phenotype et non sur le nom du diplotype."""
    ph = diplotypes(genes, sym)[1]
    return bool(ph) and all(p in NORMAUX or p.endswith("negative") for p in ph)


# -------------------------------------------------------- recommandations
def recommandations(rep, genes, mesures, les_ambigus, tr, meds):
    """Les consignes restituees, sous forme brute et triable.

    Une consigne n'est retenue que si elle est forte ou moderee, qu'elle porte
    sur un medicament du perimetre, et que tous les genes de sa cle sont rendus
    et non ambigus : un diplotype ambigu ne declenche aucune conduite.
    """
    brutes = []
    for nom, m in rep["drugs"].get("CPIC Guideline Annotation", {}).items():
        if nom not in meds:
            continue
        cites = sorted(m.get("citations") or [], key=lambda c: c.get("year") or 0)
        pmid = cites[-1].get("pmid", "") if cites else ""
        for gl in m.get("guidelines") or []:
            for a in gl.get("annotations") or []:
                cl = a.get("classification")
                if cl not in FORCE_FR:
                    continue
                cle = [(g, str(v)) for d in (a.get("lookupKey") or [])
                       for g, v in d.items()]
                if not cle or any(g not in mesures or g in les_ambigus
                                  for g, _ in cle):
                    continue
                # le resultat qui declenche la consigne : pas un gene normal,
                # pas un allele absent
                gcle = [g for g, v in cle
                        if not gene_normal(genes, g) and not v.endswith("negative")] \
                    or [g for g, _ in cle]
                src = norm(a.get("drugRecommendation", ""))
                t = tr["traductions"].get(src)
                if t:
                    fr, ty, ok = t["fr"], t["type"], True
                else:
                    # texte inconnu de la ressource : rendu tel quel, jamais tu
                    fr, ok = src, False
                    ty = "adaptation" if (a.get("alternateDrugAvailable")
                                          or a.get("dosingInformation")) else "standard"
                if ty in RESTITUES:
                    brutes.append((RESTITUES[ty], nom, fr, ty,
                                   "CPIC, %s" % FORCE_FR[cl],
                                   tuple(sorted(gcle)), pmid, ok))

    # alleles de classe 1 du RNPGx que l'interpreteur n'evalue pas
    regles = {k: v for k, v in (tr.get("regles_hla_rnpgx") or {}).items()
              if k.startswith("HLA-")}
    for cle, r in regles.items():
        g, allele = cle.split("*", 1)
        if g in mesures and g not in les_ambigus \
           and ("*" + allele) in " ".join(diplotypes(genes, g)[0]):
            for nom in r["medicaments"]:
                brutes.append((RESTITUES[r["type"]], nom, r["fr"], r["type"],
                               r["source"], (g,), "", True))
    return brutes


def par_medicament(brutes):
    """Un medicament, la consigne la plus severe. A severite egale, les textes
    distincts s'additionnent plutot que de s'ecraser : deux genes peuvent
    concourir au meme medicament pour des raisons differentes, et n'en garder
    qu'une retirerait une conduite du compte rendu.
    """
    par_med = {}
    for b in sorted(brutes):
        s, nom = b[0], b[1]
        cur = par_med.get(nom)
        if cur is None or s < cur["sev"]:
            par_med[nom] = {"sev": s, "fr": [b[2]], "type": b[3],
                            "source": [b[4]], "genes": set(b[5]),
                            "pmid": {b[6]} - {""}, "ok": b[7]}
        elif s == cur["sev"]:
            if b[2] not in cur["fr"]:
                cur["fr"].append(b[2])
            if b[4] not in cur["source"]:
                cur["source"].append(b[4])
            cur["genes"] |= set(b[5])
            cur["pmid"] |= {b[6]} - {""}
            cur["ok"] &= b[7]
    return par_med
