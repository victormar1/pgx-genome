# -*- coding: utf-8 -*-
"""Compte rendu pharmacogenetique en francais, depuis la sortie structuree de PharmCAT.

Hierarchie de restitution (cahier des charges du POC PFMG) : d'abord les medicaments pour
lesquels une recommandation s'applique, puis les resultats par gene, une conclusion, le
perimetre et les limites, et la ligne de validation par le biologiste.

Regles de securite :
 1. seuls les genes reellement sequences sont rendus. Les autres sortent en « non analyse »,
    jamais en « normal » ;
 2. ne sont restituees que les recommandations CPIC fortes ou moderees qui modifient la prise
    en charge. Un diplotype ambigu n'en declenche aucune ;
 3. aucune recommandation n'est tue faute de traduction : le texte source est alors rendu,
    signale comme tel.
"""
import argparse, html, json, os, re, datetime
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_p = argparse.ArgumentParser(description="Compte rendu pharmacogenetique, francais.")
_p.add_argument("--rapport", required=True, help="le JSON produit par PharmCAT")
_p.add_argument("--perimetre", required=True, help="perimetre.json produit par le controle qualite")
_p.add_argument("--sortie", required=True, help="chemin du PDF a ecrire")
_p.add_argument("--echantillon", default="")
_p.add_argument("--patient", default="", help="identite du patient, fournie par le SIL")
_p.add_argument("--naissance", default="", help="date de naissance, fournie par le SIL")
_p.add_argument("--contexte", default="", help="preindication PFMG, par ex. « néphropathies »")
_p.add_argument("--traductions", default=os.path.join(RACINE, "ressources", "traductions_cpic_fr.json"))
_p.add_argument("--sans-mention-prototype", action="store_true",
                help="a utiliser pour un vrai patient ; par defaut, mention prototype 1000 Genomes")
_a = _p.parse_args()
SRC, OUT = _a.rapport, _a.sortie
ECH = _a.echantillon or os.path.basename(SRC).split(".")[0]
os.makedirs(os.path.dirname(os.path.abspath(OUT)) or ".", exist_ok=True)
TR = json.load(open(_a.traductions, encoding="utf-8"))
MEDS = TR["medicaments_restitues"]


def _perimetre(rapport):
    """Croise ce que le controle qualite a mesure et ce que l'interpreteur a rendu.
      complet  toutes les positions diagnostiques lues -> rendu sans reserve
      partiel  au moins une position non lue           -> rendu avec reserve
      absent   aucune position lue, ou aucun diplotype -> non rendu
    Les genes venus d'un outil dedie (Cyrius, OptiType) ont le perimetre de cet outil."""
    meta = json.load(open(_a.perimetre, encoding="utf-8"))
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
            absents.append((sym, "aucun diplotype rendu" if v else "non appelable depuis le fichier de variants"))
        elif e.get("callSource") == "OUTSIDE" or v is None or statut == "complet":
            rendus.append(sym)
        elif statut == "partiel":
            rendus.append(sym)
            reserves[sym] = v["positions_perdues"]
        else:
            absents.append((sym, "aucune position lue"))
    for sym in sorted(clinique - set(rapport.get("genes", {}))):
        absents.append((sym, "gène absent des tables de l'interpréteur" if sym in hors_interpreteur
                        else "non rendu par l'interpréteur"))
    return sorted(rendus), sorted(absents), reserves


rep = json.load(open(SRC, encoding="utf-8"))
genes = rep["genes"]
MESURES, ABSENTS, RESERVES = _perimetre(rep)

# ------------------------------------------------------------------ libelles
PHENO_FR = {
    "Normal Metabolizer": "métaboliseur normal", "Intermediate Metabolizer": "métaboliseur intermédiaire",
    "Poor Metabolizer": "métaboliseur lent", "Rapid Metabolizer": "métaboliseur rapide",
    "Ultrarapid Metabolizer": "métaboliseur ultrarapide",
    "Likely Intermediate Metabolizer": "métaboliseur intermédiaire probable",
    "Likely Poor Metabolizer": "métaboliseur lent probable",
    "Normal Function": "fonction normale", "Decreased Function": "fonction diminuée",
    "Possible Decreased Function": "fonction possiblement diminuée", "Poor Function": "fonction faible",
    "Increased Function": "fonction augmentée", "Indeterminate": "non interprétable",
    "No Result": "non appelé", "n/a": "—",
}
CYP3A5_EXPR = {"Normal Metabolizer": "expresseur", "Intermediate Metabolizer": "expresseur",
               "Poor Metabolizer": "non-expresseur"}
FORCE_FR = {"Strong": "forte", "Moderate": "modérée"}
RESTITUES = {"contre-indication": 0, "éviter": 1, "adaptation": 2, "vigilance": 3}
STATINES = {"simvastatin", "atorvastatin", "rosuvastatin", "pravastatin", "fluvastatin"}
NORMAUX = {"Normal Metabolizer", "Normal Function"}
ordre_med = lambda n: (n == "fosphenytoin", MEDS[n])   # phenytoine avant fosphenytoine
e_ = lambda t: html.escape(str(t), quote=False)
norm = lambda t: re.sub(r"\s+", " ", html.unescape(str(t))).strip()


def diplotypes(sym):
    e = genes.get(sym) or {}
    dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
    labs = sorted({x.get("label", "") for x in dips if x.get("label")})
    phen = sorted({p for x in dips for p in (x.get("phenotypes") or [])})
    return labs, phen


def libelle(sym):
    labs, _ = diplotypes(sym)
    out = []
    for lab in labs:
        m = re.findall(r"(rs\d+) (?:reference|variant) \((\w+)\)", lab)   # ABCG2 : rsID par allele
        if m and len(m) == 2:
            lab = f"{m[0][0]} {m[0][1]}/{m[1][1]}"
        lab = lab.replace("Reference", "référence").replace("Unknown", "non déterminé")
        out.append(lab)
    return ("ambigu : " + " ou ".join(out)) if len(out) > 1 else (out[0] if out else "—")


def phenotype(sym):
    _, phen = diplotypes(sym)
    if sym.startswith("HLA-"):
        pos = [p.replace(" positive", "") for p in phen if p.endswith("positive")]
        neg = [p.replace(" negative", "") for p in phen if p.endswith("negative")]
        if pos:
            return ", ".join(f"{a} présent" for a in pos)
        return (", ".join(neg) + (" absents" if len(neg) > 1 else " absent")) if neg else "—"
    fr = []
    for p in phen:
        t = PHENO_FR.get(p, p)
        if sym == "CYP3A5" and p in CYP3A5_EXPR:
            t += f" ({CYP3A5_EXPR[p]})"
        fr.append(t)
    return ", ".join(fr) or "—"


AMBIGUS = {s for s in MESURES if len(diplotypes(s)[0]) > 1}


def gene_normal(sym):
    """Phenotype normal, ou aucun allele HLA a risque : ce gene ne declenche pas la consigne.
    (La cle de PharmCAT porte parfois un score d'activite, d'ou ce jugement sur le phenotype.)"""
    ph = diplotypes(sym)[1]
    return bool(ph) and all(p in NORMAUX or p.endswith("negative") for p in ph)

# --------------------------------------------------- recommandations restituees
brutes = []   # (severite, med, fr, type, force, genes, pmid, traduit)
for nom, m in rep["drugs"].get("CPIC Guideline Annotation", {}).items():
    if nom not in MEDS:
        continue
    cites = sorted(m.get("citations") or [], key=lambda c: c.get("year") or 0)
    pmid = cites[-1].get("pmid", "") if cites else ""
    for gl in m.get("guidelines") or []:
        for a in gl.get("annotations") or []:
            cl = a.get("classification")
            if cl not in FORCE_FR:
                continue
            cle = [(g, str(v)) for d in (a.get("lookupKey") or []) for g, v in d.items()]
            if not cle or any(g not in MESURES or g in AMBIGUS for g, _ in cle):
                continue
            # le resultat qui declenche la consigne : pas un gene normal ni un allele absent
            gcle = [g for g, v in cle if not gene_normal(g) and not v.endswith("negative")] or [g for g, _ in cle]
            src = norm(a.get("drugRecommendation", ""))
            t = TR["traductions"].get(src)
            if t:
                fr, ty, ok = t["fr"], t["type"], True
            else:   # texte inconnu de la ressource : rendu tel quel, jamais tu
                fr, ok = src, False
                ty = "adaptation" if (a.get("alternateDrugAvailable") or a.get("dosingInformation")) else "standard"
            if ty in RESTITUES:
                brutes.append((RESTITUES[ty], nom, fr, ty, f"CPIC, {FORCE_FR[cl]}", tuple(sorted(gcle)), pmid, ok))

# allèles de classe 1 du RNPGx que PharmCAT n'évalue pas
for cle, r in {k: v for k, v in TR.get("regles_hla_rnpgx", {}).items() if k.startswith("HLA-")}.items():
    g, allele = cle.split("*", 1)
    if g in MESURES and g not in AMBIGUS and f"*{allele}" in " ".join(diplotypes(g)[0]):
        for nom in r["medicaments"]:
            brutes.append((RESTITUES[r["type"]], nom, r["fr"], r["type"], r["source"], (g,), "", True))

# un medicament = la consigne la plus severe ; a severite egale, les textes distincts s'additionnent
par_med = {}
for b in sorted(brutes):
    s, nom = b[0], b[1]
    cur = par_med.get(nom)
    if cur is None or s < cur["sev"]:
        par_med[nom] = {"sev": s, "fr": [b[2]], "type": b[3], "source": [b[4]], "genes": set(b[5]),
                        "pmid": {b[6]} - {""}, "ok": b[7]}
    elif s == cur["sev"]:
        if b[2] not in cur["fr"]:
            cur["fr"].append(b[2])
        if b[4] not in cur["source"]:
            cur["source"].append(b[4])
        cur["genes"] |= set(b[5]); cur["pmid"] |= {b[6]} - {""}; cur["ok"] &= b[7]

# regroupe les medicaments a consigne identique (IPP, phenytoine/fosphenytoine) ; statines sur une ligne
lignes, vus = [], set()
for nom, d in sorted(par_med.items(), key=lambda kv: (kv[1]["sev"], MEDS[kv[0]])):
    if nom in vus:
        continue
    if nom in STATINES and sum(n in STATINES for n in par_med) > 1:
        grp = [n for n in sorted(par_med, key=lambda n: (par_med[n]["sev"], MEDS[n])) if n in STATINES]
        vus |= set(grp)
        lignes.append({"noms": ["statines"], "sev": min(par_med[n]["sev"] for n in grp),
                       "action": "<br/>".join(f"<b>{MEDS[n]}</b> : {e_(' '.join(par_med[n]['fr']))}"
                                              + ("" if par_med[n]['source'] == ["CPIC, forte"] else
                                                 f" <font color='#666666'>({par_med[n]['source'][0].split(', ')[-1]})</font>")
                                              for n in grp),
                       "genes": set().union(*(par_med[n]["genes"] for n in grp)),
                       "source": "CPIC", "pmid": set().union(*(par_med[n]["pmid"] for n in grp)),
                       "ok": all(par_med[n]["ok"] for n in grp), "meds": grp})
        continue
    freres = [n for n, x in par_med.items() if n not in vus and n not in STATINES
              and x["fr"] == d["fr"] and x["genes"] == d["genes"] and x["source"] == d["source"]]
    vus |= set(freres)
    lignes.append({"noms": [MEDS[n] for n in sorted(freres, key=ordre_med)], "sev": d["sev"],
                   "action": e_(" ".join(d["fr"])), "genes": d["genes"], "source": " ; ".join(d["source"]),
                   "pmid": set().union(*(par_med[n]["pmid"] for n in freres)), "ok": d["ok"],
                   "meds": sorted(freres, key=ordre_med)})

MED_PAR_GENE = {}
for l in lignes:
    for g in l["genes"]:
        for n in l["meds"]:
            MED_PAR_GENE.setdefault(g, []).append((par_med[n]["sev"], MEDS[n]))
NON_TRADUITS = sorted({MEDS[n] for l in lignes for n in l["meds"] if not par_med[n]["ok"]})

# ------------------------------------------------------------------ styles
NOIR, GRIS = colors.HexColor("#1a1a1a"), colors.HexColor("#666666")
TRAIT, ROUGE = colors.HexColor("#cccccc"), colors.HexColor("#a02020")
FOND_R, FOND_G, FOND_T = colors.HexColor("#fbeeee"), colors.HexColor("#f4f4f2"), colors.HexColor("#ececea")


def S(n, **k):
    base = dict(fontName="Helvetica", fontSize=8.2, leading=10.4, textColor=NOIR, alignment=TA_LEFT)
    base.update(k)
    return ParagraphStyle(n, **base)


st_titre = S("t", fontName="Helvetica-Bold", fontSize=13, leading=15)
st_section = S("s", fontName="Helvetica-Bold", fontSize=8.8, leading=10.5, spaceBefore=8, spaceAfter=3)
st_corps = S("c", fontSize=7.9, leading=9.8)
st_rouge = S("r", fontSize=7.9, leading=9.8, textColor=ROUGE, fontName="Helvetica-Bold")
st_petit = S("p", fontSize=7.2, leading=9.0, textColor=GRIS)
st_tete = S("h", fontSize=7.2, leading=9.0, textColor=GRIS, fontName="Helvetica-Bold")
st_gras = S("g", fontName="Helvetica-Bold", fontSize=7.9, leading=9.8)


def tableau(rows, widths, fonds=None):
    t = Table(rows, colWidths=[w * mm for w in widths], repeatRows=1)
    style = [("VALIGN", (0, 0), (-1, -1), "TOP"),
             ("LINEBELOW", (0, 0), (-1, 0), 0.6, NOIR),
             ("LINEBELOW", (0, 1), (-1, -1), 0.3, TRAIT),
             ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
             ("TOPPADDING", (0, 0), (-1, -1), 1.6), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6)]
    for i, f in (fonds or {}).items():
        style.append(("BACKGROUND", (0, i), (-1, i), f))
    t.setStyle(TableStyle(style))
    return t


# ------------------------------------------------------------------ document
h = [Paragraph("COMPTE RENDU DE PHARMACOGÉNÉTIQUE", st_titre), Spacer(1, 4)]
proto = not _a.sans_mention_prototype
ident = []
if _a.patient:
    ident.append(("Patient", e_(_a.patient) + (f" — né(e) le {e_(_a.naissance)}" if _a.naissance else "")))
ident.append(("Prélèvement", f"{e_(ECH)}" + (" — génome public du projet 1000 Genomes, aucun patient" if proto else "")))
ident.append(("Contexte", f"Préindication PFMG : {e_(_a.contexte)}" if _a.contexte else
              "Sous-ensemble néphrologie et épilepsie du panel socle RNPGx 2026"))
ident.append(("Analyse", "Extraction pharmacogénétique secondaire à partir du génome (lectures courtes, GRCh38) — "
              f"{len(MESURES)} gènes rendus sur 12"))
ident.append(("Date / version", f"{datetime.date.today().strftime('%d/%m/%Y')} — PharmCAT {e_(rep.get('pharmcatVersion', '?'))}, "
              f"référentiel CPIC, traductions v{e_(TR.get('version', '?'))}"))
t = Table([[Paragraph(f"<b>{k}</b>", st_petit), Paragraph(v, st_corps)] for k, v in ident],
          colWidths=[28 * mm, 142 * mm])
t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("BACKGROUND", (0, 0), (-1, -1), FOND_T),
                       ("LEFTPADDING", (0, 0), (-1, -1), 4), ("TOPPADDING", (0, 0), (-1, -1), 1.4),
                       ("BOTTOMPADDING", (0, 0), (-1, -1), 1.4)]))
h += [t]

# --- 1. medicaments
h.append(Paragraph("MÉDICAMENTS POUR LESQUELS UNE RECOMMANDATION PHARMACOGÉNÉTIQUE S'APPLIQUE", st_section))
if lignes:
    rows = [[Paragraph(x, st_tete) for x in ("Médicament", "Action à retenir", "Résultat PGx", "Source")]]
    fonds = {}
    for i, l in enumerate(lignes, 1):
        res = "<br/>".join(f"<b>{g}</b> {e_(libelle(g))} — {e_(phenotype(g))}" for g in sorted(l["genes"]))
        src = l["source"] + ("".join(f" ; PMID {p}" for p in sorted(l["pmid"])))
        ci = l["sev"] == 0
        rows.append([Paragraph(f"<b>{' / '.join(l['noms']).upper()}</b>", st_rouge if ci else st_gras),
                     Paragraph(l["action"] + ("" if l["ok"] else " <i>(texte CPIC d'origine, traduction à établir)</i>"),
                               st_rouge if ci else st_corps),
                     Paragraph(res, st_corps), Paragraph(src, st_petit)])
        if ci:
            fonds[i] = FOND_R
    h.append(tableau(rows, [27, 75, 44, 24], fonds))
else:
    h.append(Paragraph("Aucune recommandation pharmacogénétique modifiant la prise en charge n'est identifiée "
                       "dans le périmètre testé.", st_corps))

# --- 2. resultats par gene
h.append(Paragraph("RÉSULTATS PHARMACOGÉNÉTIQUES DU PÉRIMÈTRE", st_section))
rows = [[Paragraph(x, st_tete) for x in ("Gène", "Résultat", "Phénotype / fonction", "Conséquence")]]
ordre = sorted(MESURES, key=lambda g: (g not in MED_PAR_GENE, g))
for g in ordre:
    if g in AMBIGUS:
        cons = "Diplotype ambigu : non restitué, expertise manuelle nécessaire."
    elif g in MED_PAR_GENE:
        meds = sorted(set(MED_PAR_GENE[g]))
        ci = [m for s, m in meds if s == 0]
        autres = [m for s, m in meds if s != 0]
        cons = ((f"Contre-indication : {', '.join(ci)}. " if ci else "")
                + (f"Recommandation pour {', '.join(autres)}. " if autres else "") + "Voir ci-dessus.")
    else:
        cons = "Pas de modification de la prise en charge."
    if g in RESERVES:
        n = RESERVES[g]
        cons += f" Sous réserve : {n} position{'s' if n > 1 else ''} non lue{'s' if n > 1 else ''}."
    rows.append([Paragraph(g, st_gras), Paragraph(e_(libelle(g)), st_corps),
                 Paragraph(e_(phenotype(g)), st_corps), Paragraph(cons, st_corps)])
for g, motif in ABSENTS:
    rows.append([Paragraph(g, st_gras), Paragraph("non analysé", st_corps), Paragraph("—", st_corps),
                 Paragraph(f"Aucune conclusion : {e_(motif)}. Ni normal ni anormal, non mesuré.", st_corps)])
h.append(tableau(rows, [17, 32, 45, 76], {i: FOND_G for i in range(2, len(rows), 2)}))

# --- 3. conclusion
h.append(Paragraph("CONCLUSION", st_section))
tous = [MEDS[n] for l in lignes for n in l["meds"]]
ci = [MEDS[n] for l in lignes for n in l["meds"] if par_med[n]["sev"] == 0]
if tous:
    txt = (f"{len(tous)} médicament{'s' if len(tous) > 1 else ''} avec une recommandation pharmacogénétique "
           f"directement actionnable dans le périmètre testé : {', '.join(tous)}. "
           + (f"<b>Contre-indication : {', '.join(ci)}.</b> " if ci else "")
           + "Les adaptations proposées doivent être intégrées au contexte clinique, aux interactions "
             "médicamenteuses et, le cas échéant, au suivi thérapeutique pharmacologique.")
else:
    txt = ("Aucune recommandation pharmacogénétique directement actionnable n'est identifiée dans le périmètre "
           "testé. Ce résultat ne dispense pas de la surveillance clinique habituelle.")
h.append(Paragraph(txt, st_corps))

# --- 4. perimetre et limites
h.append(Paragraph("PÉRIMÈTRE ET LIMITES", st_section))
limites = [
    "Restitution limitée au sous-ensemble néphrologie et épilepsie du panel socle RNPGx 2026 : ce compte rendu "
    "n'est pas l'inventaire exhaustif des résultats calculés à partir du génome, disponibles au format structuré.",
    "Ne sont pas restitués : les recommandations optionnelles ou non graduées, les diplotypes ambigus et les "
    "divergences entre référentiels. Le référentiel appliqué est CPIC.",
    "Les recommandations sont versionnées et peuvent évoluer ; le génotype reste disponible au format structuré "
    "pour une réinterprétation ultérieure.",
    "Les remaniements de structure ne sont analysés que pour CYP2D6. Seuls les allèles répertoriés (PharmVar) sont "
    "recherchés : un résultat normal n'exclut pas un déficit.",
]
if RESERVES:
    limites.append("Sous réserve pour " + ", ".join(f"{g} ({n} position{'s' if n > 1 else ''} non lue{'s' if n > 1 else ''})"
                                                     for g, n in sorted(RESERVES.items()))
                   + " : un allèle défini par l'une de ces positions n'aurait pas été vu.")
for x in limites:
    h.append(Paragraph("• " + x, st_petit))

# --- validation, references, methode
h.append(Spacer(1, 7))
h.append(Paragraph("Compte rendu validé par : ...............................................  (biologiste médical)"
                   "      Date : ....../....../..........", st_corps))
refs = sorted({(" / ".join(l["noms"]), p) for l in lignes for p in l["pmid"]})
h.append(Spacer(1, 4))
h.append(HRFlowable(width="100%", thickness=0.5, color=TRAIT, spaceAfter=3))
if refs:
    h.append(Paragraph("<b>Références.</b> Recommandations CPIC : " + " ; ".join(f"{m} (PMID {p})" for m, p in refs)
                       + ".", st_petit))
h.append(Paragraph(
    "<b>Méthode.</b> Séquençage du génome entier, lectures courtes appariées, alignement sur GRCh38. CYP2D6 par Cyrius "
    "sur alignement complet ; HLA de classe I par OptiType ; autres gènes et interprétation par PharmCAT. "
    f"Traductions des recommandations CPIC : version {e_(TR.get('version', '?'))}"
    + (" (proposition, à valider)." if TR.get("statut") != "validée" else " (validée).")
    + (f" Texte d'origine conservé pour : {', '.join(NON_TRADUITS)}." if NON_TRADUITS else ""), st_petit))
if proto:
    h.append(Paragraph("<b>PROTOTYPE.</b> Établi sur un génome public du projet 1000 Genomes, à des fins de mise au "
                       "point. Ce document n'est pas un compte rendu de biologie médicale et ne concerne aucun patient.",
                       st_petit))

SimpleDocTemplate(OUT, pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm, topMargin=12 * mm, bottomMargin=10 * mm,
                  title=f"Compte rendu pharmacogénétique {ECH}").build(h)
print(f"{ECH} : {len(MESURES)} genes rendus, {len(tous)} medicaments restitues"
      + (f", {len(NON_TRADUITS)} texte(s) non traduit(s)" if NON_TRADUITS else "") + f" -> {os.path.basename(OUT)}")
