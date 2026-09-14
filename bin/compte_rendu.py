# -*- coding: utf-8 -*-
"""Compte rendu pharmacogenetique, une page, francais, depuis la sortie structuree de PharmCAT.

Deux regles de securite :
 1. seuls les genes REELLEMENT sequences sont rapportes. Les autres sortent en "non analyse",
    jamais en "normal", car l'option de remplissage les rendrait faussement rassurants ;
 2. les recommandations sont triees par force, les optionnelles sont comptees et non enumerees.

Usage : python generer_cr.py <echantillon>
"""
import argparse, json, os, sys, datetime
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable)

_p = argparse.ArgumentParser(description="Compte rendu pharmacogenetique, une page, francais.")
_p.add_argument("--rapport", required=True, help="le JSON produit par PharmCAT")
_p.add_argument("--perimetre", required=True, help="perimetre.json produit par le controle qualite")
_p.add_argument("--sortie", required=True, help="chemin du PDF a ecrire")
_p.add_argument("--echantillon", default="")
_a = _p.parse_args()
SRC, OUT = _a.rapport, _a.sortie
ECH = _a.echantillon or os.path.basename(SRC).split(".")[0]
os.makedirs(os.path.dirname(os.path.abspath(OUT)) or ".", exist_ok=True)


def _perimetre(rapport):
    """Croise ce que le controle qualite a mesure et ce que l interpreteur a rendu.

    Trois niveaux, parce que deux ne suffisent pas. Un gene dont une position sur
    cent cinquante n a pas ete lue n est ni sur ni ininterpretable : il est rendu,
    avec la reserve ecrite noir sur blanc.
      complet  toutes les positions diagnostiques lues -> rendu sans reserve
      partiel  au moins une position non lue           -> rendu avec reserve
      absent   aucune position lue, ou aucun diplotype -> non rendu
    Les genes venus d un outil dedie n ont pas de position dans le fichier de
    variants : leur perimetre est celui de cet outil, pas celui du controle qualite.
    """
    meta = json.load(open(_a.perimetre, encoding="utf-8"))
    qc = meta.get("genes", {})
    # Le perimetre clinique fixe ce que le compte rendu rapporte. Un gene
    # hors de ce perimetre, meme lu, n est pas un resultat rendu : il est
    # nomme a part. Sans liste (ancien perimetre.json), on garde tout.
    clinique = set(meta.get("perimetre_clinique") or [])
    hors_interpreteur = list(meta.get("perimetre_clinique_hors_interpreteur") or [])
    rendus, reserves, absents, hors_perimetre = [], {}, [], []
    for sym, e in rapport.get("genes", {}).items():
        if clinique and sym not in clinique:
            dips0 = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
            if any(x.get("label") and "Unknown" not in x.get("label") for x in dips0):
                hors_perimetre.append(sym)
            continue
        dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
        lab = ";".join(sorted({x.get("label", "") for x in dips if x.get("label")}))
        appele = bool(lab) and "Unknown" not in lab
        externe = e.get("callSource") == "OUTSIDE"
        v = qc.get(sym)
        statut = (v or {}).get("statut")
        if not appele:
            absents.append((sym, "aucun diplotype rendu" if v else "non appelable depuis un fichier de variants"))
        elif externe or v is None:
            rendus.append(sym)
        elif statut == "complet":
            rendus.append(sym)
        elif statut == "partiel":
            rendus.append(sym)
            reserves[sym] = v["positions_perdues"]
        else:
            absents.append((sym, "aucune position lue"))
    vus = set(rapport.get("genes", {}))
    for sym in sorted(clinique - vus):
        motif = ("gene absent des tables de l interpreteur"
                 if sym in hors_interpreteur else "non rendu par l interpreteur")
        absents.append((sym, motif))
    return sorted(rendus), sorted(absents), reserves, sorted(hors_perimetre)


NOIR, GRIS = colors.HexColor("#1a1a1a"), colors.HexColor("#666666")
TRAIT, ROUGE = colors.HexColor("#cccccc"), colors.HexColor("#a02020")
FOND_R, FOND_G = colors.HexColor("#fdf0f0"), colors.HexColor("#f4f4f2")

def S(n, **k):
    base = dict(fontName="Helvetica", fontSize=8.4, leading=10.8, textColor=NOIR, alignment=TA_LEFT)
    base.update(k)
    return ParagraphStyle(n, **base)

st_titre   = S("t", fontName="Helvetica-Bold", fontSize=13, leading=15)
st_ident   = S("i", fontSize=8.2, textColor=GRIS)
st_section = S("s", fontName="Helvetica-Bold", fontSize=8.8, leading=10.5, spaceBefore=7, spaceAfter=3)
st_alerte  = S("a", fontName="Helvetica-Bold", fontSize=10.5, leading=13, textColor=ROUGE)
st_corps   = S("c", fontSize=8.0, leading=9.9)
st_petit   = S("p", fontSize=7.4, leading=9.2, textColor=GRIS)
st_gene    = S("g", fontName="Helvetica-Bold", fontSize=8.0, leading=9.9)

PHENO_FR = {
    "Normal Metabolizer": "métaboliseur normal", "Intermediate Metabolizer": "métaboliseur intermédiaire",
    "Poor Metabolizer": "métaboliseur lent", "Rapid Metabolizer": "métaboliseur rapide",
    "Ultrarapid Metabolizer": "métaboliseur ultrarapide", "Normal Function": "fonction normale",
    "Decreased Function": "fonction diminuée", "Poor Function": "fonction faible",
    "Increased Function": "fonction augmentée", "Possible Intermediate Metabolizer": "métaboliseur intermédiaire possible",
    "Indeterminate": "non interprétable", "No Result": "non appelé",
    "Uncertain Susceptibility": "susceptibilité incertaine à l'hyperthermie maligne",
    "Malignant Hyperthermia Susceptibility": "susceptibilité à l'hyperthermie maligne",
    "ivacaftor non-responsive in CF patients": "pas de réponse attendue à l'ivacaftor",
    "ivacaftor responsive in CF patients": "réponse attendue à l'ivacaftor",
    "Normal": "activité normale", "Deficient": "déficit",
    "Deficient with CNSHA": "déficit avec anémie hémolytique chronique",
    "Variable": "activité variable", "Favorable Response Genotype": "génotype de réponse favorable",
    "Unfavorable Response Genotype": "génotype de réponse défavorable",
    "normal risk of aminoglycoside-induced hearing loss":
        "risque usuel de surdité sous aminoside",
    "increased risk of aminoglycoside-induced hearing loss":
        "risque accru de surdité sous aminoside",
    "n/a": "-",
}
MED_FR = {
    "carbamazepine": "carbamazépine", "oxcarbazepine": "oxcarbazépine", "phenytoin": "phénytoïne",
    "fosphenytoin": "fosphénytoïne", "lamotrigine": "lamotrigine", "clopidogrel": "clopidogrel",
    "tacrolimus": "tacrolimus", "tramadol": "tramadol", "codeine": "codéine",
    "azathioprine": "azathioprine", "mercaptopurine": "mercaptopurine", "thioguanine": "thioguanine",
    "allopurinol": "allopurinol", "voriconazole": "voriconazole", "omeprazole": "oméprazole",
    "pantoprazole": "pantoprazole", "lansoprazole": "lansoprazole", "dexlansoprazole": "dexlansoprazole",
    "esomeprazole": "ésoméprazole", "rabeprazole": "rabéprazole", "amitriptyline": "amitriptyline",
    "citalopram": "citalopram", "escitalopram": "escitalopram", "sertraline": "sertraline",
    "paroxetine": "paroxétine", "fluvoxamine": "fluvoxamine", "venlafaxine": "venlafaxine",
    "nortriptyline": "nortriptyline", "imipramine": "imipramine", "clomipramine": "clomipramine",
    "doxepin": "doxépine", "trimipramine": "trimipramine", "desipramine": "désipramine",
    "atomoxetine": "atomoxétine", "vortioxetine": "vortioxétine", "tamoxifen": "tamoxifène",
    "ondansetron": "ondansétron", "metoprolol": "métoprolol", "flecainide": "flécaïnide",
    "propafenone": "propafénone", "simvastatin": "simvastatine", "atorvastatin": "atorvastatine",
    "rosuvastatin": "rosuvastatine", "pravastatin": "pravastatine", "fluvastatin": "fluvastatine",
    "lovastatin": "lovastatine", "pitavastatin": "pitavastatine", "hydrocodone": "hydrocodone",
    "oxycodone": "oxycodone", "mavacamten": "mavacamten", "pimozide": "pimozide",
    "zuclopenthixol": "zuclopenthixol", "eliglustat": "éliglustat", "aripiprazole": "aripiprazole",
    "risperidone": "rispéridone", "brexpiprazole": "brexpiprazole", "quetiapine": "quétiapine",
}
FORCE_FR = {"Strong": "forte", "Moderate": "modérée", "Optional": "optionnelle",
            "No recommendation": "aucune", "Unspecified": "non précisée", "n/a": ""}
GENE_MED = {  # medicaments a ne retenir que si le gene correspondant est mesure
    "carbamazepine": "HLA-B", "oxcarbazepine": "HLA-B", "phenytoin": "HLA-B", "fosphenytoin": "HLA-B",
    "lamotrigine": "HLA-B", "allopurinol": "HLA-B", "abacavir": "HLA-B",
}

# traduction des conduites, par motif d'ouverture de la recommandation d'origine
CONDUITE_FR = [
 ("Avoid codeine use because of possibility of diminished analgesia",
  "Éviter la codéine : analgésie insuffisante attendue. Si un opioïde est nécessaire, en choisir un autre que le tramadol."),
 ("Avoid codeine use because of potential for serious toxicity",
  "Éviter la codéine : risque de toxicité opioïde. Si un opioïde est nécessaire, en choisir un autre que le tramadol."),
 ("Avoid tramadol use because of possibility of diminished analgesia",
  "Éviter le tramadol : analgésie insuffisante attendue. Si un opioïde est nécessaire, en choisir un autre que la codéine."),
 ("Avoid tramadol use because of potential for serious toxicity",
  "Éviter le tramadol : risque de toxicité opioïde. Si un opioïde est nécessaire, en choisir un autre que la codéine."),
 ("Avoid tricyclic use due to potential for side effects",
  "Éviter les tricycliques : risque d'effets indésirables. Préférer un médicament non métabolisé par CYP2D6. Si un tricyclique s'impose, réduire de moitié la dose initiale et suivre les concentrations plasmatiques."),
 ("Avoid tricyclic use due to potential for lack of efficacy",
  "Éviter les tricycliques : efficacité insuffisante attendue. Préférer un médicament non métabolisé par CYP2D6. Si un tricyclique s'impose, suivre les concentrations plasmatiques."),
 ("Recommend alternative hormonal therapy such as an aromatase inhibitor",
  "Préférer une autre hormonothérapie : inhibiteur de l'aromatase après la ménopause, ou inhibiteur de l'aromatase avec suppression ovarienne avant la ménopause."),
 ("Select a therapy other than tamoxifen",
  "Choisir un traitement autre que le tamoxifène."),
 ("If patient is carbamazepine-na", "Ne pas utiliser si le patient n'en a jamais reçu."),
 ("If patient is oxcarbazepine-na", "Ne pas utiliser si le patient n'en a jamais reçu."),
 ("If patient is phenytoin-naive", "Ne pas utiliser si le patient n'en a jamais reçu. Éviter aussi carbamazépine et oxcarbazépine."),
 ("Avoid standard dose (75 mg) clopidogrel", "Éviter la dose standard de 75 mg. Préférer prasugrel ou ticagrélor."),
 ("Avoid clopidogrel if possible", "Éviter le clopidogrel. Préférer prasugrel ou ticagrélor."),
 ("Consider an alternative P2Y12", "Envisager un autre inhibiteur P2Y12 à dose standard."),
 ("Initiate therapy with reduced starting doses (30-80%", "Dose initiale réduite à 30 ou 80 % si la dose usuelle dépasse 2 mg/kg/j. Surveiller l'hémogramme."),
 ("Initiate therapy with decreased starting doses (30-80%", "Dose initiale réduite à 30 ou 80 % de la dose usuelle. Surveiller l'hémogramme."),
 ("Increase starting dose 1.5 to 2 times", "Dose initiale de 1,5 à 2 fois la dose usuelle, sans dépasser 0,3 mg/kg/j. Suivi thérapeutique."),
 ("Initiate standard starting daily dose. Consider increasing dose by 50-100%",
  "Dose standard. Augmentation de 50 à 100 % à envisager pour H. pylori ou une œsophagite érosive."),
 ("Initiate standard starting daily dose. For chronic therapy",
  "Dose standard. Au-delà de douze semaines et si l'efficacité est obtenue, réduction de 50 % envisageable."),
 ("Initiate therapy with recommended starting dose. Consider a slower titration",
  "Dose initiale usuelle, titration plus lente et dose d'entretien plus basse."),
 ("Consider a lower starting dose, slower titration", "Dose initiale plus basse, titration plus lente, entretien réduit de moitié."),
 ("Consider a clinically appropriate antidepressant not predominantly metabolized by CYP2C19",
  "Préférer un antidépresseur peu métabolisé par CYP2C19."),
 ("Choose an alternative agent that is not dependent on CYP2C19",
  "Préférer un antifongique indépendant de CYP2C19."),
 ("Consider a 25% reduction of recommended starting dose", "Réduire la dose initiale d'environ 25 %, avec suivi thérapeutique."),
 ("Avoid tertiary amine use", "Éviter les amines tertiaires. Préférer un antidépresseur non métabolisé par CYP2C19."),
 ("Avoid amitriptyline use", "Éviter l'amitriptyline. Si nécessaire, réduire la dose initiale de moitié."),
 ("Consider hormonal therapy such as an aromatase inhibitor",
  "Envisager un inhibiteur de l'aromatase. Éviter les inhibiteurs du CYP2D6."),
 ("Initiate with a dose of 40 mg/day", "Débuter à 40 mg/j, puis adapter selon la réponse et la concentration plasmatique."),
 ("Initiate with a dose of 0.5 mg/kg/day", "Débuter à 0,5 mg/kg/j, puis adapter selon la réponse et la concentration plasmatique."),
 ("Prescribe an alternative statin", "Préférer une autre statine, selon la puissance recherchée."),
 ("Prescribe ", "Plafonner la dose initiale et adapter selon les recommandations de la pathologie."),
 ("Use tramadol label recommended", "Dose selon le résumé des caractéristiques du produit. Si inefficace, changer d'opioïde."),
 ("Use codeine label recommended", "Dose selon le résumé des caractéristiques du produit. Si inefficace, changer d'opioïde."),
]
def conduite_fr(txt):
    for motif, trad in CONDUITE_FR:
        if txt.startswith(motif):
            return trad
    return None

fr_ph  = lambda p: PHENO_FR.get(p, p)
fr_med = lambda m: MED_FR.get(m, m)

rep = json.load(open(SRC, encoding="utf-8"))
genes = rep["genes"]
MESURES, _NA, RESERVES, HORS_PERIMETRE = _perimetre(rep)
NON_ANALYSES = [x[0] for x in _NA]
MOTIF = {g: m for g, m in _NA}

def appel(sym):
    e = genes.get(sym) or {}
    dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
    lab = ";".join(sorted({x.get("label", "") for x in dips if x.get("label")}))
    phen = sorted({p for x in dips for p in (x.get("phenotypes") or [])})
    return lab, phen, e.get("callSource", "")

NORMAL = {"Normal Metabolizer", "Normal Function", "*31:01 negative", "*15:02 negative",
          "*57:01 negative", "*58:01 negative", "*15:11 negative", "*13:01 negative"}

# --------------------------------------------------------- recommandations retenues
retenues = []       # (force, gene, medicament, conduite, source)
positifs_hla = []
for sym in MESURES:
    lab, phen, _ = appel(sym)
    for p in phen:
        if "positive" in p:
            positifs_hla.append((sym, p.replace(" positive", "")))

for section, meds in rep["drugs"].items():
    if not section.startswith(("CPIC", "DPWG")):
        continue
    for nom, m in meds.items():
        g_requis = GENE_MED.get(nom)
        if g_requis and g_requis not in MESURES:
            continue
        for gl in (m.get("guidelines") or []):
            for a in (gl.get("annotations") or []):
                if not (a.get("alternateDrugAvailable") or a.get("dosingInformation")):
                    continue
                impl = a.get("implications") or []
                if isinstance(impl, list):
                    impl = " ; ".join(str(x) for x in impl)
                # ne garder que si l'implication cite un gene mesure
                if not any(g in impl for g in MESURES):
                    continue
                cl = a.get("classification", "")
                rec = str(a.get("drugRecommendation", "")).replace("&gt;", ">").replace("&quot;", '"')
                retenues.append((cl, nom, rec, section[:4], impl))

ORDRE = {"Strong": 0, "Moderate": 1, "Unspecified": 2, "Optional": 3, "No recommendation": 4}
retenues.sort(key=lambda r: (ORDRE.get(r[0], 5), r[1]))
fortes = [r for r in retenues if r[0] == "Strong"]
moderees = sorted({r[1] for r in retenues if r[0] == "Moderate"} - {r[1] for r in fortes})
autres = sorted({r[1] for r in retenues if r[0] not in ("Strong", "Moderate")}
                - {r[1] for r in fortes} - set(moderees))

# ------------------------------------------------------------------ document
h = []
h.append(Paragraph("COMPTE RENDU PHARMACOGÉNÉTIQUE", st_titre))
h.append(Paragraph(
    f"Échantillon <b>{ECH}</b> &nbsp;·&nbsp; analyse du {datetime.date.today().strftime('%d/%m/%Y')}"
    f" &nbsp;·&nbsp; génome entier, lectures courtes &nbsp;·&nbsp; {len(MESURES)} gènes rendus"
    + (f", {len(NON_ANALYSES)} non analysé" + ("s" if len(NON_ANALYSES) > 1 else "")
       if NON_ANALYSES else ""), st_ident))
h.append(Spacer(1, 3))
h.append(HRFlowable(width="100%", thickness=1.1, color=NOIR, spaceAfter=6))

# Chaque allele a risque porte sa propre contre-indication : les nommer toutes
# ensemble ferait dire au compte rendu ce qui ne vaut que pour l'une d'elles.
RISQUE_HLA = {
    "*15:02": ("carbamazépine, oxcarbazépine, phénytoïne et fosphénytoïne",
               "syndrome de Stevens-Johnson et nécrolyse épidermique toxique"),
    "*15:11": ("carbamazépine et oxcarbazépine",
               "syndrome de Stevens-Johnson et nécrolyse épidermique toxique"),
    "*31:01": ("carbamazépine",
               "réaction cutanée grave, dont le syndrome DRESS"),
    "*57:01": ("abacavir",
               "syndrome d'hypersensibilité à l'abacavir"),
    "*58:01": ("allopurinol",
               "syndrome d'hypersensibilité sévère, dont le syndrome DRESS"),
    "*13:01": ("dapsone",
               "syndrome d'hypersensibilité à la dapsone"),
}
if positifs_hla:
    phrases = []
    for g, a in positifs_hla:
        med, risque = RISQUE_HLA.get(a, ("", ""))
        if med:
            phrases.append(f"<b>{g}{a} présent.</b> Ne pas prescrire {med}. Risque de {risque}.")
        else:
            phrases.append(f"<b>{g}{a} présent.</b> Allèle à risque signalé par le référentiel ; "
                           "se reporter à la recommandation correspondante.")
    corps = " ".join(phrases) + (
        " Chez un patient traité depuis plus de trois mois sans réaction, le risque "
        "d'apparition tardive est faible." if any(a in ("*15:02", "*31:01") for _, a in positifs_hla) else "")
    t = Table([[Paragraph("CONTRE-INDICATION", st_alerte)], [Paragraph(corps, st_corps)]],
              colWidths=[170 * mm])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), FOND_R),
                           ("BOX", (0, 0), (-1, -1), 1.2, ROUGE),
                           ("LEFTPADDING", (0, 0), (-1, -1), 8), ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                           ("TOPPADDING", (0, 0), (-1, 0), 6), ("BOTTOMPADDING", (0, -1), (-1, -1), 6)]))
    h.append(t); h.append(Spacer(1, 7))

# --- genotypes
h.append(Paragraph("GÉNOTYPES", st_section))
lignes = [[Paragraph("<b>Gène</b>", st_petit), Paragraph("<b>Diplotype</b>", st_petit),
           Paragraph("<b>Interprétation</b>", st_petit)]]
for sym in MESURES:
    lab, phen, src = appel(sym)
    if not lab:
        continue
    # Un libelle long n'est jamais resume par une affirmation de contenu : il est
    # dit ambigu, ou coupe. Ecrire "genotype de reference" sur un diplotype variant
    # inverserait le sens du compte rendu.
    _ref = {"Reference/Reference", "B (reference)/B (reference)", "*1/*1"}
    if lab in _ref:
        court = "génotype de référence"
    elif ";" in lab:
        court = "plusieurs diplotypes compatibles : " + lab.split(";")[0].strip() + ", …"
    elif len(lab) > 46:
        court = lab[:44] + "…"
    else:
        court = lab
    ph = ", ".join(fr_ph(p) for p in phen if "negative" not in p) or \
         ("allèles à risque absents" if phen else "-")
    if "positive" in " ".join(phen):
        ph = ", ".join(p.replace(" positive", " présent").replace(" negative", " absent") for p in phen)
    if sym in RESERVES:
        court += " <font color='#a02020'>·</font>"
    lignes.append([Paragraph(sym, st_gene), Paragraph(court, st_corps), Paragraph(ph, st_corps)])
t = Table(lignes, colWidths=[24 * mm, 62 * mm, 84 * mm])
t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                       ("LINEBELOW", (0, 0), (-1, 0), 0.6, TRAIT),
                       ("LINEBELOW", (0, -1), (-1, -1), 0.6, TRAIT),
                       ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, FOND_G]),
                       ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                       ("TOPPADDING", (0, 0), (-1, -1), 0.8), ("BOTTOMPADDING", (0, 0), (-1, -1), 0.8)]))
h.append(t)

# --- conduites
if fortes:
    h.append(Paragraph("CONDUITES À TENIR", st_section))
    vus, lignes, non_traduites = {}, [], []
    for cl, nom, rec, src, impl in fortes:
        if vus.get(nom) == rec:
            continue          # meme medicament, meme texte : un seul affichage
        t = conduite_fr(rec)
        if t is None:
            # jamais de silence sur une recommandation forte : on rend le texte source
            t = (rec[:400] + "…") if len(rec) > 400 else rec
            non_traduites.append(nom)
        vus[nom] = rec
        lignes.append([Paragraph(f"<b>{fr_med(nom)}</b>", st_corps),
                       Paragraph(t, st_corps),
                       Paragraph(f"{src} {FORCE_FR.get(cl, cl)}", st_petit)])
    t = Table(lignes, colWidths=[30 * mm, 112 * mm, 28 * mm])
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                           ("LINEBELOW", (0, 0), (-1, -2), 0.4, TRAIT),
                           ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                           ("TOPPADDING", (0, 0), (-1, -1), 1.4), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.4)]))
    h.append(t)
    if non_traduites:
        h.append(Paragraph(
            "Texte source conservé, faute de traduction validée, pour : "
            + ", ".join(fr_med(x) for x in sorted(set(non_traduites))) + ".", st_petit))

if moderees:
    h.append(Spacer(1, 3))
    h.append(Paragraph(
        "<b>Recommandation modérée</b> du référentiel sur ces mêmes gènes, à lire dans le fichier "
        "structuré pour le détail de chacune : "
        + ", ".join(fr_med(m) for m in moderees) + ".", st_corps))
if autres:
    h.append(Paragraph(
        f"<b>Recommandation optionnelle ou non graduée</b> pour {len(autres)} autres médicaments : "
        + ", ".join(fr_med(m) for m in autres[:16]) + (", et d'autres." if len(autres) > 16 else "."), st_petit))

# --- limites
h.append(Paragraph("CE QUI N'A PAS ÉTÉ ANALYSÉ", st_section))
_liste = ", ".join(f"{g} ({MOTIF.get(g, '')})" for g in NON_ANALYSES) if NON_ANALYSES else "aucun gène du panel PharmCAT"
h.append(Paragraph(
    f"Gènes du panel restés sans résultat sur cet échantillon : <b>{_liste}</b>. "
    + ("Aucune conclusion ne peut en être tirée : il n'est ni normal ni anormal, il n'a pas "
       "été mesuré. " if len(NON_ANALYSES) == 1 else
       "Aucune conclusion ne peut en être tirée : ils ne sont ni normaux ni anormaux, ils "
       "n'ont pas été mesurés. ") + "Les remaniements de structure des gènes autres que CYP2D6 ne "
    "sont pas analysés. Seuls les allèles répertoriés au catalogue PharmVar sont recherchés : "
    "un variant délétère non répertorié donne un résultat de métaboliseur normal. Un résultat "
    "normal ne vaut donc pas exclusion d'un déficit."
    + ((" <b>Sous réserve</b> pour "
        + ", ".join(f"{g} ({n} position{'s' if n > 1 else ''} non lue{'s' if n > 1 else ''})"
                    for g, n in sorted(RESERVES.items()))
        + " : un allèle défini par l'une de ces positions n'aurait pas été vu.") if RESERVES else ""),
    st_corps))
if HORS_PERIMETRE:
    # L'interpréteur voit des gènes que ce compte rendu ne rapporte pas : ils
    # sortent du périmètre RNPGx retenu. Les nommer, sans les interpréter, évite
    # de laisser croire qu'ils ont été omis par erreur.
    h.append(Paragraph(
        "<b>Hors périmètre.</b> L'interpréteur produit aussi un résultat pour "
        + ", ".join(HORS_PERIMETRE) + ". Ces gènes ne font pas partie du sous-ensemble "
        "néphrologie et épilepsie retenu ici et ne sont pas rapportés : leur validation "
        "clinique n'a pas été conduite sur ce banc.", st_petit))
h.append(Spacer(1, 6))
h.append(HRFlowable(width="100%", thickness=0.6, color=TRAIT, spaceAfter=4))
h.append(Paragraph(
    "<b>Méthode.</b> Séquençage du génome entier, lectures courtes appariées, alignement sur GRCh38. "
    "CYP2D6 par Cyrius 1.1.1 sur alignement complet. HLA de classe I par OptiType 1.3.5. "
    "Autres gènes et interprétation par PharmCAT 3.4.0, référentiels CPIC et DPWG. "
    "Périmètre conforme au panel socle du RNPGx 2026, sous-ensemble néphrologie et épilepsie.", st_petit))
h.append(Spacer(1, 2))
h.append(Paragraph(
    "<b>PROTOTYPE.</b> Établi sur un génome public du projet 1000 Genomes, à des fins de mise au point. "
    "Ce document n'est pas un compte rendu de biologie médicale et ne concerne aucun patient.", st_petit))

SimpleDocTemplate(OUT, pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm,
                  topMargin=12 * mm, bottomMargin=10 * mm,
                  title=f"Compte rendu pharmacogénétique {ECH}").build(h)
print(f"{ECH} : {len(MESURES)} genes rendus, {len(fortes)} recommandations fortes -> {os.path.basename(OUT)}")
