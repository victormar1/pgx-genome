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
import argparse, csv, html, json, os, re, sys, datetime
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cr_logique as L   # la logique, testable, sans mise en page
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
_p.add_argument("--typage-complementaire", default=None,
                help="pypgx.tsv produit par l'etage 4b ; par defaut cherche a cote du perimetre")
_a = _p.parse_args()
SRC, OUT = _a.rapport, _a.sortie
ECH = _a.echantillon or os.path.basename(SRC).split(".")[0]
os.makedirs(os.path.dirname(os.path.abspath(OUT)) or ".", exist_ok=True)
TR = json.load(open(_a.traductions, encoding="utf-8"))
MEDS = TR["medicaments_restitues"]


def _perimetre(rapport):
    global META
    META = json.load(open(_a.perimetre, encoding="utf-8"))
    return L.perimetre(rapport, META)


META = {}
rep = json.load(open(SRC, encoding="utf-8"))
genes = rep["genes"]
MESURES, ABSENTS, RESERVES = _perimetre(rep)
COMPLEMENT = META.get("complement_rnpgx") or []

# Typage complementaire (etage 4b). MT-RNR1 est deja restitue plus haut : son
# allele est passe par l'interpreteur, qui lui attache la recommandation CPIC.
# Les trois autres genes n'ont pas de table clinique, et ne sortent donc qu'avec
# un allele nomme, dit comme tel.
_tc = _a.typage_complementaire or os.path.join(os.path.dirname(os.path.abspath(_a.perimetre)), "pypgx.tsv")
TYPAGE = []
if os.path.exists(_tc):
    with open(_tc, encoding="utf-8") as _fh:
        TYPAGE = list(csv.DictReader(_fh, delimiter="\t"))

# ---------------------------------------------- libelles et logique
# Tables et regles vivent dans cr_logique ; on les reprend sous leurs noms
# d'origine pour que la mise en page ci-dessous reste inchangee.
PHENO_FR, CYP3A5_EXPR = L.PHENO_FR, L.CYP3A5_EXPR
FORCE_FR, RESTITUES = L.FORCE_FR, L.RESTITUES
STATINES, AMINOSIDES, NORMAUX = L.STATINES, L.AMINOSIDES, L.NORMAUX
e_, norm = L.e_, L.norm
ordre_med = lambda n: (n == "fosphenytoin", MEDS[n])   # phenytoine avant fosphenytoine


diplotypes = lambda sym: L.diplotypes(genes, sym)
libelle = lambda sym: L.libelle(genes, sym)
phenotype = lambda sym: L.phenotype(genes, sym)
gene_normal = lambda sym: L.gene_normal(genes, sym)
AMBIGUS = L.ambigus(genes, MESURES)


# --------------------------------------------------- recommandations restituees
brutes = L.recommandations(rep, genes, MESURES, AMBIGUS, TR, MEDS)
par_med = L.par_medicament(brutes)

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
    if nom in AMINOSIDES and sum(n in AMINOSIDES for n in par_med) > 1:
        grp = sorted((n for n in par_med if n in AMINOSIDES), key=lambda n: MEDS[n])
        vus |= set(grp)
        lignes.append({"noms": ["aminosides"], "sev": min(par_med[n]["sev"] for n in grp),
                       "action": e_(" ".join(d["fr"]))
                                 + "<br/><font color='#666666'>Concerne : "
                                 + e_(", ".join(MEDS[n] for n in grp)) + ".</font>",
                       "genes": set().union(*(par_med[n]["genes"] for n in grp)),
                       "source": " ; ".join(d["source"]),
                       "pmid": set().union(*(par_med[n]["pmid"] for n in grp)),
                       "ok": all(par_med[n]["ok"] for n in grp), "meds": grp})
        continue
    freres = [n for n, x in par_med.items() if n not in vus and n not in STATINES
              and n not in AMINOSIDES
              and x["fr"] == d["fr"] and x["genes"] == d["genes"] and x["source"] == d["source"]]
    vus |= set(freres)
    lignes.append({"noms": [MEDS[n] for n in sorted(freres, key=ordre_med)], "sev": d["sev"],
                   "action": e_(" ".join(d["fr"])), "genes": d["genes"], "source": " ; ".join(d["source"]),
                   "pmid": set().union(*(par_med[n]["pmid"] for n in freres)), "ok": d["ok"],
                   "meds": sorted(freres, key=ordre_med)})

MED_PAR_GENE = {}
for l in lignes:
    for g in l["genes"]:
        # Une ligne groupee est citee par son nom de classe dans le tableau des
        # genes : y reciter les onze molecules noierait la conclusion.
        if len(l["noms"]) == 1 and len(l["meds"]) > 1:
            MED_PAR_GENE.setdefault(g, []).append((l["sev"], l["noms"][0]))
            continue
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
_TYPE_RENDU = {r["gene"]: r for r in TYPAGE if r.get("statut") == "rendu"}
for g, motif in ABSENTS:
    # Un gene absent des tables de l'interpreteur mais type par l'etage
    # complementaire a bien ete mesure : le declarer « non mesure » ici et
    # publier son genotype deux sections plus bas serait se contredire.
    if g in _TYPE_RENDU:
        rows.append([Paragraph(g, st_gras), Paragraph(e_(_TYPE_RENDU[g]["genotype"]), st_corps),
                     Paragraph("—", st_corps),
                     Paragraph(f"Génotypé hors interpréteur : {e_(motif)}. Aucun phénotype ni "
                               "recommandation ne lui est attaché ; voir le typage complémentaire.",
                               st_corps)])
        continue
    rows.append([Paragraph(g, st_gras), Paragraph("non analysé", st_corps), Paragraph("—", st_corps),
                 Paragraph(f"Aucune conclusion : {e_(motif)}. Ni normal ni anormal, non mesuré.", st_corps)])
h.append(tableau(rows, [17, 32, 45, 76], {i: FOND_G for i in range(2, len(rows), 2)}))

# --- 2 bis. positions du core panel RNPGx hors definitions de l'interpreteur
# Elles sont mesurees pour repondre a l'exigence d'extraction des ROI, mais aucune
# n'est interpretee : on rend le genotype, jamais un phenotype.
# Seul MT-RNR1 voit ses positions retirees d'ici : son allele est restitue plus
# haut avec son phenotype et sa recommandation, le citer deux fois n'apporte
# rien.
#
# Pour les autres, les positions restent. L'etage 4b lit desormais les deux
# haplotypes et ne tait plus rien, mais le genotype mesure position par position
# ne depend d'aucun typeur : il reste la source qui fait foi, et c'est lui qui
# avait revele le silence du typeur sur rs1799807.
_TYPES_OK = {r["gene"] for r in TYPAGE
             if r.get("statut") == "rendu" and r.get("phenotype")}
COMPLEMENT = [x for x in COMPLEMENT if x["gene"] not in _TYPES_OK]
if COMPLEMENT:
    ref = [x for x in COMPLEMENT if x["genotype"] == "référence"]
    var = [x for x in COMPLEMENT if x not in ref and x["statut"] in ("mesuree", "couverture seule")]
    perdu = [x for x in COMPLEMENT if x["statut"] not in ("mesuree", "reference, lue", "couverture seule")]
    h.append(Paragraph("POSITIONS DU CORE PANEL RNPGx HORS PÉRIMÈTRE INTERPRÉTÉ", st_section))
    if var:
        rows = [[Paragraph(x, st_tete) for x in ("Gène", "Variant", "Classe", "Génotype", "Commentaire")]]
        for x in sorted(var, key=lambda x: (x["classe"] or 9, x["gene"])):
            rows.append([Paragraph(f"<b>{e_(x['gene'])}</b>", st_gras), Paragraph(e_(x["rsid"]), st_corps),
                         Paragraph(str(x["classe"] or "—"), st_corps), Paragraph(e_(x["genotype"]), st_gras),
                         Paragraph(e_(x["note"]), st_petit)])
        h.append(tableau(rows, [20, 24, 14, 22, 90]))
    resume = (f"{len(ref)} position{'s' if len(ref) > 1 else ''} à l'état de référence"
              if ref else "") + (" ; " if ref and perdu else "") + (
              f"{len(perdu)} non mesurée{'s' if len(perdu) > 1 else ''} "
              + "(" + ", ".join(sorted({x['rsid'] for x in perdu})) + ")" if perdu else "")
    h.append(Paragraph(
        (resume + ". " if resume else "")
        + "Ces positions sont extraites et contrôlées au titre du core panel RNPGx 2026, mais "
          "<b>aucune n'est interprétée</b> : le génotype est rendu, sans phénotype ni recommandation. "
        + "Le VNTR de TYMS n'est pas génotypable en lectures courtes : seule sa couverture est mesurée.", st_petit))

# --- 2 ter. typage complementaire : genes hors de l'interpreteur
# MT-RNR1 n'apparait pas ici : son allele a ete rendu a l'interpreteur, qui lui
# attache la recommandation CPIC, et il figure donc avec les autres genes. Les
# trois autres n'ont pas de table clinique : un allele nomme, et rien de plus.
_TYP = [r for r in TYPAGE if r.get("gene") != "MT-RNR1"]
if _TYP:
    h.append(Paragraph("TYPAGE COMPLÉMENTAIRE : ALLÈLE NOMMÉ, SANS INTERPRÉTATION", st_section))
    rows = [[Paragraph(x, st_tete) for x in ("Gène", "Classe", "Génotype", "Positions lues", "Commentaire")]]
    for r in sorted(_TYP, key=lambda r: (r.get("classe_rnpgx") or "9", r["gene"])):
        rendu = r.get("statut") == "rendu"
        rows.append([
            Paragraph(f"<b>{e_(r['gene'])}</b>", st_gras),
            Paragraph(e_(r.get("classe_rnpgx") or "—"), st_corps),
            Paragraph(e_(r.get("genotype") or "—"), st_gras if rendu else st_corps),
            Paragraph("%s/%s" % (r.get("positions_lues", "?"), r.get("positions_attendues", "?")), st_corps),
            Paragraph("Aucune interprétation clinique disponible pour ce gène." if rendu
                      else f"Non conclusif : {e_(r.get('motif') or 'non précisé')}.", st_petit)])
    h.append(tableau(rows, [20, 14, 34, 22, 80]))
    h.append(Paragraph(
        "Ces gènes appartiennent au core panel RNPGx mais sont absents des tables de l'interpréteur. "
        "L'allèle est nommé à partir des positions qui le définissent, <b>toutes vérifiées lues</b> sur "
        "l'alignement ; aucun phénotype ni aucune recommandation ne lui est attaché, faute de source "
        "clinique opposable. Un gène non conclusif n'est ni normal ni anormal : il n'a pas été mesuré.", st_petit))

_MT = next((r for r in TYPAGE if r.get("gene") == "MT-RNR1"), None)
if _MT and _MT.get("heteroplasmie") == "True":
    h.append(Paragraph(
        f"<b>MT-RNR1 : hétéroplasmie suspectée</b> — les deux haplotypes diffèrent"
        + (f", fraction allélique observée {_MT['fraction']}" if _MT.get("fraction") else "")
        + f". L'allèle non référence ({e_(_MT.get('allele_externe') or '')}) a été retenu : un porteur "
          "hétéroplasmique reste exposé. Le taux d'hétéroplasmie n'est pas quantifié par ce module et "
          "ne doit pas être interprété comme tel.", st_petit))

# --- 3. conclusion
h.append(Paragraph("CONCLUSION", st_section))
_cite = lambda l: ([l["noms"][0]] if len(l["noms"]) == 1 and len(l["meds"]) > 1
                   else [MEDS[n] for n in l["meds"]])
tous = [n for l in lignes for n in _cite(l)]
ci = [n for l in lignes if l["sev"] == 0 for n in _cite(l)]
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
