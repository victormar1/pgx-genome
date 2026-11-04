# -*- coding: utf-8 -*-
"""Comparaison d'un diplotype rendu a la verite GeT-RM.

Un diplotype est une paire non ordonnee. Les suffixes de sous-allele ne comptent
pas (*2.001 vaut *2). GeT-RM ecrit parfois plusieurs hypotheses separees par
"or", ou entoure de parentheses un allele incertain : on accepte la concordance
avec l'une quelconque des hypotheses.
"""
import collections, itertools, re

# equivalences de nomenclature documentees
# Un ancien nom peut avoir plusieurs lectures actuelles : GeT-RM ecrit *1B la ou
# PharmVar a tantot *37, tantot *1 depuis la fusion de *1A et *1B.
REN = {
    "SLCO1B1": {"*1A": ["*1"], "*1B": ["*37", "*1"]},
    "CYP3A4":  {"*1B": ["*1"]},
    # *80 est en desequilibre de liaison avec *28 : PharmCAT ecrit l'haplotype
    # *80+*28 la ou GeT-RM ecrit *28 ou *28 + *60. *60 seul n'a pas d'effet connu.
    "UGT1A1":  {"*28": ["*80+*28", "*28"], "*28+*60": ["*80+*28"], "*60": ["*1"],
                "*80+*28": ["*80+*28", "*28"], "*36": ["*80+*36", "*36"]},
    "CYP4F2":  {},
}

# GeT-RM ecrit *1 pour tout haplotype SLCO1B1 de fonction normale, sans distinguer
# *1, *14, *20 ni *37. CPIC leur donne la meme fonction, la difference est de
# resolution et non de resultat.
EQUIVALENTS_VERITE = {
    "SLCO1B1": {"*1": ["*1", "*14", "*20", "*37"], "*1A": ["*1"], "*1B": ["*37", "*1"]},
    # PharmVar a fait de *38 la reference de CYP2C19 ; le *1 de GeT-RM couvre les deux
    "CYP2C19": {"*1": ["*1", "*38"]},
}

# genes dont GeT-RM utilise une nomenclature differente : la comparaison n'a pas de sens
NOMENCLATURE_DIVERGENTE = {"NAT2", "CYP1A2", "GSTP1", "SLCO2B1", "UGT2B15", "CYP2E1"}

# GeT-RM ecrit VKORC1 en genotype de rs9923231, brin cDNA
VKORC1 = {"GG": ("ref", "ref"), "GA": ("ref", "var"), "AG": ("ref", "var"), "AA": ("var", "var")}

def _vkorc1(d):
    d = (d or "").strip()
    if d.upper() in VKORC1:
        return [tuple(sorted(VKORC1[d.upper()]))]
    t = ("var" if "variant" in d else "ref", "var" if d.count("variant") == 2 else "ref")
    if "/" in d:
        a, b = d.split("/", 1)
        t = ("var" if "variant" in a else "ref", "var" if "variant" in b else "ref")
    return [tuple(sorted(t))]

def _sans_appel(paires):
    return [p for p in paires if not any("Unknown" in a or "No " in a for a in p)]


def _alleles(d):
    """Decoupe un diplotype en deux haplotypes, tandems conserves."""
    d = d.strip()
    if "/" not in d:
        return None
    a, b = d.split("/", 1)
    return [a.strip(), b.strip()]

# Le typage de reference est ecrit en trois champs (*02:01:01), le module en
# deux (*02:01). Le troisieme champ designe une substitution synonyme : il ne
# change ni la proteine ni la recommandation, et CPIC ecrit ses conduites en deux
# champs. Les deux ecritures se comparent donc a deux champs.
GENES_HLA = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1")


def _deux_champs(hap):
    m = re.match(r"^(\*?\d+:\d+)", hap.strip())
    return m.group(1) if m else hap.strip()


# UGT1A1 : *60 et *80 accompagnent *28, *36 et *37 sans fonction propre. La
# reference les ecrit, le module les ecrit ailleurs, et aucun des deux n'en tire
# de conduite. On compare donc les alleles porteurs.
ACCOMPAGNATEURS = {"UGT1A1": ("*60", "*80")}


def _sans_accompagnateurs(gene, hap):
    acc = ACCOMPAGNATEURS.get(gene)
    if not acc:
        return hap
    morceaux = [m for m in hap.split("+") if m.strip() not in acc]
    return "+".join(morceaux) if morceaux else "*1"


def canon(gene, hap):
    """Retourne la liste des ecritures actuelles possibles d'un haplotype."""
    if gene in GENES_HLA:
        return [_deux_champs(hap)]
    h = hap.strip().replace(" ", "").strip("()")
    h = re.sub(r"\.\d+$", "", h)            # sous-allele Aldy
    h = re.sub(r"\.\d+(?=\+)", "", h)
    h = _sans_accompagnateurs(gene, h)
    lectures = REN.get(gene, {}).get(h, [h])
    return sorted({_sans_accompagnateurs(gene, x) for x in lectures})

def variantes(gene, d):
    """Toutes les lectures possibles d'une chaine GeT-RM."""
    if not d: return []
    d = d.replace("’", "'").strip()
    if d.lower() in ("n/a", "nd", "none", "-", "unknown"): return []
    if gene == "VKORC1": return _vkorc1(d)
    if gene == "CYP2D6":
        # idiomes GeT-RM : *2(*35) designe l'allele *35 precise, (*68)+*4 le tandem
        # *68+*4, et *36x2+*10 deux copies de *36 suivies de *10
        d = re.sub(r"\*\d+\s*\(\s*(\*\d+[A-Za-z]*)\s*\)", lambda m: m.group(1), d)
        d = re.sub(r"\((\*\d+[A-Za-z]*)\)", lambda m: m.group(1), d)
        d = re.sub(r"(\*\d+[A-Za-z]*)x(\d+)",
                   lambda m: "+".join([m.group(1)] * int(m.group(2))), d)
    # un allele declare indetermine par la reference n'est pas un desaccord
    if re.search(r"undefined", d, re.I):
        return []
    # "(A;B)/C" : la reference hesite entre A et B pour ce seul allele. Le
    # point-virgule sert ailleurs a separer deux diplotypes entiers ; ici il est
    # a l'interieur des parentheses et ne separe qu'un haplotype. Sans ce
    # developpement, le decoupage sur le point-virgule perd l'hypothese A.
    m = re.search(r"\(([^()]*;[^()]*)\)", d)
    if m:
        alt = [x.strip() for x in m.group(1).split(";") if x.strip()]
        d = " or ".join(d[:m.start()] + "(" + x + ")" + d[m.end():] for x in alt)
    # GeT-RM ecrit "X/Y (Z)" quand un laboratoire a lu Z a la place de X ou de Y
    m = re.match(r"^\s*([^()]+?)\s*\(([^()]+)\)\s*$", d)
    if m and "/" in m.group(1):
        a, b = [x.strip() for x in m.group(1).split("/", 1)]
        z = m.group(2).strip()
        d = " or ".join([f"{a}/{b}", f"{a}/{z}", f"{z}/{b}"])
    # PharmCAT ecrit "Reference" la ou GeT-RM ecrit *1
    d = re.sub(r"(Reference|B \(reference\)|WT)", "*1", d)
    out = []
    for part in re.split(r"\s+or\s+|\s*;\s*", d):
        a = _alleles(part)
        if not a: continue
        for x in canon(gene, a[0]):
            for y in canon(gene, a[1]):
                out.append(tuple(sorted((x, y))))
    return sorted(set(out))

# Sans phasage, deux variants heterozygotes se lisent en cis ou en trans : les
# deux ecritures sont indiscernables. GeT-RM tranche par analyse de trio, pas nous.
AMBIGU_SANS_PHASE = {
    "CYP4F2": [({"*1", "*4"}, {"*2", "*3"})],
}

def _ambigu(gene, a, b):
    """Deux ecritures portent-elles les memes variants, autrement reparties ?

    Se deduit des definitions d'alleles plutot que d'une liste tenue a la main :
    si l'ensemble des variants est le meme des deux cotes, les deux diplotypes
    sont indiscernables sur un fichier non phase. La liste ne sert plus que de
    secours quand les definitions manquent pour le gene.
    """
    va, vb = _variants(gene, a), _variants(gene, b)
    if va is not None and vb is not None and va and va == vb:
        return True
    for x, y in AMBIGU_SANS_PHASE.get(gene, []):
        if (set(a) == x and set(b) == y) or (set(a) == y and set(b) == x):
            return True
    return False


# Univers des alleles que la methode de reference sait distinguer. Rempli par
# l'appelant depuis la table GeT-RM : un allele hors de cet univers ne peut pas
# etre contredit par elle, c'est une zone que son test ne couvre pas.
UNIVERS = {}

# Alleles que PharmCAT sait nommer, extraits de ses propres tables de traduction.
# Si la verite cite un allele absent de cette liste, l'outil ne pouvait pas le rendre :
# c'est une limite de perimetre, pas une erreur d'appel.
ALLELES_OUTIL = {}

# Variants qui definissent chaque allele, extraits des tables de traduction de
# PharmCAT. Sert a trancher le sens d'un desaccord : avons-nous vu un variant de
# plus que la methode de reference, ou un de moins ?
DEFINITIONS = {}


def _variants(gene, paire):
    d = DEFINITIONS.get(gene)
    if not d: return None
    s = set()
    for hap in paire:
        for morceau in hap.split("+"):
            v = d.get(morceau)
            if v is None: return None
            s.update(v)
    return s

def _elargir(gene, paires):
    """Etend les lectures de la verite aux alleles que sa nomenclature ne separe pas."""
    eq = EQUIVALENTS_VERITE.get(gene)
    if not eq: return paires
    out = set()
    for x, y in paires:
        for a in eq.get(x, [x]):
            for b in eq.get(y, [y]):
                out.add(tuple(sorted((a, b))))
    return sorted(out)


# ---------------------------------------------------------------------------
# Genes rendus en liste de variants plutot qu'en paire d'alleles.
#
# PharmCAT ecrit DPYD "c.85T>C (*9A)/c.1627A>G (*5)" ou "c.2846A>T
# (heterozygous)" ; GeT-RM tantot "*1/*9", tantot "c.1601G>A / c.2846A>T". Les
# deux ecritures designent le meme objet : le compte de chaque variant, 1 pour
# un heterozygote, 2 pour un homozygote. C'est cet objet qu'on compare, pas le
# nom du diplotype, qu'aucun des deux cotes n'ecrit de la meme facon.
RENDU_EN_VARIANTS = {"DPYD"}

_RIEN = re.compile(r"no\s+variants?\s+detected|^none$|^negative$|^wild\s*type$|^wt$"
                   r"|^\*1$|^reference$|^ref$", re.I)
_STRUCTUREL = re.compile(r"deletion|duplication|\bexon\b", re.I)
_INDEX = {}


def _index_hgvs(gene):
    """alias d'ecriture -> lectures possibles, chacune un tuple de positions."""
    if gene in _INDEX:
        return _INDEX[gene]
    d = DEFINITIONS.get(gene)
    if not d:
        return None
    idx = {}
    def pose(cle, t):
        idx.setdefault(cle, set()).add(t)
    for cle, pos in d.items():
        t = tuple(sorted(pos))
        pose(cle, t)
        m = re.match(r"^(.*?)\s*\(\*(\w+)\)\s*$", cle)
        if m:
            pose(m.group(1).strip(), t)
            pose("*" + m.group(2), t)
            base = re.match(r"^(\d+)", m.group(2))
            if base and base.group(1) != m.group(2):
                pose("*" + base.group(1), t)     # *9 couvre *9A et *9B
    _INDEX[gene] = {k: sorted(v) for k, v in idx.items()}
    return _INDEX[gene]


def _terme(mot, idx):
    """(lectures, etat). etat vaut "ok", "structurel" ou "inconnu"."""
    mot = re.sub(r"\s+", " ", mot.strip().strip(".").strip())
    mot = re.sub(r"\s*\((?:het|hom)\w*\)\s*$", "", mot, flags=re.I)
    mot = re.sub(r"\s+(?:het|hom|hemizygous)\.?$", "", mot, flags=re.I)
    if not mot or _RIEN.search(mot):
        return [()], "ok"
    if mot in idx:
        return idx[mot], "ok"
    sans = re.sub(r"\s*\(\*\w+\)\s*$", "", mot).strip()
    if sans in idx:
        return idx[sans], "ok"
    return None, ("structurel" if _STRUCTUREL.search(mot) else "inconnu")


def _lectures(texte, idx):
    """(liste de compteurs de positions, etat)."""
    if not texte or not texte.strip():
        return None, "inconnu"
    t = texte.replace("\u2019", "'").strip()
    if re.search(r"unknown|indeterminate|not\s+tested", t, re.I):
        return None, "inconnu"
    zyg = 2 if re.search(r"\(homozygous\)", t, re.I) else 1
    t = re.sub(r"\s*\((?:heterozygous|homozygous)\)\s*", " ", t, flags=re.I).strip()
    cotes = [c.strip() for c in t.split("/", 1)] if "/" in t else None

    def morceaux(c):
        out = []
        for m in re.split(r",", c):
            m = m.strip()
            if not m:
                continue
            # un terme entierement entre parentheses est declare incertain par
            # GeT-RM ; "(*4)" en est un, "c.85T>C (*9A)" n'en est pas un
            incertain = bool(re.match(r"^\([^()]*\)$", m))
            out.append((m.strip("()") if incertain else m, incertain))
        return out or [("", False)]

    groupes, etat = [], "ok"
    for c in (cotes if cotes is not None else [t]):
        for mot, incertain in morceaux(c):
            lec, e = _terme(mot, idx)
            if lec is None:
                return None, e
            groupes.append(lec + [()] if incertain else lec)
    sorties = []
    for choix in itertools.product(*groupes) if groupes else [()]:
        cpt = collections.Counter()
        for positions in choix:
            for p in positions:
                cpt[p] += (zyg if cotes is None else 1)
        sorties.append(cpt)
    return sorties, etat


def _compare_variants(gene, rendu, verite):
    idx = _index_hgvs(gene)
    if idx is None:
        return None
    lr, er = _lectures(rendu, idx)
    lv, ev = _lectures(verite, idx)
    if lr is None:
        return "non appele" if er == "inconnu" else "non representable"
    if lv is None:
        return "allele hors des tables de l outil" if ev == "inconnu" else "non representable"
    for a in lr:
        for b in lv:
            if a == b:
                return "concordant"
    # Meme regle directionnelle que pour les genes rendus en paires d'alleles :
    # porter un variant que la verite ne cite pas n'est pas une erreur d'appel.
    # Un "*1/*1" de GeT-RM veut dire "aucun variant parmi ceux que ce
    # laboratoire teste", et son panel n'est pas ecrit dans la table. Ce sens du
    # desaccord est donc sorti du denominateur, mais la fiche le nomme.
    for a in lr:
        for b in lv:
            if all(a.get(k, 0) >= v for k, v in b.items()):
                return "hors perimetre de la verite"
    return "DISCORDANT"


# Un allele entre parentheses : la reference l'a recu d'un laboratoire et pas des
# autres. Elle ne tranche pas, donc elle ne peut ni confirmer ni contredire.
# CYP2D6 est exclu : ses parentheses precisent un sous-allele, elles ne doutent
# pas.
INCERTAIN = re.compile(r"\(\s*\*")


def _incertitude_explique(gene, rendu, verite):
    """L'ecart tient-il au seul allele que la reference met entre parentheses ?

    Il faut deux choses : que la reference en marque un, et que le module
    retrouve celui qu'elle donne pour certain. Sinon le desaccord porte aussi sur
    l'allele ferme, et l'incertitude n'explique rien.
    """
    if gene == "CYP2D6" or not verite or not INCERTAIN.search(verite):
        return False
    cotes = _alleles(verite.replace(";", "/"))
    if not cotes:
        return False
    certains = [c for c in cotes if not re.match(r"^\s*\(", c)]
    if not certains:
        return False
    lus = set()
    for paire in variantes(gene, rendu):
        lus.update(paire)
    for c in certains:
        if set(canon(gene, c)) & lus:
            return True
    return False


def compare(gene, rendu, verite):
    """(verdict, forme comparee du rendu, formes acceptees)"""
    if gene in RENDU_EN_VARIANTS:
        v = _compare_variants(gene, rendu, verite)
        if v is not None:
            return v, [rendu], [verite]
    vr = _sans_appel(variantes(gene, rendu))
    vv = _elargir(gene, _sans_appel(variantes(gene, verite)))
    if not vv: return "pas de verite", vr, vv
    if not vr: return "non appele", vr, vv
    if any(x in vv for x in vr): return "concordant", vr, vv
    if gene in NOMENCLATURE_DIVERGENTE: return "nomenclature differente", vr, vv
    if any(_ambigu(gene, x, y) for x in vr for y in vv): return "ambigu sans phasage", vr, vv
    o = ALLELES_OUTIL.get(gene)
    if o and all(any(a.split("+")[0] not in o for a in x) for x in vv):
        return "allele hors des tables de l outil", vr, vv
    # notre appel porte-t-il strictement plus de variants que la verite ?
    for x in vr:
        a = _variants(gene, x)
        if a is None: continue
        for y in vv:
            b = _variants(gene, y)
            if b is None: continue
            if b < a:
                return "hors perimetre de la verite", vr, vv
    # Un allele manquant ne conclut pas tout de suite : la reference met peut-etre
    # entre parentheses justement celui qu'on ne retrouve pas, et alors elle n'a
    # rien tranche. On epuise donc la boucle avant de decider.
    if _incertitude_explique(gene, rendu, verite):
        return "reference incertaine", vr, vv
    return "DISCORDANT", vr, vv
