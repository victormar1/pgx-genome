# -*- coding: utf-8 -*-
"""Fiche de performance du module, calculee sur les sorties du module lui-meme.

Aucun chiffre n est repris d une execution anterieure : tout vient des fichiers
que pgx_lot.sh vient de produire, confrontes a la table de reference du CDC et
au typage HLA Sanger. C'est elle qui produit les chiffres du § 3.1 de
doc/VALIDATION.md.

Usage :
  python bench/fiche_performance.py --lot DOSSIER --sortie DOSSIER
"""
import argparse, collections, csv, glob, json, os, statistics, sys

SEP = chr(9)      # tabulation
FL = chr(10)      # fin de ligne

ICI = os.path.dirname(os.path.abspath(__file__))
RACINE = os.path.dirname(ICI)
VERITES = os.path.join(ICI, "verites")
sys.path.insert(0, ICI)
sys.path.insert(1, os.path.join(RACINE, "bin"))   # sensibilite_seuils.py
import comparer
from comparer import compare


def charger_tables(res):
    """Les tables d'alleles de l'interpreteur, qui servent a reconnaitre deux
    noms differents pour le meme diplotype."""
    for f, d in (("alleles_pharmcat.json", comparer.ALLELES_OUTIL),
                 ("definitions_alleles.json", comparer.DEFINITIONS)):
        p = os.path.join(res, f)
        if not os.path.exists(p):
            sys.exit("table introuvable : %s" % p)
        d.update(json.load(open(p, encoding="utf-8")))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lot", required=True)
    p.add_argument("--sortie", required=True)
    p.add_argument("--verite", default=os.path.join(VERITES, "verite_getrm.tsv"))
    p.add_argument("--tables", default=VERITES,
                   help="dossier des tables d'alleles de l'interpreteur")
    p.add_argument("--ressources", default=os.path.join(RACINE, "ressources"))
    p.add_argument("--hla-verite", dest="hla_verite",
                   default=os.path.join(VERITES, "hla_sanger.tsv"),
                   help="typage Sanger HLA de reference (Gourraud 2014)")
    p.add_argument("--gq-defaut", type=int, default=20, dest="gq_defaut")
    p.add_argument("--profondeur-defaut", type=int, default=10, dest="prof_defaut")
    a = p.parse_args()
    os.makedirs(a.sortie, exist_ok=True)
    charger_tables(a.tables)
    # Le perimetre clinique fixe les douze genes que le produit revendique.
    # Le reste est mesure mais rapporte a part : le taux qui engage le produit
    # est celui du perimetre.
    pc = os.path.join(a.ressources, "perimetre_rnpgx.json")
    CLINIQUE = set(json.load(open(pc, encoding="utf-8"))["genes"]) if os.path.exists(pc) else set()

    # ------------------------------------------------------------- la verite
    VER = {}
    for r in csv.DictReader(open(a.verite, encoding="utf-8"), delimiter="\t"):
        for g, v in r.items():
            if g != "sample" and v:
                VER[(r["sample"], g)] = v
    u = collections.defaultdict(set)
    for (s, g), v in VER.items():
        for pair in comparer.variantes(g, v):
            u[g].update(pair)
    comparer.UNIVERS.update(u)

    # ---------------------------------------------------- ce que le lot a rendu
    echantillons, appels, temps, perimetres = {}, [], collections.defaultdict(list), {}
    for prov in sorted(glob.glob(os.path.join(a.lot, "*", "sortie", "provenance.json"))):
        d = json.load(open(prov, encoding="utf-8"))
        s = d["echantillon"]
        echantillons[s] = d
        for nom, v in d.get("etages", {}).items():
            temps[nom].append(v["secondes"])
        pp = os.path.join(os.path.dirname(os.path.dirname(prov)), "travail", "perimetre.json")
        if os.path.exists(pp):
            perimetres[s] = json.load(open(pp, encoding="utf-8"))
        rap = os.path.join(os.path.dirname(prov), s + ".report.json")
        if not os.path.exists(rap):
            # Deux situations, qui ne se comptent pas pareil. Un module en echec
            # a accepte l entree et n a pas abouti : ses couples restent au
            # denominateur. Un genome refuse a l entree n a jamais ete traite, et
            # le module a dit pourquoi : il sort du denominateur, et la fiche le
            # nomme dans sa propre section.
            refuse = d.get("etages", {}).get("recevabilite", {}).get("etat") not in ("OK", None)
            verdict = "refuse a l entree" if refuse else "module en echec"
            for (ech, g), v in VER.items():
                if ech == s:
                    appels.append({"echantillon": s, "gene": g, "rendu": "", "source": "",
                                   "phenotype": "", "verite": v, "statut_qc": "",
                                   "verdict": verdict})
            continue
        j = json.load(open(rap, encoding="utf-8"))
        for g, e in j["genes"].items():
            dips = e.get("recommendationDiplotypes") or e.get("sourceDiplotypes") or []
            lab = "; ".join(sorted({x.get("label") for x in dips if x.get("label")}))
            phen = "; ".join(sorted({q for x in dips for q in (x.get("phenotypes") or [])}))
            appels.append({"echantillon": s, "gene": g, "rendu": lab,
                           "source": e.get("callSource", ""), "phenotype": phen,
                           "verite": VER.get((s, g), ""),
                           "statut_qc": (perimetres.get(s, {}).get("genes", {}).get(g) or {}).get("statut", "")})

    for x in appels:
        if x.get("verdict") in ("module en echec", "refuse a l entree"):
            continue
        x["verdict"] = compare(x["gene"], x["rendu"], x["verite"])[0] if x["verite"] else "sans verite"

    # ------------------------------------------------- HLA : verite Sanger
    # OptiType rend deux champs ; le typage Sanger porte des ambiguites
    # (allele/allele). On compare a deux champs, paire non ordonnee.
    def _2c(a):
        a = a.strip().lstrip("A*").lstrip("B*").lstrip("C*").lstrip("*")
        p = a.split(":")
        return ":".join(p[:2]) if len(p) >= 2 else a
    def _opts(champ):
        return {_2c(x) for x in (champ or "").split("/") if x.strip()}
    hla = {}
    if not os.path.exists(a.hla_verite):
        sys.exit("typage HLA de reference introuvable : %s" % a.hla_verite)
    for r in csv.DictReader(open(a.hla_verite, encoding="utf-8"), delimiter="\t"):
        hla[r["sample"]] = r
    CH = {"HLA-A": ("A", "A.1"), "HLA-B": ("B", "B.1")}
    for x in appels:
        if x["gene"] not in CH or x["echantillon"] not in hla:
            continue
        k1, k2 = CH[x["gene"]]
        vA, vB = _opts(hla[x["echantillon"]].get(k1)), _opts(hla[x["echantillon"]].get(k2))
        rendu = [_2c(y) for y in x["rendu"].replace(" ", "").split("/") if y.strip()]
        x["verite"] = "%s / %s" % ("/".join(sorted(vA)) or "?", "/".join(sorted(vB)) or "?")
        if len(rendu) != 2 or not (vA and vB):
            x["verdict"] = "sans verite" if not (vA and vB) else "non appele"
            continue
        mA, mB = rendu
        if (mA in vA and mB in vB) or (mA in vB and mB in vA):
            x["verdict"] = "concordant"
        elif mA in vA or mA in vB or mB in vA or mB in vB:
            x["verdict"] = "DISCORDANT"; x["nuance"] = "un allele sur deux"
        else:
            x["verdict"] = "DISCORDANT"

    # ------------------------------------------------- chaque ecart, explique
    # Un taux ne se defend pas tout seul. Pour chaque couple qui n'est pas
    # concordant, on redit ce que le controle qualite a mesure sur ce gene chez
    # cet echantillon : les positions qu'il a ecartees, et leur qualite. C'est
    # la difference entre "l'outil s'est trompe" et "la position n'etait pas
    # lisible", et elle ne se lit pas dans le diplotype.
    ecartees, fragiles = {}, {}
    for s_ech in echantillons:
        f = os.path.join(a.lot, s_ech, "travail", "qc_positions.tsv")
        if not os.path.exists(f):
            continue
        for r in csv.DictReader(open(f, encoding="utf-8"), delimiter=SEP):
            genes = [g.strip() for g in (r.get("genes") or "").split(",") if g.strip()]
            if r.get("retenue") != "oui":
                for g in genes:
                    ecartees.setdefault((s_ech, g), []).append(
                        "%s:%s %s GQ=%s prof=%s (%s)"
                        % (r["chrom"], r["pos"], r.get("rsid") or "", r.get("GQ"),
                           r.get("profondeur_alignement"), r.get("statut")))
                continue
            # Position retenue, mais l appeleur n a garde qu une fraction des
            # lectures qui s y trouvent : le placement y est incertain, et le
            # seuil de profondeur seul ne le voit pas. C est ce qui separe le
            # seul desaccord SLCO1B1 du lot des seize autres lectures de la
            # meme position.
            try:
                pa, pv = int(r["profondeur_alignement"]), int(r["profondeur_VCF"])
            except (ValueError, KeyError, TypeError):
                continue
            if pa >= 10 and pv * 2 < pa:
                for g in genes:
                    fragiles.setdefault((s_ech, g), []).append(
                        "%s:%s %s %s lectures retenues sur %s, GQ=%s"
                        % (r["chrom"], r["pos"], r.get("rsid") or "", pv, pa, r.get("GQ")))
    ecarts = []
    for x in appels:
        if x["verdict"] in ("concordant", "sans verite", "pas de verite"):
            continue
        e = dict(x)
        perdues = ecartees.get((x["echantillon"], x["gene"]), [])
        e["positions_ecartees"] = "; ".join(perdues)
        e["nb_positions_ecartees"] = len(perdues)
        faibles = fragiles.get((x["echantillon"], x["gene"]), [])
        e["positions_fragiles"] = "; ".join(faibles)
        e["nb_positions_fragiles"] = len(faibles)
        ecarts.append(e)
    ecarts.sort(key=lambda e: (e["verdict"] != "DISCORDANT", e["gene"], e["echantillon"]))

    # ------------------------------------------------------------- agregation
    par_gene = collections.defaultdict(collections.Counter)
    total = collections.Counter()
    for x in appels:
        if not x["verite"]:
            continue
        par_gene[x["gene"]][x["verdict"]] += 1
        total[x["verdict"]] += 1
    comparables = total["concordant"] + total["DISCORDANT"]
    # deux taux : celui des appels rendus, et celui rapporte a tout ce que la
    # verite couvre, echecs du module compris. Le second est le taux de service.
    couvert = comparables + total["non appele"] + total["module en echec"]
    # "refuse a l entree" n y figure pas : le module n a pas pretendu traiter ces
    # genomes, il les a ecartes en nommant le defaut de l alignement.

    rendu_par_gene = collections.Counter()
    for x in appels:
        if x["rendu"] and "Unknown" not in x["rendu"]:
            rendu_par_gene[x["gene"]] += 1

    complet = [s for s, d in echantillons.items() if d["reussite_complete"]]
    # Un genome refuse a l entree n est pas un genome en echec : le module a
    # constate que l alignement n etait pas dans son domaine et s est arrete
    # avant de rien produire. C est le garde-fou qui a joue.
    refuses = [{"echantillon": s,
                "motif": d.get("etages", {}).get("recevabilite", {}).get("detail", "")}
               for s, d in sorted(echantillons.items())
               if d.get("etages", {}).get("recevabilite", {}).get("etat") not in ("OK", None)]
    noms_refuses = {r["echantillon"] for r in refuses}
    recevables = [s for s in echantillons if s not in noms_refuses]
    en_echec = [s for s in recevables if not echantillons[s]["reussite_complete"]]
    etages = {}
    for nom, v in temps.items():
        v = [x for x in v if x]
        if v:
            # SANS_RESULTAT : l etage a tourne jusqu au bout et n a rien trouve.
            # Un CYP2D6 que l outil ne resout pas n est pas une panne du module,
            # et le compte rendu rend alors le gene non analyse.
            # Les genomes refuses a l entree n ont fait tourner aucun etage : les
            # compter ici ferait porter a chaque etage un echec qui n est pas le
            # sien.
            etats_nom = [d.get("etages", {}).get(nom, {}).get("etat")
                         for s_, d in echantillons.items() if s_ not in noms_refuses]
            etages[nom] = {"mediane_s": int(statistics.median(v)), "min_s": min(v), "max_s": max(v),
                           "echecs": sum(1 for e in etats_nom
                                         if e not in ("OK", "SANS_RESULTAT", None)),
                           "sans_resultat": sum(1 for e in etats_nom if e == "SANS_RESULTAT")}

    # ------------------------------------------- sensibilite aux deux seuils
    # Rejoue le controle qualite a d autres seuils sur les memes fichiers. Sert
    # a montrer ou se situe le defaut : loin d une falaise, ou juste au bord.
    grille = []
    try:
        import sensibilite_seuils as sens
        brut = sens.charge(a.lot)
        for gq in (10, 20, 30, 40):
            for pr in (10, 20):
                c = p_ = z = 0
                for lignes in brut.values():
                    att, ret = sens.rejoue(lignes, gq, pr)
                    for g in att:
                        if ret[g] == att[g]:
                            c += 1
                        elif ret[g] == 0:
                            z += 1
                        else:
                            p_ += 1
                t = c + p_ + z
                grille.append({"gq": gq, "profondeur": pr, "genes": t, "complets": c,
                               "partiels": p_, "absents": z,
                               "part_complets": round(100 * c / t, 1) if t else None,
                               "defaut": (gq, pr) == (a.gq_defaut, a.prof_defaut)})
    except Exception as e:                       # l analyse est un complement
        print("sensibilite non calculee :", e)

    couverture = []
    for s, d in perimetres.items():
        couverture.append(d["positions_retenues"] / d["positions_visees"])
    genes_complets = collections.Counter()
    for s, d in perimetres.items():
        for g, v in d["genes"].items():
            genes_complets[v["statut"]] += 1

    # ----------------------------------------------- vue du perimetre clinique
    clin = collections.Counter()
    clin_gene = collections.defaultdict(collections.Counter)
    for x in appels:
        if not x["verite"] or (CLINIQUE and x["gene"] not in CLINIQUE):
            continue
        clin[x["verdict"]] += 1
        clin_gene[x["gene"]][x["verdict"]] += 1
    clin_comp = clin["concordant"] + clin["DISCORDANT"]
    clin_couvert = clin_comp + clin["non appele"] + clin["module en echec"]
    perimetre_clinique = {
        "genes": sorted(CLINIQUE),
        "comparaisons": clin_comp,
        "concordants": clin["concordant"],
        "discordants": clin["DISCORDANT"],
        "taux": round(100 * clin["concordant"] / clin_comp, 1) if clin_comp else None,
        "taux_de_service": round(100 * clin["concordant"] / clin_couvert, 1) if clin_couvert else None,
        "autres_verdicts": {k: v for k, v in clin.items() if k not in ("concordant", "DISCORDANT")},
        "par_gene": {g: {"comparables": c["concordant"] + c["DISCORDANT"],
                          "concordants": c["concordant"], "discordants": c["DISCORDANT"],
                          "taux": round(100 * c["concordant"] / (c["concordant"] + c["DISCORDANT"]), 1)
                                  if (c["concordant"] + c["DISCORDANT"]) else None,
                          "autres": {k: v for k, v in c.items() if k not in ("concordant", "DISCORDANT")}}
                     for g, c in sorted(clin_gene.items())},
        "genes_sans_reference": sorted(g for g in CLINIQUE if not clin_gene.get(g)),
    }
    # --------------------------------------- alleles HLA a risque (RNPGx)
    # Le referentiel ne demande pas un typage complet mais des alleles nommes.
    # On mesure la sensibilite du module sur chacun : porteurs selon le Sanger,
    # retrouves par OptiType.
    RISQUE = {"HLA-A": ["31:01", "02:01"],
              "HLA-B": ["15:02", "15:11", "58:01", "57:01"]}
    MED_RISQUE = {"15:02": "carbamazépine (SJS/NET)", "15:11": "carbamazépine",
                  "58:01": "allopurinol", "57:01": "abacavir", "31:01": "carbamazépine",
                  "02:01": "—"}
    mod_hla = {}
    for pth in glob.glob(os.path.join(a.lot, "*", "travail", "hla.tsv")):
        e = os.path.basename(os.path.dirname(os.path.dirname(pth)))
        rs = list(csv.DictReader(open(pth, encoding="utf-8"), delimiter="\t"))
        if rs:
            mod_hla[e] = rs[0]
    alleles_risque = []
    for gene, lst in RISQUE.items():
        col = {"HLA-A": ("A", "A.1", "A1", "A2"), "HLA-B": ("B", "B.1", "B1", "B2")}[gene]
        for al in lst:
            porteurs, trouves, manques = 0, 0, []
            for e in hla:
                if e not in mod_hla or e not in echantillons:
                    continue
                vset = _opts(hla[e].get(col[0])) | _opts(hla[e].get(col[1]))
                if al not in vset:
                    continue
                porteurs += 1
                mset = {_2c(mod_hla[e].get(col[2], "")), _2c(mod_hla[e].get(col[3], ""))}
                if al in mset:
                    trouves += 1
                else:
                    manques.append(e)
            if porteurs:
                alleles_risque.append({"allele": gene[-1] + "*" + al, "medicament": MED_RISQUE.get(al, ""),
                                       "porteurs": porteurs, "retrouves": trouves, "manques": manques})
    fiche = {
        "echantillons_traites": len(echantillons),
        "alleles_a_risque_hla": alleles_risque,
        "perimetre_clinique": perimetre_clinique,
        "reussite_complete": len(complet),
        "echantillons_recevables": len(recevables),
        "refuses_a_l_entree": refuses,
        "recevables_en_echec": en_echec,
        "echantillons_avec_verite": len({x["echantillon"] for x in appels if x["verite"]}),
        "comparaisons": comparables,
        "concordants": total["concordant"],
        "discordants": total["DISCORDANT"],
        "taux": round(100 * total["concordant"] / comparables, 1) if comparables else None,
        "couples_couverts_par_la_verite": couvert,
        "taux_de_service": round(100 * total["concordant"] / couvert, 1) if couvert else None,
        "non_appeles": total["non appele"],
        "module_en_echec": total["module en echec"],
        "autres_verdicts": {k: v for k, v in total.items()
                            if k not in ("concordant", "DISCORDANT")},
        "par_gene": {g: {"comparables": c["concordant"] + c["DISCORDANT"],
                         "concordants": c["concordant"],
                         "discordants": c["DISCORDANT"],
                         "taux": round(100 * c["concordant"] / (c["concordant"] + c["DISCORDANT"]), 1)
                                 if (c["concordant"] + c["DISCORDANT"]) else None,
                         "autres": {k: v for k, v in c.items()
                                    if k not in ("concordant", "DISCORDANT")}}
                    for g, c in sorted(par_gene.items())},
        "genes_rendus_par_echantillon": (round(sum(rendu_par_gene.values()) / len(echantillons), 1)
                                         if echantillons else 0),
        "etages": etages,
        "couverture_positions": {
            "mediane": round(100 * statistics.median(couverture), 1) if couverture else None,
            "minimum": round(100 * min(couverture), 1) if couverture else None},
        "statut_des_genes": dict(genes_complets),
        "sensibilite_aux_seuils": grille,
        "ecarts": ecarts,
        "phenotype_non_normal_hors_verite": [
            {"echantillon": e["echantillon"], "gene": e["gene"], "rendu": e["rendu"],
             "phenotype": e["phenotype"], "verite": e["verite"]}
            for e in ecarts
            if e["verdict"] == "hors perimetre de la verite" and e["phenotype"]
            and "Normal" not in e["phenotype"] and "Reference" not in e["phenotype"]],
    }
    json.dump(fiche, open(os.path.join(a.sortie, "fiche_performance.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open(os.path.join(a.sortie, "ecarts.tsv"), "w", newline=FL, encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["echantillon", "gene", "rendu", "phenotype", "verite",
                                           "verdict", "statut_qc", "nb_positions_ecartees",
                                           "positions_ecartees", "nb_positions_fragiles",
                                           "positions_fragiles"], delimiter=SEP,
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(ecarts)
    with open(os.path.join(a.sortie, "appels.tsv"), "w", newline="\n", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["echantillon", "gene", "rendu", "source", "phenotype",
                                           "statut_qc", "verite", "verdict", "nuance"], delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(appels)

    print("echantillons traites %d, recevables %d, reussite complete %d"
          % (fiche["echantillons_traites"], len(recevables), fiche["reussite_complete"]))
    for r in refuses:
        print("  refuse a l entree : %-10s %s" % (r["echantillon"], r["motif"]))
    for e in en_echec:
        print("  recevable en echec : %s" % e)
    print("comparaisons %d, concordants %d, taux %s %%" % (comparables, total["concordant"], fiche["taux"]))
    print("couples couverts par la verite %d, taux de service %s %% (non appeles %d, module en echec %d)"
          % (couvert, fiche["taux_de_service"], total["non appele"], total["module en echec"]))
    pcx = fiche["perimetre_clinique"]
    print("--- perimetre clinique (12 genes) : %d comparaisons, %d concordants, taux %s %%, service %s %%"
          % (pcx["comparaisons"], pcx["concordants"], pcx["taux"], pcx["taux_de_service"]))
    if pcx["genes_sans_reference"]:
        print("    genes du perimetre sans reference GeT-RM : %s" % ", ".join(pcx["genes_sans_reference"]))
    if alleles_risque:
        print("--- alleles HLA a risque (verite Sanger) :")
        for x in alleles_risque:
            print("      %-10s %d/%d retrouves%s" % (x["allele"], x["retrouves"], x["porteurs"],
                  "  MANQUES: " + ",".join(x["manques"]) if x["manques"] else ""))
    print("%-10s %10s %11s %7s" % ("gene", "comparables", "concordants", "taux"))
    for g, v in fiche["par_gene"].items():
        if v["comparables"]:
            print("%-10s %10d %11d %6s" % (g, v["comparables"], v["concordants"], v["taux"]))
    for nom, v in etages.items():
        print("etage %-9s mediane %4d s, %d echec(s)" % (nom, v["mediane_s"], v["echecs"]))
    for g in grille:
        print("seuils GQ %-3d profondeur %-3d : %d genes complets sur %d (%s %%)%s"
              % (g["gq"], g["profondeur"], g["complets"], g["genes"], g["part_complets"],
                 "  <- defaut" if g["defaut"] else ""))
    print("ecarts consignes : %d, dont %d discordants"
          % (len(ecarts), sum(1 for e in ecarts if e["verdict"] == "DISCORDANT")))


if __name__ == "__main__":
    main()
