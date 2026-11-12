# -*- coding: utf-8 -*-
"""Construit ressources/genes_hgnc.json : l'identifiant HGNC de chaque gene.

Un symbole de gene peut changer ; l'identifiant HGNC, non. Les sorties
structurees portent donc les deux, et l'identifiant sert de cle de jointure.

La liste des genes n'est pas saisie a la main : elle est lue dans les
ressources du module, de sorte qu'un gene qui y entre sans identifiant fait
echouer la construction. Chaque identifiant est lu aupres du HGNC, et n'est
retenu que si le symbole demande est le symbole approuve.

Prerequis : acces au reseau.
Usage : python3 outils/construire_hgnc.py [dossier_ressources]
"""
import datetime
import json
import os
import re
import sys
import urllib.request

RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ressources")
SERVICE = "https://rest.genenames.org/fetch/symbol/"
# Etiquettes des ressources qui ne designent pas un gene : le complement RNPGx
# range sous ce nom deux positions d'un haplotype qui couvre plusieurs genes.
HORS_HGNC = {"CYP2C_cluster"}


def genes_des_ressources(res):
    """Tout symbole que le module peut ecrire dans une sortie."""
    out = set()
    for nom in ("pharmcat_positions.vcf", "rnpgx_complement.vcf"):
        with open(os.path.join(res, nom), encoding="utf-8") as fh:
            for ligne in fh:
                if not ligne.startswith("#"):
                    out |= set(re.findall(r"PX=([^;\t\n]+)", ligne))
    for nom in ("pypgx_genes.json", "perimetre_rnpgx.json"):
        with open(os.path.join(res, nom), encoding="utf-8") as fh:
            out |= set(json.load(fh)["genes"])
    return sorted(out - HORS_HGNC)


def identifiant(symbole):
    req = urllib.request.Request(SERVICE + symbole,
                                 headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        docs = json.load(r)["response"]["docs"]
    approuves = [d for d in docs
                 if d.get("symbol") == symbole and d.get("status") == "Approved"]
    if len(approuves) != 1:
        sys.exit("%s : pas de symbole approuve unique au HGNC" % symbole)
    return approuves[0]["hgnc_id"]


def main():
    genes = {s: identifiant(s) for s in genes_des_ressources(RES)}
    with open(os.path.join(RES, "genes_hgnc.json"), "w", encoding="utf-8",
              newline="\n") as fh:
        json.dump({"source": "HGNC, rest.genenames.org, symboles approuves",
                   "releve_le": datetime.date.today().isoformat(),
                   "genes": genes}, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print("%d genes -> genes_hgnc.json" % len(genes))


if __name__ == "__main__":
    main()
