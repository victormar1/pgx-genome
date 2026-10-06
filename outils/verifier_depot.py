# -*- coding: utf-8 -*-
"""Controles d'integrite du depot. Rend 0 si tout passe, 1 sinon.

Chaque controle correspond a une faute deja publiee : un SVG non conforme au
XML qu'aucune plateforme n'affichait, un nombre de genes qui differait entre le
perimetre et la documentation, des chemins d'une arborescence de travail restes
dans la doc, et une ressource generee qui aurait pu dater d'une version
anterieure du typeur.

Usage : python outils/verifier_depot.py [--regenerer]
"""
import io
import json
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
echecs = []


def dire(ok, titre, detail=""):
    print("  [%s] %-46s %s" % ("ok " if ok else "ECHEC", titre, detail))
    if not ok:
        echecs.append(titre)


def suivis():
    r = subprocess.run(["git", "-C", RACINE, "ls-files"], capture_output=True)
    return [l for l in r.stdout.decode("utf-8", "replace").splitlines() if l]


# --------------------------------------------------------------- 1. les SVG
def controle_svg():
    print("\nFigures")
    for f in sorted(os.listdir(os.path.join(RACINE, "doc", "figures"))):
        if not f.endswith(".svg"):
            continue
        p = os.path.join(RACINE, "doc", "figures", f)
        texte = io.open(p, encoding="utf-8", errors="replace").read()
        # Une figure du depot n'a aucune raison de declarer un DOCTYPE ni une
        # entite : les refuser ferme la porte aux entites externes et a
        # l'expansion en cascade, sans ajouter de dependance a l'analyseur.
        if "<!DOCTYPE" in texte or "<!ENTITY" in texte:
            dire(False, f, "declaration DOCTYPE ou ENTITY")
            continue
        try:
            ET.fromstring(texte)
            dire(True, f, "XML conforme")
        except ET.ParseError as e:
            dire(False, f, str(e))


# ------------------------------------------- 2. coherence du nombre de genes
NOMBRES = {12: "douze", 13: "treize", 14: "quatorze"}


def controle_perimetre():
    print("\nPerimetre")
    p = os.path.join(RACINE, "ressources", "perimetre_rnpgx.json")
    d = json.load(io.open(p, encoding="utf-8"))
    genes = d.get("genes") or {}
    n = len(genes) if isinstance(genes, (dict, list)) else 0
    if not n:
        # le fichier decrit le compte dans sa note quand il n'indexe pas par gene
        m = re.search(r"(Douze|Treize|Quatorze)\s+g[eè]nes", d.get("note", ""), re.I)
        n = {v: k for k, v in NOMBRES.items()}.get((m.group(1).lower() if m else ""), 0)
    dire(bool(n), "nombre de genes lu dans le perimetre", "%d" % n)
    if not n:
        return
    attendu = NOMBRES.get(n, "")
    for rel in ("README.md", "doc/MODE_EMPLOI.md", "doc/FLUX_FICHIERS.md"):
        texte = io.open(os.path.join(RACINE, rel), encoding="utf-8").read()
        mauvais = set()
        for k, mot in NOMBRES.items():
            if k == n:
                continue
            # « douze genes du perimetre » ailleurs que dans un texte de figure
            for m in re.finditer(r"%s g[eè]nes du p[eé]rim[eè]tre" % mot, texte, re.I):
                debut = texte.rfind("\n", 0, m.start()) + 1
                ligne = texte[debut:texte.find("\n", m.start())]
                if "<img" in ligne or "confrontab" in ligne:
                    continue      # alt-text d'une figure, ou sous-ensemble a verite
                mauvais.add(mot)
            if re.search(r"revendique \*\*%s g[eè]nes\*\*" % mot, texte, re.I):
                mauvais.add(mot)
        dire(not mauvais, rel,
             "coherent avec %s" % attendu if not mauvais
             else "annonce %s" % ", ".join(sorted(mauvais)))


# --------------------------------------------------- 3. motifs interdits
INTERDITS = [
    (r"[A-Z]:\\\\Users", "chemin Windows d'un poste"),
    (r"/mnt/[a-z]/", "chemin monte d'un poste"),
    (r"\b0[0-9]_[a-z_]+/", "arborescence de travail"),
    (r"\blot[0-9]+/", "dossier de lot"),
]


def controle_motifs():
    print("\nMotifs interdits dans les fichiers suivis")
    trouves = []
    for rel in suivis():
        if rel.startswith("bench/") or rel.endswith((".png", ".svg")):
            continue
        p = os.path.join(RACINE, rel)
        try:
            texte = io.open(p, encoding="utf-8").read()
        except (UnicodeDecodeError, OSError):
            continue
        for motif, quoi in INTERDITS:
            for m in re.finditer(motif, texte):
                if rel in ("outils/verifier_depot.py", "CONTRIBUTING.md"):
                    continue      # c'est ici qu'on les enumere
                trouves.append("%s : %s (%s)" % (rel, m.group(0), quoi))
    dire(not trouves, "aucun chemin de poste ni de lot",
         "" if not trouves else trouves[0])
    for t in trouves[1:6]:
        print("        %s" % t)


# ------------------------------- 4. les ressources generees se regenerent
def controle_ressources(regenerer):
    print("\nReproductibilite des ressources generees")
    cible = os.path.join(RACINE, "ressources", "pypgx_genes.json")
    try:
        import pypgx  # noqa: F401
    except Exception:
        print("  [passe] pypgx absent : regeneration non verifiee")
        return
    import tempfile
    tmp = tempfile.mkdtemp()
    for f in os.listdir(os.path.join(RACINE, "ressources")):
        if f.startswith("rnpgx_complement"):
            import shutil
            shutil.copy2(os.path.join(RACINE, "ressources", f),
                         os.path.join(tmp, f))
    r = subprocess.run([sys.executable,
                        os.path.join(RACINE, "outils", "construire_pypgx.py"),
                        tmp], capture_output=True)
    if r.returncode != 0:
        dire(False, "construire_pypgx.py",
             r.stderr.decode("utf-8", "replace").strip().splitlines()[-1:][0]
             if r.stderr.strip() else "code %d" % r.returncode)
        return
    neuf = json.load(io.open(os.path.join(tmp, "pypgx_genes.json"),
                             encoding="utf-8"))
    vieux = json.load(io.open(cible, encoding="utf-8"))
    pareil = neuf.get("genes") == vieux.get("genes")
    dire(pareil, "pypgx_genes.json se regenere a l'identique",
         "" if pareil else "le catalogue livre differe de ce que PyPGx rend")
    if not pareil and regenerer:
        import shutil
        shutil.copy2(os.path.join(tmp, "pypgx_genes.json"), cible)
        print("        regenere sur place")


def main():
    regenerer = "--regenerer" in sys.argv
    controle_svg()
    controle_perimetre()
    controle_motifs()
    controle_ressources(regenerer)
    print()
    if echecs:
        print("%d controle(s) en echec : %s" % (len(echecs), ", ".join(echecs)))
        return 1
    print("tous les controles passent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
