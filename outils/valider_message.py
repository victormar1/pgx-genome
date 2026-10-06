# -*- coding: utf-8 -*-
"""Valide un message de commit contre les conventions de CONTRIBUTING.md.

Appele par le crochet commit-msg, il refuse le commit plutot que de laisser
passer une faute qu'il faudra reecrire ensuite : les seize premiers commits de
ce depot ont ete reecrits pour une ligne d'attribution, et six autres pour des
corps replies a soixante-dix-huit colonnes.

Les longueurs sont comptees en caracteres et non en octets : un accent en UTF-8
pese deux octets, et un controle qui compte des octets declare trop longs des
sujets qui ne le sont pas.
"""
import io
import re
import sys

LARGEUR = 72
SUJET_VISE = 50

DOMAINES = ("pipeline", "typage", "qc", "rendu", "appels", "provenance",
            "ressources", "outils", "doc", "figures", "env", "bench", "tests")

INTERDITS = [
    (re.compile(r"[A-Z]:\\"), "chemin Windows d'un poste"),
    (re.compile(r"/mnt/[a-z]/"), "chemin monte d'un poste"),
    (re.compile(r"\b0[0-9]_[a-z_]+/"), "arborescence de travail"),
    (re.compile(r"\b(HG|NA)[0-9]{5}\b"), "identifiant d'echantillon"),
    (re.compile(r"\b(finalement|apr[eè]s relecture|on a pr[eé]f[eé]r[eé]|"
                r"d[eé]cid[eé] de|v[0-9]+ du script|version rat[eé]e)\b",
                re.I), "recit du travail ou deliberation"),
    (re.compile(r"\b([aà] valider par|reste [aà] trancher|TODO|FIXME)\b", re.I),
     "question ouverte"),
]


def verifie(texte):
    # Les commentaires de git et les conflits de fusion ne sont pas du message.
    lignes = [l for l in texte.splitlines() if not l.startswith("#")]
    while lignes and not lignes[0].strip():
        lignes.pop(0)
    if not lignes:
        return ["message vide"]

    fautes = []
    sujet = lignes[0].rstrip()

    if len(sujet) > LARGEUR:
        fautes.append("sujet de %d caracteres, %d au maximum"
                      % (len(sujet), LARGEUR))
    if sujet.endswith("."):
        fautes.append("sujet termine par un point")

    m = re.match(r"^([a-z][a-z0-9-]*): (.+)$", sujet)
    if not m:
        fautes.append("sujet hors format « domaine: description » "
                      "(domaine en minuscules, deux-points, espace)")
    else:
        domaine, description = m.group(1), m.group(2)
        if domaine not in DOMAINES:
            fautes.append("domaine « %s » inconnu ; attendus : %s"
                          % (domaine, ", ".join(DOMAINES)))
        if description[:1].isupper():
            fautes.append("majuscule apres le deux-points : « %s »"
                          % description[:24])

    if len(lignes) > 1 and lignes[1].strip():
        fautes.append("pas de ligne vide entre le sujet et le corps")

    for i, l in enumerate(lignes[2:], start=3):
        # Une ligne de fin, une URL ou un chemin ne se replient pas.
        if re.match(r"^[A-Z][A-Za-z-]+: ", l) or "://" in l:
            continue
        if len(l.rstrip()) > LARGEUR:
            fautes.append("ligne %d du corps : %d caracteres, %d au maximum"
                          % (i, len(l.rstrip()), LARGEUR))

    # Casse des lignes de fin : le projet Git ne met en majuscule que
    # l'initiale du jeton, et « Co-Authored-By: » est la forme fautive la plus
    # repandue parce que c'est celle qu'insere l'outillage.
    for l in lignes:
        m = re.match(r"^([A-Za-z][A-Za-z-]*)(: .+)$", l.rstrip())
        if not m:
            continue
        jeton = m.group(1)
        attendu = jeton[:1].upper() + jeton[1:].lower()
        if jeton != attendu and "-" in jeton:
            fautes.append("ligne de fin « %s: » — seule l'initiale se met en "
                          "majuscule, donc « %s: »" % (jeton, attendu))

    entier = "\n".join(lignes)
    for motif, quoi in INTERDITS:
        trouve = motif.search(entier)
        if trouve:
            fautes.append("%s : « %s »" % (quoi, trouve.group(0)))

    return fautes


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("usage : valider_message.py <fichier>\n")
        return 2
    texte = io.open(sys.argv[1], encoding="utf-8", errors="replace").read()
    fautes = verifie(texte)
    if not fautes:
        sujet = texte.splitlines()[0] if texte.splitlines() else ""
        if len(sujet) > SUJET_VISE:
            sys.stderr.write(
                "message accepte ; sujet de %d caracteres, %d vises\n"
                % (len(sujet), SUJET_VISE))
        return 0
    sys.stderr.write("\nMessage de commit refuse — voir CONTRIBUTING.md\n")
    for f in fautes:
        sys.stderr.write("  - %s\n" % f)
    sys.stderr.write("\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
