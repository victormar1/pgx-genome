# -*- coding: utf-8 -*-
"""Etage 1 et 4b : la colonne d'un individu, tiree d'un fichier a plusieurs.

Entree standard : le fichier de variants reduit a la colonne de l'individu
(bcftools view -s NOM). Sortie standard : le meme, debarrasse de ce qu'une
fusion de fichiers individuels y a ajoute sans rien mesurer de lui.

Un fichier a plusieurs echantillons porte une ligne a chaque site ou l'un
d'eux varie. Pour les autres, deux cas :

  - un appel conjoint les a lus : la colonne porte des lectures, et un genotype
    de reference ou un refus d'appeler. C'est une mesure : la ligne passe telle
    quelle, et le controle qualite la juge comme toute ligne d'un fichier
    individuel ;
  - une fusion de fichiers individuels ne sait rien d'eux. Elle laisse une
    colonne vide, ou ramenee a la reference sans une seule valeur ; et, sur la
    ligne d'un individu qui varie, les alleles des autres, sans lectures.

Deux retraits, donc, et rien d'autre :

  - la ligne sans aucune valeur. Gardee, elle ferait ecarter la position, et le
    gene sortirait avec une reserve que le meme individu, lu dans son propre
    fichier, n'a pas. Retiree, la position se juge comme partout ou le fichier
    n'a pas de ligne, sur la profondeur de l'alignement ;
  - l'allele que l'individu ne porte pas et pour lequel sa colonne ne donne
    aucune lecture. Garde, il laisse des champs que ni le controle qualite ni
    le typeur ne savent lire.

Un allele alternatif affirme sans valeur reste, lui : ce n'est pas une absence,
et le controle qualite doit l'ecarter. Un allele non porte dont les lectures
sont comptees, fut-ce zero, reste aussi : un appel conjoint l'a mesure.

Les deux comptes sont ecrits sur la sortie d'erreur.
"""
import io
import re
import sys

VIDES = (".", "")
FORMAT = re.compile(r"^##FORMAT=<ID=([^,>]+),Number=([^,>]+)")


def sans_valeur(fmt, valeur):
    """Vrai si la colonne ne porte ni allele alternatif ni aucune valeur."""
    cles = fmt.split(":")
    if not cles or cles[0] != "GT":
        return False
    champs = valeur.split(":")
    alleles = champs[0].replace("|", "/").split("/")
    if any(a not in (".", "0") for a in alleles):
        return False
    return all(c in VIDES for c in champs[1:])


def retirer_alleles_sans_lectures(c, nombres):
    """Retire de la ligne `c` les alleles non portes et sans lectures.

    `nombres` donne, par champ, le nombre de valeurs que declare l'en-tete :
    R une par allele, A une par allele alternatif, G une par genotype. AD est
    lu par allele quoi que dise l'en-tete : d'anciens appelants le declarent
    de taille libre. Rend vrai si la ligne a change ; elle est alors modifiee
    en place, et son champ INFO, qui decrivait tous les alleles, est vide.
    """
    cles, vals = c[8].split(":"), c[9].split(":")
    if cles[0] != "GT" or "AD" not in cles or cles.index("AD") >= len(vals):
        return False
    ad = vals[cles.index("AD")].split(",")
    alts = c[4].split(",")
    n = len(ad)
    if n != len(alts) + 1:
        return False
    gt = re.split(r"([/|])", vals[0])
    portes = {int(x) for x in gt[::2] if x.isdigit()}
    garder = [0] + [i for i in range(1, n) if i in portes or ad[i] != "."]
    if len(garder) == n:
        return False
    rang = {ancien: nouveau for nouveau, ancien in enumerate(garder)}
    c[4] = ",".join(alts[i - 1] for i in garder[1:]) or "."
    vals[0] = "".join(str(rang[int(x)]) if x.isdigit() else x for x in gt)
    for k in range(1, min(len(cles), len(vals))):
        if vals[k] in VIDES:
            continue
        nombre = "R" if cles[k] == "AD" else nombres.get(cles[k], ".")
        v = vals[k].split(",")
        if nombre == "R" and len(v) == n:
            v = [v[i] for i in garder]
        elif nombre == "A" and len(v) == n - 1:
            v = [v[i - 1] for i in garder[1:]] or ["."]
        elif nombre == "G" and len(v) == n * (n + 1) // 2:
            # diploide : le genotype (a, b), a <= b, est au rang b(b+1)/2 + a
            v = [v[b * (b + 1) // 2 + a] for b in garder for a in garder if a <= b]
        elif nombre == "G" and len(v) == n:
            v = [v[i] for i in garder]
        else:
            continue
        vals[k] = ",".join(v)
    c[9] = ":".join(vals)
    c[7] = "."
    return True


def filtrer(entree, sortie):
    """Recopie le flux ; rend le nombre de lignes retirees et de lignes allegees."""
    retirees, allegees, colonnes, nombres = 0, 0, None, {}
    for l in entree:
        if l.startswith("##"):
            m = FORMAT.match(l)
            if m:
                nombres[m.group(1)] = m.group(2)
            sortie.write(l)
            continue
        c = l.rstrip("\n").split("\t")
        if l.startswith("#CHROM"):
            colonnes = len(c) - 9
            if colonnes != 1:
                raise ValueError("%d colonnes d'echantillon, une seule attendue" % colonnes)
            sortie.write(l)
            continue
        if len(c) < 10:
            sortie.write(l)
        elif sans_valeur(c[8], c[9]):
            retirees += 1
        elif retirer_alleles_sans_lectures(c, nombres):
            allegees += 1
            sortie.write("\t".join(c) + "\n")
        else:
            sortie.write(l)
    if colonnes is None:
        raise ValueError("en-tete des colonnes absent")
    return retirees, allegees


def main():
    # Un champ INFO peut porter des octets qui ne sont pas de l'UTF-8 : ils
    # traversent tels quels, la ligne n'est pas a nous.
    entree = io.TextIOWrapper(sys.stdin.buffer, encoding="utf-8",
                              errors="surrogateescape", newline="")
    sortie = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="surrogateescape", newline="")
    try:
        retirees, allegees = filtrer(entree, sortie)
        sortie.flush()
    except ValueError as e:
        sys.stderr.write("colonne_vcf : %s\n" % e)
        return 1
    sys.stderr.write("colonne_vcf : %d ligne(s) sans valeur pour l'individu, retiree(s) ; "
                     "%d ligne(s) allegee(s) d'alleles sans lectures\n" % (retirees, allegees))
    return 0


if __name__ == "__main__":
    sys.exit(main())
