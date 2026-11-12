# -*- coding: utf-8 -*-
"""Le lot : recevabilite du manifeste et tableau de bord.

    lot.py manifeste FICHIER   les lignes recevables sur la sortie standard,
                               les refus sur la sortie d erreur
    lot.py tableau DOSSIER     ecrit tableau_de_bord.tsv et imprime le resume

Le tableau porte une ligne par echantillon du manifeste, y compris ceux qui
n ont produit aucune sortie : un echantillon qui manque au tableau est un
echantillon dont personne ne remarque l absence.
"""
import collections, csv, json, os, re, sys

from provenance import ATTENDUS

# Le typage complementaire est optionnel : il figure au tableau de bord mais ne
# compte pas dans la reussite, et son absence n est pas un manque.
OPTIONNELS = ["pypgx"]


def verifier(manifeste):
    """(recevables, refus). Un identifiant sert de nom de dossier : il ne doit
    contenir ni separateur ni espace, et ne doit pas apparaitre deux fois."""
    vus, recevables, refus = set(), [], []
    with open(manifeste, encoding="utf-8") as fh:
        for i, l in enumerate(fh, 1):
            l = l.rstrip("\n")
            if not l.strip() or l.startswith("#"):
                continue
            c = l.split("\t")
            if len(c) < 3:
                refus.append("ligne %d : moins de trois colonnes" % i)
                continue
            ech, cram, vcf = c[0].strip(), c[1].strip(), c[2].strip()
            if not re.fullmatch(r"[A-Za-z0-9._-]{1,64}", ech):
                refus.append("ligne %d : identifiant inutilisable comme nom de dossier : %r"
                             % (i, ech))
            elif ech in vus:
                refus.append("ligne %d : identifiant en double : %s" % (i, ech))
            elif not os.path.exists(cram):
                refus.append("ligne %d : alignement introuvable : %s" % (i, cram))
            elif not os.path.exists(vcf):
                refus.append("ligne %d : fichier de variants introuvable : %s" % (i, vcf))
            else:
                vus.add(ech)
                recevables.append((ech, cram, vcf))
    return recevables, refus


def tableau(sortie):
    """Une ligne par echantillon du manifeste verifie, lue dans sa trace."""
    lignes = []
    with open(os.path.join(sortie, "manifeste_verifie.tsv"), encoding="utf-8") as fh:
        echantillons = [l.split("\t")[0].strip() for l in fh]
    for ech in echantillons:
        if not ech:
            continue
        p = os.path.join(sortie, ech, "sortie", "provenance.json")
        ligne = {"echantillon": ech}
        if not os.path.exists(p):
            ligne["reussite_complete"] = False
            ligne["date"] = ""
            for e in ATTENDUS + OPTIONNELS:
                ligne[e] = "AUCUNE SORTIE"
                ligne[e + "_s"] = ""
        else:
            with open(p, encoding="utf-8") as fh:
                d = json.load(fh)
            ligne["reussite_complete"] = d["reussite_complete"]
            ligne["date"] = d["date"]
            for e in ATTENDUS + OPTIONNELS:
                v = d.get("etages", {}).get(e)
                ligne[e] = v["etat"] if v else "ABSENT"
                ligne[e + "_s"] = v["secondes"] if v else ""
        lignes.append(ligne)
    return lignes


def ecrire(sortie, lignes):
    cols = (["echantillon", "reussite_complete", "date"]
            + [x for e in ATTENDUS + OPTIONNELS for x in (e, e + "_s")])
    p = os.path.join(sortie, "tableau_de_bord.tsv")
    with open(p, "w", newline="\n", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(lignes)
    return p


def resume(lignes):
    """Ce que l exploitant lit en fin de lot : qui a failli, et a quel etage."""
    def cites(noms):
        return "%s%s" % (" ".join(noms[:8]), " …" if len(noms) > 8 else "")

    ok = sum(1 for l in lignes if l["reussite_complete"])
    out = ["echantillons du manifeste : %d, reussite complete : %d" % (len(lignes), ok)]
    for e in ATTENDUS:
        ko = [l["echantillon"] for l in lignes if l.get(e) not in ("OK", "SANS_RESULTAT")]
        sr = [l["echantillon"] for l in lignes if l.get(e) == "SANS_RESULTAT"]
        if ko:
            out.append("  etage %-13s en echec chez %d : %s" % (e, len(ko), cites(ko)))
        if sr:
            # l etage a tourne jusqu au bout sans rien conclure : le gene est
            # rendu non analyse, le module n a pas failli
            out.append("  etage %-13s sans resultat chez %d : %s" % (e, len(sr), cites(sr)))
    for e in OPTIONNELS:
        c = collections.Counter(l.get(e) for l in lignes)
        out.append("  etage %-13s (optionnel) : %s"
                   % (e, ", ".join("%s chez %d" % (k, v) for k, v in sorted(c.items()))))
    return out


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ("manifeste", "tableau"):
        sys.stderr.write("usage : lot.py manifeste FICHIER | lot.py tableau DOSSIER\n")
        return 2
    if sys.argv[1] == "manifeste":
        recevables, refus = verifier(sys.argv[2])
        for r in recevables:
            print("\t".join(r))
        for r in refus:
            print(r, file=sys.stderr)
        return 0
    lignes = tableau(sys.argv[2])
    chemin = ecrire(sys.argv[2], lignes)
    for l in resume(lignes):
        print(l)
    print("tableau de bord :", chemin)
    return 0


if __name__ == "__main__":
    sys.exit(main())
