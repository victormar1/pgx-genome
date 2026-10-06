# -*- coding: utf-8 -*-
"""Provenance d une execution : de quoi la rejouer a l identique.

Enregistre les entrees et leur empreinte, la reference, les ressources qui
definissent le perimetre, les versions des outils telles qu ils les declarent
eux-memes, les seuils, et l etat de chaque etage.

`reussite_complete` est faux des qu un etage manque ou echoue. Un dictionnaire
d etages vide donne faux, jamais vrai.
"""
import argparse, datetime, hashlib, json, os, subprocess, sys


def sha256(chemin, limite=None):
    if not chemin or not os.path.exists(chemin):
        return ""
    h = hashlib.sha256()
    lu = 0
    with open(chemin, "rb") as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            h.update(b)
            lu += len(b)
            if limite and lu >= limite:
                return h.hexdigest() + "+partiel:%d" % limite
    return h.hexdigest()


def version_image(image):
    """Identifiant immuable de l image, pas seulement son etiquette."""
    try:
        out = subprocess.run(["docker", "image", "inspect", image, "--format", "{{index .RepoDigests 0}}"],
                             capture_output=True, text=True, timeout=60)
        d = out.stdout.strip()
        if d:
            return d
        out = subprocess.run(["docker", "image", "inspect", image, "--format", "{{.Id}}"],
                             capture_output=True, text=True, timeout=60)
        return out.stdout.strip() or image
    except Exception:
        return image


def version_pypgx(cmd):
    """Le typeur complementaire ne tourne pas en conteneur : sa version est celle
    de l installation, et elle doit figurer au meme titre que les images."""
    if not cmd:
        return {}
    try:
        out = subprocess.run(cmd.split() + ["--version"], capture_output=True, text=True, timeout=120)
        v = (out.stdout or out.stderr).strip().splitlines()
        return {"commande": cmd, "version": v[-1] if v else ""}
    except Exception:
        return {"commande": cmd, "version": ""}


def version_module(racine):
    """La version du module lui-meme. Le fichier VERSION fait foi : une copie
    deployee n'est pas forcement un depot git. Le commit s'y ajoute quand il
    est disponible, ce qui leve l'ambiguite entre deux etats portant la meme
    version pendant le developpement."""
    d = {"version": "", "commit": ""}
    f = os.path.join(racine, "VERSION")
    if os.path.exists(f):
        with open(f, encoding="utf-8") as fh:
            d["version"] = fh.read().strip()
    d["commit"] = version_git(racine)
    return d


def version_git(chemin):
    if not chemin or not os.path.isdir(chemin):
        return ""
    try:
        out = subprocess.run(["git", "-C", chemin, "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=30)
        return out.stdout.strip()
    except Exception:
        return ""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sortie", required=True)
    p.add_argument("--echantillon", required=True)
    p.add_argument("--cram", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--fasta", default="")
    p.add_argument("--ressources", required=True)
    p.add_argument("--gq", type=int, required=True)
    p.add_argument("--profondeur", type=int, required=True)
    p.add_argument("--images", default="")
    p.add_argument("--cyrius", default="")
    p.add_argument("--pypgx", default="")
    p.add_argument("--equilibre", type=float, default=0.25)
    p.add_argument("--filtre", default="respecte",
                   help="respecte ou ignore : traitement du filtrage de l appelant")
    p.add_argument("--perimetre-clinique", default="", dest="perimetre_clinique")
    a = p.parse_args()

    etats, ordre = {}, []
    pe = os.path.join(a.sortie, "travail", "etats.tsv")
    if os.path.exists(pe):
        for i, l in enumerate(open(pe, encoding="utf-8")):
            if i == 0:
                continue
            c = l.rstrip("\n").split("\t")
            if len(c) >= 4:
                etats[c[0]] = {"etat": c[1], "secondes": int(c[2] or 0), "detail": c[3]}
                ordre.append(c[0])

    ATTENDUS = ["recevabilite", "filtre", "tranche", "qc", "cyp2d6", "hla",
                "appels", "pharmcat", "rendu"]
    manquants = [e for e in ATTENDUS if e not in etats]
    # SANS_RESULTAT : l etage a tourne jusqu au bout et n a rien conclu. Un
    # CYP2D6 que l outil ne resout pas n est pas une defaillance du module, et
    # le compte rendu rend alors le gene non analyse.
    # IGNORE : l etage est optionnel et n etait pas arme sur cette plateforme.
    # Le typage complementaire est dans ce cas quand PGX_PYPGX n est pas defini.
    echecs = [e for e, v in etats.items()
              if v["etat"] not in ("OK", "SANS_RESULTAT", "IGNORE")]
    sans = [e for e, v in etats.items() if v["etat"] == "SANS_RESULTAT"]
    ignores = [e for e, v in etats.items() if v["etat"] == "IGNORE"]
    reussite = bool(etats) and not manquants and not echecs

    ressources = {}
    for f in sorted(os.listdir(a.ressources)):
        chemin = os.path.join(a.ressources, f)
        if os.path.isfile(chemin):
            ressources[f] = sha256(chemin)[:32]

    d = {
        "echantillon": a.echantillon,
        "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "entrees": {
            "alignement": os.path.abspath(a.cram),
            "empreinte_alignement": sha256(a.cram, limite=512 << 20),
            "variants": os.path.abspath(a.vcf),
            "empreinte_variants": sha256(a.vcf),
            "reference": os.path.abspath(a.fasta) if a.fasta else "",
            "empreinte_reference": sha256(a.fasta, limite=64 << 20) if a.fasta else "",
        },
        "module": version_module(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "ressources": ressources,
        "outils": {img: version_image(img) for img in a.images.split(",") if img},
        "cyrius": {"chemin": a.cyrius, "version": version_git(a.cyrius)} if a.cyrius else {},
        "pypgx": version_pypgx(a.pypgx),
        "seuils": {"GQ": a.gq, "profondeur": a.profondeur,
                   "equilibre_allelique_min": a.equilibre, "filtre_appelant": a.filtre},
        "perimetre_clinique": a.perimetre_clinique,
        "etages": etats,
        "ordre_des_etages": ordre,
        "etages_manquants": manquants,
        "etages_en_echec": echecs,
        "etages_sans_resultat": sans,
        "etages_ignores": ignores,
        "reussite_complete": reussite,
    }
    os.makedirs(os.path.join(a.sortie, "sortie"), exist_ok=True)
    json.dump(d, open(os.path.join(a.sortie, "sortie", "provenance.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("provenance ecrite, reussite complete :", reussite,
          ("| manquants : " + ",".join(manquants)) if manquants else "",
          ("| echecs : " + ",".join(echecs)) if echecs else "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
