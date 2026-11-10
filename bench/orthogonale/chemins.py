# -*- coding: utf-8 -*-
"""Chemins du banc orthogonal.

Aucun chemin de poste n'est ecrit ici. Le dossier de travail, qui porte les
haplotypes extraits et les sorties des outils, vient de l'appel ou de
PGX_VALIDATION ; le dossier des rapports du module vient de l'appel ou de
PGX_LOT ; le depot se deduit de l'emplacement de ce fichier.
"""
import os

ICI = os.path.dirname(os.path.abspath(__file__))
RACINE = os.path.dirname(os.path.dirname(ICI))
VERITES = os.path.join(RACINE, "bench", "verites")
RESSOURCES = os.path.join(RACINE, "ressources")


def _requis(valeur, variable, option):
    if not valeur:
        raise SystemExit(f"dossier absent : donner {option} ou {variable}")
    return valeur.rstrip("/\\")


def travail(arg=None):
    """Dossier de travail : haplotypes/, verite/, decypher/, immuannot/, outils/."""
    return _requis(arg or os.environ.get("PGX_VALIDATION"), "PGX_VALIDATION", "--travail")


def lot(arg=None):
    """Dossier des sorties du module, une par echantillon."""
    return _requis(arg or os.environ.get("PGX_LOT"), "PGX_LOT", "--lot")
