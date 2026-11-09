# AGENTS.md

Consignes pour les assistants de code qui travaillent dans ce dépôt. Les
conventions complètes sont dans `CONTRIBUTING.md`, le vocabulaire dans
`doc/GLOSSAIRE.md`.

## Mise en place

```bash
git config core.hooksPath outils/crochets
pip install -r env/requirements.txt
```

## Avant de livrer

```bash
python -m unittest discover -s tests
python outils/verifier_depot.py
```

## Règles

- Un commit, un changement logique, au format de `CONTRIBUTING.md` ; le hook
  `commit-msg` le vérifie.
- Un commit rédigé avec un assistant porte le trailer
  `Co-authored-by: Claude <noreply@anthropic.com>`, sans nom de modèle.
- Un garde-fou ajouté vient avec son test, et le test se vérifie par mutation
  testing : réintroduire le défaut doit le faire échouer.
- Un changement qui modifie un résultat rendu (seuil, périmètre, version d'un
  outil, table clinique) incrémente au moins la version mineure et déclenche la
  revalidation de `doc/VALIDATION.md`, § 5.
- Aucun chemin de poste dans un fichier suivi ; ni chemin ni identifiant
  d'échantillon dans un message de commit.
