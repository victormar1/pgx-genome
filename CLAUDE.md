# Règles locales

Les conventions de ce dépôt sont dans `CONTRIBUTING.md` : format des messages de
commit, ce qu'un message ne contient jamais, concision de la documentation. Les
lire avant tout commit et toute rédaction.

Les crochets les vérifient, encore faut-il qu'ils soient actifs :

```bash
git config core.hooksPath outils/crochets
```

Avant de livrer un changement : `python -m unittest discover -s tests` et
`python outils/verifier_depot.py`. Un test ajouté se justifie par le défaut
qu'il empêche de revenir, et se vérifie en réintroduisant ce défaut.
