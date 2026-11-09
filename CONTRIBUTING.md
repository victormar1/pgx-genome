# Conventions du dépôt

## Contrôles

Une règle écrite n'est pas un garde-fou : ces conventions sont vérifiées par des
scripts, et c'est la vérification qui fait foi.

```bash
git config core.hooksPath outils/crochets   # une fois par copie de travail
```

| Contrôle | Quand | Ce qu'il refuse |
|---|---|---|
| `outils/valider_message.py` | hook `commit-msg` ; CI, sur tout l'historique | sujet hors format, corps au-delà de 72 caractères, motif interdit, casse ou forme d'un trailer |
| `python -m unittest discover -s tests` | hook `pre-push` ; CI | régression sur un garde-fou : contrôle qualité, typage PyPGx, outside calls, rendu, trace, découpage en étages, outils |
| `outils/verifier_depot.py` | hook `pre-push` ; CI | SVG non conforme XML, nombre de gènes ou de tests incohérent avec la documentation, chemin de poste, script shell invalide, ressource générée périmée |

La CI (`.github/workflows/ci.yml`) rejoue ces contrôles à chaque push et à
chaque pull request. `git push --no-verify` contourne le hook local, pas la CI.

Les tests valent ce que vaut leur capacité à échouer. Chaque garde-fou ajouté
est vérifié par mutation testing : le défaut qu'il empêche est réintroduit, et
sa suite doit échouer. Cumul sur les versions : 89 mutants, 89 tués, aucun
survivant (`doc/VALIDATION.md`, § 3.6).

## Messages de commit

Le style suit celui du projet Git (`Documentation/SubmittingPatches`), et non
Conventional Commits : ses types anglais jureraient dans un dépôt rédigé en
français, là où le préfixe de domaine est neutre en langue.

**Sujet** — `domaine: description`

| Règle | |
|---|---|
| Longueur | 50 caractères visés, 72 au maximum |
| Mode | impératif : « corrige », non « correction » ni « corrigé » |
| Contenu | ce qui change, nommé : un fichier, une option, un outil, un gène, un seuil |
| Casse | premier mot après le `:` en minuscule, sauf nom propre |
| Ponctuation | pas de point final |

**Domaines** : `pipeline`, `typage`, `qc`, `rendu`, `appels`, `provenance`,
`ressources`, `outils`, `tests`, `doc`, `figures`, `env`, `bench`, `ci`.

**Corps** — séparé du sujet par une ligne vide, replié à 72 colonnes, de deux à
six lignes. Il dit ce qui n'allait pas, ce qui change, et le chiffre qui le
vérifie, avec son dénominateur. Le diff dit déjà comment le code s'y prend.

**Trailers** — en fin de message, au format `Jeton: valeur`, seule l'initiale
en majuscule.

**Un commit, un changement logique.** Une correction et une mise à jour de
documentation sans lien ne voyagent pas ensemble.

## Assistance par IA

Une partie du code et de la documentation est rédigée avec un assistant,
Claude, d'Anthropic. Ces commits le déclarent par un trailer unique, sans nom de
modèle, que le validateur impose :

```
Co-authored-by: Claude <noreply@anthropic.com>
```

Ils passent les mêmes contrôles que les autres, et l'auteur du commit en répond.
Les consignes données aux assistants sont dans `AGENTS.md`.

## Ce qu'un message ne contient jamais

| Interdit | Pourquoi |
|---|---|
| Chemin de machine ou d'arborescence de travail | n'existe pas chez le lecteur, et expose le poste |
| Identifiant d'échantillon ou donnée de patient | aucun n'a sa place dans un historique public |
| Récit du travail — « v2 », « finalement », « après relecture » | l'historique présente un état, pas son chantier |
| Délibération — « on a préféré », « décidé de » | la justification technique suffit |
| Question ouverte — « à valider par », « reste à trancher » | elle appartient au suivi, pas au commit |
| Métrique de mise au point — « 3 versions ratées » | ne documente pas le code |

## Documentation

Tableaux plutôt que prose. Deux ou trois phrases par point, pas davantage. Un
chiffre vient toujours avec son dénominateur et la cohorte sur laquelle il est
mesuré. Le vocabulaire suit `doc/GLOSSAIRE.md`.
