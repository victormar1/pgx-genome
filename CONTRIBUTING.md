# Conventions du dépôt

## Messages de commit

Le style suit celui du projet Git (`Documentation/SubmittingPatches`), et non
Conventional Commits : ses types anglais jureraient dans un dépôt rédigé en
français, là où le préfixe de domaine est neutre en langue.

**Sujet** — `domaine: description`

| Règle | |
|---|---|
| Longueur | 50 caractères visés, 72 au maximum |
| Mode | impératif : « corrige », non « correction » ni « corrigé » |
| Casse | premier mot après le `:` en minuscule, sauf nom propre |
| Ponctuation | pas de point final |

**Domaines** : `pipeline`, `typage`, `qc`, `rendu`, `appels`, `provenance`,
`ressources`, `outils`, `doc`, `figures`, `env`, `bench`.

**Corps** — séparé du sujet par une ligne vide, replié à 72 colonnes. Il énonce
le problème que le changement résout et pourquoi cette solution, pas comment le
code s'y prend : le diff le dit déjà. Les mesures qui justifient le changement y
figurent avec leur dénominateur.

**Lignes de fin** (*trailers*) — en fin de message, après le corps, au format
`Jeton: valeur` avec la seule initiale en majuscule. `Co-authored-by:` est
accepté, sauf sur le commit initial.

**Un commit, un changement logique.** Une correction et une mise à jour de
documentation sans lien ne voyagent pas ensemble.

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
mesuré.
