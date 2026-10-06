# Journal des versions

Les versions suivent `MAJEUR.MINEUR.CORRECTIF`. Un changement qui modifie un
résultat rendu — seuil, périmètre, version d'un outil, table clinique —
incrémente au moins le rang mineur et déclenche la revalidation décrite dans
`doc/VALIDATION.md`.

La version du module figure dans `provenance.json` et au pied du compte rendu.
Elle ne vaut aucune déclaration de conformité.

## 0.1.0

Premier état validé sur banc.

**Périmètre.** Treize gènes du sous-ensemble néphrologie et épilepsie du panel
socle RNPGx 2026. CYP2D6 par Cyrius, HLA de classe I par OptiType, MT-RNR1,
BCHE, MTHFR et POR par PyPGx, autres gènes et interprétation par PharmCAT.

**Validation.** 243 génomes publics ; 99,1 % de concordance sur les douze gènes
pour lesquels une vérité existe ; 81 porteurs d'allèle HLA à risque retrouvés
sans faux positif ; typage complémentaire concordant 243/243 avec la mesure
position par position. Limites déclarées dans `doc/VALIDATION.md`, dont
l'absence de porteur de MT-RNR1 sur le banc.

**Vérification.** 159 cas de test, audités par mutation — 47 défauts
réintroduits, aucune mutation survivante. Contrôles d'intégrité du dépôt et
crochets de validation des messages et de la poussée.

**Intégration.** Toutes les dépendances en conteneur, y compris le typeur
complémentaire. Moteur de conteneurs interchangeable, `docker` ou `apptainer`.
Recette d'ordonnanceur dans `exemples/`.
