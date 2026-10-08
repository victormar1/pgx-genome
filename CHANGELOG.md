# Journal des versions

Les versions suivent `MAJEUR.MINEUR.CORRECTIF`. Un changement qui modifie un
résultat rendu — seuil, périmètre, version d'un outil, table clinique —
incrémente au moins le rang mineur et déclenche la revalidation décrite dans
`doc/VALIDATION.md`.

La version du module figure dans `provenance.json` et au pied du compte rendu.
Elle ne vaut aucune déclaration de conformité.

## 0.2.0

**Domaine de validité opposable.** Un génome dont la couverture médiane sur les
positions du périmètre clinique est inférieure à 18 × est refusé : arrêt au
sortir du contrôle qualité, code de sortie 2, aucun compte rendu. Le seuil est
placé par la titration — 36 exécutions, 9 génomes : toute exécution sans faux
appel mesure 18 × ou plus, toute exécution qui en porte neuf mesure 13 × ou
moins. `--couverture-min 0` lève la porte. La couverture mesurée et le seuil
appliqué figurent dans la trace et au pied du compte rendu.

**Indépendance à l'appeleur de variants, mesurée sur la base entière.** 243
génomes, 103 683 positions comparées en bases et non en indices d'allèles :
99,924 % de concordance, et aucun écart sur un gène rendu depuis le fichier de
variants. La transférabilité du seuil GQ, jusqu'ici mesurée sur 374 positions,
l'est désormais sur toutes : le GQ se déduit des PL, et la déduction est
vérifiée contre le GQ écrit avant d'être utilisée. 0,735 % des positions
changent de décision selon l'appeleur ; aucune n'atteint le compte rendu.

**Vérification.** 178 cas de test ; la porte auditée par mutation, 11 défauts
réintroduits, aucun survivant.

## 0.1.0

Premier état validé sur banc.

**Périmètre.** Treize gènes du sous-ensemble néphrologie et épilepsie du panel
socle RNPGx 2026. CYP2D6 par Cyrius, HLA de classe I par OptiType, MT-RNR1,
BCHE, MTHFR et POR par PyPGx, autres gènes et interprétation par PharmCAT.

**Validation.** 243 génomes publics ; 99,1 % de concordance sur les douze gènes
pour lesquels une vérité existe ; 81 porteurs d'allèle HLA à risque retrouvés
sans faux positif ; typage complémentaire concordant 243/243 avec la mesure
position par position. Répétabilité : 3 exécutions d'un même génome, 29 éléments
comparés, 0 écart. Indépendance à l'appeleur de variants : 8 711 génotypes
concordants sur 8 720, les 9 écarts en région paralogue déjà traitée à part.
Limites déclarées dans `doc/VALIDATION.md`, dont l'absence de porteur de
MT-RNR1 sur le banc — son appel positif est exercé par un cas construit.

**Vérification.** 159 cas de test, audités par mutation — 47 défauts
réintroduits, aucune mutation survivante. Contrôles d'intégrité du dépôt et
crochets de validation des messages et de la poussée. Le catalogue généré se
régénère à l'identique.

**Intégration.** Toutes les dépendances en conteneur, y compris le typeur
complémentaire. Moteur de conteneurs interchangeable, `docker` ou `apptainer`,
vérifié avant tout étage. Contrat d'exécution documenté : codes de sortie,
verrou de concurrence, champ à surveiller. Recette d'ordonnanceur dans
`exemples/`.

**Traçabilité et données personnelles.** La version du module figure dans la
trace d'exécution et au pied du compte rendu. L'identité du patient se lit dans
un fichier et non en argument de ligne de commande. La signature de reprise
couvre la référence, les ressources et l'image du typeur HLA.
