# Journal des versions

Les versions suivent `MAJEUR.MINEUR.CORRECTIF`. Un changement qui modifie un
résultat rendu — seuil, périmètre, version d'un outil, table clinique —
incrémente au moins le rang mineur et déclenche la revalidation décrite dans
`doc/VALIDATION.md`.

La version du module figure dans `provenance.json` et au pied du compte rendu.
Elle ne vaut aucune déclaration de conformité.

## 0.3.0

**Spécificité à un mélange d'échantillons, mesurée.** Un contaminant injecté de
2 à 30 % dans les lectures d'un hôte : sans conséquence jusqu'à 10 %, le
document change à partir de 20 %, et à 30 % un faux appel unique change un
phénotype. Un indicateur d'ensemble est ajouté — part des hétérozygotes dont
l'allèle mineur est trop peu soutenu — et l'échantillon est refusé au-delà de
35 %, seuil posé sur 229 échantillons purs dont le plus haut est à 31,9 %. La
porte arrête un mélange de 30 % et **ne détecte pas celui de 20 %** : un
contrôle de contamination en amont reste nécessaire, et la limite est déclarée.

**Non-régression établie.** 5 génomes repassés sur les mêmes entrées : 60
couples gène × échantillon, **aucun écart** de diplotype, de phénotype ni de
recommandation. Un alignement réduit à 1 % de ses lectures est refusé.
`outils/comparer_rapports.py` rend cette comparaison exécutable, avec un code de
sortie qui en fait une porte.

**Relecture des traductions outillée.** `outils/fiche_traductions.py` produit
la fiche de relecture des 66 traductions, classée par le nombre de comptes
rendus qui les emploient, avec une colonne de visa. Sur le banc, 51 sont rendues
et **aucun texte rendu n'est laissé sans traduction** ; les 15 entrées
inemployées sont listées à part, et leurs causes attendues nommées.

**Ce qui tient lieu de vérité pour ABCG2 et POR.** Faute de matériau de
référence, les deux positions sans vérité sont caractérisées : lues sur 243
génomes sur 243, à 33 × et 53 × de profondeur médiane, et pour ABCG2 — seul des
deux appelé depuis le fichier de variants — **243 génotypes sur 243 concordants
entre deux appeleurs indépendants**, avec des hétérozygotes équilibrés. Le
dossier dit aussi ce que ces mesures ne prouvent pas.

**Vérification.** 219 cas de test ; l'indicateur de mélange audité par
mutation, 8 défauts réintroduits, aucun survivant.

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
Vérifié jusqu'au rendu sur 29 génomes repassés en entier des deux côtés : 232
couples gène × génome, 6 écarts de diplotype, tous une ambiguïté d'un côté
tranchée en l'un de ses propres membres de l'autre, et **aucun écart de
phénotype ni de recommandation**.

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
