# Journal des versions

Les versions suivent `MAJEUR.MINEUR.CORRECTIF`. Un changement qui modifie un
résultat rendu — seuil, périmètre, version d'un outil, table clinique —
incrémente au moins le rang mineur et déclenche la revalidation décrite dans
`doc/VALIDATION.md`.

La version du module figure dans `provenance.json` et au pied du compte rendu.
Elle ne vaut aucune déclaration de conformité.

## 0.4.2

**Compte rendu de patient.** La mention « prototype » tombe dès que l'identité
est fournie ; la préindication se lit dans le fichier d'identité ou se passe
par `--contexte`. Sans cette correction, aucun document de patient ne pouvait
sortir du module.

**Portabilité sur une plateforme.** Le typeur complémentaire tourne dans son
image (`PGX_IMG_PYPGX`), lancée comme les trois autres. Sous Docker, les
conteneurs tournent sous l'identité de l'appelant. Sous Apptainer, un fichier
`.sif` est accepté et son empreinte tracée. Un verrou posé depuis un autre hôte
n'est jamais repris. L'appeleur CYP2D6 n'écrit plus dans `/tmp`. Les empreintes
de `env/conteneurs.txt`, tronquées, sont complètes.

**Qualification d'une installation.** `bench/qualification/` porte un génome
public, son fichier de variants et le rapport de référence ; le mode d'emploi
décrit la mise en service, étape par étape.

**Documentation.** La liste des prélèvements du banc est déduite du lot
(`outils/construire_banc.py`) : 256 soumis, 242 analysés. Le README est
resserré sur l'emploi ; les preuves sont dans `doc/VALIDATION.md`, seule source.
Les affirmations périmées sur l'ancien lot de 105 génomes, le profil de
ressources et les comptes de tests sont retirées ou corrigées.

Aucun résultat rendu ne change : chaque correction du module est vérifiée sur
3 génomes contre le lot de validation, 0 écart. 247 cas de test.

## 0.4.1

**OptiType sur `--fils` fils.** La configuration de son image en fixait seize,
que le typage n'occupe pas : à quatre, il rend le même typage dans le même
temps. Seize fils par génome effondraient un lot à six ; il rend désormais 38,8
génomes par heure à froid.

**Recette Nextflow déclarée sur mesure.** Elle réservait seize cœurs au HLA, qui
en occupe un ; l'exécuteur faisait passer les typages un par un. Chaque étage
déclare désormais sa consommation mesurée, transmise au module comme nombre de
fils : 51,4 génomes par heure à froid au lieu de 16,8.

Aucun résultat rendu ne change : 0 écart sur 720 couples gène × génome, trente
génomes distincts. Débit et consommation par étage : `doc/VALIDATION.md` § 3.10.

## 0.4.0

**Un processus par étage.** `--etage` fait tourner les seuls étages nommés.
L'orchestration reste dans le script : ce qui a été validé d'un seul tenant est
ce qui tourne étage par étage. Chaque étage retrouve ses prérequis dans le
dossier de travail — la recevabilité de ces mêmes entrées, un refus antérieur,
le verdict du contrôle qualité — et refuse s'ils manquent. Deux instances d'un
même étage s'excluent, des étages différents tournent ensemble. La recette
Nextflow passe à un processus par étage et a été exécutée.

Sans `--etage`, rien ne change, et c'est mesuré : rendus identiques au lot. Le
gain de temps est d'une minute par génome ; le bénéfice tient à l'allocation par
étage et à la relance ciblée.

**Contrat d'exécution corrigé.** Le code 2 était décrit comme « aucun étage n'a
démarré », ce qui est faux depuis les refus après contrôle qualité.

**Vérification.** 233 cas de test ; le découpage audité par mutation, 11 défauts
réintroduits, aucun survivant. Le contrôle du dépôt vérifie désormais la
syntaxe des scripts shell, qu'aucun test n'exécute.

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
