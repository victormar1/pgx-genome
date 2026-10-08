# Validation — état des preuves et protocole de revalidation

Ce document rassemble ce qui a été mesuré, ce qui ne l'a pas été, et ce qui
déclenche une revalidation. Il est destiné à être versé dans un dossier
technique.

## 1. Usage prévu

Le module rend un avis pharmacogénétique sur treize gènes du sous-ensemble
néphrologie et épilepsie du panel socle RNPGx 2026, à partir d'un alignement et
d'un fichier de variants **produits en amont par la plateforme**.

| | |
|---|---|
| Entrées | alignement CRAM ou BAM sur GRCh38, fichier de variants mono-échantillon, référence indexée |
| Sorties | avis structuré (`report.json`), compte rendu en français (`CR_<ECH>.pdf`), trace d'exécution (`provenance.json`) |
| Ne fait pas | ni séquençage, ni alignement, ni appel de variants, ni validation biologique |
| Décision finale | un biologiste médical ; le compte rendu porte sa ligne de validation |

## 2. Conditions d'emploi mesurées

| Condition | Valeur | Origine |
|---|---|---|
| Assemblage | GRCh38, la référence ayant servi à l'alignement | contrôlé à l'étage 0 |
| Profondeur | ≥ 10× par position diagnostique | seuil par défaut |
| Couverture | médiane ≥ 18 × sur les positions du périmètre clinique | 36 exécutions ; **sous le seuil le génome est refusé**, voir 3.4 |
| Qualité de génotype | GQ ≥ 20 | seuil par défaut, **voir la limite 4.3** |
| Équilibre allélique | ≥ 0,25 pour un hétérozygote, hors CYP2D6 | mesuré sur 105 génomes |
| Parallélisme en lot | 2 génomes | OptiType lance `razers3` sur seize threads |
| Dossier de sortie | disque local | 1 361 Mo/s contre 99 sur un disque monté |
| Moteur de conteneurs | `docker`, `apptainer` ou `singularity` | vérifié avant tout étage ; un moteur absent est refusé, non diagnostiqué comme une donnée fautive |
| Identité du patient | fichier, par `--identite` | un argument de ligne de commande est lisible dans la table des processus |

## 3. Preuves disponibles

### 3.1 Exactitude contre vérité indépendante

246 génomes publics du 1000 Genomes (jeu NYGC, GRCh38), vérité GeT-RM (CDC) et
typage HLA Sanger.

| Mesure | Résultat |
|---|---|
| Concordance sur les 12 gènes à vérité | **1000 / 1009 = 99,1 %** |
| Gènes à allèles étoile | 100 % sur les huit |
| HLA-A / HLA-B, quatre chiffres | 190 / 194 et 189 / 194 |
| **Allèles HLA à risque** | **81 / 81 porteurs retrouvés, 0 faux positif** |

Les neuf écarts portent tous sur le second champ d'un allèle HLA, sans
conséquence de prescription.

### 3.2 Robustesse de la chaîne

| Mesure | Résultat |
|---|---|
| Génomes analysables menés au bout | **242 / 243** |
| Refus à l'entrée sur alignement non conforme | 13 / 13, aucun résultat produit |
| Non-régression entre deux versions | **2 673 diplotypes, 0 écart** |
| Garde-fous de l'étage 4b exercés sur donnée réelle | 5 / 5, chacun refusant le seul gène concerné |
| **Répétabilité** | 3 exécutions indépendantes d'un même génome, **29 éléments comparés, 0 écart** |
| Reprise d'une exécution interrompue | résultat identique ; le typage HLA est réutilisé, 201 s → 2 s |

### 3.3 Typage complémentaire

243 génomes. Aucune vérité externe n'existe pour ces quatre gènes ; la mesure de
référence est le génotype relevé position par position par le contrôle qualité,
sur un chemin technique indépendant.

| Gène | Rendus | Concordance avec la mesure |
|---|---|---|
| BCHE | 243 / 243 | 243 / 243 |
| MTHFR | 243 / 243 | 243 / 243 |
| POR | 243 / 243 | 243 / 243 |
| MT-RNR1 | 242 / 243 | **aucun porteur sur ce banc** |

**L'appel positif de MT-RNR1 est exercé par un cas construit.** Un fichier de
variants portant `m.1555A>G` en homoplasmie traverse la chaîne entière :
l'étage le type, l'allèle part en appel externe, l'interpréteur rend « risque
augmenté de surdité sous aminoside », **onze recommandations CPIC fortes** sont
générées, et la consigne « éviter les aminosides » figure au compte rendu avec
les onze molécules. Cela valide la chaîne logicielle, **non la détection** :
l'alignement ne porte pas le variant (`outils/cas_construit.py`).

### 3.4 Domaine de validité en couverture

L'alignement est sous-échantillonné, les variants rappelés depuis l'alignement
réduit — sans quoi seule la mesure de profondeur changerait, pas les génotypes —
puis le module repasse en entier. **9 génomes, 36 exécutions.** Les niveaux sont
des fractions ; la profondeur indiquée est celle que le module mesure lui-même.

| Couverture | Profondeur médiane | Gènes complets | partiels | absents |
|---|---|---|---|---|
| pleine | 34 × | 9,0 | 0,0 | 0,0 |
| 67 % | 24 × | 8,4 | 0,6 | 0,0 |
| 50 % | 17 × | 5,2 | 3,7 | 0,1 |
| 33 % | 11 × | 1,0 | 7,6 | 0,4 |

Un changement de diplotype n'a pas la même portée selon sa nature. Une
abstention et une ambiguïté sont rendues comme telles et ne déclenchent aucune
recommandation ; un appel unique et différent, lui, affirme.

**Dans le périmètre clinique — les treize gènes que le compte rendu restitue :**

| Couverture | Faux | Abstentions | Ambiguïtés |
|---|---|---|---|
| 67 % | **1** | 0 | 0 |
| 50 % | **1** | 1 | 1 |
| 33 % | **9** | 5 | 5 |

**Le seul faux appel à 24 × et à 17 × est le même génome et le même gène** :
HLA-B rendu `*45:04/*57:01` au lieu de `*45:01/*57:01`. C'est un second champ,
et **l'allèle à risque `*57:01` est préservé** — la faiblesse déjà déclarée en
4.4, devenue visible quand la couverture baisse, et sans conséquence de
prescription.

**À 11 ×, neuf faux appels apparaissent**, dont huit de la même forme : un
allèle variant lu comme référence.

| Génome | Gène | Pleine couverture | À 11 × |
|---|---|---|---|
| HG00101, HG00142 | CYP3A5 | `*3/*3` | `*1/*1` |
| HG00131 | CYP3A5 | `*1/*3` | `*1/*1` |
| HG00111 | CYP3A4 | `*1/*22` | `*1/*1` |
| HG00122 | CYP3A4 | `*1/*10` | `*1/*1` |
| HG00122 | SLCO1B1 | `*1/*20` | `*1/*1` |
| HG00131 | SLCO1B1 | `*1/*14` | `*1/*1` |
| HG00118 | CYP2D6 | `*1/*41` | `*34/*119` |

La conséquence est vérifiée jusqu'à la conduite à tenir, et non supposée :

| Couverture | CYP3A5 | Phénotype | Conduite pour le tacrolimus |
|---|---|---|---|
| 34 × | `*3/*3` | métaboliseur lent | dose standard |
| 11 × | `*1/*1` | métaboliseur normal | **augmenter la dose initiale de 1,5 à 2 fois** |

À 11 ×, le compte rendu recommanderait donc de doubler la dose de tacrolimus
chez un patient qui est en réalité métaboliseur lent. Deux génomes sur neuf.

**À la couverture du banc, en revanche, le masquage n'a aucune conséquence
clinique.** Mesure directe : pour chaque gène du périmètre portant une position
masquée, l'interpréteur a été relancé sur le fichier non masqué, les appels
externes fournis des deux côtés.

| Mesure | n |
|---|---|
| Couples gène × génome avec une position masquée | 35 |
| Diplotype inchangé | 25 |
| Devenu une ambiguïté **contenant** l'appel non masqué | 9 |
| Devenu un appel unique différent | **1** |

Les neuf ambiguïtés sont le garde-fou qui fonctionne : la position masquée
aurait discriminé, et le module refuse de choisir — un gène ambigu ne déclenche
aucune recommandation. Le seul appel différent, SLCO1B1 `*1/*1` au lieu de
`*1/*37`, porte **le même phénotype, fonction normale, et les sept mêmes
recommandations de statines**.

**Conséquence : le domaine de validité s'arrête au niveau 50 %**, qui mesure
18 × sur les positions du périmètre clinique et 17 × sur le panel entier.
Jusque-là, aucun allèle à conséquence clinique n'est perdu.

**Le module refuse un génome hors de ce domaine.** Le contrôle qualité mesure la
couverture médiane des positions du périmètre sur l'alignement ; sous le seuil,
le génome est refusé — code de sortie 2, aucun document produit. Le refus plutôt
que l'avertissement : rien, dans le compte rendu, ne distinguerait un appel
fondé sur une position lue d'un appel fondé sur une position supposée.

Le seuil est placé par la mesure, sur les mêmes 36 exécutions :

| Niveau | Couverture médiane du périmètre | Faux appels | Porte |
|---|---|---|---|
| pleine | 35 à 40 × | 0 | accepté |
| 67 % | 24 à 28 × | 1, sans conséquence de prescription | accepté |
| 50 % | 18 à 21 × | 1, le même | accepté |
| 33 % | 11 à 13 × | **9** | **refusé, 9 sur 9** |

Tout seuil de 14 à 18 × sépare identiquement les deux groupes. **18 × est
retenu** : c'est le point le plus bas réellement mesuré, et le seuil
n'extrapole donc pas sous lui. `--couverture-min 0` lève la porte ; la
couverture mesurée et le seuil appliqué figurent dans la trace et sur le compte
rendu.

Hors périmètre, l'interpréteur rend dix autres gènes, mesurés mais non
rapportés : ils accusent un faux appel à 17 × et treize à 11 ×. Ces changements
n'atteignent pas le compte rendu, et c'est une raison de plus de ne pas les
rapporter.

**Ce que ces faux appels enseignent sur le mécanisme.** Les positions fautives
avaient bien été écartées par le contrôle qualité, qui écrit `./.` pour dire
qu'il ne leur fait pas confiance. Mais le préprocesseur de l'interpréteur
**supprime l'enregistrement**, et l'interpréteur suppose alors la référence.
Réinjecter le `./.` après le prétraitement ne change rien : l'interpréteur
l'ignore. Et `--absent-to-ref` ne peut pas être retiré — sans lui, 15 gènes sur
23 sortent indéterminés, parce que l'interpréteur a besoin de bien plus de
positions que les 1 226 mesurées par le module.

**Il n'existe donc aucun canal, par cette interface, pour signaler à
l'interpréteur qu'une position précise n'est pas fiable.** Le module ne peut
agir que dans son propre rendu : il déclare le gène « partiel » et affiche
désormais la proportion de positions perdues. Un gène dont une part notable des
positions manque ne doit pas fonder une décision de posologie.

### 3.5 Profil de ressources

| Mesure | Valeur |
|---|---|
| Mémoire de pointe du module, conteneurs exclus | 668 à 733 Mio |
| Mémoire de pointe du typage HLA, dans son conteneur | **1 235 Mio** |
| Durée d'un génome, somme des médianes par étage | 295 s |

L'interpréteur n'a pas été capté : son étage dure sept secondes, trop peu pour
être échantillonné. Une réservation de 4 Gio par génome couvre les valeurs
mesurées avec une marge.

### 3.6 Indépendance à l'appeleur de variants

**243 génomes**, mêmes alignements, mêmes positions, appel par l'appeleur de
la plateforme (GATK) et par bcftools. Les génotypes sont comparés en bases et
non en indices d'allèles : les deux fichiers n'ont pas la même liste d'ALT, et
comparer les indices fabriquerait des désaccords qui n'en sont pas.

| Mesure | Résultat |
|---|---|
| Concordance des génotypes | **103 604 / 103 683 = 99,924 %** |
| Localisation des 79 écarts | CYP2D6 (59) ; CYP4F2 (18) ; NAT2 (1) ; CYP2B6 (1) |
| Écarts sur un gène rendu depuis le fichier de variants | **0** |

Les régions concernées sont les régions paralogues du périmètre, déjà
traitées à part : CYP2D6 est appelé par un outil dédié sur l'alignement
complet, les autres sont hors périmètre clinique. Les écarts ne sont pas
dispersés : `rs3915951` (CYP2D6) chez 26 génomes ; `rs1058172` (CYP2D6) chez 16 génomes ; `rs4020346` (CYP4F2) chez 15 génomes.

**Le seuil GQ ≥ 20 prend-il la même décision ?** La mesure ne portait
jusqu'ici que sur les positions où bcftools écrit un GQ, soit une quinzaine par
génome. Or le GQ est l'écart entre la deuxième meilleure vraisemblance et la
meilleure, et `mpileup` produit des PL sur **toutes** les positions. La
déduction est vérifiée avant d'être utilisée :

| Côté | Positions | GQ déduit égal au GQ écrit | Même décision au seuil |
|---|---|---|---|
| plateforme, ses propres PL | 106 358 | **99,99 %** | — |
| bcftools, PL de `mpileup` contre GQ de l'appel | 4 269 | 71,98 % | **99,91 %** |

Le GQ déduit n'est pas numériquement identique au GQ écrit par bcftools —
l'appeleur y ajoute son a priori — mais il prend la même décision au seuil. La
mesure passe ainsi de 374 à **103 683 positions**.

| Plateforme | bcftools | n | part |
|---|---|---|---|
| retenue | retenue | 102 872 | 99,218 % |
| écartée | retenue | 719 | 0,693 % |
| écartée | écartée | 49 | 0,047 % |
| retenue | écartée | 43 | 0,041 % |

**762 positions sur 103 683 changent de décision, soit 0,735 %.** Elles se
concentrent sur CYP2B6 (477) ; SLCO1B1 (151) ; CYP2D6 (86).

Une bascule ne change le rendu que si l'un des deux appeleurs lit un **variant**
à cette position : une position écartée est supprimée par le prétraitement et
l'interpréteur suppose alors la référence, donc une référence retenue d'un
côté et écartée de l'autre donne le même résultat. Le critère porte sur les
deux génotypes : dans un sens le module verrait un allèle qu'il ne voyait pas,
dans l'autre il en perdrait un.

**42 bascules sur 762 portent un variant, et toutes sur
CYP2D6 et CYP4F2 — appelés hors du fichier de variants ou hors
périmètre. Aucune n'atteint donc le compte rendu.** Sur SLCO1B1, où se
concentrent 151 bascules, aucune ne porte un variant : `rs71581941` est une
position où l'appeleur de la plateforme n'a aucune confiance et où l'autre lit
une référence, ce que l'interpréteur suppose de toute façon.

Vérification sur le rendu, et non sur la seule position : les 9 génomes
passés en entier des deux côtés donnent **72 / 72 couples
gène × génome identiques** en diplotype, en phénotype et en recommandations,
sur les huit gènes que l'appeleur peut atteindre. Quatre d'entre eux portent une
bascule sur `rs71581941` : le diplotype SLCO1B1 est inchangé dans les quatre.

**Ce que la comparaison ne couvre pas.** 2 190 positions, environ 9 par génome, portent une référence différente d'un fichier à l'autre — des indels ancrés sur une autre base — et sont hors comparaison : c'est là que deux appeleurs divergent le plus. 350 positions portent plusieurs enregistrements, et c'est l'enregistrement appelé qui est retenu, comme à l'étage de contrôle.

### 3.7 Vérification du code

| Mesure | Résultat |
|---|---|
| Cas de test | 159 |
| Audit par mutation | **47 défauts réintroduits, 47 suites en échec, 0 mutation survivante** |
| Contrôles d'intégrité du dépôt | 4 fautes réintroduites, 4 attrapées |
| Reproductibilité du catalogue généré | 3 régénérations, **même empreinte**, égale au fichier livré |

**Deux défauts ont été trouvés par ces contrôles eux-mêmes**, et non par
relecture. Le générateur de catalogue n'était pas déterministe : deux positions
de même coordonnée changeaient de place selon l'ordre d'itération d'un
ensemble. Et la réutilisation du typage HLA lors d'une reprise ne vérifiait pas
l'image qui l'avait produit — une montée de version du typeur suivie d'une
reprise aurait rendu l'ancien résultat sans rien signaler. Les deux sont
corrigés.

## 4. Limites connues, à déclarer

**4.1 MT-RNR1 n'a aucun porteur sur le banc.** Ses trois positions du panel sont
de référence sur les 243 génomes. Le gène porte pourtant toute la valeur
clinique du typage complémentaire, puisqu'il est le seul à rendre un phénotype
et une recommandation. L'appel négatif est vérifié 242 fois ; l'appel positif
l'est par un cas construit (`outils/cas_construit.py`), qui valide la chaîne
logicielle et non la détection. **Un matériau de référence porteur reste
nécessaire.**

**4.2 POR et ABCG2 sont rendus sans vérité.** Aucun matériau de référence
disponible sur ce banc.

**4.3 Le seuil GQ ≥ 20 est mesuré sur deux appeleurs, pas sur tous.** Le GQ n'a
pas la même échelle d'un appeleur à l'autre : plafond 99 pour GATK, 127 pour
bcftools. Sur 103 683 positions et 243 génomes, aucun écart de décision
n'atteint le compte rendu (§ 3.6). **La réserve porte désormais sur les
familles non mesurées** — DRAGEN et DeepVariant — et sur les positions dont la
référence s'écrit autrement, hors comparaison.

Trois positions méritent une vérification sur l'appeleur de la plateforme
d'accueil, parce que leur décision change déjà entre les deux appeleurs
mesurés : `rs186335453` et `rs36060847` (CYP2B6, hors périmètre clinique) et
`rs71581941` (SLCO1B1, 62 % des génomes). Aucune ne porte de variant sur ce
banc, et c'est pourquoi aucune n'y a de conséquence ; sur un autre appeleur,
cela se vérifie plutôt que de se supposer.

**4.4 Le second champ des allèles HLA vacille sur A\*02 et B\*27.** Documenté,
sans conséquence de prescription, le premier champ restant résolu.

**4.5 Les traductions des recommandations CPIC n'ont pas de validation
pharmacologique.**

**4.6 Aucun contrôle de qualité externe.** Une participation à un programme
(EMQN, GenQA) reste à organiser.

## 5. Protocole de revalidation

Toute version des outils et des tables est figée dans `provenance.json` par
empreinte immuable, non par étiquette : un compte rendu ancien est rejouable à
l'identique, tables cliniques comprises.

### 5.1 Déclencheurs

| Changement | Revalidation |
|---|---|
| Version d'une image (interpréteur, HLA, bcftools) | complète |
| Tables cliniques CPIC embarquées | complète |
| Dépôt de l'appeleur CYP2D6, version du typeur complémentaire | complète |
| Une ressource de `ressources/` | complète |
| Un seuil par défaut | complète |
| Code d'un étage | suite de tests, puis non-régression |
| Appeleur de variants en amont, ou plateforme d'accueil | **concordance inter-appeleur et seuil de qualité**, section 3.4 |

### 5.2 Procédure

1. `python -m unittest discover -s tests` et `python outils/verifier_depot.py`
   doivent passer. Le crochet `pre-push` les impose.
2. Rejouer le banc complet et comparer les diplotypes gène par gène à
   l'exécution précédente. Un écart n'est pas un échec : c'est un élément à
   instruire, dont la matérialité clinique se juge au cas par cas.
3. Consigner dans le dossier technique : les deux empreintes d'images, le
   nombre de diplotypes comparés, le nombre et la nature des écarts.
4. Si un écart touche une recommandation forte ou modérée, la décision de
   diffuser la nouvelle version appartient au biologiste responsable.

### 5.3 Ce que le dépôt fournit pour l'exécuter

| Outil | Rôle |
|---|---|
| `tests/` | 139 cas, audités par mutation |
| `outils/verifier_depot.py` | intégrité des figures, cohérence du périmètre, reproductibilité des ressources générées |
| `outils/cas_construit.py` | exercice d'un allèle absent de la cohorte |
| `provenance.json` | empreintes des entrées, des images, des ressources, et seuils employés |
