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

### 3.4 Indépendance à l'appeleur de variants

20 génomes, mêmes alignements, mêmes positions, appel par l'appeleur de la
plateforme (GATK) et par bcftools.

| Mesure | Résultat |
|---|---|
| Concordance des génotypes | **8 711 / 8 720 = 99,90 %** |
| Localisation des 9 écarts | 8 dans CYP2D6, 1 dans CYP4F2 |

Les deux régions concernées sont les seules régions paralogues du périmètre, et
toutes deux sont déjà traitées à part : CYP2D6 est appelé par un outil dédié sur
l'alignement complet, CYP4F2 est hors périmètre clinique. **Aucun gène rendu
depuis le fichier de variants n'est touché.**

### 3.5 Vérification du code

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

**4.3 Le seuil GQ ≥ 20 suppose une échelle comparable à celle de GATK.** Le GQ
n'a pas la même définition d'un appeleur à l'autre : médiane 99 et plafond 99
pour GATK, médiane 127 et plafond 127 pour bcftools. Sur les positions
comparables, aucune décision ne change — mais la mesure ne porte que sur 374
positions, les sites variants, bcftools n'émettant pas de GQ sur un site
homozygote de référence. **Le seuil doit être revérifié sur l'appeleur de la
plateforme d'accueil**, en particulier DRAGEN ou DeepVariant, non mesurés.

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
