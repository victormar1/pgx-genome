# Module pharmacogénétique pour génome entier

Ce module lit un génome déjà séquencé et déjà aligné par la plateforme, et rend
un avis pharmacogénétique. Il ne demande ni nouveau prélèvement, ni nouvelle
course de séquençage, ni réalignement, ni second appel de variants.

Règle de conception : **aucun résultat faux ne doit pouvoir sortir sans être
signalé**. Chaque étage teste son code de retour, inscrit son état, et refuse de
se rabattre sur une entrée de secours. Le code de retour du module ne vaut zéro
que si les neuf étages ont abouti.

---

## 1. Ce qu'il faut lui donner

| Entrée | Format | Contrainte |
|---|---|---|
| Alignement | CRAM ou BAM indexé | GRCh38, le même que celui du rendu diagnostique |
| Variants | VCF d'un seul échantillon | GRCh38, avec les champs `GT`, `GQ` et `DP` |
| Référence | FASTA indexé | celle qui a servi à l'alignement, obligatoire pour un CRAM |

Le module **refuse** à l'entrée, sans lancer aucun étage :

- un fichier de variants portant plus d'un échantillon ;
- un identifiant demandé qui ne concorde ni avec le nom porté par le fichier de
  variants, ni avec le champ `SM` de l'alignement ;
- un alignement tronqué, c'est-à-dire dont le marqueur de fin manque ;
- un alignement dont les contigs visés sont absents, ou dont l'un d'eux n'a pas
  la longueur que lui donne la référence fournie ;
- un alignement sur lequel l'accès à une région ciblée échoue, faute d'index
  lisible.

Les deux dernières vérifications ne sont pas de confort. Un alignement sur un
autre assemblage ne lève aucune erreur : les régions visées n'y existent pas, la
tranche revient vide, et tout ce qui suit conclut à une absence de couverture.
Le contrôle ne nomme aucun assemblage, il confronte les contigs et leurs
longueurs à ceux de la référence passée en argument.

**Un refus arrête tout.** Aucun étage ne démarre, le dossier de travail est
vidé de ses intermédiaires, et la seule sortie est `provenance.json`, qui porte
le motif. La purge n'est pas une précaution de confort : avec `--reprise`, la
tranche d'une exécution antérieure serait retrouvée par l'étage 2a, qui se
déclarerait abouti, et le module rendrait un compte rendu pour un génome qu'il
vient de refuser.

La présence de `GQ` n'est pas un confort. Une position dont la qualité est
absente ou illisible est **écartée**, pas acceptée : le module préfère déclarer
qu'il n'a pas mesuré plutôt que de laisser passer une mesure qu'il ne peut pas
juger. Sans `GQ` nulle part, le périmètre se réduit aux positions identiques à la
référence, et le compte rendu le dit.

**Sur les blocs de référence d'un gVCF.** Le module ne les lit pas. Il n'en a pas
besoin : il va chercher la profondeur dans l'alignement, ce qui est une mesure
directe plutôt qu'une déclaration. Un gVCF est donc accepté, ses blocs sont
ignorés, et le résultat est le même qu'avec un VCF de variants.

## 2. Ce qu'il rend

| Sortie | Contenu |
|---|---|
| `sortie/<ECH>.report.json` | l'avis complet, structuré, exploitable par une machine |
| `sortie/CR_<ECH>.pdf` | compte rendu d'une page en français, provisoire |
| `sortie/provenance.json` | entrées, empreintes, versions d'images, seuils, état de chaque étage |
| `travail/perimetre.json` | par gène : positions attendues, retenues, perdues, statut |
| `travail/qc_positions.tsv` | une ligne par position diagnostique, ce qui a été lu |
| `travail/etats.tsv` | l'état et la durée de chaque étage |

`provenance.json` porte le champ `reussite_complete`. Il vaut `true` seulement si
les neuf étages attendus sont présents et qu'aucun n'est en échec. Un dictionnaire
d'étages vide donne `false`. C'est ce champ qu'une plateforme surveille.

Un étage peut porter l'état **sans résultat** : il a tourné jusqu'au bout et n'a
rien conclu. Ce n'est pas une défaillance du module, le gène concerné est rendu
non analysé, et `etages_sans_resultat` les nomme. Confondre les deux ferait
passer pour une panne ce qui est une limite de l'outil sur ce génome, et
inversement.

L'empreinte enregistrée pour l'alignement porte sur ses **512 premiers
mégaoctets**, pas sur la totalité : lire trente gigaoctets pour tracer une
exécution coûterait plus cher que l'analyse elle-même. Celle du fichier de
variants et celle des ressources portent sur la totalité. Les conteneurs sont
tracés par leur empreinte immuable, pas par leur étiquette.

## 3. Les neuf étages

**0. Recevabilité.** Les cinq contrôles de la section 1. Aucun étage ne démarre
si l'un d'eux échoue.

**1. Filtre.** Restreint le fichier de variants aux régions qui portent les 1 226
positions interrogées. Ces régions font 42 kb au total, en 258 intervalles, parce
que chaque position est élargie de cinquante bases de part et d'autre pour que
les insertions et délétions décalées soient récupérées. Le fichier réduit porte
donc quelques enregistrements hors cible ; l'étage suivant ne retient que les
1 226 positions nommées.

**2a. Tranche.** Extrait de l'alignement les régions utiles : les mêmes 42 kb, le
complexe majeur d'histocompatibilité, et les contigs HLA alternatifs de la
référence, qui sont lus dans l'en-tête de l'alignement et non figés dans le
module. Une seule traversée du fichier. Tout ce qui suit travaille sur cette
tranche.

**2b. Contrôle qualité.** Confronte le fichier réduit et la tranche, position par
position. Produit le périmètre : la liste des gènes dont toutes les positions
diagnostiques ont été lues à une profondeur et une qualité suffisantes. Écarte
les génotypes qui ne tiennent pas le seuil en les ramenant à « non appelé »,
jamais à « référence ». Échoue si la profondeur n'a pas pu être mesurée sur
l'alignement, si le fichier porte plus d'un échantillon, ou si un enregistrement
se perd à la réécriture.

Cinq motifs écartent une position, tous mesurés sur le lot de cent cinq génomes
avant d'être retenus : qualité de génotype sous le seuil, profondeur sous le
seuil, **site rejeté par le filtrage de l'appelant**, **génotype calculé sur
moins de lectures que le seuil de profondeur** même si l'alignement en porte
davantage, et **déséquilibre allélique**, quand un allèle d'un hétérozygote porte
moins du quart des lectures. Ce dernier motif ne s'applique pas aux gènes dont le
nombre de copies varie, où une duplication déplace légitimement l'équilibre :
seul CYP2D6 est dans ce cas, et son diplotype vient d'un outil dédié. Le filtrage
de l'appelant peut être ignoré par `--ignorer-filtre`, à n'utiliser que si le
fichier de variants a déjà été refiltré en amont.

Le contrôle marque aussi chaque gène comme appartenant ou non au **périmètre
clinique**, les douze gènes du sous-ensemble RNPGx. Cette distinction commande ce
que le compte rendu rapporte.

**3. CYP2D6.** Ce gène ne se lit pas dans un fichier de variants : il porte des
duplications, des délétions et des hybrides avec son pseudogène voisin. Un outil
dédié lit l'alignement complet et rend un diplotype. L'étage échoue si l'outil
échoue. S'il tourne jusqu'au bout sans conclure, l'étage porte l'état **sans
résultat**, qui n'est pas un échec : le gène est rendu non analysé dans le compte
rendu, et le module reste en réussite complète.

Deux situations donnent ce résultat. L'outil peut ne rendre aucun diplotype, sur
un locus trop remanié : sur ce banc, un génome porte onze copies du locus avec
sept hybrides en série. Il peut aussi rendre plusieurs lectures possibles, qu'il
sépare par un point-virgule. Elles ne sont **pas** transmises à l'interpréteur.
En choisir une serait affirmer ce que l'outil n'a pas dit, et les phénotypes
diffèrent d'une lecture à l'autre. L'interpréteur refuse d'ailleurs une telle
ligne et s'arrête, ce qui ferait perdre tout le génome.

**4. HLA.** Même situation. Les lectures du complexe majeur sont extraites de la
tranche et typées séparément. L'étage échoue en dessous de deux cents paires.

**5. Appels externes.** Met les deux résultats précédents au format attendu par
l'interpréteur. Les tandems de CYP2D6 sont écrits avec des espaces autour du
plus, faute de quoi l'interpréteur rend un résultat indéterminé sans lever
d'erreur. L'étage échoue si un étage marqué `OK` n'a pas déposé sa ligne.

**6. Interprétation.** L'interpréteur reçoit le fichier qualifié et les appels
externes. **Aucun repli** : si le fichier qualifié est absent, l'étage échoue.
Se rabattre sur le fichier filtré remettrait en service les génotypes écartés sur
leur qualité, et la garantie du module disparaîtrait sans une ligne de trace.

**7. Compte rendu.** Le document d'une page. Le périmètre de l'étage 2b commande
ce qu'il a le droit d'écrire.

Les étages 3 et 4 ne dépendent pas de l'étage 2b.

## 4. Les garde-fous, et pourquoi ils existent

**Le périmètre se déduit des données, jamais d'une liste écrite à la main.**
Une version antérieure du banc excluait du compte rendu une liste de gènes tenue
à la main. Un gène absent du fichier de variants passait au travers et ressortait
de type référence, avec la mention « appelé » : quatre faux négatifs sur six pour
un seul gène, sans le moindre avertissement.

**L'option de remplissage par la référence est limitée aux positions absentes.**
L'interpréteur propose de considérer comme identiques à la référence à la fois
les positions absentes du fichier et les génotypes explicitement non appelés. La
seconde moitié annulerait le contrôle qualité.

**Les génotypes sont filtrés sur leur qualité avant interprétation.**
Une position d'un transporteur hépatique a une qualité médiane de 6 sur
soixante-treize génomes, quand les trente-quatre autres positions du même gène
sont au-dessus de 81. Elle se trouve dans une zone où les lectures se placent
mal, et elle passe pourtant les filtres de cohorte.

**Une qualité absente écarte la position.** Un seuil qu'on ne peut pas appliquer
n'est pas un seuil satisfait.

**Aucun fichier intermédiaire n'est choisi par ordre alphabétique.** Les noms
attendus sont construits, jamais devinés par un motif. Un reste d'exécution
antérieure ne peut pas être repris à la place du fichier courant.

**Les intermédiaires sont liés à une empreinte des entrées.** Sans `--reprise`,
le dossier de travail est purgé. Avec `--reprise`, les intermediaires ne sont
réutilisés que si l'alignement, le fichier de variants et les seuils sont
identiques.

**Une recommandation forte non traduite sort en texte source, jamais omise**, et
les médicaments concernés sont nommés en bas du compte rendu. Une version
antérieure mettait de côté sans les afficher celles qu'elle ne savait pas
traduire : sur un métaboliseur lent, quatre recommandations sur neuf
disparaissaient, dont deux contre-indications d'antalgiques.

**Un diplotype long n'est jamais résumé par une affirmation de contenu.** Il est
dit ambigu, ou coupé. Écrire « génotype de référence » sur un diplotype variant
inverserait le sens du compte rendu.

## 5. Le périmètre clinique

Le produit revendique **douze gènes**, le sous-ensemble néphrologie et épilepsie
du panel socle du RNPGx 2026, décrit dans `ressources/perimetre_rnpgx.json` :
CYP2C9, HLA-A et HLA-B pour l'épilepsie ; CYP3A5, CYP3A4, POR, CYP2C19, CYP2D6,
SLCO1B1, ABCG2, TPMT et NUDT15 pour la néphrologie et la transplantation.

Le compte rendu ne rapporte que ces gènes. L'interpréteur en analyse vingt-deux :
les dix autres, mesurés au passage, sont **nommés mais non rapportés**, sous le
titre « Hors périmètre ». Leur validation clinique n'a pas été conduite.

Deux réserves sur le périmètre lui-même. **POR** n'a pas de table dans
l'interpréteur : il est toujours déclaré non analysé. **ABCG2** n'a aucun
échantillon de référence dans la table du CDC, et NUDT15 seulement deux : la
fiche ne prouve donc rien sur l'allopurinol par ce gène, et peu sur l'azathioprine.

**Le HLA, lui, est bien couvert**, à condition d'utiliser la bonne référence. La
table du CDC ne contient un typage HLA que pour un seul génome ; mais le typage
Sanger de Gourraud 2014 en couvre quatre-vingt-huit sur cent cinq. Confronté à
lui, le module concorde à quatre-vingt-seize pour cent à la résolution de deux
champs, et surtout retrouve **les vingt-deux porteurs des allèles à risque** du
référentiel : les douze B\*15:02 et le B\*15:11 qui contre-indiquent la
carbamazépine, les quatre B\*58:01 liés à l'allopurinol, le B\*57:01, et les
quatre A\*31:01. Le référentiel demande justement des allèles nommés, pas un
typage complet : c'est exactement ce que la cohorte permet de vérifier.

## 5 bis. Ce que le module ne fait pas

- **MT-RNR1** n'est pas appelable depuis un fichier de variants avec cet
  interpréteur, qui ne déclare aucune position pour ce gène.
- **Les remaniements de structure** ne sont analysés que pour CYP2D6.
- **Seuls les allèles répertoriés** au catalogue de référence sont recherchés. Un
  variant délétère non répertorié donne un résultat de métaboliseur normal.
- **La profondeur est comptée sans filtre de qualité d'alignement.** Une position
  couverte uniquement par des lectures mal placées est comptée comme lue, à moins
  que la profondeur de l'appelant, elle, ne tombe sous le seuil.
- **Le module ne valide pas biologiquement.** Il prépare, il ne signe pas.

## 6. Prérequis

Sur la machine qui exécute le module :

| Outil | Rôle |
|---|---|
| `bash` 4 ou plus | l'orchestrateur |
| `docker` | les trois conteneurs ci-dessous |
| `samtools` 1.13 ou plus | tranche, profondeur, extraction du complexe majeur |
| `python3` 3.8 ou plus, avec `reportlab` | contrôle qualité, provenance, compte rendu |
| GNU `coreutils` et `findutils` | `stat -c`, `xargs -d` dans le traitement de lot |
| dépôt Cyrius, avec ses dépendances Python | l'étage CYP2D6 |

Conteneurs, versions figées. Un changement de version change les résultats et
doit être remesuré.

```
quay.io/biocontainers/bcftools:1.24--h118bc1c_2
quay.io/biocontainers/optitype:1.3.5--hdfd78af_3
pgkb/pharmcat:3.4.0
```

## 7. Appel

```
export PGX_FASTA=/chemin/GRCh38.fa
export PGX_CYRIUS=/chemin/Cyrius

pgx_genome.sh \
  --cram    /chemin/echantillon.cram \
  --vcf     /chemin/echantillon.vcf.gz \
  --sortie  /chemin/resultats/echantillon \
  --echantillon ECH001
```

| Option | Défaut | Rôle |
|---|---|---|
| `--fasta` | `$PGX_FASTA` | référence d'alignement |
| `--gq` | 20 | seuil de qualité de génotype |
| `--profondeur` | 10 | seuil de profondeur par position |
| `--fils` | 4 | fils pour samtools et l'étage CYP2D6 |
| `--reprise` | absente | réutilise les intermédiaires si les entrées sont identiques |
| `--forcer` | absente | écrit dans un dossier non vide sans purger |

Variables équivalentes : `PGX_GQ`, `PGX_PROFONDEUR`, `PGX_FILS`, et
`PGX_IMG_BCFTOOLS`, `PGX_IMG_PHARMCAT`, `PGX_IMG_OPTITYPE` pour les conteneurs.

## 8. Traitement d'un lot

```
pgx_lot.sh --manifeste liste.tsv --sortie /chemin/resultats --parallele 4
```

Le manifeste porte trois colonnes séparées par des tabulations : identifiant,
chemin de l'alignement, chemin du fichier de variants. Un échantillon déjà
complet est sauté. Le lot produit `tableau_de_bord.tsv` : une ligne par
échantillon, l'état et la durée de chaque étage.

## 9. Choisir les seuils

Les deux seuils sont des choix, pas des constantes. Avant d'en hériter, une
plateforme peut rejouer le contrôle qualité à d'autres valeurs sur ses propres
résultats, sans relire un seul alignement :

```
python3 bin/sensibilite_seuils.py --lot /chemin/resultats --sortie /chemin/analyse
```

Le script relit les `qc_positions.tsv` du lot et recompte, pour chaque couple de
seuils, combien de gènes sortent complets, sous réserve ou sans résultat. Il
nomme aussi les positions dont la qualité est juste au-dessus du seuil : ce sont
elles qui basculent en premier.

Mesuré sur cent cinq génomes de ce banc, la profondeur pèse environ huit fois
plus que la qualité de génotype. Passer la profondeur de dix à vingt coûte treize
points de gènes complets, quand passer la qualité de vingt à trente n'en coûte
qu'un et demi. Le seuil de qualité par défaut n'est donc pas au bord d'une
falaise, celui de profondeur en est proche.

## 10. Coût mesuré

Mesuré sur ce banc : vingt-quatre cœurs, soixante-deux gigaoctets de mémoire,
données sur disque local, cinq échantillons en parallèle.

| Étage | Médiane |
|---|---|
| 0 · recevabilité | 43 s |
| 1 · filtre | 151 s |
| 2a · tranche | 246 s |
| 2b · contrôle qualité | 164 s |
| 3 · CYP2D6 | 697 s |
| 4 · HLA | 265 s |
| 5 · appels externes | 1 s |
| 6 · interprétation | 118 s |
| 7 · compte rendu | 1 s |

Vingt-neuf minutes par génome, dont près de la moitié pour le seul étage
CYP2D6 : l'outil lit trois mille régions de normalisation dispersées dans tout
le génome. Le banc était limité par le disque, jamais par le processeur, qui
restait à quatre-vingt-huit pour cent d'inactivité. Sur un stockage plus rapide,
la marge est donc importante.

Les durées et les volumes par étage sont dans `tableau_de_bord.tsv` du lot, et
résumés dans la fiche de performance.

## 11. Ce que le module a rendu sur cent cinq génomes publics

Cent huit génomes du projet 1000 Genomes. Trois refusés à l'entrée : deux
alignements tronqués, un aligné sur un autre assemblage. Les cent cinq recevables
ont tous abouti.

**Sur les douze gènes du périmètre**, confrontés à la table du CDC :

| | |
|---|---|
| génomes recevables | 105 |
| exécutions complètes | 105 sur 105 |
| couples comparables dans le périmètre | 522 |
| concordance dans le périmètre | 100 % |
| taux de service dans le périmètre | 99,8 % |

La concordance ne porte que sur les couples que la référence sait trancher. Le
taux de service, plus bas, rapporte les concordants à tout ce que la référence
couvre, un couple CYP2D6 non résolu compris.

**Ce que ce chiffre recouvre.** Le HLA est validé à part, contre le typage
Sanger de Gourraud 2014, qui couvre 88 des 105 génomes : concordance de 96 % à
deux champs, et 22 sur 22 des porteurs d'allèles à risque retrouvés. Restent deux
gènes sans preuve sur ce banc : ABCG2 et POR n'ont aucun échantillon de
référence. NUDT15 n'en a que deux. Le taux CPIC/PharmCAT repose donc surtout sur
CYP2C9, CYP2C19, CYP2D6, CYP3A4, CYP3A5, SLCO1B1 et TPMT.

**Les écarts nommés.** Dix couples SLCO1B1 sortent du dénominateur parce que le
module lit un variant à fonction diminuée ou inconnue que le panel de référence,
plus ancien, ne couvre pas. Plusieurs changent le phénotype rendu au clinicien,
vers une fonction diminuée ou un résultat indéterminé sur les statines. Ce ne
sont pas des erreurs de lecture, mais ils ne sont pas confirmés par la référence.
Un couple CYP2D6 reste sans résultat, l'outil dédié n'ayant pas su trancher entre
deux diplotypes. La fiche de performance liste chacun de ces écarts, avec ce que
le contrôle qualité a mesuré.

**Hors périmètre.** L'interpréteur rend dix autres gènes, mesurés mais non
rapportés. Un audit position par position y a trouvé une dizaine d'appels
douteux du module, en particulier de faux CYP4F2\*17 dans une région où trois
gènes se ressemblent. Ils ont motivé les garde-fous de qualité de l'étage 2b, et
ne touchent aucun gène du périmètre sur ce lot.
