# Validation orthogonale — protocole

Préenregistré le 9 novembre 2026, avant tout résultat du module sur les 211
génomes. Les amendements sont datés et motivés plus bas.

## Objet

Mesurer le module sur des génomes qui n'ont jamais servi à le régler, contre une
vérité tirée d'une autre technologie : les assemblages diploïdes du HPRC
Release 2, en HiFi et ONT ultra-long.

## Ce qui est figé

| Élément | Valeur |
|---|---|
| Module | 0.5.5, commit `103797c`, images épinglées par digest (`env/conteneurs.txt`) |
| Seuils | ceux de la 0.5.5, sans retouche |
| Méthode de vérité | les scripts et ressources de ce dossier, empreintes dans `methode.sha256` |
| Sources de vérité | assemblages HPRC R2 v1.0.1 ; alignements winnowmap sur GRCh38_no_alt v1.1 |

## Génomes

Les 211 génomes du HPRC R2 qui ont un alignement NYGC 30× et ne sont pas dans la
cohorte de réglage (`meta/hprc_hors_cohorte.txt`). Les 12 génomes présents dans
les deux servent à éprouver la méthode et ne comptent pas dans les mesures
(`meta/pont_12.txt`).

## Bras

| Bras | Alignement | Variants |
|---|---|---|
| A | CRAM NYGC 30× | VCF NYGC, GATK, appel conjoint de 2020 |
| B | CRAM NYGC 30× | VCF DRAGEN 3.7.6, réanalyse Illumina |

Le bras A prépare ses entrées par `bench/preparer_banc.sh`, celui de la cohorte,
sans modification : les deux bras ne diffèrent que par la source des variants.

## Vérité par gène

| Gènes | Méthode |
|---|---|
| CYP2C9, CYP2C19, CYP3A5, CYP3A4, TPMT, NUDT15, SLCO1B1, ABCG2 | variants des deux haplotypes, phasés, puis PharmCAT 3.4.0 |
| POR | génotype de rs1057868 lu dans les deux haplotypes |
| CYP2D6 | deCYPher (commit `54a0f9d`) sur la séquence du locus de chaque haplotype, plus la règle ci-dessous |
| HLA-A, HLA-B | Immuannot (commit `31362c3`, données IPD du 2 février 2024), consensus à deux champs, une copie par haplotype exigée |

Règle de spécificité puis de sévérité. Quand plusieurs allèles étoile ont toutes
leurs variations définissantes présentes dans une copie :
1. on écarte ceux dont la définition est contenue dans celle d'un autre
   candidat, car l'allèle le plus spécifique l'emporte, comme dans PharmVar
   (\*49 contient \*10) ;
2. entre les candidats restants, on retient le moins actif selon la table CPIC de
   PharmCAT 3.4.0 (\*56 contre \*10, \*4 contre \*10). À égalité, l'appel de
   deCYPher est gardé.

Un gène qu'un des deux haplotypes ne couvre pas d'un seul alignement est écarté
de la mesure et compté comme tel. Il n'est jamais lu comme référence.

POR n'a pas de table PharmCAT : le module le rend par son typage complémentaire,
et la confrontation le lit dans `travail/pypgx.tsv`.

## Mesures

- Concordance par gène du périmètre, en diplotype, avec l'intervalle de confiance
  à 95 % de Wilson, dans chaque bras.
- Allèles HLA à risque (B\*15:02, B\*58:01, A\*31:01, B\*57:01) : sensibilité et
  spécificité.
- Génomes refusés (code 2) et en échec (code 1), dans chaque bras.
- Différence entre les bras A et B, sur les mêmes génomes.
- Analyse de sensibilité : les mêmes mesures sans les loci à chaîne douteuse.

Tout couple gène × génome écarté l'est avec son motif, et les motifs sont
publiés avec les taux.

## Classement d'un écart

Chaque écart est classé avant toute retouche : erreur du module, erreur de la
vérité (assemblage ou méthode), limite de nomenclature, non conclusif. Les
preuves sont la séquence de l'assemblage, les lectures courtes à la position et
GeT-RM quand il existe. Le classement ne change aucun seuil.

## Critères d'acceptation

Arrêtés avant tout résultat sur les 211 :
- borne basse de l'intervalle à 95 % de la concordance du périmètre au moins
  égale à 98 % dans le bras A ;
- aucun allèle HLA à risque manqué ;
- aucun allèle sans fonction rendu en fonction normale ;
- bras B : mêmes critères. Un écart au bras A se décrit, il ne se corrige pas.

## Amendements

| Date | Amendement | Motif |
|---|---|---|
| 9 novembre 2026, avant tout résultat du module | règle de spécificité ajoutée devant la règle de sévérité | la sévérité seule remplaçait huit appels de deCYPher par un allèle moins spécifique, dont trois \*49 en \*10 |
| 11 et 12 novembre 2026, avant toute mesure | alignements rapatriés par plages, recollés sur la plage servie, gardés sur la taille puis le MD5 de l'index | le transport seul change ; l'alignement retenu est identique à l'octet |
| 11 novembre 2026, avant toute mesure | les deux bras d'un génome sont tentés séparément | la perte d'un bras retirait l'autre de la mesure |

## Épreuve de la méthode sur les 12 génomes du pont

| Confrontation | Résultat |
|---|---|
| Assemblage contre GeT-RM, gènes à variants du périmètre | 31 sur 32 |
| Assemblage contre puce, ABCG2 et POR | 20 sur 20 |
| Assemblage contre GeT-RM, CYP2D6 | 6 sur 7 ; le 7e a les mêmes copies dans une autre phase |
| Assemblage contre Sanger 2014, HLA-A et HLA-B | 16 sur 16 |
| Module 0.5.5 contre assemblage, périmètre | 118 sur 118 |

Le pont CYP2D6 couvre deux délétions \*5, deux duplications et les tandems
\*36+\*10. Le pont HLA couvre quatre porteurs de A\*31:01, un de B\*58:01, un de
B\*57:01, et deux B\*15:01 distingués de B\*15:02.

Le seul désaccord du périmètre est SLCO1B1 chez un génome : l'assemblage et le
module rendent \*37/\*41, GeT-RM \*1/\*1, faute d'avoir testé \*41.

Deux écarts relèvent de la nomenclature : PharmCAT 3.4.0 n'a pas de NUDT15\*2,
que son \*3 englobe ; et une séquence CYP2C19 sans variant est \*1 pour deCYPher
et \*38 pour PharmCAT.

## Vérité des 211 génomes

| Gène | Génomes conclusifs | Contenu |
|---|---|---|
| CYP2C9, CYP2C19, CYP3A5, CYP3A4, TPMT, NUDT15, SLCO1B1, ABCG2 | 211 chacun | un appel PharmCAT unique, aucune exclusion |
| CYP2D6 | 211 | 18 délétions \*5, 13 duplications, 48 tandems ou hybrides sur 422 haplotypes |
| HLA-A | 211 | 13 porteurs de A\*31:01 |
| HLA-B | 211 | 6 porteurs de B\*15:02, 10 de B\*58:01, 10 de B\*57:01 |

Exclusions : MT-RNR1 (pas de chromosome mitochondrial dans les assemblages) et,
hors périmètre, G6PD chez les hommes (104).

ABCG2 et POR n'ont pas de référence GeT-RM. Confrontés à la puce Omni 2.5, ils
donnent 147 accords sur 147 et 142 sur 142 ; les autres génomes ne sont pas sur
la puce.

Loci à chaîne douteuse, sur 422 haplotypes :
- 8 CYP2D6, tous des tandems (7 \*36+\*10, 1 \*68+\*4) ;
- 18 HLA-A : 7 chaînes de 84,5 kb pour une fenêtre de 26 kb, 11 non colinéaires.

Ils restent dans la mesure principale, et l'analyse de sensibilité les retire.
