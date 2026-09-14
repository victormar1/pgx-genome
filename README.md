# pgx-genome

Module de pharmacogénétique pour génome entier. Il lit un génome déjà séquencé et
aligné, et rend un avis pharmacogénétique sur douze gènes. Il ne demande ni
nouveau prélèvement, ni nouvelle course de séquençage, ni réalignement, ni second
appel de variants : il se branche en aval d'un pipeline de séquençage existant.

Règle de conception unique : **aucun résultat faux ne sort sans être signalé**.
Chaque étage teste son code de retour, inscrit son état, et refuse de se rabattre
sur une entrée de secours. Le code de retour du module ne vaut zéro que si les
neuf étages ont abouti.

## Ce qu'il fait

À partir d'un alignement (CRAM/BAM) et d'un fichier de variants (VCF), il produit
pour un génome :

- un avis structuré exploitable par une machine (`report.json`) ;
- un compte rendu d'une page en français (`CR_<ECH>.pdf`) ;
- une trace d'exécution complète (`provenance.json`).

Il couvre le sous-ensemble néphrologie et épilepsie du panel socle RNPGx 2026 :
CYP2C9, CYP2C19, CYP2D6, CYP3A4, CYP3A5, SLCO1B1, TPMT, NUDT15, ABCG2, et HLA-A,
HLA-B. POR est dans le périmètre mais l'interpréteur ne le couvre pas ; il sort
« non analysé ».

Trois outils font le travail biologique, orchestrés par un script bash :

| Étape | Outil |
|---|---|
| CYP2D6 (duplications, hybrides) | Cyrius, sur l'alignement complet |
| HLA de classe I | OptiType, sur les lectures du complexe majeur |
| Les autres gènes et l'interprétation | PharmCAT |

Un étage de contrôle qualité maison, en amont de PharmCAT, distingue « position
lue, identique à la référence » de « position non séquencée ». C'est ce qui évite
que PharmCAT rende un génotype de référence sur une position simplement absente du
fichier.

## Prérequis

Sur la machine qui exécute le module :

| Dépendance | Version | Rôle |
|---|---|---|
| Linux ou WSL2 | — | l'environnement |
| `bash` | ≥ 4 | l'orchestrateur |
| `docker` | ≥ 20 | les trois conteneurs |
| `samtools` | ≥ 1.13 | tranche, profondeur, extraction du complexe majeur |
| `python3` | ≥ 3.8, avec `reportlab` | contrôle qualité, provenance, compte rendu |
| GNU coreutils / findutils | — | `stat -c`, `xargs -d` pour le traitement de lot |
| dépôt Cyrius | commit `7c060db` | l'étage CYP2D6 (apporte ses propres dépendances Python) |

Deux ressources externes, volumineuses, non versionnées dans ce dépôt :

- **La référence GRCh38.** Celle qui a servi à l'alignement, indexée (`.fai`).
  Ce banc utilise `GRCh38_full_analysis_set_plus_decoy_hla.fa` (assemblage NYGC du
  1000 Genomes). Un CRAM en exige une.
- **Le dépôt Cyrius.** Cloné et épinglé au commit ci-dessus.

Les trois images Docker sont épinglées par étiquette **et** par empreinte immuable
dans `env/conteneurs.txt`. Un changement de version change les résultats et doit
être remesuré.

## Installation

```bash
git clone <url-de-ce-depot> pgx-genome
cd pgx-genome
pip install -r env/requirements.txt

# Cyrius, épinglé
git clone https://github.com/Illumina/Cyrius.git /opt/Cyrius
git -C /opt/Cyrius checkout 7c060db
pip install -r /opt/Cyrius/requirements.txt

# les trois conteneurs (voir env/conteneurs.txt pour les empreintes)
docker pull quay.io/biocontainers/bcftools:1.24--h118bc1c_2
docker pull pgkb/pharmcat:3.4.0
docker pull quay.io/biocontainers/optitype:1.3.5--hdfd78af_3
```

Deux variables d'environnement pointent la référence et Cyrius :

```bash
export PGX_FASTA=/chemin/GRCh38_full_analysis_set_plus_decoy_hla.fa
export PGX_CYRIUS=/opt/Cyrius
```

## Utilisation

Un génome :

```bash
bin/pgx_genome.sh \
  --cram    /chemin/echantillon.cram \
  --vcf     /chemin/echantillon.vcf.gz \
  --sortie  /chemin/resultats/echantillon \
  --echantillon ECH001
```

Un lot, avec tableau de bord :

```bash
bin/pgx_lot.sh --manifeste liste.tsv --sortie /chemin/resultats --parallele 4
```

Le manifeste porte trois colonnes séparées par des tabulations : identifiant,
chemin de l'alignement, chemin du fichier de variants. Un exemple est dans
`exemples/manifeste.exemple.tsv`. Un échantillon déjà complet est sauté.

Les options complètes (`--gq`, `--profondeur`, `--fils`, `--reprise`, `--forcer`)
sont décrites dans `doc/MODE_EMPLOI.md`.

## Entrées et sorties

| Entrée | Contrainte |
|---|---|
| Alignement | CRAM ou BAM indexé, GRCh38, le même que celui du rendu diagnostique |
| Variants | VCF d'un seul échantillon, GRCh38, avec `GT`, `GQ`, `DP` |
| Référence | FASTA indexé, celle de l'alignement, obligatoire pour un CRAM |

Le module **refuse à l'entrée**, sans lancer aucun étage : un VCF multi-échantillon,
une identité qui ne concorde pas, un alignement tronqué, un alignement sur un
autre assemblage, un alignement sans index lisible.

| Sortie | Contenu |
|---|---|
| `sortie/<ECH>.report.json` | l'avis complet, structuré |
| `sortie/CR_<ECH>.pdf` | le compte rendu d'une page |
| `sortie/provenance.json` | entrées, empreintes, versions d'images, seuils, état de chaque étage |
| `travail/perimetre.json` | par gène : positions attendues, retenues, perdues, statut |
| `travail/qc_positions.tsv` | une ligne par position diagnostique, ce qui a été lu |

`provenance.json` porte le champ `reussite_complete`, vrai seulement si les neuf
étages sont présents et qu'aucun n'est en échec. C'est le champ qu'une plateforme
surveille.

## Traçabilité

Chaque exécution enregistre dans `provenance.json` : l'empreinte des entrées, le
digest immuable des trois images Docker, le commit de Cyrius, les empreintes des
ressources, les seuils, et l'état de chaque étage. Deux exécutions à provenance
identique produisent le même résultat.

## Ce que le module ne fait pas

- **POR** n'a pas de table dans l'interpréteur : il sort toujours « non analysé ».
- **Les remaniements de structure** ne sont analysés que pour CYP2D6.
- **Seuls les allèles répertoriés** au catalogue de référence sont recherchés ;
  un variant délétère non répertorié donne un métaboliseur normal.
- Le module **prépare, il ne signe pas** : il ne remplace pas la validation
  biologique.

## Structure du dépôt

```
bin/            les sept scripts du pipeline
ressources/     les quatre ressources figées, lues au runtime
doc/            le mode d'emploi complet et le schéma du flux de fichiers
env/            requirements Python et empreintes des conteneurs
exemples/       un manifeste type
```

`doc/MODE_EMPLOI.md` est le manuel de référence, étage par étage, avec la
justification de chaque garde-fou. `doc/FLUX_FICHIERS.md` décrit le cheminement
des fichiers d'un bout à l'autre.

## Licence

À définir. Voir `LICENSE`.
