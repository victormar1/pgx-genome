# Banc de validation

De quoi refaire, à partir des données publiques, les mesures qui fondent les
chiffres de `doc/VALIDATION.md` § 3.1 et limite 4.2. C'est la revalidation
complète du protocole (§ 5).

| Fichier | Rôle |
|---|---|
| `echantillons.csv` | les 256 prélèvements soumis, leur statut et leurs vérités disponibles ; produit par `outils/construire_banc.py` |
| `preparer_banc.sh` | télécharge les alignements, construit les fichiers de variants, fait tourner le lot, calcule la fiche |
| `fiche_performance.py`, `comparer.py` | confrontent le lot à GeT-RM et au typage HLA Sanger, gène par gène |
| `confronter_puce.sh`, `confronter_puce.py` | confrontent ABCG2 et POR à la puce Omni 2.5 |
| `verites/` | les tables de vérité et les tables d'allèles de l'interpréteur, figées |
| `construire_verite_getrm.py` | reconstruit `verites/verite_getrm.tsv` depuis le classeur du CDC |
| `resultats_bench.csv` | la fiche condensée, une ligne par gène |
| `qualification/` | le génome de qualification d'une installation, voir son `LISEZMOI.md` |

## Rejouer le banc

```bash
export PGX_CYRIUS=/opt/Cyrius PGX_IMG_PYPGX=pgx-genome/pypgx:0.27.0 PYPGX_BUNDLE=/opt/pypgx-bundle
bench/preparer_banc.sh --donnees /disque/banc --fasta /chemin/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

| Étape | Volume, durée |
|---|---|
| `alignements` | 256 alignements de 15 à 20 Go, environ 4 To ; plusieurs jours selon le lien |
| `variants` | quelques Mo par génome ; une heure environ pour l'extraction par chromosome |
| `lot` | environ 40 génomes par heure à six de front sur 24 cœurs |
| `fiche` | quelques minutes |

Chaque étape reprend où la précédente s'est arrêtée. `--genomes` restreint le
banc à une liste ; `--etapes` choisit les étapes.

## Ce qui doit en sortir

| Mesure | Attendu, version 0.4.3 |
|---|---|
| Prélèvements soumis, refusés à l'entrée, menés au bout | 256, 13, 242 |
| Concordance sur les dix gènes à vérité | 1000 / 1009 |
| Allèles HLA à risque | 81 / 81, 0 faux positif |
| ABCG2 et POR contre la puce | 228 / 228 et 204 / 204 |

Un écart à ces valeurs après un changement est un élément à instruire, gène par
gène, avec `fiche/ecarts.tsv`.

**Reproductibilité éprouvée** sur deux génomes du banc, des fichiers de
variants au rendu : 24 couples gène × génome identiques au lot de validation,
typage complémentaire identique, ABCG2 et POR concordants avec la puce. Les
fichiers de variants reconstruits diffèrent de ceux de la validation par trois à
sept indels dans des répétitions (polyT, VNTR), que bcftools appelle selon la
référence employée ; aucun n'atteint un rendu.

## Les vérités, et d'où elles viennent

| Table | Source |
|---|---|
| `verites/verite_getrm.tsv` | CDC GeT-RM, *Consolidated PGx-HLA table*, version du 22 janvier 2025 |
| `verites/hla_sanger.tsv` | typage HLA par séquençage Sanger des génomes du 1000 Genomes (Gourraud *et al.*, PLoS One 2014) |
| `verites/alleles_pharmcat.json`, `verites/definitions_alleles.json` | allèles et définitions de PharmCAT 3.4.0, pour reconnaître deux noms d'un même diplotype |
| puce Omni 2.5 | lue à l'exécution sur le serveur du 1000 Genomes (`release/20130502/supporting/hd_genotype_chip`) |

La table GeT-RM s'écrit dans des idiomes que la comparaison doit lire —
`*2(*35)`, `(*68)+*4`, `*36x2+*10`, `*1B`, `undefined` — et DPYD s'y compare en
variants, pas en noms : `comparer.py` porte ces règles.

**Porteurs MT-RNR1.** Les deux porteurs de m.1095T>C de la limite 4.1 viennent
des appels mitochondriaux de la phase 3 (`release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz`),
confirmés par comptage des lectures sur leur alignement 30×.
