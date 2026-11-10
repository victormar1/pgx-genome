# Banc orthogonal — validation sur des génomes assemblés

Ce benchmark mesure le module contre une vérité tirée d'une **autre technologie** :
les assemblages diploïdes du HPRC Release 2, faits en lecture longue HiFi et ONT
ultra-long. Il répond à la première limite écrite dans la validation : les seuils
du module ont été posés sur les génomes qui les mesurent, et GeT-RM comme le
typage Sanger ont servi à choisir les outils.

Le protocole, préenregistré, est dans [protocole.md](protocole.md). Il fixe les
génomes, les bras, la vérité par gène, les mesures et les critères d'acceptation
**avant** tout résultat. Les empreintes de la méthode sont dans
[methode.sha256](methode.sha256).

## Ce qui sert de vérité

| Gènes | Vérité | Pourquoi |
|---|---|---|
| CYP2C9, CYP2C19, CYP3A5, CYP3A4, TPMT, NUDT15, SLCO1B1, ABCG2 | PharmCAT sur un VCF phasé, tiré des deux haplotypes | les variants ponctuels suffisent à décider |
| POR | génotype de rs1057868 lu sur chaque haplotype | PharmCAT n'a pas de table POR ; le module le rend par son typage complémentaire, et la confrontation l'y lit |
| CYP2D6 | deCYPher sur la séquence du locus | le gène se décide par sa structure, qu'un VCF ne porte pas |
| HLA-A, HLA-B | Immuannot, à deux champs | idem, et le référentiel demande des allèles nommés |
| MT-RNR1 | aucune | les assemblages n'ont pas de chromosome mitochondrial |

Une région qu'un seul des deux haplotypes couvre ne devient jamais une référence :
le gène est écarté de la mesure et compté comme tel.

## Enchaînement

```bash
export PGX_VALIDATION=/chemin/travail PGX_FASTA=/chemin/GRCh38.fa
bash installer.sh --travail "$PGX_VALIDATION"
bash construire_verite.sh --travail "$PGX_VALIDATION" --fasta "$PGX_FASTA" \
     --echantillons meta/hprc_hors_cohorte.txt
```

`construire_verite.sh` enchaîne l'extraction des haplotypes, le VCF de vérité
puis PharmCAT, et le typage par deCYPher et Immuannot. Chaque étape saute ce qui
est déjà fait : une exécution interrompue se reprend par la même commande. Les
alignements et les assemblages sont lus **à distance**, par tranches ; rien n'est
rapatrié en entier.

Les entrées du module se préparent ensuite, un bras à la fois :

```bash
bash recuperer_cram.sh --donnees DOSSIER --echantillons LISTE     # alignements
bash ../preparer_banc.sh --donnees DOSSIER --fasta "$PGX_FASTA" \
     --genomes LISTE --etapes variants                            # bras A, GATK
bash vcf_dragen.sh --donnees DOSSIER --echantillons LISTE         # bras B, DRAGEN
```

Le module tourne ensuite une fois par bras, avec les mêmes seuils et le même
alignement, le fichier de variants seul changeant :

```bash
bash ../../bin/pgx_genome.sh --cram DOSSIER/cram/ECH.cram \
     --vcf DOSSIER/vcf/ECH.vcf.gz --sortie LOT_A/ECH --echantillon ECH
bash ../../bin/pgx_genome.sh --cram DOSSIER/cram/ECH.cram \
     --vcf DOSSIER/vcf_dragen/ECH.vcf.gz --sortie LOT_B/ECH --echantillon ECH
```

Puis la confrontation et les mesures :

```bash
python confronter.py --travail "$PGX_VALIDATION" --lot LOT_A \
       --echantillons LISTE --sortie bras_a.tsv --bras a
python confronter.py --travail "$PGX_VALIDATION" --lot LOT_B \
       --echantillons LISTE --sortie bras_b.tsv --bras b
python mesurer.py bras_a.tsv bras_b.tsv --sensibilite
```

ABCG2 et POR n'ont pas de référence GeT-RM. Ils se confrontent à la puce du
1000 Genomes, que `../confronter_puce.sh` lit à distance :

```bash
python puce_assemblage.py --travail "$PGX_VALIDATION" --puce DOSSIER/puce.tsv \
       --echantillons-puce DOSSIER/puce_echantillons.txt --genomes LISTE
```

## Fichiers

| Fichier | Rôle |
|---|---|
| `installer.sh` | deCYPher et Immuannot, épinglés au commit, dans le dossier de travail |
| `construire_verite.sh` | enchaîne les trois étapes ci-dessous |
| `lancer_extraction.sh`, `extraire_haplotype.py`, `manifeste.py` | statut de couverture, variants et séquences de loci, par haplotype |
| `pharmcat_verite.sh`, `vcf_verite.py` | VCF phasé des deux haplotypes, puis PharmCAT |
| `typer_haplotypes.sh`, `extraire_hla.py`, `verite_cyp2d6.py` | deCYPher et Immuannot, et la règle de choix d'allèle CYP2D6 |
| `recuperer_cram.sh`, `vcf_dragen.sh` | entrées du module : alignement, et variants du bras B |
| `recuperer_cram_plages.sh` | le même alignement en plages simultanées, pour la source qui bride la connexion |
| `lancer_bras.sh`, `lancer_bras_recouvert.sh` | pilotes des deux bras en flux ; le second recouvre le réseau et le calcul |
| `confronter.py` | une ligne par gène et par génome : vérité, référence, rendu, verdicts |
| `puce_assemblage.py` | ABCG2 et POR : la vérité contre la puce Omni 2.5, hors de toute lecture courte |
| `mesurer.py` | concordance et intervalle de Wilson, allèles HLA à risque, écart entre bras |
| `bilan_construction.py` | complétude et contenu de la vérité construite, avant toute mesure |
| `ressources/loci.bed` | fenêtres des loci dont la séquence est extraite |
| `meta/` | les listes de génomes et les deux tables d'index du HPRC |

Le périmètre des positions n'est pas un fichier de ce dossier : il se refait à
chaque exécution à partir des deux BED du module, pour qu'il ne puisse pas s'en
écarter en silence.

## Ce que ce benchmark ne prouve pas

- **Les tandems et duplications CYP2D6 sont la vérité de deCYPher**, elle-même
  non confrontée à une méthode longue-lecture indépendante ailleurs que sur les
  sept génomes du pont que GeT-RM documente.
- **MT-RNR1 n'est pas mesuré du tout**, faute de chromosome mitochondrial dans
  les assemblages.
- La vérité et le module partagent PharmCAT 3.4.0 pour les gènes à variants
  ponctuels. Ce qui est mis à l'épreuve est donc la chaîne qui mène au VCF, pas
  la table d'interprétation.
