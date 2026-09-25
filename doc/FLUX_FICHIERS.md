# Flux de fichiers — module `pgx_genome.sh`

**Source de vérité :** `05_produit/bin/pgx_genome.sh`
**Échantillon inspecté :** `lot2/NA06991/` (exécution complète, 9 étages OK)

## Conventions et repères

- **RES** = `05_produit/ressources/` (ressources figées, montées en lecture seule dans les conteneurs).
- **T** = `<SORTIE>/travail/` (fichiers intermédiaires).
- **sortie** = `<SORTIE>/sortie/` (livrables).
- Le module ne modifie **jamais** ses deux entrées externes (le CRAM et le VCF), montées `:ro`.
- Deux entrées externes de la plateforme : **le CRAM** (`--cram`, + `--fasta` obligatoire pour un CRAM) et **le VCF mono-échantillon** (`--vcf`). Le FASTA de référence est une troisième entrée externe.

### Ressources figées — usage réel

| Ressource | Réellement lue par | Rôle |
|---|---|---|
| `pharmcat_positions.bed` | étage 0, 1, 2a | régions cibles (258 intervalles, ±50 pb), test d'assemblage, base du BED de tranche |
| `positions_exactes.bed` | étage 2b | positions exactes où mesurer la profondeur (`samtools depth -b`) |
| `pharmcat_positions.vcf` | étage 2b | les 1 226 positions diagnostiques attendues + tag `PX=` gène |
| `perimetre_rnpgx.json` | étage 2b + provenance | les 12 gènes du périmètre clinique RNPGx |
| `traductions_cpic_fr.json` | étage 7 | traduction de chaque texte CPIC émis, type de consigne, médicaments restitués, règle HLA-B*15:11 hors PharmCAT |
| `tranche.bed` | non embarquée | le script construit son propre `travail/tranche.bed` dynamiquement (contigs HLA lus dans l'en-tête) |

Note : `alleles_pharmcat.json` et `definitions_alleles.json` ne sont pas utilisées par le pipeline de production (uniquement par l'outillage de validation, hors de ce dépôt).

---

## Table étape par étape

| Étape | Entrées | Traitement | Type de code / outil | Sorties | Échec si |
|---|---|---|---|---|---|
| **0. recevabilité** | CRAM + index, VCF, FASTA+`.fai`, `pharmcat_positions.bed` | Cinq contrôles : VCF mono-échantillon ; concordance d'identité ; CRAM non tronqué ; contigs de même longueur que la référence ; accès aléatoire possible | Bash + `samtools` (hôte) + `bcftools` (conteneur) + awk | `travail/etats.tsv`, `entete_alignement.txt`, `signature.txt`. Si refus : purge de `travail/`, 8 étages `ECHEC`, `provenance.json` seul, `exit 2` | l'un des 5 contrôles échoue |
| **1. filtre** | VCF externe, `pharmcat_positions.bed` | Restreint le VCF aux 258 régions cibles puis indexe | Docker **bcftools** 1.24 | `travail/filtre.vcf.gz` + `.tbi` | `bcftools view` ou `index` échoue |
| **2a. tranche** | CRAM + FASTA, `pharmcat_positions.bed`, en-tête CRAM | Construit `tranche.bed` (cible + MHC + contigs HLA), une seule traversée `view -M -L` → sort → index | Bash + `samtools` (hôte) | `travail/tranche.bam` + `.bai`, `tranche.bed` | `view`, `sort` ou `index` échoue |
| **2b. contrôle qualité** | `filtre.vcf.gz`, `tranche.bam`, `pharmcat_positions.vcf`, `perimetre_rnpgx.json` | Profondeur sur la tranche → `profondeur.txt` ; écarte les positions (5 motifs), établit le périmètre, marque le périmètre clinique ; réécrit le VCF GT seul | `samtools` + **Python `qc_perimetre.py`** + Docker bcftools (bgzip) | `profondeur.txt`, `qc_positions.tsv`, `perimetre.json`, `qualifie.vcf.gz` + `.tbi` | profondeur non mesurée, script échoue, >1 échantillon, perte d'enregistrement |
| **3. CYP2D6** | CRAM complet + FASTA, dépôt Cyrius | `star_caller.py -g 38` sur le CRAM entier ; diplotype non-unitaire → SANS_RESULTAT | **Dépôt Cyrius** (`python3`, hôte) | `cyp2d6.tsv`, `cyp2d6.json` | Cyrius échoue ou sortie vide |
| **4. HLA** | `tranche.bam`, contigs HLA | Extrait les lectures MHC → paires ; si ≥ 200 paires, OptiType type HLA I | `samtools` + Docker **OptiType** 1.3.5 | `mhc_1.fq`, `mhc_2.fq`, `hla/<ts>/<ts>_result.tsv`, `hla.tsv` | < 200 paires ou OptiType échoue |
| **5. appels externes** | `cyp2d6.tsv`, `hla.tsv` | Traduit au format PharmCAT `-po` ; espaces autour du `+` des tandems | **Python `appels_externes.py`** | `appels_externes.tsv` | un étage OK n'a pas déposé sa ligne |
| **6. interprétation** | `qualifie.vcf.gz` (aucun repli), `appels_externes.tsv` | Préprocesseur `--absent-to-ref` puis `pharmcat.jar -reporterJson` ; 22 gènes | Docker **PharmCAT** 3.4.0 | `qualifie.vcf.preprocessed.vcf.bgz`, `sortie/<ECH>.report.json`, `match_warnings.txt` | `qualifie.vcf.gz` absent ; préproc. ou jar échoue |
| **7. compte rendu** | `<ECH>.report.json`, `perimetre.json`, `traductions_cpic_fr.json` | CR PDF : médicaments concernés d'abord, puis 12 gènes, conclusion, limites ; recommandations CPIC fortes ou modérées qui modifient la prise en charge | **Python `compte_rendu.py`** + reportlab | `sortie/CR_<ECH>.pdf` | `report.json`, `perimetre.json` ou traductions absents |
| **provenance** | `etats.tsv`, empreintes CRAM/VCF/FASTA, hash RES, digests images, git Cyrius | Agrège tout ; `reussite_complete` vrai ssi 9 étages présents et aucun échec (SANS_RESULTAT n'est pas un échec) | **Python `provenance.py`** | `sortie/provenance.json` | toujours exécuté ; code de sortie du module |

---

## Pièges pour le bioinformaticien

1. **Le VCF qualifié ne garde que le `GT`.** Les autres champs de format font échouer le préprocesseur ; les génotypes hors seuils sont ramenés à `./.` (jamais « référence »). PharmCAT ne voit plus ni `GQ` ni `DP` ni `AD`.
2. **La profondeur vient de l'alignement, pas du VCF.** Mesurée par `samtools depth` sur `tranche.bam` ; le `DP` du VCF n'est qu'un repli. Un gVCF est accepté, ses blocs ignorés.
3. **Cyrius lit le CRAM entier, pas la tranche.** Seule branche à recevoir le CRAM complet ; ~3 000 régions de normalisation dispersées, d'où sa durée (~la moitié du total).
4. **L'appel externe CYP2D6 exige des espaces autour du `+`.** `*36+*10` deviendrait « indéterminé » en silence ; `*36 + *10` est interprété.
5. **Le périmètre clinique de 12 gènes ne restreint pas ce que PharmCAT reçoit, mais ce que le CR rapporte.** PharmCAT interprète 22 gènes ; le filtrage clinique est en aval, à l'étage 7.
6. **Aucun repli, purge liée à l'empreinte.** L'étage 6 échoue si `qualifie.vcf.gz` manque (pas de rabattement sur `filtre.vcf.gz`). `travail/` n'est réutilisé qu'avec `--reprise` et signature d'entrées identique.
