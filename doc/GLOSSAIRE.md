# Glossaire

Les noms d'étage sont ceux de l'option `--etage` ; le reste suit le vocabulaire
du domaine.

## Pipeline

| Terme | Sens | En anglais |
|---|---|---|
| Étage | étape du pipeline, lançable seule par `--etage` ; neuf sont requis, le 4b est optionnel | stage |
| Recevabilité (`recevabilite`) | contrôle des entrées avant tout calcul : VCF mono-échantillon déclarant `AD`, identité concordante, CRAM non tronqué, référence conforme ; au sens de la recevabilité d'un échantillon en biologie médicale | input validation |
| Tranche (`tranche`, `tranche.bam`) | extraction de l'alignement restreinte aux régions du panel et au complexe HLA, en une seule lecture du CRAM | targeted extraction, sub-BAM |
| Contrôle qualité (`qc`) | mesure profondeur et qualité position par position, écarte ce qui ne passe pas les seuils et établit le périmètre | QC |
| Typage complémentaire (`pypgx`, étage 4b) | typage par PyPGx de MT-RNR1, BCHE, MTHFR et POR ; nom de la section correspondante du compte rendu | PyPGx typing |
| Activé | se dit de l'étage 4b quand une image (`PGX_IMG_PYPGX`) ou une installation (`PGX_PYPGX`) de PyPGx est fournie | enabled |
| Outside calls (`appels`) | génotypes établis hors du VCF (CYP2D6 par Cyrius, HLA par OptiType, MT-RNR1 par PyPGx), transmis à PharmCAT par son option `-po` | outside calls |
| Interpréteur | PharmCAT : traduit les génotypes en diplotypes, phénotypes et recommandations CPIC | interpreter |
| Trace d'exécution (`provenance.json`) | hash des entrées, digests des images, versions, seuils et état de chaque étage | provenance |
| Ordonnanceur | ce qui répartit les étages sur les machines : Nextflow, gestionnaire de grappe | scheduler |

## Résultats

| Terme | Sens | En anglais |
|---|---|---|
| Périmètre clinique | les gènes rendus au compte rendu, fixés par `ressources/perimetre_rnpgx.json` | clinical scope |
| Complément RNPGx | positions du panel RNPGx absentes de la table de PharmCAT, mesurées à part et rendues sans interprétation | — |
| Non conclusif | gène dont une position définissante n'est pas lue ou pas qualifiée : aucun diplotype n'est rendu | no-call |
| Partiel | gène rendu alors qu'une partie de ses positions manque ; la proportion lue figure au compte rendu | partial call |
| Critère d'acceptation | seuil qui refuse le génome au sortir du contrôle qualité, sans compte rendu : couverture médiane, part d'hétérozygotes déséquilibrés | acceptance criterion, QC gate |
| `reussite_complete` | vrai si les neuf étages requis ont abouti et qu'aucun étage, 4b compris s'il est activé, n'est en échec | — |

## Validation

| Terme | Sens | En anglais |
|---|---|---|
| Benchmark | génomes publics du 1000 Genomes confrontés à GeT-RM, au typage HLA Sanger et à une puce de génotypage (`bench/`) | benchmark |
| Génome de qualification | génome public livré avec son rapport de référence ; une installation est qualifiée quand elle le reproduit sans écart | installation qualification |
| Variant caller | outil qui produit le VCF en amont du module, par exemple GATK ou bcftools | variant caller |
| Downsampling | réduction aléatoire des lectures d'un alignement, pour mesurer le comportement du module à couverture réduite | downsampling |
| Contrôle positif synthétique | VCF construit pour porter un allèle absent du benchmark (`outils/cas_construit.py`) ; il vérifie la chaîne logicielle, pas la détection | synthetic positive control |
| Couple gène × génome | unité de comparaison entre deux exécutions : diplotype, phénotype et recommandations d'un gène pour un génome | gene–sample pair |
| Lot | génomes traités ensemble par `bin/pgx_lot.sh` | batch |

## Développement

| Terme | Sens | En anglais |
|---|---|---|
| Mutation testing | réintroduire volontairement un défaut, le mutant, pour vérifier qu'un test échoue ; un mutant détecté est dit tué | mutation testing |
| Hook | script Git lancé au commit (`commit-msg`) ou au push (`pre-push`), rangé dans `outils/crochets/` | Git hook |
| CI | GitHub Actions rejoue tests, contrôles du dépôt et validation des messages à chaque push | continuous integration |
| Digest | identifiant sha256 d'une image de conteneur ; contrairement à un tag, il ne peut pas être repointé | image digest |
| Trailer | ligne de fin d'un message de commit, au format `Jeton: valeur` | trailer |
