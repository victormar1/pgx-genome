# Génome de qualification

Un génome public, son fichier de variants restreint aux régions que le module
lit, et le rapport que ce module en a produit. Une installation est qualifiée
quand elle rend le même rapport.

| | |
|---|---|
| Prélèvement | HG00699, 1000 Genomes Project, séquençage 30× NYGC, GRCh38 |
| Alignement | `ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241747/HG00699.final.cram` et son `.crai`, MD5 `0bd9d2f5f91bdbecda1f19968c59a4e6` |
| Variants | `HG00699.vcf.gz`, appels de la plateforme NYGC (GATK) sur les régions de l'interpréteur, complétés par des appels bcftools sur les régions du typage complémentaire (`doc/VALIDATION.md`, § 3.3), restreints aux régions utiles et au contig mitochondrial ; rend le même rapport que le fichier complet |
| Référence | `GRCh38_full_analysis_set_plus_decoy_hla.fa`, celle de l'alignement |
| Rapport de référence | `HG00699/sortie/HG00699.report.json`, produit par le module en version 0.4.2 avec les images de `env/conteneurs.txt` |

```bash
bin/pgx_genome.sh --cram /chemin/HG00699.final.cram --vcf bench/qualification/HG00699.vcf.gz \
  --fasta /chemin/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --sortie /chemin/qualification/HG00699 --echantillon HG00699
python3 outils/comparer_rapports.py --avant bench/qualification --apres /chemin/qualification
```

Attendu : code de retour `0` aux deux commandes, `0 écart` sur 12 couples gène ×
génome. Ce génome porte HLA-B\*58:01 et CYP3A5 \*1/\*3 : le compte rendu doit
contre-indiquer l'allopurinol et majorer la dose initiale de tacrolimus.

Un écart désigne une différence d'installation — version d'une image, référence,
typeur complémentaire — à résoudre avant tout usage sur un patient.
