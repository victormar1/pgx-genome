#!/usr/bin/env nextflow
/*
 * Recette d'integration pour une plateforme, a adapter.
 *
 * Le module orchestre lui-meme ses conteneurs : il appelle l'interpreteur, le
 * typeur HLA et bcftools dans des images epinglees par empreinte. Ce processus
 * ne doit donc PAS etre execute dans un conteneur — il doit tourner sur un
 * noeud ou le moteur de conteneurs est disponible, et c'est lui qui lance les
 * images.
 *
 * Ce choix se discute. Il a un avantage : la chaine reste identique en ligne
 * de commande et sous ordonnanceur, donc ce qui a ete valide est ce qui tourne.
 * Il a un cout : l'ordonnanceur ne voit qu'une tache par genome, et ne peut pas
 * allouer les ressources etage par etage. Une decomposition en un processus par
 * etage est possible et souhaitable a terme ; elle demande de sortir
 * l'orchestration du script bash, et elle devra etre revalidee.
 *
 * Usage :
 *   nextflow run exemples/pgx_genome.nf \
 *     --manifeste manifeste.tsv --sortie resultats --racine /chemin/du/depot \
 *     --fasta /chemin/GRCh38.fa -profile slurm
 *
 * Le manifeste est celui du module : echantillon, alignement, variants,
 * separes par des tabulations, sans en-tete.
 */

nextflow.enable.dsl = 2

params.manifeste = null
params.sortie    = 'resultats'
params.racine    = projectDir.parent          // la racine du depot
params.fasta     = null
params.cyrius    = null                       // depot de l'appeleur CYP2D6
params.pypgx     = ''                         // vide : etage 4b non arme
params.moteur    = 'docker'                   // ou apptainer
params.identites = null                       // dossier de fichiers d'identite

process PGX_GENOME {
    tag "${echantillon}"

    // Un genome occupe un coeur pour l'essentiel de sa duree, et seize le temps
    // du typage HLA. Les durees viennent de la mesure sur deux cent quarante
    // trois genomes : mediane 295 s, maximum observe 306 s pour le seul etage
    // HLA. La memoire n'a pas ete mesuree par etage — a regler sur le site
    // apres une premiere mesure, et non a deviner ici.
    cpus   16
    time   '2h'
    errorStrategy 'retry'
    maxRetries 1

    publishDir "${params.sortie}", mode: 'copy', pattern: "${echantillon}/sortie/*"

    input:
    tuple val(echantillon), path(alignement), path(variants), path(index_al), path(index_va)

    output:
    tuple val(echantillon), path("${echantillon}/sortie/provenance.json"), emit: provenance
    path "${echantillon}/sortie/*", emit: sortie
    path "${echantillon}/journal.txt", emit: journal, optional: true

    script:
    def ident = params.identites ? "--identite ${params.identites}/${echantillon}.json" : ''
    """
    export PGX_MOTEUR='${params.moteur}'
    export PGX_FASTA='${params.fasta}'
    ${params.cyrius ? "export PGX_CYRIUS='${params.cyrius}'" : ''}
    ${params.pypgx  ? "export PGX_PYPGX='${params.pypgx}'"   : ''}

    bash '${params.racine}/bin/pgx_genome.sh' \\
      --cram '${alignement}' --vcf '${variants}' --fasta '${params.fasta}' \\
      --sortie '${echantillon}' --echantillon '${echantillon}' ${ident}
    """
}

process BILAN {
    publishDir "${params.sortie}", mode: 'copy'

    input:
    path provenances

    output:
    path 'bilan.tsv'

    script:
    """
    printf 'echantillon\\treussite_complete\\tetages_en_echec\\n' > bilan.tsv
    for f in ${provenances}; do
      python3 -c "
import json, sys
d = json.load(open(sys.argv[1], encoding='utf-8'))
print('%s\\t%s\\t%s' % (d.get('echantillon', '?'),
                        d.get('reussite_complete'),
                        ','.join(d.get('etages_en_echec') or [])))
" "\$f" >> bilan.tsv
    done
    """
}

workflow {
    if (!params.manifeste) { exit 1, 'il manque --manifeste' }
    if (!params.fasta)     { exit 1, 'il manque --fasta : un CRAM ne se decode pas sans sa reference' }

    genomes = Channel
        .fromPath(params.manifeste)
        .splitCsv(sep: '\t')
        .map { ligne ->
            def al = file(ligne[1]), va = file(ligne[2])
            // Les index accompagnent obligatoirement leurs fichiers : sans eux
            // l'etage de recevabilite refuse, et l'acces par region est
            // impossible.
            def ial = file(al.toString() + (al.name.endsWith('.cram') ? '.crai' : '.bai'))
            def iva = file(va.toString() + '.tbi')
            tuple(ligne[0], al, va, ial, iva)
        }

    PGX_GENOME(genomes)
    BILAN(PGX_GENOME.out.provenance.map { it[1] }.collect())
}
