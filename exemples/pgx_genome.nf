#!/usr/bin/env nextflow
/*
 * Recette d'integration pour une plateforme, a adapter.
 *
 * Un processus par etage. Le module garde son orchestration : chaque processus
 * appelle le meme script avec --etage, si bien que ce qui a ete valide d'un
 * seul tenant est ce qui tourne ici, etage par etage. Ce que l'ordonnanceur
 * gagne : reserver a chaque etage ce qu'il consomme, faire tourner ensemble
 * CYP2D6, le HLA et le typage complementaire, et relancer un etage sans rejouer
 * les autres.
 *
 * Les coeurs declares sont mesures, etage par etage, sur genomes reels : un
 * genome occupe en moyenne 1,3 coeur. Le module recoit la reservation comme
 * nombre de fils : samtools, Cyrius et OptiType n'en lancent pas davantage.
 *
 * Le module orchestre lui-meme ses conteneurs. Ces processus ne doivent donc PAS
 * etre executes dans un conteneur : ils tournent sur un noeud ou le moteur de
 * conteneurs est disponible, et ce sont eux qui lancent les images.
 *
 * Les etages se transmettent leur travail par un dossier par genome, sous
 * --sortie, et non par le repertoire de travail de Nextflow : ce dossier doit
 * etre sur un systeme de fichiers partage par les noeuds. Les entrees passent
 * de meme en chemins absolus, et non en fichiers mis en scene : Nextflow les
 * mettrait en liens symboliques dans le repertoire de la tache, que le module
 * monte dans ses conteneurs, ou le lien pointerait vers un chemin non monte.
 *
 * Codes de retour d'un etage lance seul :
 *   0  l'etage a abouti, ou n'a rien conclu, ou etait optionnel et non arme
 *   1  l'etage a echoue ; l'echec est inscrit, la chaine continue, et le
 *      genome sortira sans reussite complete
 *   2  le genome est refuse, a l'entree ou apres le controle qualite ; aucun
 *      etage suivant ne tourne pour lui
 *   3  un autre processus tient deja cet etage pour ce genome
 *
 * Executee sous Nextflow 25.04, executeur local, vingt-quatre coeurs. Sur deux
 * genomes et un alignement tronque : rendus identiques a l'execution complete,
 * refus du tronque sans arret des autres. Sur des lots de dix genomes : rendus
 * identiques au lot, 51 genomes par heure a froid, contre 39 au lanceur de lot
 * a six. Sous un ordonnanceur de grappe, elle reste a eprouver sur le site.
 *
 * Usage :
 *   nextflow run exemples/pgx_genome.nf \
 *     --manifeste manifeste.tsv --sortie /partage/resultats \
 *     --racine /chemin/du/depot --fasta /chemin/GRCh38.fa -profile slurm
 */

nextflow.enable.dsl = 2

params.manifeste = null
params.sortie    = null                       // chemin partage, absolu
params.racine    = projectDir.parent          // la racine du depot
params.fasta     = null
params.cyrius    = null                       // depot de l'appeleur CYP2D6
params.pypgx     = ''                         // vide : etage 4b non arme
params.moteur    = 'docker'                   // ou apptainer
params.identites = null                       // dossier de fichiers d'identite

// L'appel commun. Le code 1 est un echec inscrit : la chaine continue, comme en
// execution complete. Les codes 2 et 3 remontent a l'ordonnanceur.
def etage(String noms, String e, al, va, fils, String extra = '') {
    """
    export PGX_MOTEUR='${params.moteur}'
    export PGX_FASTA='${params.fasta}'
    export PGX_FILS='${fils}'
    ${params.cyrius ? "export PGX_CYRIUS='${params.cyrius}'" : ''}
    ${params.pypgx  ? "export PGX_PYPGX='${params.pypgx}'"   : ''}

    rc=0
    bash '${params.racine}/bin/pgx_genome.sh' \\
      --cram '${al}' --vcf '${va}' --fasta '${params.fasta}' \\
      --sortie '${params.sortie}/${e}' --echantillon '${e}' \\
      --etage '${noms}' ${extra} || rc=\$?
    [ "\$rc" -eq 1 ] && rc=0
    exit "\$rc"
    """
}

// Un genome refuse sort de la chaine ; un etage tenu ailleurs est relance plus
// tard ; toute autre erreur arrete tout, parce qu'elle n'est pas prevue.
def strategie(code) {
    code == 2 ? 'ignore' : (code == 3 ? 'retry' : 'terminate')
}

// Coeurs occupes, mesures etage par etage : moyenne, puis pointe d'une seconde.
// La reservation couvre la moyenne ; une pointe breve deborde sans dommage.

// Recevabilite, filtre et tranche : 1,2 en moyenne, 3,4 en pointe au tri.
process OUVRIR {
    tag "${e}"
    cpus 2
    errorStrategy { strategie(task.exitStatus) }
    maxRetries 3
    input:  tuple val(e), val(al), val(va)
    output: tuple val(e), val(al), val(va)
    script: etage('recevabilite,filtre,tranche', e, al, va, task.cpus)
}

// Un seul fil : 1,0 en moyenne.
process CONTROLER {
    tag "${e}"
    cpus 1
    errorStrategy { strategie(task.exitStatus) }
    maxRetries 3
    input:  tuple val(e), val(al), val(va)
    output: tuple val(e), val(al), val(va)
    script: etage('qc', e, al, va, task.cpus)
}

// Les trois etages independants.
// CYP2D6 : 1,9 en moyenne, 3,8 en pointe, sur quatre fils.
process CYP2D6 {
    tag "${e}"
    cpus 4
    errorStrategy { strategie(task.exitStatus) }
    maxRetries 3
    input:  tuple val(e), val(al), val(va)
    output: tuple val(e), val('cyp2d6')
    script: etage('cyp2d6', e, al, va, task.cpus)
}

// HLA : 1,1 en moyenne, 4,4 en pointe au placement des lectures, sur quatre
// fils. Le typage est identique a quatre et a seize fils, dans le meme temps.
process HLA {
    tag "${e}"
    cpus 4
    errorStrategy { strategie(task.exitStatus) }
    maxRetries 3
    input:  tuple val(e), val(al), val(va)
    output: tuple val(e), val('hla')
    script: etage('hla', e, al, va, task.cpus)
}

// Typage complementaire : 3,7 en moyenne, 10,9 en pointe ; il ne prend pas de
// nombre de fils.
process COMPLEMENT {
    tag "${e}"
    cpus 4
    errorStrategy { strategie(task.exitStatus) }
    maxRetries 3
    input:  tuple val(e), val(al), val(va)
    output: tuple val(e), val('pypgx')
    script: etage('pypgx', e, al, va, task.cpus)
}

// Appels, interpretation, compte rendu et trace : 2,0 en moyenne.
process INTERPRETER {
    tag "${e}"
    cpus 2
    errorStrategy { strategie(task.exitStatus) }
    maxRetries 3
    input:  tuple val(e), val(al), val(va)
    output: val(e)
    script:
    def ident = params.identites ? "--identite ${params.identites}/${e}.json" : ''
    etage('appels,pharmcat,rendu,provenance', e, al, va, task.cpus, ident)
}

process BILAN {
    publishDir "${params.sortie}", mode: 'copy'
    input:
    val(echantillons)
    val(attente)
    output:
    path 'bilan.tsv'
    // La trace se lit dans le dossier partage : un genome refuse en amont n'a
    // pas traverse l'interpretation, mais il a sa trace, et elle figure ici.
    script:
    """
    printf 'echantillon\\treussite_complete\\tetages_en_echec\\n' > bilan.tsv
    for e in ${echantillons.join(' ')}; do
      f='${params.sortie}'/"\$e"/sortie/provenance.json
      if [ ! -s "\$f" ]; then printf '%s\\tabsente\\t\\n' "\$e" >> bilan.tsv; continue; fi
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
    if (!params.sortie)    { exit 1, 'il manque --sortie : un chemin partage par les noeuds' }
    if (!params.fasta)     { exit 1, 'il manque --fasta : un CRAM ne se decode pas sans sa reference' }

    // Chemins absolus : le module les lit lui-meme. Leurs index doivent les
    // accompagner, faute de quoi l'etage de recevabilite refuse le genome.
    genomes = Channel
        .fromPath(params.manifeste)
        .splitCsv(sep: '\t')
        .map { ligne ->
            tuple(ligne[0],
                  file(ligne[1]).toAbsolutePath().toString(),
                  file(ligne[2]).toAbsolutePath().toString())
        }

    controles = CONTROLER(OUVRIR(genomes))

    // Les trois etages independants partent ensemble ; l'interpretation attend
    // les trois, pour les seuls genomes qu'aucun etage n'a refuses.
    trois = CYP2D6(controles)
        .join(HLA(controles))
        .join(COMPLEMENT(controles))
    INTERPRETER(controles.join(trois).map { it[0..2] })

    // Le bilan attend la fin de l'interpretation, et porte sur tous les genomes
    // du manifeste, refuses compris. La valeur fictive evite qu'un lot ou tout
    // est refuse ne laisse le bilan attendre indefiniment.
    attente = INTERPRETER.out.mix(Channel.of('-')).collect()
    BILAN(genomes.map { it[0] }.collect(), attente)
}
