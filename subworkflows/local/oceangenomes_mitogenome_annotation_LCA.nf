/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//mitogenome
include { DOWNLOAD_BLAST_DB } from '../../modules/local/download_blast_db/main'
include { DOWNLOAD_TAXONKIT_DB } from '../../modules/local/download_taxonkit_db/main'
include { EMMA } from '../../modules/local/EMMA/main'
// include { MITOZ } from '../../modules/local/mitoz/main'
include { BLAST_BLASTN } from '../../modules/nf-core/blast/blastn/main'
// include { BLAST_BLASTP } from '../../modules/nf-core/blast/blastp/main'
include { LCA } from '../../modules/local/LCA/main'

// MultiQC and helper functions
include { MULTIQC                   } from '../../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../../subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_ANNOTATION {

    take:
    mito_assembly //  tuple val(meta), path(fasta)
    curated_blast_db
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Download taxonomy database
    DOWNLOAD_BLAST_DB(Channel.value("taxdb"))
    // Download taxonkit database
    DOWNLOAD_TAXONKIT_DB(Channel.value("taxdump"))

    // Extracts the assembly name from the fasta file and creates new tuple
    fasta_with_assembly_prefix = mito_assembly 
        .map { meta, fasta ->
            def assembly_prefix = fasta.baseName  // Gets "OG898.ilmn.250131.getorg1770"
            [meta, fasta, assembly_prefix]
        }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ANNOTATION, EMMA FOR MAIN ANNOTATION AND MITOZ FOR COMPARISON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
               
    EMMA (
        fasta_with_assembly_prefix // tuple val(meta), path(fasta), val(assembly_prefix)
    )

    // // Function to extract assembly name // This i can do within the process using fast.getbasename()
    def getAnnotationName = { filename ->
        def name = filename.toString().replaceAll(/\.fa$/, '')
        def parts = name.split('\\.', 2)
        return parts.size() > 1 ? parts[1] : name
    }

    // Use mix() to process Co1, 12s and 16s sequences through blast
    combined_sequences = EMMA.out.co1_sequences
        .map { meta, assembly_prefix, files -> 
            def annotationName = getAnnotationName(files.name)
            [meta, assembly_prefix, files, 'CO1', annotationName] 
        }
        .mix(
            EMMA.out.s12_sequences.map { meta, assembly_prefix, files -> 
                def annotationName = getAnnotationName(files.name)
                [meta, assembly_prefix, files, '12s', annotationName] 
            },
            EMMA.out.s16_sequences.map { meta, assembly_prefix, files -> 
                def annotationName = getAnnotationName(files.name)
                [meta, assembly_prefix, files, '16s', annotationName] 
            }
        )
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USING CO1,12s and 16s RUN BLAST TO DETERMINE LCA FOR SPECIES VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    BLAST_BLASTN (
        combined_sequences, // tuple val(meta), val(assembly_name), path(fasta), val(gene_type), val(annotation_name)
        curated_blast_db,
        DOWNLOAD_BLAST_DB.out.db_files // path(db)
    )

    LCA (
        BLAST_BLASTN.out.filtered, // tuple val(meta), path(blast_filtered), val(gene_type), val(assembly_name), val(annotation_name)
        DOWNLOAD_TAXONKIT_DB.out.db_files // path(db)
    )



    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oceangenomes_draftgenomes_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    // //
    // // MODULE: MultiQC
    // //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
  


}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
