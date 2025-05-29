/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//mitogenome
include { GETORGANELLE_CONFIG } from '../../modules/nf-core/getorganelle/config/main'
include { GETORGANELLE_FROMREADS } from '../../modules/nf-core/getorganelle/fromreads/main'
// include { EMMA } from '../../modules/local/emma/main'
// include { MITOZ } from '../../modules/local/mitoz/main'
include { BLAST_BLASTN } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTP } from '../../modules/nf-core/blast/blastp/main'
// include { LCA } from '../../modules/local/lca/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOMES {

    take:
    fastp_reads
    organelle_type
    
    
    main:
    // ch_versions = Channel.empty()
    // ch_multiqc_files = Channel.empty()
    // Input channels that are set up to run from different stages. If you provide the assembled genome path then it will skip the getorganelle
    // Input channels
    // If params.fastq is provided to the filtered and trimmed read the subworkflow can be run seperate from the main workflow
    // ch_fastq = params.fastq ? Channel.fromPath(params.fastq) : Channel.empty()
    // // If params.mito is provided to the assembled mitogenome then the workflow will start from annotation
    // ch_mito = params.mito ? Channel.fromPath(params.mito) : Channel.empty()
    // // if params.anno is provided to the annotated mitogenome fasta file then the workflow will start from the LCA part of the workflow
    // ch_anno = params.anno ? Channel.fromPath(params.anno) : Channel.empty()
    
    // take:
    // ch_samplesheet // channel: samplesheet read in from --input
    // main:

    // ch_versions = Channel.empty()
    // ch_multiqc_files = Channel.empty()

 
    
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ASSEMBLY USING GET ORGANELLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    GETORGANELLE_CONFIG (
        organelle_type // val(organelle_type)
    )

    GETORGANELLE_FROMREADS (
        fastp_reads, // tuple val(meta), path(fastq)
        GETORGANELLE_CONFIG.out.db // tuple val(organelle_type), path(db)  // getOrganelle has a database and config file
    )

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     RUN ANNOTATION, EMMA FOR MAIN ANNOTATION AND MITOZ FOR COMPARISON
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

//     EMMA (
//         // tuple val(og_num), val(prefix), path(fasta)
//     )

//     MITOZ (
        
//     )

//     // maybe a python module to run the code to extract the cds and translate them
//     EXTRACT_SEQUENCE (

//     )

//     TRANSLATE (

//     )

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     USING CO1,12s and 16s RUN BLAST TO THE DETERMINE LCA FOR SPECIES VALIDATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

//     BLAST_BLASTN (
//         // tuple val(meta) , path(fasta)
//         // tuple val(meta2), path(db)
//     )

//     LCA (
//         // tuple val(og_num), val(prefix), val(emma_prefix), path(lca)
//     )



    //
    // Collate and save software versions
    //
    // softwareVersionsToYAML(ch_versions)
    //     .collectFile(
    //         storeDir: "${params.outdir}/pipeline_info",
    //         name: 'nf_core_'  +  'oceangenomesdraftgenomes_software_'  + 'mqc_'  + 'versions.yml',
    //         sort: true,
    //         newLine: true
    //     ).set { ch_collated_versions }


    // //
    // // MODULE: MultiQC
    // //
    // ch_multiqc_config        = Channel.fromPath(
    //     "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config = params.multiqc_config ?
    //     Channel.fromPath(params.multiqc_config, checkIfExists: true) :
    //     Channel.empty()
    // ch_multiqc_logo          = params.multiqc_logo ?
    //     Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
    //     Channel.empty()

    // summary_params      = paramsSummaryMap(
    //     workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
    //     file(params.multiqc_methods_description, checkIfExists: true) :
    //     file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(
    //     methodsDescriptionText(ch_multiqc_custom_methods_description))

    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_methods_description.collectFile(
    //         name: 'methods_description_mqc.yaml',
    //         sort: true
    //     )
    // )

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList(),
    //     [],
    //     []
    // )

    // emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    // versions       = ch_versions                 // channel: [ path(versions.yml) ]



}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
