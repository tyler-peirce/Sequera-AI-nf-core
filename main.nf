#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oceangenomes_draftgenomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/oceangenomes_draftgenomes
    Website: https://nf-co.re/oceangenomes_draftgenomes
    Slack  : https://nfcore.slack.com/channels/oceangenomes_draftgenomes
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DRAFTGENOMES  } from './workflows/draftgenomes'
include { MITOGENOMES  } from './subworkflows/local/oceangenomes_mitogenomes'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow OCEANGENOMES_DRAFTGENOMES {

    take:
    //samplesheet // channel: samplesheet read in from --input
        run_id
        bs_config

    main:

    //
    // WORKFLOW: Run pipeline
    //
    DRAFTGENOMES (
        //samplesheet
        run_id,
        bs_config,
    )

    println "DRAFTGENOMES emits: ${DRAFTGENOMES.out.fastp_reads}"


    MITOGENOMES (
        fastp_reads    = DRAFTGENOMES.out.fastp_reads,
        organelle_type = "animal_mt"  // << pass it in here
    )
    emit:
    multiqc_report = DRAFTGENOMES.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        // params.input // this is now inncliuded in the samplesheetHybrid wf
    )

    //
    // WORKFLOW: Run main workflow
    //
    OCEANGENOMES_DRAFTGENOMES (
        //PIPELINE_INITIALISATION.out.samplesheet
        params.run,
        params.bs_config
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        OCEANGENOMES_DRAFTGENOMES.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
