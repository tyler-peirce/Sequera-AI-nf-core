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

// Helper functions
include { softwareVersionsToYAML    } from '../../subworkflows/nf-core/utils_nfcore_pipeline'


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

    // Collect MultiQC files
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.blast.collect{it[1]})
    ch_versions = ch_versions.mix(EMMA.out.versions.first())
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
    ch_versions = ch_versions.mix(LCA.out.versions.first())



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




    emit:
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]
  


}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
