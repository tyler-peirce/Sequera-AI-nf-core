process DOWNLOAD_TAXONKIT_DB {
    tag "${db_name}"
    label 'process_low'
    storeDir params.taxonkit_db_dir ?: "${launchDir}/taxonkit_dbs"  // Cache the database
    
    input:
    val db_name
    
    output:
    tuple path("*.dmp"), path("gc.prt"), path("readme.txt"), emit: db_files
   
    script:
    """
    echo "Downloading taxonomy database..."
    
    # Check if files already exist
    if [ -f "${launchDir}/taxonkit_dbs/citations.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/delnodes.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/division.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/gencode.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/images.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/merged.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/names.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/nodes.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/gc.prt" ]; then
        echo "Taxonomy files already exist, skipping download"
    else
        echo "Downloading fresh taxonomy database..."
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        
        wait

        tar -xzf taxdump.tar.gz
        rm taxdump.tar.gz
    fi
    
    # Verify files
    ls -la *.dmp *.prt
    """
    
    stub:
    """
    touch citations.dmp delnodes.dmp division.dmp gencode.dmp images.dmp merged.dmp names.dmp nodes.dmp gc.prt
    """
}