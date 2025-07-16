process EMMA {
    tag "meta.id" // tag is is being used for publish dir
    label 'process_medium'
    // No native Singularity container available, using Docker image for both engines
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/emma:0.8' :
        'tylerpeirce/emma:0.8' }"

    input:
    tuple val(meta), path(fasta), val(assembly_prefix)

    output:
    tuple val(meta), val(assembly_prefix), path("${prefix}/emma/cds/*CO1*.fa"),  emit: co1_sequences
    tuple val(meta), val(assembly_prefix), path("${prefix}/emma/cds/*RNR1*.fa"),  emit: s12_sequences  
    tuple val(meta), val(assembly_prefix), path("${prefix}/emma/cds/*RNR2*.fa"),  emit: s16_sequences
    tuple val(meta), path("${assembly_prefix}") 
    path "emma/versions_emma.yml"

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: assembly_prefix

        """
            emma_version=\$(cat /opt/Emma/Project.toml | grep version | sed -E 's/version = "([0-9]+)\\.([0-9]+)\\.([0-9]+)"/\\1\\2\\3/')
                     
            emma_prefix="${prefix}.emma"\${emma_version}""
            
            mkdir -p tempdir cds emma
            
            /opt/julia-1.10.5/bin/julia \\
                /opt/Emma/src/command.jl \\
                $args \\
                --fa emma/\${emma_prefix}.fa \\
                --gff emma/\${emma_prefix}.gff \\
                --tbl emma/\${emma_prefix}.tbl \\
                --svg emma/\${emma_prefix}.svg \\
                --tempdir tempdir/ \\
                --loglevel debug \\
                ${fasta} 
            
            /opt/julia-1.10.5/bin/julia /opt/extract_proteins.jl \\
                emma/ \\
                emma/
            
            mv cds emma/  
            
            # Add in the sample to the file names in cds
            for file in emma/cds/*; do
                new_file="\${file%.fa}.\${emma_prefix}.fa"     
                mv "\$file" "\$new_file"
            done
            # Creat file structure and move results
            mkdir -p ${prefix}
            mv emma ${prefix}/

            cat <<-END_VERSIONS > emma/versions_emma.yml
            "${task.process}":
                julia: \$(julia -v | sed 's/^julia version //g' )
                emma: \$(cat /opt/Emma/Project.toml | grep version)
            END_VERSIONS
    
        """

 stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${assembly_prefix}"
    
        """
        # Create mock emma version for stub
        emma_version=\$(cat /opt/Emma/Project.toml | grep version | sed -E 's/version = "([0-9]+)\\.([0-9]+)\\.([0-9]+)"/\\1\\2\\3/')
        emma_prefix="${prefix}.emma\${emma_version}"
        
        # Create directory structure
        mkdir -p emma/cds
        
        # Create mock output files
        touch emma/\${emma_prefix}.fa
        touch emma/\${emma_prefix}.gff
        touch emma/\${emma_prefix}.tbl
        touch emma/\${emma_prefix}.svg
        
        # Create mock CDS files
        touch emma/cds/gene1.\${emma_prefix}.fa
        touch emma/cds/gene2.\${emma_prefix}.fa
        touch emma/cds/gene3.\${emma_prefix}.fa
        
        # Create mock versions file
        cat <<-END_VERSIONS > emma/versions_emma.yml
        "${task.process}":
            julia: \$(julia -v | sed 's/^julia version //g' )
            emma: \$(cat /opt/Emma/Project.toml | grep version)
        END_VERSIONS

        """
}