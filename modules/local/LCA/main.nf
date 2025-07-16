process LCA {
    tag "$run_id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/lca:0.1' :
        'tylerpeirce/lca:0.1' }"

    input:
    tuple val(meta), path(blast_filtered), val(gene_type), val(assembly_name), val(annotation_name)
    path(taxdump)

    output:
    path("lca"), emit: fastqs

    script:
    """          
    export TAXONKIT_DB=$taxdump

    python /computeLCA.py \\
        $blast_filtered \\
        > lca.${gene_type}.${annotation_name}.tsv
            
    sed -i "s/\$/\\t\$(date +%y%m%d)/" lca.${gene_type}.${annotation_name}.tsv
    wait
    sed -i "s/\$/\\t${gene_type}/" lca.${gene_type}.${annotation_name}.tsv

    mv lca.${gene_type}.${annotation_name}.tsv ${assembly_name}/lca/

    cat <<-END_VERSIONS > versions_LCA.yml
    "${task.process}":
        Python: \$(python -V | sed 's/Python //g')
        TaxonKit: \$(/taxonkit version | sed 's/taxonkit //g')
    END_VERSIONS
    
    """
}