module load nextflow/24.04.3
module load singularity/4.1.0

nextflow run main.nf -profile singularity,pawsey_setonix \
    --input assets/samplesheet.csv \
    --outdir /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR \
    --gxdb \
    --binddir /scratch \
    --tempdir /scratch/pawsey0964/tpeirce/tmp \
    -resume \
    -with-report
    
