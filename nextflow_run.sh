module load nextflow/24.04.3
module load singularity/4.1.0-nompi

nextflow run main.nf -profile singularity \
    --outdir /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR \
    --curated_blast_db /scratch/pawsey0964/pbayer/OceanGenomes.CuratedNT.NBDLTranche1and2.CuratedBOLD.fasta \
    --bs_config ~/.basespace/default.cfg \
    --gxdb \
    --binddir /scratch \
    --tempdir /scratch/pawsey0964/tpeirce/tmp \
    -c pawsey_profile.config \
    -resume \
    --refresh-modules \
    -with-report
    #--input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    
