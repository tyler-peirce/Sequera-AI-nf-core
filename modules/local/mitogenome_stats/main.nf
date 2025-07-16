process BASESPACE {
    tag "$run_id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2' :
        'tylerpeirce/psycopg2' }"

    input:
    val run_id
    path config

    output:
    path "*fastq.gz", emit: fastqs
    path "*json"    , emit: jsons

    script:
    """
# this code pulls all of the stats for the mitogenomes on Acaia

# Create output file
echo -e "OG_id\ttech\tseq_date\tcode\tannotation\tstats\tlength\tlength_emma\tseqlength_12s\tseqlegth_16s\tseqlength_CO1\tcds_no\ttrna_no\trrna_no\tstatus\tgenbank\trrna12s\trrna16s\tatp6\tatp8\tcox1\tcox2\tcox3\tcytb\tnad1\tnad2\tnad3\tnad4\tnad4l\tnad5\tmad6\ttRNA_Phe\ttRNA_Val\ttRNA_LeuUAG\ttRNA_LeuUAA\ttRNA_Ile\ttRNA_Met\ttRNA_Thr\ttRNA_Pro\ttRNA_Lys\ttRNA_Asp\ttRNA_Glu\ttRNA_SerGCU\ttRNA_SerUGA\ttRNA_Tyr\ttRNA_Cys\ttRNA_Trp\ttRNA_Ala\ttRNA_Asn\ttRNA_Gly\ttRNA_Arg\ttRNA_His\ttRNA_Gln"  > mtdnastat.${meta.id}.tsv

#create list of all the samples that have mitogenomes
rundir=/scratch/pawsey0964/tpeirce/_MITOGENOMES/NOVA_250606_AD

for og in $rundir/OG*/OG*; do
    echo "og = $og"

        # Pull out the stats of the mitogenome 
        if [[ $og == *.hifi.* ]]; then

            # Details on whether the mitogenome is circular if hifi
            stats=$(rclone cat "$og"/mtdna/contigs_stats.tsv | grep "final_mitogenome" | awk -F "\t" '{print $6}' | tr -d '[:space:]') # Remove whitespace
            if [ "$stats" = True ]; then
                stats="circular genome"
            else
                stats="not circular"
            fi
        
        else
            # Details on whether the mitogenome is circular if getorg
            stats=$(rclone cat "$og"/mtdna/ --include "*get_org.log.txt" --max-depth 1 | grep "Result status of animal_mt: " | tail -n 1 | awk -F 'animal_mt: ' '{print $NF}')
        fi                
           
        # Check if there is a genbank file
        genbank=$(rclone cat "$og"/genbank/ --include "*.{fa,fasta}" --max-depth 1 | grep ">" | sed 's/>//g')

        # Get the annoataion code
        annotation=$(rclone lsf "$og"/emma/ --include "*.{fa,fasta}" --max-depth 1 | head -n 1 | awk -F'.' '{print $5}')

        # Count the length of the sequence
        length=$(rclone cat "$og"/mtdna/ --include "*.{fa,fasta}" --max-depth 1 | grep -v "^>" | tr -d '\n' | awk '{print length}')
        length_emma=$(rclone cat "$og"/emma/ --include "*.{fa,fasta}" --max-depth 1 | grep -v "^>" | tr -d '\n' | awk '{print length}')

        # Count the length of the CO1, 12s and 16s sequences
        seq_12s=$(rclone cat "$og"/emma/cds/ --include "MT-RNR1.*" --include "*12s*" | grep -v "^>" | tr -d '\n' | awk '{print length}')
        seq_16s=$(rclone cat "$og"/emma/cds/ --include "MT-RNR2.*" --include "*16s*" | grep -v "^>" | tr -d '\n' | awk '{print length}')
        seq_CO1=$(rclone cat "$og"/emma/cds/ --include "MT-CO1.*" --include "*COX1*" | grep -v "^>" | tr -d '\n' | awk '{print length}')
        
        ###need to change this as it doesnt do it for the scaffolds### check this for scaffolds... maybe there is just no gff
        gff=""$og"/emma/ --include "*.gff""
        cds_no=$(rclone cat $gff | awk '$3 == "CDS" { count++ } END { print count }') # count the number of cds in mitogenome
        trna_no=$(rclone cat $gff | awk '$3 == "tRNA" { count++ } END { print count }') # count the number of tRNAs in mitogenome
        rrna_no=$(rclone cat $gff | awk '$3 == "rRNA" { count++ } END { print count }') # count the number of rRNAs in mitogenome

        # Length of genome feature from gff
        rrna12s=$(rclone cat $gff | grep -i 12S | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        rrna16s=$(rclone cat $gff | grep -i 16S | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        atp6=$(rclone cat $gff | grep ATP6 | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        atp8=$(rclone cat $gff | grep ATP8 | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        cox1=$(rclone cat $gff | grep -E "(CO1|COX1)" | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        cox2=$(rclone cat $gff | grep -E "(CO2|COX2)" | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        cox3=$(rclone cat $gff | grep -E "(CO3|COX3)" | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        cytb=$(rclone cat $gff | grep CYTB | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        nad1=$(rclone cat $gff | grep ND1 | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        nad2=$(rclone cat $gff | grep ND2 | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        nad3=$(rclone cat $gff | grep ND3 | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        nad4=$(rclone cat $gff | grep -w "ND4" | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        nad4l=$(rclone cat $gff | grep ND4L | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        nad5=$(rclone cat $gff | grep ND5 | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        mad6=$(rclone cat $gff | grep ND6 | grep CDS | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')

        tRNA_Phe=$(rclone cat $gff | grep GAA | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Val=$(rclone cat $gff | grep UAC | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_LeuUAG=$(rclone cat $gff | grep UAG | grep tRNA | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Ile=$(rclone cat $gff | grep GAU | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Met=$(rclone cat $gff | grep CAU | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Thr=$(rclone cat $gff | grep UGU | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Pro=$(rclone cat $gff | grep -w UGG | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Lys=$(rclone cat $gff | grep UUU | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Asp=$(rclone cat $gff | grep GUC | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Glu=$(rclone cat $gff | grep UUC | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_SerGCU=$(rclone cat $gff | grep GCU | grep tRNA | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Tyr=$(rclone cat $gff | grep GUA | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Cys=$(rclone cat $gff | grep GCA | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Trp=$(rclone cat $gff | grep UCA | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Ala=$(rclone cat $gff | grep UGC | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Asn=$(rclone cat $gff | grep GUU | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_LeuUAA=$(rclone cat $gff | grep UAA | grep tRNA | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_SerUGA=$(rclone cat $gff | grep UGA | grep tRNA | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Gly=$(rclone cat $gff | grep UCC | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Arg=$(rclone cat $gff | grep UCG | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_His=$(rclone cat $gff | grep GUG | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')
        tRNA_Gln=$(rclone cat $gff | grep UUG | awk -F '\t' '$4 != "" && $5 != "" {len = $5 - $4 + 1; if (len > max) max = len} END {print max}')

        # Test to make sure that all the features are present
        # Define all variables to check
        variables=(
            "$rrna12s" "$rrna16s" "$atp6" "$atp8" "$cox1" "$cox2" "$cox3" "$cytb"
            "$nad1" "$nad2" "$nad3" "$nad4" "$nad4l" "$nad5" "$mad6"
            "$tRNA_Phe" "$tRNA_Val" "$tRNA_LeuUAG" "$tRNA_Ile" "$tRNA_Met" "$tRNA_Thr"
            "$tRNA_Pro" "$tRNA_Lys" "$tRNA_Asp" "$tRNA_Glu" "$tRNA_SerGCU" "$tRNA_Tyr"
            "$tRNA_Cys" "$tRNA_Trp" "$tRNA_Ala" "$tRNA_Asn" "$tRNA_LeuUAA" "$tRNA_SerUGA"
            "$tRNA_Gly" "$tRNA_Arg" "$tRNA_His" "$tRNA_Gln"
        )

        # Initialize status
        status="all"

        # Check if any variable is missing
        for var in "${variables[@]}"; do
            if [ -z "$var" ]; then
                status="missing"
                break
            fi
        done
            
        # Split up the assmbly name for the database
        key=$(echo $(basename $og) | tr '.' '\t')
        echo -e "$key\t$annotation\t$stats\t$length\t$length_emma\t$seq_12s\t$seq_16s\t$seq_CO1\t$cds_no\t$trna_no\t$rrna_no\t$status\t$genbank\t$rrna12s\t$rrna16s\t$atp6\t$atp8\t$cox1\t$cox2\t$cox3\t$cytb\t$nad1\t$nad2\t$nad3\t$nad4\t$nad4l\t$nad5\t$mad6\t$tRNA_Phe\t$tRNA_Val\t$tRNA_LeuUAG\t$tRNA_LeuUAA\t$tRNA_Ile\t$tRNA_Met\t$tRNA_Thr\t$tRNA_Pro\t$tRNA_Lys\t$tRNA_Asp\t$tRNA_Glu\t$tRNA_SerGCU\t$tRNA_SerUGA\t$tRNA_Tyr\t$tRNA_Cys\t$tRNA_Trp\t$tRNA_Ala\t$tRNA_Asn\t$tRNA_Gly\t$tRNA_Arg\t$tRNA_His\t$tRNA_Gln" >> mtdnastat."$(date '+%y%m%d')".tsv
    done

done




"""
}