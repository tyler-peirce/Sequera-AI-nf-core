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
    import psycopg2
import pandas as pd
import os
import glob
import numpy as np  # Required for handling infinity values

# run using: singularity run $SING/psycopg2:0.1.sif python 05_push_mitogenome_results_to_sqldb.py
# PostgreSQL connection parameters
db_params = {
    'dbname': 'oceanomics',
    'user': 'postgres',
    'password': 'oceanomics',
    'host': '203.101.227.69',
    'port': 5432
}


# Access file paths


# File containing mitogenome data
mito_path = "/scratch/pawsey0964/tpeirce/_MITOGENOMES/OceanOmics-Mitogenome-Nextflow/mtdnastat.250714.tsv"  # Update this with the correct file path

# Import Mitogenome Data
print(f"Importing data from {mito_path}")

# Load and preprocess mito data
mito = pd.read_csv(mito_path, sep="\t")

# Replace NaN with None
mito = mito.replace({np.nan: None})

# Normalize column names (remove spaces & make lowercase)
mito.columns = mito.columns.str.strip().str.lower()

try:
    # Connect to PostgreSQL
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    row_count = 0  # Track number of processed rows

    for index, row in mito.iterrows():
        row_dict = row.to_dict()
        # Extract primary key values
        og_id, tech, seq_date, code = row['og_id'], row['tech'], row['seq_date'], row['code']

        # UPSERT: Insert if not exists, otherwise update the row
        upsert_query = """
       INSERT INTO mitogenome_data (
            og_id, tech, seq_date, code, annotation, stats, length, length_emma, seqlength_12s,
            seqlength_16s, seqlength_co1, cds_no, trna_no, rrna_no, status, genbank, rrna12s,
            rrna16s, atp6, atp8, cox1, cox2, cox3, cytb, nad1, nad2, nad3, nad4, nad4l, 
            nad5, mad6, trna_phe, trna_val, trna_leuuag, trna_leuuaa, trna_ile, trna_met, 
            trna_thr, trna_pro, trna_lys, trna_asp, trna_glu, trna_sergcu, trna_seruga, 
            trna_tyr, trna_cys, trna_trp, trna_ala, trna_asn, trna_gly, trna_arg, trna_his, 
            trna_gln
        )
        VALUES (
            %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s, %(stats)s, %(length)s, %(length_emma)s, %(seqlength_12s)s,
            %(seqlength_16s)s, %(seqlength_co1)s, %(cds_no)s, %(trna_no)s, %(rrna_no)s, %(status)s, %(genbank)s, %(rrna12s)s,
            %(rrna16s)s, %(atp6)s, %(atp8)s, %(cox1)s, %(cox2)s, %(cox3)s, %(cytb)s, %(nad1)s, %(nad2)s, %(nad3)s, %(nad4)s, %(nad4l)s,
            %(nad5)s, %(mad6)s, %(trna_phe)s, %(trna_val)s, %(trna_leuuag)s, %(trna_leuuaa)s, %(trna_ile)s, %(trna_met)s,
            %(trna_thr)s, %(trna_pro)s, %(trna_lys)s, %(trna_asp)s, %(trna_glu)s, %(trna_sergcu)s, %(trna_seruga)s,
            %(trna_tyr)s, %(trna_cys)s, %(trna_trp)s, %(trna_ala)s, %(trna_asn)s, %(trna_gly)s, %(trna_arg)s, %(trna_his)s,
            %(trna_gln)s
        )
        ON CONFLICT (og_id, tech, seq_date, code) 
        DO UPDATE SET
            annotation = EXCLUDED.annotation,
            stats = EXCLUDED.stats,
            length = EXCLUDED.length,
            length_emma = EXCLUDED.length_emma,
            seqlength_12s = EXCLUDED.seqlength_12s,
            seqlength_16s = EXCLUDED.seqlength_16s,
            seqlength_co1 = EXCLUDED.seqlength_co1,
            cds_no = EXCLUDED.cds_no,
            trna_no = EXCLUDED.trna_no,
            rrna_no = EXCLUDED.rrna_no,
            status = EXCLUDED.status,
            genbank = EXCLUDED.genbank,
            rrna12s = EXCLUDED.rrna12s,
            rrna16s = EXCLUDED.rrna16s,
            atp6 = EXCLUDED.atp6,
            atp8 = EXCLUDED.atp8,
            cox1 = EXCLUDED.cox1,
            cox2 = EXCLUDED.cox2,
            cox3 = EXCLUDED.cox3,
            cytb = EXCLUDED.cytb,
            nad1 = EXCLUDED.nad1,
            nad2 = EXCLUDED.nad2,
            nad3 = EXCLUDED.nad3,
            nad4 = EXCLUDED.nad4,
            nad4l = EXCLUDED.nad4l,
            nad5 = EXCLUDED.nad5,
            mad6 = EXCLUDED.mad6,
            trna_phe = EXCLUDED.trna_phe,
            trna_val = EXCLUDED.trna_val,
            trna_leuuag = EXCLUDED.trna_leuuag,
            trna_leuuaa = EXCLUDED.trna_leuuaa,
            trna_ile = EXCLUDED.trna_ile,
            trna_met = EXCLUDED.trna_met,
            trna_thr = EXCLUDED.trna_thr,
            trna_pro = EXCLUDED.trna_pro,
            trna_lys = EXCLUDED.trna_lys,
            trna_asp = EXCLUDED.trna_asp,
            trna_glu = EXCLUDED.trna_glu,
            trna_sergcu = EXCLUDED.trna_sergcu,
            trna_seruga = EXCLUDED.trna_seruga,
            trna_tyr = EXCLUDED.trna_tyr,
            trna_cys = EXCLUDED.trna_cys,
            trna_trp = EXCLUDED.trna_trp,
            trna_ala = EXCLUDED.trna_ala,
            trna_asn = EXCLUDED.trna_asn,
            trna_gly = EXCLUDED.trna_gly,
            trna_arg = EXCLUDED.trna_arg,
            trna_his = EXCLUDED.trna_his,
            trna_gln = EXCLUDED.trna_gln;
        """

        params = {
            "og_id": row_dict["og_id"],
            "tech": row_dict["tech"],
            "seq_date": row_dict["seq_date"],
            "code": row_dict["code"],
            "annotation": row_dict["annotation"],
            "stats": row_dict.get("stats"),
            "length": row_dict.get("length"),
            "length_emma": row_dict.get("length_emma"),
            "seqlength_12s": row_dict.get("seqlength_12s"),
            "seqlength_16s": row_dict.get("seqlength_16s"),
            "seqlength_co1": row_dict.get("seqlength_co1"),
            "cds_no": row_dict.get("cds_no"),
            "trna_no": row_dict.get("trna_no"),
            "rrna_no": row_dict.get("rrna_no"),
            "status": row_dict.get("status"),
            "genbank": row_dict.get("genbank"),
            "rrna12s": row_dict.get("rrna12s"),
            "rrna16s": row_dict.get("rrna16s"),
            "atp6": row_dict.get("atp6"),
            "atp8": row_dict.get("atp8"),
            "cox1": row_dict.get("cox1"),
            "cox2": row_dict.get("cox2"),
            "cox3": row_dict.get("cox3"),
            "cytb": row_dict.get("cytb"),
            "nad1": row_dict.get("nad1"),
            "nad2": row_dict.get("nad2"),
            "nad3": row_dict.get("nad3"),
            "nad4": row_dict.get("nad4"),
            "nad4l": row_dict.get("nad4l"),
            "nad5": row_dict.get("nad5"),
            "mad6": row_dict.get("mad6"),
            "trna_phe": row_dict.get("trna_phe"),
            "trna_val": row_dict.get("trna_val"),
            "trna_leuuag": row_dict.get("trna_leuuag"),
            "trna_leuuaa": row_dict.get("trna_leuuaa"),
            "trna_ile": row_dict.get("trna_ile"),
            "trna_met": row_dict.get("trna_met"),
            "trna_thr": row_dict.get("trna_thr"),
            "trna_pro": row_dict.get("trna_pro"),
            "trna_lys": row_dict.get("trna_lys"),
            "trna_asp": row_dict.get("trna_asp"),
            "trna_glu": row_dict.get("trna_glu"),
            "trna_sergcu": row_dict.get("trna_sergcu"),
            "trna_seruga": row_dict.get("trna_seruga"),
            "trna_tyr": row_dict.get("trna_tyr"),
            "trna_cys": row_dict.get("trna_cys"),
            "trna_trp": row_dict.get("trna_trp"),
            "trna_ala": row_dict.get("trna_ala"),
            "trna_asn": row_dict.get("trna_asn"),
            "trna_gly": row_dict.get("trna_gly"),
            "trna_arg": row_dict.get("trna_arg"),
            "trna_his": row_dict.get("trna_his"),
            "trna_gln": row_dict.get("trna_gln"),
        }


        # Debugging Check
        print(f"Number of rows being passed: {len(row)}")
        print("row:", row_dict)
        print("params:", params)

        cursor.execute(upsert_query, params)
        row_count += 1  

        conn.commit()
        print(f"✅ Successfully processed {row_count} rows!")

except Exception as e:
    conn.rollback()
    print(f"❌ Error: {e}")

finally:
    cursor.close()
    conn.close()


"""
}