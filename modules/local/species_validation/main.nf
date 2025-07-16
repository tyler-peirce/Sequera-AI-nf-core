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
    need to take the 3rd tab seperated section of LCA file and then use that to compare to the NomID

    If there isnt a match then check the filtered results

    give a table that is gene \t nom_id \t LCA \t nom_id in blast results
    
    
import psycopg2
import csv

# --- CONFIGURATION ---
DB_CONFIG = {
    'host': 'your_host',
    'port': 5432,
    'dbname': 'your_db',
    'user': 'your_user',
    'password': 'your_password'
}

def concatenate_files(file_list, output_file):
    with open(output_file, 'w') as outfile:
        for filename in file_list:
            with open(filename, 'r') as infile:
                for line in infile:
                    outfile.write(line)

lca_files = ["lca_part1.tsv", "lca_part2.tsv", "lca_part3.tsv"]
blast_files = ["blast1.txt", "blast2.txt", "blast3.txt"]

concatenate_files(lca_files, "lca_combined.tsv")
concatenate_files(blast_files, "blast_combined.txt")

INPUT_TSV = "lca_combined.tsv"
BLAST_FILE = "blast_combined.txt"

# --- HELPER FUNCTIONS ---
def normalise_name(name):
    return name.strip().lower().replace('_', ' ')

def get_species_map(og_id):
    query = """
        SELECT og_id, nominal_species_id
        FROM sample s
        WHERE s.og_id = ANY(%s)
    """
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        with conn.cursor() as cur:
            cur.execute(query, (og_id,))
            return {og_id: species for og_id, species in cur.fetchall()}
    except Exception as e:
        print(f"Database error: {e}")
        return {}
    finally:
        if conn:
            conn.close()

def load_blast_species_set(filepath):
    with open(filepath, 'r') as f:
        text = f.read().lower()
    return text

# --- MAIN SCRIPT ---
def compare_lca_and_blast():
    input_rows = []
    og_id = $meta.id

    # Load input TSV
    with open(INPUT_TSV, newline='') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            if len(row) < 3:
                continue
            og_id = row[0].strip()
            lca = row[2].strip()
            og_id.add(og_id)
            input_rows.append((og_id, lca))

    # Get DB species map
    species_map = get_species_map(og_id)

    # Load BLAST results as big lowercase blob
    blast_blob = load_blast_species_set(BLAST_FILE)

    # Compare and write results
    with open(OUTPUT_TSV, "w", newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["og_id", "LCA_result", "DB_species", "Match_YN", "Found_in_blast_YN"])

        for og_id, lca in input_rows:
            db_species = species_map.get(og_id)
            if db_species is None:
                continue

            match = "Yes" if normalise_name(lca) == normalise_name(db_species) else "No"
            in_blast = "Yes" if normalise_name(db_species) in blast_blob else "No"

            writer.writerow([og_id, lca, db_species, match, in_blast])

    print(f"Done. Output written to {OUTPUT_TSV}")

# --- RUN ---
if __name__ == "__main__":
    compare_lca_and_blast()

"""
}