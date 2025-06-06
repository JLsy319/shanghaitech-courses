from Bio import SeqIO

ALIGNED_FASTA_FILE = "aligned_mpro.fasta"
REFERENCE_ID = "6Y2F_A"
RESIDUE_POSITIONS = [26, 41, 49, 140, 141, 142, 143, 145, 163, 164, 165, 166, 167, 168, 172, 187, 189]
OUTPUT_FASTA_FILE = "pocket_signatures.fasta"

# Find reference sequence and determine column indices
column_indices = []
with open(ALIGNED_FASTA_FILE, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id.startswith(REFERENCE_ID):
            ref_seq = str(record.seq)
            res_count = 0
            
            for i, char in enumerate(ref_seq):
                if char != '-':  # if not a gap, increment residue counter
                    res_count += 1
                if res_count in RESIDUE_POSITIONS:
                    column_indices.append(i)
                    RESIDUE_POSITIONS.remove(res_count)
            
            column_indices.sort()  # sort indices for sequential extraction
            break

# Extract pocket signature sequences based on column indices
pocket_signatures = []
with open(ALIGNED_FASTA_FILE, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        original_seq = str(record.seq)
        pocket_seq = "".join([original_seq[i] for i in column_indices])
        record.seq = type(record.seq)(pocket_seq)
        pocket_signatures.append(record)

# Write new sequences to output file
with open(OUTPUT_FASTA_FILE, "w") as output_handle:
    SeqIO.write(pocket_signatures, output_handle, "fasta")

print(f"Created pocket signature file: {OUTPUT_FASTA_FILE}")
print(f"New sequence length: {len(column_indices)}")