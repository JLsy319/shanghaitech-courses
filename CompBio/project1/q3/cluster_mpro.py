import subprocess
import os

def run_cdhit_clustering(input_fasta="q2_mpro_sequences.fasta", identity_threshold=0.85):
    """
    Run CD-HIT clustering on protein sequences.
    
    Args:
        input_fasta: Path to input FASTA file (default: q2_mpro_sequences.fasta)
        identity_threshold: Sequence identity threshold (default: 0.85 or 85%)
    Returns:
        Path to the output cluster file
    """
    output_prefix = "temp_clustering"
    output_clstr = "q3_mpro_clusters.clstr"
    # Run CD-HIT
    cmd = f"cd-hit -i {input_fasta} -o {output_prefix} -c {identity_threshold} -n 5 -M 16000 -d 0"
    print(f"Running CD-HIT with command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    
    os.rename(f"{output_prefix}.clstr", output_clstr)
    
    # Clean up temporary files
    if os.path.exists(output_prefix):
        os.remove(output_prefix)
    
    return output_clstr

def count_clusters(cluster_file):
    """
    Count the number of clusters in CD-HIT output file.
    
    Args:
        cluster_file: Path to CD-HIT .clstr file
    Returns:
        Number of clusters
    """
    cluster_count = 0
    with open(cluster_file) as f:
        for line in f:
            if line.startswith('>Cluster'):
                cluster_count += 1
    return cluster_count

def main():
    print("Starting CD-HIT clustering...")
    cluster_file = run_cdhit_clustering()
    num_clusters = count_clusters(cluster_file)
    
    print(f"\nClustering Results Summary:")
    print(f"Total number of clusters: {num_clusters}")
    print(f"Results saved to: {cluster_file}")

if __name__ == "__main__":
    main()