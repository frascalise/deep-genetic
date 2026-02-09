"""
Deepsomatic with Docker
BAM: Aligned reads of DNA sample
FASTA: Reference genome
VCF: Final output, only positions with variants
"""

import argparse
import subprocess
from pathlib import Path

def run_deepsomatic(verbose):
    base_dir = Path(__file__).parent
    test_dir = base_dir / "deepsomatic_test"
    output_dir = base_dir / "deepsomatic_output"
    output_dir.mkdir(exist_ok=True)
    
    # Docker mount point
    docker_mount = "/data"
    
    docker_ref = f"{docker_mount}/deepsomatic_test/GRCh38_no_alt_analysis_set.chr20.fa.gz"
    docker_tumor_bam = f"{docker_mount}/deepsomatic_test/somatic_quickstart_tumor.bam"
    docker_normal_bam = f"{docker_mount}/deepsomatic_test/somatic_quickstart_normal.bam"
    
    # Output files (Docker paths)
    docker_output_vcf = f"{docker_mount}/deepsomatic_output/output.vcf.gz"
    docker_output_gvcf = f"{docker_mount}/deepsomatic_output/output.g.vcf.gz"
    docker_logs_dir = f"{docker_mount}/deepsomatic_output/logs"
    docker_intermediate_dir = f"{docker_mount}/deepsomatic_output/intermediate"
    
    # Model and sample configuration
    model_type = "WGS" 
    sample_name_tumor = "TUMOR"
    sample_name_normal = "NORMAL"
    num_shards = "1"
    
    if verbose:
        print("Verbose mode enabled")
        print(">>> Setup completed")
        print(f">>> Normal BAM: {docker_normal_bam}")
        print(f">>> Tumor BAM:  {docker_tumor_bam}")
        print(">>> Running DeepSomatic (Paired Mode)")

    #** Run DeepSomatic
    # Note: Paths must be relative to the Docker mount point (/data)
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{base_dir}:{docker_mount}",
        "google/deepsomatic:1.8.0",
        "run_deepsomatic",
        "--model_type", model_type,
        "--ref", docker_ref,
        "--reads_tumor", docker_tumor_bam,
        "--reads_normal", docker_normal_bam,
        "--output_vcf", docker_output_vcf,
        "--output_gvcf", docker_output_gvcf,
        "--sample_name_tumor", sample_name_tumor,
        "--sample_name_normal", sample_name_normal,
        "--num_shards", num_shards,
        "--logging_dir", docker_logs_dir,
        "--intermediate_results_dir", docker_intermediate_dir
    ]
    
    #* Execute the command in the terminal
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if verbose:
        print(">>> DeepSomatic completed")
        print(">>> Visualizing results")
        subprocess.run(["python", "visualize_variants.py"], check=True)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Deepsomatic runner")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()
    run_deepsomatic(args.verbose)
        