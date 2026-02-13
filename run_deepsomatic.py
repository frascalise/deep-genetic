import argparse
import subprocess
import yaml
from pathlib import Path


# Load parameters from YAML
with open("params.yaml", "r") as f:
    params = yaml.safe_load(f)

TUMOR_BAM = f"/input/{params['tumor_bam']}" # File with tumor reads (with docker path)
NORMAL_BAM = f"/input/{params['normal_bam']}" # File with normal reads (with docker path)
REFERENCE = f"/input/{params['reference']}" # Reference genome (with docker path)
OUTPUT_VCF = f"/output/{params['output_vcf']}" # Output VCF file (with docker path)
OUTPUT_GVCF = f"/output/{params['output_gvcf']}" if params.get('output_gvcf') else None # Output gVCF file (with docker path)
NUM_SHARDS = params["num_shards"] # Number of shards (the work is divided in N parts)
REGIONS = params.get("regions") # Regions to analyze (optional)
MODEL_TYPE = "WGS" # Model type (WGS for genome or WES for exome)

verbose = False # Can be set to True with the --verbose flag


def run_deepsomatic(): 
    cmd = [
        "docker", "run", 
        "--rm",
        "--gpus", "all",
        "-v", f"{Path(__file__).parent}:/input",
        "-v", f"{Path(__file__).parent / params['output_dir']}:/output",
        "google/deepsomatic:1.8.0",
        "run_deepsomatic",
        "--model_type", MODEL_TYPE,
        "--ref", REFERENCE,
        "--reads_tumor", TUMOR_BAM,
        "--reads_normal", NORMAL_BAM,
        "--output_vcf", OUTPUT_VCF,
        "--num_shards", str(NUM_SHARDS)
    ]
    if REGIONS:
        cmd.extend(["--regions", REGIONS])
    if OUTPUT_GVCF:
        cmd.extend(["--output_gvcf", OUTPUT_GVCF])
    
    if verbose:
        subprocess.run(cmd, check=True)
    else:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Deepsomatic runner")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()

    print("Tumor BAM: \t", TUMOR_BAM)
    print("Normal BAM: \t", NORMAL_BAM)
    print("Reference: \t", REFERENCE)
    print("Output VCF: \t", OUTPUT_VCF)
    print("Output GVCF: \t", OUTPUT_GVCF)
    print("Num shards: \t", NUM_SHARDS)
    print("Regions: \t", REGIONS)
    print("Model type: \t", MODEL_TYPE)
    input("Press Enter to continue... (may take a while)")
    print("\nDeepSomatic is running...\n")

    run_deepsomatic()