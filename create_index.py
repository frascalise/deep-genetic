import os
import subprocess

def index_bams():
    cwd = os.getcwd()
    print(f"Creating index for BAM files in {cwd}/downsampled")

    # Using staphb/samtools image to create the index
    docker_cmd = f"docker run --rm -v \"{cwd}\":/data staphb/samtools:latest"

    # Index Tumor
    print("  -> Indexing Tumor...")
    subprocess.run(f'{docker_cmd} samtools index /data/downsampled/tumor.bam', shell=True, check=True)

    # Index Normal
    print("  -> Indexing Normal...")
    subprocess.run(f'{docker_cmd} samtools index /data/downsampled/normal.bam', shell=True, check=True)

    print("Index created, you can now run DeepSomatic.")

if __name__ == "__main__":
    index_bams()