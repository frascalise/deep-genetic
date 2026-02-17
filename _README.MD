# Deep-Genetic: Somatic Variant Calling using DeepSomatic

![Project Status](https://img.shields.io/badge/Status-Development-yellow) ![Python](https://img.shields.io/badge/Python-3.8%2B-blue) ![Docker](https://img.shields.io/badge/Docker-Required-blue)

A pipeline for somatic variant calling using Google's DeepSomatic deep learning model, with comprehensive analysis and comparison tools.

## Overview

This project provides a streamlined workflow for:
- Running DeepSomatic for tumor-normal somatic variant calling
- Analyzing variant calling results against ground truth data
- Computing precision, recall, and F1 scores with VAF filtering
- Visualizing variant calling performance

## Project Structure

```
deep-genetic/
├── run_deepsomatic.py      # Main script to run DeepSomatic via Docker
├── somatic_analysis.py     # Analyze and compare VCF outputs with ground truth
├── create_index.py         # Utility to create BAM/FASTA index files
├── params.yaml             # Configuration file for all pipeline parameters
├── requirements.txt        # Python dependencies
├── deepsomatic_data/       # Input data (reference genome, ground truth VCF)
├── downsampled/            # Downsampled BAM files (tumor/normal)
└── deepsomatic_output/     # Output directory for VCF results
```

## Prerequisites

- **Docker**: Required to run DeepSomatic
- **Python 3.8+**: For running analysis scripts
- **Input Data**:
  - Tumor BAM file with index (.bam + .bai)
  - Normal BAM file with index (.bam + .bai)
  - Reference genome FASTA with index (.fa + .fai)
  - Ground truth VCF (optional, for validation)

## Installation

1. **Clone the repository**:
```bash
git clone <repository-url>
cd deep-genetic
```

2. **Install Python dependencies**:
```bash
pip install -r requirements.txt
```

3. **Verify Docker installation**:
```bash
docker --version
```

## Configuration

Edit `params.yaml` to configure your analysis:

```yaml
# Output directory
output_dir: "deepsomatic_output/"

# Input files (relative paths from project root)
tumor_bam: "downsampled/tumor.bam"
normal_bam: "downsampled/normal.bam"
reference: "deepsomatic_data/reference.fa"
ground_truth_vcf: "downsampled/truth_tumor_only.vcf.gz"

# Optional parameters
num_shards: 1                        # Parallel processing shards
regions: "chr17:5000000-8000000"     # Specific genomic region (optional)
batch_size: 256                      # Batch size for variant calling
```

## Usage

### 1. Run DeepSomatic Variant Calling

```bash
python3 run_deepsomatic.py
```

This script will:
- Load configuration from `params.yaml`
- Run DeepSomatic via Docker with the specified parameters
- Generate output VCF and gVCF files in the output directory

**Output files**:
- `deepsomatic_output/output.vcf.gz` - Variant calls
- `deepsomatic_output/output.g.vcf.gz` - Genomic VCF with all sites

### 2. Analyze Results

```bash
python3 somatic_analysis.py
```

This script will:
- Read the DeepSomatic output VCF
- Compare against ground truth VCF
- Filter variants by VAF (Variant Allele Frequency) threshold
- Calculate and display:
  - True Positives (TP)
  - False Positives (FP)
  - False Negatives (FN)
  - Precision, Recall, F1 Score

**Configuration in script**:
- `VAF_THRESHOLD`: Default 0.15 (15% minimum variant frequency)
- `CHROMOSOMES`: List of chromosomes to analyze
- `MIN_REGION` / `MAX_REGION`: Region boundaries for filtering

### 3. Create Index Files (if needed)

```bash
python3 create_index.py
```

Use this utility to generate index files for BAM or FASTA files.

## Docker Configuration

The pipeline uses the official DeepSomatic Docker image:
- **Image**: `google/deepsomatic:1.8.0`
- **Model**: WGS (Whole Genome Sequencing)
- **GPU Support**: Currently disabled (uncomment `--gpus all` in `run_deepsomatic.py` to enable)

## Output Interpretation

### VCF Filters
- **PASS**: High-quality somatic variant
- **RefCall**: Reference call (filtered out)
- **GERMLINE**: Germline variant (filtered out)

### Analysis Metrics

The analysis provides two sets of metrics:

1. **Unfiltered**: All PASS variants
2. **VAF Filtered**: Only variants with VAF ≥ threshold (default 0.15)

**Metrics**:
- **Precision**: TP / (TP + FP) - Proportion of called variants that are true
- **Recall**: TP / (TP + FN) - Proportion of true variants that were called
- **F1 Score**: Harmonic mean of precision and recall

## Troubleshooting

### Common Issues

1. **Missing index files**: Run `create_index.py` to generate .bai or .fai files
2. **Docker permission errors**: Add user to docker group or run with sudo
3. **Segmentation faults**: Reduce `num_shards` to 1 for small regions
4. **GPU not detected**: Ensure NVIDIA Docker runtime is installed

### File Paths

- All paths in `params.yaml` are relative to the project root
- Docker mounts the project directory as `/input` and output directory as `/output`
- Ensure all input files exist before running

## Data Requirements

### Test Data Structure

```
deepsomatic_data/
├── reference.fa          # Reference genome
├── reference.fa.fai      # FASTA index
└── truth_big.vcf.gz      # Ground truth variants

downsampled/
├── tumor.bam             # Tumor sample
├── tumor.bam.bai         # BAM index
├── normal.bam            # Normal sample
├── normal.bam.bai        # BAM index
└── truth_tumor_only.vcf.gz  # Truth set for region
```

## Advanced Configuration

### Region-Specific Analysis

To analyze specific genomic regions, set in `params.yaml`:
```yaml
regions: "chr17:7670000-7680000"
```

### Performance Tuning

- **num_shards**: Increase for larger datasets (e.g., 8-16 for WGS)
- **batch_size**: Adjust based on available memory (default: 256)

## References

- [DeepSomatic GitHub](https://github.com/google/deepsomatic)
- [DeepSomatic Paper](https://www.nature.com/articles/s41587-023-01681-4)

## License

This project is for research purposes. DeepSomatic is licensed under the BSD 3-Clause License.
