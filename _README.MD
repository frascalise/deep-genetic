# Somatic Variant Calling using Deep Learning (WIP)

![Project Status](https://img.shields.io/badge/Status-Development-yellow) ![Python](https://img.shields.io/badge/Python-3.8%2B-blue) ![Framework](https://img.shields.io/badge/Framework-PyTorch%20%7C%20Scikit--Learn-orange)

## Test Data (`deepsomatic_test/`)

### Tumor Sample
- **`NA12878_S1.chr20.10_10p1mb.bam`** - Tumor sample BAM file (chromosome 20, 10-10.1Mb region)
- **`NA12878_S1.chr20.10_10p1mb.bam.bai`** - BAM index file

### Healthy/Normal Sample
- ⚠️ **Not currently included** - DeepSomatic is running in tumor-only mode

### Reference Genome
- **`ucsc.hg19.chr20.unittest.fasta`** - Reference genome (hg19, chromosome 20)
- **`ucsc.hg19.chr20.unittest.fasta.fai`** - FASTA index
- **`ucsc.hg19.chr20.unittest.fasta.gz`** - Compressed reference
- **`ucsc.hg19.chr20.unittest.fasta.gz.fai`** - Compressed FASTA index
- **`ucsc.hg19.chr20.unittest.fasta.gz.gzi`** - GZI index

### Validation Data
- **`test_nist.b37_chr20_100kbp_at_10mb.vcf.gz`** - NIST truth set VCF
- **`test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi`** - VCF index
- **`test_nist.b37_chr20_100kbp_at_10mb.bam.bai`** - BAM index for validation

## Usage

### 1. Run DeepSomatic

```bash
python3 run_deepsomatic.py
```

Optional flags:
- `--verbose` - Enable verbose output and automatic visualization

### 2. Compare with Ground Truth

After running DeepSomatic, compare the output VCF with the NIST ground truth:

```bash
python3 compare_vcf.py \
  --truth deepsomatic_test/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz \
  --query deepsomatic_output/output.vcf.gz
```

This will generate:
- **Console output**: Precision, Recall, F1 scores
- **`comparison_results.csv`**: Detailed variant-by-variant results
- **`comparison_report.png`**: Visualization of metrics and variant counts

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```
