# Somatic Variant Calling using Deep Learning

![Project Status](https://img.shields.io/badge/Status-Development-yellow) ![Python](https://img.shields.io/badge/Python-3.8%2B-blue) ![Framework](https://img.shields.io/badge/Framework-PyTorch%20%7C%20Scikit--Learn-orange)

## 🧬 Project Overview
This project focuses on **Somatic Variant Calling**, the process of identifying genetic alterations that occur in somatic cells (tumors) but are not present in the germline (healthy tissue).

The goal is to develop an AI pipeline that processes DNA sequencing data (BAM files) to distinguish true somatic variants from sequencing errors and artifacts. The project implements two distinct approaches:
1.  **Feature-based MLP:** A Multi-Layer Perceptron using extracted genomic features.
2.  **Image-based CNN:** A Deep Learning model (CNN) that classifies variants based on **Pileup Images**, treating the genomic alignment as a visual recognition problem.

This repository contains the source code, documentation, and demo scripts required for the **Bioinformatics AI** exam.

## 🎯 Objectives
As defined in the project specifications:
1.  **Feature Extraction & MLP:** Create a Multi-Layer Perceptron to distinguish true variants from false positives using tabular features (e.g., mapping quality, depth, strand bias).
2.  **Pileup Image Generation:** Develop a script to convert genomic read alignments into RGB images (Pileups) representing Normal vs. Tumor samples.
3.  **Deep Learning Classification:** Train a Convolutional Neural Network (CNN) to classify variants using the generated pileup images.

## 📂 Dataset
The project relies on paired **Normal/Tumor** sequencing data.
* **Format:** BAM (Binary Alignment Map) for reads, VCF (Variant Call Format) for ground truth labels.
* **Source:** [Insert Source Here, e.g., SEQC2 Data / Synthetic Dream Challenge].
* **Genome Reference:** [e.g., hg19 or hg38].

> **Note:** Due to size constraints, raw BAM files are not included in this repository. Please refer to `data/download_instructions.txt` to access the dataset.

## 🛠️ Architecture & Methods

### 1. Preprocessing & Pileup Generation
We use `pysam` to parse BAM files. For every candidate variant:
* **Tabular Data:** We extract metrics like Read Depth (DP), Allele Frequency (VAF), and Base Quality.
* **Visual Data:** We generate a "Pileup Image" where:
    * **Rows:** Individual reads aligned to the reference.
    * **Columns:** Genomic positions (centered on the candidate variant).
    * **Channels:** Bases are color-coded (A=Red, C=Blue, G=Green, T=Yellow).

### 2. Models
* **Baseline (MLP):** A simple fully connected network trained on the tabular features.
* **DeepVariant-like (CNN):** A Convolutional Neural Network (e.g., ResNet or Custom CNN) that takes the (Height x Width x 3) pileup image as input and outputs the probability of the variant being true.

## 🚀 Installation

1.  Clone the repository:
    ```bash
    git clone [https://github.com/frascalise/deep-genetic.git](https://github.com/frascalise/deep-genetic.git)
    cd deep-genetic
    ```

2.  Create a virtual environment (recommended):
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate
    ```

3.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## 💻 Usage

### A. Run the Demo
To test the pipeline on a sample region (required for the exam):
```bash
python main.py --mode demo --bam_tumor data/sample_tumor.bam --bam_normal data/sample_normal.bam --pos chr22:1000500
```
