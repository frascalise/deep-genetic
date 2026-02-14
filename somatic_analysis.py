import pandas as pd
import gzip
import yaml
import os
import io
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np


DEEPSOMATIC_OUTPUT = None # Contains the path to the output VCF file
DEEPSOMATIC_OUTPUT_GVCF = None # Contains the path to the output GVCF file
GROUND_TRUTH_VCF = None # Contains the path to the ground truth VCF file


def setup_path():
    with open("params.yaml", "r") as f:
        params = yaml.safe_load(f)

    global DEEPSOMATIC_OUTPUT, DEEPSOMATIC_OUTPUT_GVCF, GROUND_TRUTH_VCF
    DEEPSOMATIC_OUTPUT = params["output_vcf"]
    DEEPSOMATIC_OUTPUT_GVCF = params['output_gvcf'] if params.get('output_gvcf') else None
    GROUND_TRUTH_VCF = params["ground_truth_vcf"]

    # Check if the files exists
    if not os.path.exists(DEEPSOMATIC_OUTPUT):
        print("ERROR: FILE", params["output_vcf"], "NOT FOUND")
        exit()

    if DEEPSOMATIC_OUTPUT_GVCF and not os.path.exists(DEEPSOMATIC_OUTPUT_GVCF):
        print("ERROR: FILE", params["output_gvcf"], "NOT FOUND")
        exit()

    if not os.path.exists(GROUND_TRUTH_VCF):
        print("ERROR: FILE", params["ground_truth_vcf"], "NOT FOUND")
        exit()
    

def read_vcf(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            lines = [l for l in f if not l.startswith('##')]
    else:
        with open(path, "r") as f:
            lines = [l for l in f if not l.startswith('##')]
    
    # lines = [                         lines is a list of strings
    # "#CHROM\tPOS\tID\tREF...\n",      Header
    # "chr1\t12345\t.\tA\tT...\n",      Data row 1
    # "chr1\t12346\t.\tG\tC...\n"       Data row 2
    # ]
    # vcf_string is a giant string that contains all the row of lines
    # pd.read_csv(io.StringIO(vcf_string), sep="\t") reads the VCF file and returns a DataFrame
    vcf_string = ''.join(lines)
    df = pd.read_csv(io.StringIO(vcf_string), sep="\t")

    # The result is a DataFrame like this:
    # CHROM      POS ID REF ALT  QUAL   FILTER INFO              FORMAT                               TUMOR
    # 0  chr17  7671133  .   A   C   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL    0/0:32:98:95,3:0.0306122:0,32,41
    # 1  chr17  7671570  .   G   T   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL  0/0:49:122:118,4:0.0327869:0,52,51
    # 2  chr17  7671587  .   A   T   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL  0/0:46:136:129,6:0.0441176:0,46,52
    # 3  chr17  7671613  .   A   C   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL  0/0:43:122:115,6:0.0491803:0,52,43
    # 4  chr17  7671630  .   A   C   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL  0/0:47:128:119,5:0.0390625:0,48,54
    # print(df.head(10), '\n\n')
    
    return df


def deepsomatic_scores(vcf_df, gt_df):
    GT_POS = gt_df['POS']

    for i in GT_POS:
        if vcf_df.loc[vcf_df['POS'] == i].empty:
            print("Record not found in VCF")
            continue
        # Check if the REF, ALT, POS and CHROM are the same between VCF and GT
        same_chrom = vcf_df.loc[vcf_df['POS'] == i, '#CHROM'].values[0] == gt_df.loc[gt_df['POS'] == i, '#CHROM'].values[0]
        same_pos   = vcf_df.loc[vcf_df['POS'] == i, 'POS'].values[0]   == gt_df.loc[gt_df['POS'] == i, 'POS'].values[0]
        same_ref   = vcf_df.loc[vcf_df['POS'] == i, 'REF'].values[0]   == gt_df.loc[gt_df['POS'] == i, 'REF'].values[0]
        same_alt   = vcf_df.loc[vcf_df['POS'] == i, 'ALT'].values[0]   == gt_df.loc[gt_df['POS'] == i, 'ALT'].values[0]
        
        if same_ref and same_alt and same_pos and same_chrom:
            print("VARIANT NUMBER: ", gt_df.loc[gt_df['POS'] == i, 'ID'].index[0])
            print("VCF CHROM: ", vcf_df.loc[vcf_df['POS'] == i, '#CHROM'].values[0], "\t| GTCHROM: ", gt_df.loc[gt_df['POS'] == i, '#CHROM'].values[0])
            print("VCF POS: ", vcf_df.loc[vcf_df['POS'] == i, 'POS'].values[0], "\t| GT POS: ", gt_df.loc[gt_df['POS'] == i, 'POS'].values[0])
            print("VCF REF: ", vcf_df.loc[vcf_df['POS'] == i, 'REF'].values[0], "\t\t| GT REF: ", gt_df.loc[gt_df['POS'] == i, 'REF'].values[0])
            print("VCF ALT: ", vcf_df.loc[vcf_df['POS'] == i, 'ALT'].values[0], "\t\t| GT ALT: ", gt_df.loc[gt_df['POS'] == i, 'ALT'].values[0])
            print("\n")
        else:
            print("No match")   


def visualize_vcf_data(df_input):
    pass


if __name__ == "__main__":
    setup_path()

    #print("DEEPSOMATIC_OUTPUT: ", DEEPSOMATIC_OUTPUT)
    #print("DEEPSOMATIC_OUTPUT_GVCF: ", DEEPSOMATIC_OUTPUT_GVCF)

    # Read the output VCF file
    vcf_df = read_vcf(DEEPSOMATIC_OUTPUT)
    # Read the ground truth VCF file
    gt_df = read_vcf(GROUND_TRUTH_VCF)

    #print(vcf_df[vcf_df['POS'] == 7675217])
    #print(vcf_df.head())
    #print(gt_df.head())
    deepsomatic_scores(vcf_df, gt_df)

    # visualize_vcf_data(vcf_df)
    # visualize_vcf_data(gt_df)