import pandas as pd
import yaml
import gzip
import os
import io
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np


DEEPSOMATIC_OUTPUT = None # Contains the path to the output VCF file
DEEPSOMATIC_OUTPUT_GVCF = None # Contains the path to the output GVCF file
GROUND_TRUTH_VCF = None # Contains the path to the ground truth VCF file
LIMIT = False # Limit the analysis to a specific region
MIN_REGION = 7670000 # Contains the minimum region size (found in the yaml file in the regions section)
MAX_REGION = 7680000 # Contains the maximum region size (found in the yaml file in the regions section)


def setup_path():
    global DEEPSOMATIC_OUTPUT, DEEPSOMATIC_OUTPUT_GVCF, GROUND_TRUTH_VCF
    DEEPSOMATIC_OUTPUT = "deepsomatic_output/output.vcf.gz"
    DEEPSOMATIC_OUTPUT_GVCF = "deepsomatic_output/output.vcf.gz.tbi"
    GROUND_TRUTH_VCF = "deepsomatic_data/truth_big.vcf.gz"

    # Check if the files exists
    if not os.path.exists(DEEPSOMATIC_OUTPUT):
        print("ERROR: FILE", DEEPSOMATIC_OUTPUT, "NOT FOUND")
        exit()

    if DEEPSOMATIC_OUTPUT_GVCF and not os.path.exists(DEEPSOMATIC_OUTPUT_GVCF):
        print("ERROR: FILE", DEEPSOMATIC_OUTPUT_GVCF, "NOT FOUND")
        exit()

    if not os.path.exists(GROUND_TRUTH_VCF):
        print("ERROR: FILE", GROUND_TRUTH_VCF, "NOT FOUND")
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
        if LIMIT:
            if i < MIN_REGION or i > MAX_REGION:
                continue
        if vcf_df.loc[vcf_df['POS'] == i].empty:
            print("Record from GT at position", i, "not found in TUMOR VCF")
            continue
        
        # Check if the REF, ALT, POS and CHROM are the same between VCF and GT
        same_chrom = vcf_df.loc[vcf_df['POS'] == i, '#CHROM'].values[0] == gt_df.loc[gt_df['POS'] == i, '#CHROM'].values[0]
        same_pos   = vcf_df.loc[vcf_df['POS'] == i, 'POS'].values[0]   == gt_df.loc[gt_df['POS'] == i, 'POS'].values[0]
        same_ref   = vcf_df.loc[vcf_df['POS'] == i, 'REF'].values[0]   == gt_df.loc[gt_df['POS'] == i, 'REF'].values[0]
        same_alt   = vcf_df.loc[vcf_df['POS'] == i, 'ALT'].values[0]   == gt_df.loc[gt_df['POS'] == i, 'ALT'].values[0]
        
        if same_ref and same_alt and same_pos and same_chrom:
            variant_num = gt_df.loc[gt_df['POS'] == i, 'ID'].index[0]
            
            # Extract values
            vcf_chrom = str(vcf_df.loc[vcf_df['POS'] == i, '#CHROM'].values[0])
            vcf_pos = str(vcf_df.loc[vcf_df['POS'] == i, 'POS'].values[0])
            vcf_ref = str(vcf_df.loc[vcf_df['POS'] == i, 'REF'].values[0])
            vcf_alt = str(vcf_df.loc[vcf_df['POS'] == i, 'ALT'].values[0])
            
            gt_chrom = str(gt_df.loc[gt_df['POS'] == i, '#CHROM'].values[0])
            gt_pos = str(gt_df.loc[gt_df['POS'] == i, 'POS'].values[0])
            gt_ref = str(gt_df.loc[gt_df['POS'] == i, 'REF'].values[0])
            gt_alt = str(gt_df.loc[gt_df['POS'] == i, 'ALT'].values[0])
            
            # Custom table with Unicode box-drawing characters
            print(f"\n┌{'─'*70}┐")
            print(f"│ VARIANT NUMBER: {variant_num:<53}│")
            print(f"├{'─'*20}┬{'─'*24}┬{'─'*24}┤")
            print(f"│ {'Field':<18} │ {'VCF':<22} │ {'Ground Truth':<22} │")
            print(f"├{'─'*20}┼{'─'*24}┼{'─'*24}┤")
            print(f"│ {'CHROM':<18} │ {vcf_chrom:<22} │ {gt_chrom:<22} │")
            print(f"│ {'POS':<18} │ {vcf_pos:<22} │ {gt_pos:<22} │")
            print(f"│ {'REF':<18} │ {vcf_ref:<22} │ {gt_ref:<22} │")
            print(f"│ {'ALT':<18} │ {vcf_alt:<22} │ {gt_alt:<22} │")
            print(f"└{'─'*20}┴{'─'*24}┴{'─'*24}┘\n")
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