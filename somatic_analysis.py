import pandas as pd
import yaml
import gzip
import os
import io
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)


DEEPSOMATIC_OUTPUT = None # Contains the path to the output VCF file
DEEPSOMATIC_OUTPUT_GVCF = None # Contains the path to the output GVCF file
GROUND_TRUTH_VCF = None # Contains the path to the ground truth VCF file
MIN_REGION = 5000000 # Contains the minimum region size (found in the yaml file in the regions section)
MAX_REGION = 8000000 # Contains the maximum region size (found in the yaml file in the regions section)

CHROMOSOMES = ['chr17'] # Contains the chromosomes to analyze
                        # 1. Check the fix_output_vcf function to see the uniques chromosomes in the output VCF (uncomment the print in the function)
                        # 2. Add them in the list above

VAF_THRESHOLD = 0.15    # Contains the VAF threshold to filter the variants
                        # VAF tell the percentage of the variant in the tumor
                        # 0.15 means that the variant must be present in at least 15% of the tumor reads

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
    

def fix_output_vcf(df):
    # Remove the rows where the FILTER column is equal to "RefCall" and "GERMLINE" (we are not interested in GERMLINE variants)
    df = df[df['FILTER'] != 'RefCall']
    df = df[df['FILTER'] != 'GERMLINE']
    # print(df['FILTER'].unique())
    
    # Check the unique CHROM 
    # print(df['#CHROM'].unique())
    return df


def fix_ground_truth_vcf(df):
    # Remove the rows where the POS column is out the MIN-MAX REGION range
    df = df[(df['POS'] >= MIN_REGION) & (df['POS'] <= MAX_REGION)]

    # Filter the rows where the #CHROM column is equal to "chr17"
    df = df[df['#CHROM'].isin(CHROMOSOMES)]
    # print(df['#CHROM'].unique())
    return df


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

    # The result is a DataFrame like this: (considering the OUTPUT VCF)
    # if Filter = RefCall, it means that the variant is not a somatic mutation (we can delete this)
    # CHROM      POS ID REF ALT  QUAL   FILTER INFO              FORMAT                               TUMOR
    # 0  chr17  7671133  .   A   C   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL    0/0:32:98:95,3:0.0306122:0,32,41
    # 1  chr17  7671570  .   G   T   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL  0/0:49:122:118,4:0.0327869:0,52,51
    # 2  chr17  7671587  .   A   T   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL  0/0:46:136:129,6:0.0441176:0,46,52
    # 3  chr17  7671613  .   A   C   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL  0/0:43:122:115,6:0.0491803:0,52,43
    # 4  chr17  7671630  .   A   C   0.0  RefCall    .  GT:GQ:DP:AD:VAF:PL  0/0:47:128:119,5:0.0390625:0,48,54
    # --------------------------------------------------------------------------------------------------------------------------------------------
    # The result is a DataFrame like this: (considering the GROUND TRUTH VCF)
    # CHROM      POS ID REF ALT QUAL FILTER                                               INFO    FORMAT      HG008-N           HG008-T
    # 0  chr12    94894  .   T   C    .   PASS  CSQ=440073|protein_coding||intron_variant|||||...  AD:DP:AF  125,0:125:0  74,90:164:0.5488
    # 1  chr12   216351  .   C   T    .   PASS  CSQ=6540|protein_coding||downstream_gene_varia...  AD:DP:AF    57,0:57:0   42,56:98:0.5714
    # 2  chr12  1029855  .   C   G    .   PASS  CSQ=23085|protein_coding||intron_variant||||||...  AD:DP:AF  104,0:104:0  125,9:136:0.0662 
    # print(df.head(5), '\n\n')
    # print(df.tail(5), '\n\n')
    return df


def deepsomatic_scores(vcf_df, gt_df):
    GT_POS = gt_df['POS']

    # True Positives (present in both VCF and GT)
    tp = pd.merge(vcf_df, gt_df, on=['#CHROM', 'POS', 'REF', 'ALT'], how='inner')
    # False Negatives (present in GT but not in VCF)
    fn = pd.merge(gt_df, vcf_df, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left', indicator=True).query("_merge == 'left_only'").drop(columns=['_merge'])
    # False Positives (present in VCF but not in GT)
    fp = pd.merge(vcf_df, gt_df, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left', indicator=True).query("_merge == 'left_only'").drop(columns=['_merge'])

    # Creating VAF Value column
    # In tp VAF Value is the 5th value after the ':' in 'HG008-T_x' (the index is in 'FORMAT_x')
    tp['VAF_VALUE'] = tp['HG008-T_x'].apply(lambda x: float(str(x).split(':')[4]) if pd.notnull(x) else 0.0)
    # In fn VAF Value is the last value after the ':' in 'HG008-T_y' (the index is in 'FORMAT_y')
    fn['VAF_VALUE'] = fn['HG008-T_y'].apply(lambda x: float(str(x).split(':')[-1]) if pd.notnull(x) else 0.0)
    # In fp VAF Value is the 5th value after the ':' in 'HG008-T_x' (the index is in 'FORMAT_x') 
    fp['VAF_VALUE'] = fp['HG008-T_x'].apply(lambda x: float(str(x).split(':')[4]) if pd.notnull(x) else 0.0)
    
    # Filter the variants where the VAF is greater than the threshold
    tp_f = tp[tp['VAF_VALUE'] >= VAF_THRESHOLD]
    fn_f = fn[fn['VAF_VALUE'] >= VAF_THRESHOLD]
    fp_f = fp[fp['VAF_VALUE'] >= VAF_THRESHOLD]
    
    print("\n")
    print("===== Unfiltered =====")
    print("True Positives:", len(tp))
    print("False Negatives:", len(fn))
    print("False Positives:", len(fp))
    print("Precision:", len(tp) / (len(tp) + len(fp)))
    print("Recall:", len(tp) / (len(tp) + len(fn)))
    print("F1 Score:", 2 * (len(tp) / (len(tp) + len(fp))) * (len(tp) / (len(tp) + len(fn))) / ((len(tp) / (len(tp) + len(fp))) + (len(tp) / (len(tp) + len(fn)))))
    print("\n")
    print(f"===== Filtered VAF (>{VAF_THRESHOLD}) =====")
    print("True Positives:", len(tp_f))
    print("False Negatives:", len(fn_f))
    print("False Positives:", len(fp_f))
    print("Precision:", len(tp_f) / (len(tp_f) + len(fp_f)))
    print("Recall:", len(tp_f) / (len(tp_f) + len(fn_f)))
    print("F1 Score:", 2 * (len(tp_f) / (len(tp_f) + len(fp_f))) * (len(tp_f) / (len(tp_f) + len(fn_f))) / ((len(tp_f) / (len(tp_f) + len(fp_f))) + (len(tp_f) / (len(tp_f) + len(fn_f)))))
    print("\n")
    

def visualize_vcf_data(df_input):
    pass


if __name__ == "__main__":
    setup_path()

    # Read the Deepsomatic's output VCF file and fix it
    vcf_df = read_vcf(DEEPSOMATIC_OUTPUT)
    vcf_df = fix_output_vcf(vcf_df)

    # Read the ground truth VCF file
    gt_df = read_vcf(GROUND_TRUTH_VCF)
    gt_df = fix_ground_truth_vcf(gt_df)

    #print(vcf_df[vcf_df['POS'] == 7675217])
    #print(vcf_df.head())
    #print(gt_df.head())
    deepsomatic_scores(vcf_df, gt_df)

    # visualize_vcf_data(vcf_df)
    # visualize_vcf_data(gt_df)