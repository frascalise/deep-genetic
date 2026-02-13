import pandas as pd
import gzip
import yaml
import os
import io


DEEPSOMATIC_OUTPUT = None # Contains the path to the output VCF file
DEEPSOMATIC_OUTPUT_GVCF = None # Contains the path to the output GVCF file


def setup_path():
    with open("params.yaml", "r") as f:
        params = yaml.safe_load(f)

    global DEEPSOMATIC_OUTPUT, DEEPSOMATIC_OUTPUT_GVCF
    DEEPSOMATIC_OUTPUT = params["output_vcf"]
    DEEPSOMATIC_OUTPUT_GVCF = params['output_gvcf'] if params.get('output_gvcf') else None

    # Check if the files exists
    if not os.path.exists(DEEPSOMATIC_OUTPUT):
        print("ERROR: FILE", params["output_vcf"], "NOT FOUND")
        exit()

    if DEEPSOMATIC_OUTPUT_GVCF and not os.path.exists(DEEPSOMATIC_OUTPUT_GVCF):
        print("ERROR: FILE", params["output_gvcf"], "NOT FOUND")
        exit()


def read_vcf():
    with gzip.open(DEEPSOMATIC_OUTPUT, "rt") as f:
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
    # print(df.head(10))
    return df


if __name__ == "__main__":
    setup_path()

    print("DEEPSOMATIC_OUTPUT: ", DEEPSOMATIC_OUTPUT)
    print("DEEPSOMATIC_OUTPUT_GVCF: ", DEEPSOMATIC_OUTPUT_GVCF)

    # Read the VCF file
    df = read_vcf()