import pysam

bam_file = "deepsomatic_test/NA12878_S1.chr20.10_10p1mb.bam"
samfile = pysam.AlignmentFile(bam_file, "rb")

for read in samfile:
    print(read)
    break

print("Hello World")