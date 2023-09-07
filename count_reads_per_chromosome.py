import pysam

save = pysam.set_verbosity(0)
samFile = "data/output.bam"
bam = pysam.AlignmentFile(samFile, "rb")
pysam.set_verbosity(save)

chromosome_counts = {}  

for read in bam:
    chromosome = read.reference_name
    if chromosome in chromosome_counts:
        chromosome_counts[chromosome] += 1
    else:
        chromosome_counts[chromosome] = 1

bam.close()

for chromosome, count in chromosome_counts.items():
    print(f"Chromosome {chromosome}: {count} reads")

total_reads = sum(chromosome_counts.values())
print(f"Total reads: {total_reads}")