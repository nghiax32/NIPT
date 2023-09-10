import pysam


# Read file.
samFile = "data/output.sorted.bam"
bam = pysam.AlignmentFile(samFile, "rb")

# See what's in a read
# print(next(bam).__dir__())

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