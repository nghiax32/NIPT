import os
import subprocess
# import rpy2.robjects as robjects
# from rpy2.robjects.packages import importr
import pysam
import numpy as np

def data_preprocessing(input_folder, output_folder):
    # create output folder if not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # list bam files
    bam_files = []
    for file in os.listdir(input_folder):
        if file.endswith(".bam"):
            bam_files.append(os.path.join(input_folder, file))
    # NGS data preprocessing
    for bam_file in bam_files:
        # remove duplicate reads
        input_bam = bam_file
        output_bam = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.remove_dup.bam'))
        marked_dup_metrics = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.marked_dup_metrics.txt'))
        command = ['java', '-jar', '/home/tori/Tori/picard/picard.jar', 'MarkDuplicates', f'I={input_bam}', f'O={output_bam}', f'M={marked_dup_metrics}', 'REMOVE_DUPLICATES=true']
        subprocess.run(command, check=True)
        os.remove(marked_dup_metrics)
        # remove mapping quality below 30
        input_bam = output_bam
        output_bam = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.above_30.bam'))
        command = ['samtools', 'view', '-b', '-h', '-q', '30', f'{input_bam}', '-o', f'{output_bam}']
        subprocess.run(command, check=True)
        os.remove(input_bam)
        # TODO: adjust the GC content and mappability
        # index bam file
        input_bam = output_bam
        command = ['samtools', 'index', f'{input_bam}']
        subprocess.run(command, check=True)

def calculate_fd(input_folder):
    # list bam files
    bam_files = []
    for file in os.listdir(input_folder):
        if file.endswith(".bam"):
            bam_files.append(os.path.join(input_folder, file))
    # process
    for bam_file in bam_files:
        # non-overlapping binning at 1 Mbp
        bam = pysam.AlignmentFile(str(bam_file), "rb")
        chr13 = bam.fetch("chr13")
        positions = []
        for read in chr13:
            # read_attributes = dir(read)
            # for read_attribute in read_attributes:
            #     print(read_attribute, getattr(read, read_attribute))
            # break
            forward = not(read.flag & (1 << 4))
            if forward: 
                positions.append(read.reference_start + 80)
            else:
                positions.append(read.reference_end - 80)
        positions = sorted(positions)
        bin_size = 1e6
        fragment_distance = [[] for _ in range(int(max(positions) // bin_size) + 1)]
        # print(positions)
        for i in range(len(positions) - 1):
            position = positions[i]
            bin_position = int(position // bin_size)
            position_next = positions[i+1]
            bin_position_next = int(position_next // bin_size)
            position_distance = position_next - position
            bin_seperate = bin_position_next * bin_size
            if bin_position == bin_position_next:
                fragment_distance[bin_position].append(position_distance)
            elif bin_seperate - position < position_next - bin_seperate:
                fragment_distance[bin_position].append(position_distance)
            else:
                fragment_distance[bin_position_next].append(position_distance)
        fd_mean = []
        fd_median = []
        fd_iqr = []
        for i in range(len(fragment_distance)):
            bin = fragment_distance[i]
            if len(bin) == 0:
                continue
            fd_mean.append(np.mean(bin))
            fd_median.append(np.median(bin))
            q1 = np.percentile(bin, 25)
            q3 = np.percentile(bin, 75)
            iqr = q3 - q1
            fd_iqr.append(iqr)
        
        fd_mean_median = np.median(fd_mean)
        fd_median_median = np.median(fd_median)
        fd_iqr_median = np.median(fd_iqr)
        print("before", fd_mean, fd_median, fd_iqr, sep='\n')
        fd_mean = [(fd_mean_value - fd_mean_median) / fd_mean_median for fd_mean_value in fd_mean]
        fd_median = [(fd_median_value - fd_median_median) / fd_median_median for fd_median_value in fd_median]
        fd_iqr = [(fd_iqr_value - fd_iqr_median) / fd_iqr_median for fd_iqr_value in fd_iqr]
        print("normalized", fd_mean, fd_median, fd_iqr, sep='\n')


def main():
    negatives_input_folder = "data/negatives"
    negatives_output_folder = "data/negatives_processed"
    # data_preprocessing(negatives_input_folder, negatives_output_folder)
    calculate_fd(negatives_output_folder)

if __name__ == "__main__":
    main()