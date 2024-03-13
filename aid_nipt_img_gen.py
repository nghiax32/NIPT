import os
import subprocess
import pysam
import numpy as np
import matplotlib.pyplot as plt

# NGS data preprocessing.
def data_preprocessing(input_folder, output_folder):
    # Create output folder if not exist.
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # List bam files.
    bam_files = []
    for file in os.listdir(input_folder):
        if file.endswith(".bam"):
            bam_files.append(os.path.join(input_folder, file))
    # Data preprocessing.
    for bam_file in bam_files:
        # Remove duplicate reads.
        input_bam = bam_file
        output_bam = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.remove_dup.bam'))
        marked_dup_metrics = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.marked_dup_metrics.txt'))
        command = ['java', '-jar', '/home/tori/Tori/picard/picard.jar', 'MarkDuplicates', f'I={input_bam}', f'O={output_bam}', f'M={marked_dup_metrics}', 'REMOVE_DUPLICATES=true']
        subprocess.run(command, check=True)
        os.remove(marked_dup_metrics)
        # Remove mapping quality below 30.
        input_bam = output_bam
        output_bam = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.above_30.bam'))
        command = ['samtools', 'view', '-b', '-h', '-q', '30', f'{input_bam}', '-o', f'{output_bam}']
        subprocess.run(command, check=True)
        os.remove(input_bam)
        # TODO: Adjust the GC content and mappability.
        # Index bam file.
        input_bam = output_bam
        command = ['samtools', 'index', f'{input_bam}']
        subprocess.run(command, check=True)

# Removal of lower 10% and upper 10% bins of the normalized FD representative values.
def remove_outliers(values):
    sorted_values = sorted(values)
    lower_index = int(len(values) * 0.1)
    upper_index = int(len(values) * 0.9)
    trimmed_values = sorted_values[lower_index:upper_index]
    result = [value for value in values if value in trimmed_values]
    return result 

# Calculation of FD.
def fd_calculatation(bam, chr):
    # Non-overlapping binning at 1 Mbp.
    chr_bam = bam.fetch(chr)
    positions = []
    for read in chr_bam:
        forward = not(read.flag & (1 << 4))
        if forward: 
            positions.append(read.reference_start + 80)
        else:
            positions.append(read.reference_end - 80)
    positions = sorted(positions)
    bin_size = 1e6
    fragment_distance = [[] for _ in range(int(max(positions) // bin_size) + 1)]
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
    fd_mean = [fd_mean_value / fd_mean_median for fd_mean_value in fd_mean]
    fd_mean = remove_outliers(fd_mean)
    fd_median = [fd_median_value / fd_median_median for fd_median_value in fd_median]
    fd_median = remove_outliers(fd_median)
    fd_iqr = [fd_iqr_value/ fd_iqr_median for fd_iqr_value in fd_iqr]
    fd_iqr = remove_outliers(fd_iqr)
    return fd_mean, fd_median, fd_iqr

# Image generation.
def image_generation(output_image, fd_tc, fd_icc):
    fig, axs = plt.subplots(6, 1, figsize=(4, 2), dpi=200)
    for i in range(3):
        axs[i*2].plot(fd_tc, color='black')
        axs[i*2].axis('off')
        axs[i*2].set_xlim(0, len(fd_tc))
        axs[i*2].set_ylim(min(fd_tc), max(fd_tc))
        axs[i*2].fill_between(range(len(fd_tc)), fd_tc, color = "black")

        axs[i*2+1].plot(fd_icc[i], color='black')
        axs[i*2+1].axis('off')
        axs[i*2+1].set_xlim(0, len(fd_icc[i]))
        axs[i*2+1].set_ylim(min(fd_icc[i]), max(fd_icc[i]))
        axs[i*2+1].fill_between(range(len(fd_icc[i])), fd_icc[i], color = "black")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_image, dpi=200, bbox_inches='tight', pad_inches=0)

    plt.cla()
    plt.clf()
    plt.close('all')
    plt.close(fig)

# TRS image generation.
def trs_image_generation(input_folder, output_folder, tc, icc):
    # Create output folder if not exist.
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # List bam files.
    bam_files = []
    for file in os.listdir(input_folder):
        if file.endswith(".bam"):
            bam_files.append(os.path.join(input_folder, file))
    # Process.
    for bam_file in bam_files:
        bam = pysam.AlignmentFile(str(bam_file), "rb")
        fd_mean_tc, fd_median_tc, fd_iqr_tc = fd_calculatation(bam, tc)
        fd_mean_icc = [0] * 3
        fd_median_icc = [0] * 3
        fd_iqr_icc = [0] * 3
        fd_mean_icc[0], fd_median_icc[0], fd_iqr_icc[0] = fd_calculatation(bam, icc[0])
        fd_mean_icc[1], fd_median_icc[1], fd_iqr_icc[1] = fd_calculatation(bam, icc[1])
        fd_mean_icc[2], fd_median_icc[2], fd_iqr_icc[2] = fd_calculatation(bam, icc[2])
        output_image_mean = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.mean.png'))
        output_image_median = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.median.png'))
        output_image_iqr = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.iqr.png'))
        image_generation(output_image_mean, fd_mean_tc, fd_mean_icc)
        image_generation(output_image_median, fd_median_tc, fd_median_icc)
        image_generation(output_image_iqr, fd_iqr_tc, fd_iqr_icc)

def main():
    negatives_input_folder = "data/negatives"
    negatives_output_folder = "data/negatives_processed"
    negatives_image_folder = "data/negatives_image"
    # data_preprocessing(negatives_input_folder, negatives_output_folder)
    trs_image_generation(negatives_output_folder, negatives_image_folder, "chr13", ["chr4", "chr5", "chr6"])

    positives_input_folder = "data/positives"
    positives_output_folder = "data/positives_processed"
    positives_image_folder = "data/positives_image"
    # data_preprocessing(positives_input_folder, positives_output_folder)
    trs_image_generation(positives_output_folder, positives_image_folder, "chr13", ["chr4", "chr5", "chr6"])

if __name__ == "__main__":
    main()