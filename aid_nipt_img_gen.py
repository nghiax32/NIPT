import os
import subprocess
import pysam
import numpy as np
import matplotlib.pyplot as plt

# NGS data preprocessing.
def data_preprocessing(output_folder, bam_file):
    # Remove duplicate reads.
    input_bam = bam_file
    output_bam = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.remove_dup.bam'))
    marked_dup_metrics = os.path.join(output_folder, os.path.basename(bam_file).replace('.bam', '.marked_dup_metrics.txt'))
    # PICARD = '/home/tinhnh/aidNIPT/picard/picard.jar'
    # PICARD = '/home/tori/Tori/picard/picard.jar'
    command = ['java', '-jar', PICARD, 'MarkDuplicates', f'I={input_bam}', f'O={output_bam}', f'M={marked_dup_metrics}', 'REMOVE_DUPLICATES=true']
    subprocess.run(command, check=True)
    os.remove(marked_dup_metrics)

    # Remove mapping quality below 30.
    input_bam = output_bam
    output_bam = os.path.join(output_folder, os.path.basename(input_bam).replace('.bam', '.above_30.bam'))
    command = ['samtools', 'view', '-b', '-h', '-q', '30', f'{input_bam}', '-o', f'{output_bam}']
    subprocess.run(command, check=True)
    os.remove(input_bam)

    # TODO: Adjust the GC content and mappability.

    # Index bam file.
    input_bam = output_bam
    output_bai = os.path.join(output_folder, os.path.basename(input_bam).replace('.bam', '.bam.bai'))
    command = ['samtools', 'index', f'{input_bam}']
    subprocess.run(command, check=True)
    return output_bam, output_bai

# Removal of lower 5% and upper 5% bins of the normalized FD representative values.
def remove_outliers(values):
    sorted_values = sorted(values)
    lower_index = int(len(values) * 0.05)
    upper_index = int(len(values) * 0.95)
    trimmed_values = sorted_values[lower_index:upper_index]
    result = [value for value in values if value in trimmed_values]
    return result 

# Calculate FD of chromosome.
def calculate_fd_of_chr(bam, chr):
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

# Save 1D array to text file.
def save_1d_array_to_file(output_file, arr, header):
    with open(output_file, 'a') as file:
        file.write(header + '\n')
        for item in arr:
            file.write(str(item) + '\n')

# Save 2D array to text file.
def save_2d_array_to_file(output_file, arr, header):
    with open(output_file, 'a') as file:
        file.write(header + '\n')
        for row in arr:
            file.write(' '.join(map(str, row)) + '\n')

# Read array from text file.
def read_array_from_file(input_file):
    with open(input_file, 'r') as file:
        data = {}
        current_list = None
        arr_type = None 
        for line in file:
            line = line.strip()
            if line.startswith("#"):
                current_list = line[1:]
                data[current_list] = []
            else:
                line = line.split(' ')
                if len(line) > 1:
                    line = [float(x) for x in line]
                else:
                    line = ''.join(line)
                    line = float(line)
                data[current_list].append(line)
        return data

# Calculation of FD.
def fd_calculatation(output_file, bam, tc, icc):
    # Calculate values.
    fd_mean_tc, fd_median_tc, fd_iqr_tc = calculate_fd_of_chr(bam, tc)
    fd_mean_icc = [0] * 3; fd_median_icc = [0] * 3; fd_iqr_icc = [0] * 3
    fd_mean_icc[0], fd_median_icc[0], fd_iqr_icc[0] = calculate_fd_of_chr(bam, icc[0])
    fd_mean_icc[1], fd_median_icc[1], fd_iqr_icc[1] = calculate_fd_of_chr(bam, icc[1])
    fd_mean_icc[2], fd_median_icc[2], fd_iqr_icc[2] = calculate_fd_of_chr(bam, icc[2])

    # Save data to txt file.
    save_1d_array_to_file(output_file, fd_mean_tc, '#fd_mean_tc')
    save_1d_array_to_file(output_file, fd_median_tc, '#fd_median_tc')
    save_1d_array_to_file(output_file, fd_iqr_tc, '#fd_iqr_tc')
    save_2d_array_to_file(output_file, fd_mean_icc, '#fd_mean_icc')
    save_2d_array_to_file(output_file, fd_median_icc, '#fd_median_icc')
    save_2d_array_to_file(output_file, fd_iqr_icc, '#fd_iqr_icc')

# TRS image generation.
def trs_image_generation(output_image, fd_tc, fd_icc):
    fig, axs = plt.subplots(6, 1, figsize=(4, 2), dpi=200)
    for i in range(3):
        axs[i*2].plot(fd_tc, color='black')
        axs[i*2].axis('off')
        axs[i*2].set_xlim(0, len(fd_tc))
        axs[i*2].set_ylim(min(fd_tc), max(fd_tc))
        axs[i*2].fill_between(range(len(fd_tc)), fd_tc, color = 'black')

        axs[i*2+1].plot(fd_icc[i], color='black')
        axs[i*2+1].axis('off')
        axs[i*2+1].set_xlim(0, len(fd_icc[i]))
        axs[i*2+1].set_ylim(min(fd_icc[i]), max(fd_icc[i]))
        axs[i*2+1].fill_between(range(len(fd_icc[i])), fd_icc[i], color = 'black')
        
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_image, bbox_inches='tight', pad_inches=0)

    plt.cla()
    plt.clf()
    plt.close('all')
    plt.close(fig)        

def handle(data_folder, value_folder, image_folder, tc, icc):
    # Create value folder and image folder if not exist.
    if not os.path.exists(value_folder):
        os.makedirs(value_folder)
    if not os.path.exists(image_folder):
        os.makedirs(image_folder)

    # List data files with absolute paths.
    bam_files = []
    for file in os.listdir(data_folder):
        if file.endswith('.bam'):
            bam_files.append(os.path.join(data_folder, file))

    # List value files with absolute paths.
    value_files = []
    for file in os.listdir(value_folder):
        if file.endswith('.txt'):
            sample_name = file.split('.txt')[0]
            if sample_name not in value_files:
                value_files.append(sample_name)

    # List image files with absolute paths.
    image_files = []
    for file in os.listdir(image_folder):
        if file.endswith('.mean.png'):
            sample_name = file.split('.mean.png')[0]
            if sample_name not in image_files:
                image_files.append(sample_name)

    # Loop through the list and process it.
    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split('.bam')[0]
        if sample_name not in value_files:
            # preprocessed_bam, preprocessed_bai = data_preprocessing(value_folder, bam_file)
            bam = pysam.AlignmentFile(str(bam_file), 'rb')
            value_file = os.path.join(value_folder, os.path.basename(bam_file).replace('.bam', '.txt'))
            fd_calculatation(value_file, bam, tc, icc)
            # os.remove(preprocessed_bam)
            # os.remove(preprocessed_bai)
        if sample_name not in image_files:
            # Get data from file. 
            value_file = os.path.join(value_folder, os.path.basename(bam_file).replace('.bam', '.txt'))
            data = read_array_from_file(value_file)
            fd_mean_tc = data['fd_mean_tc']
            fd_median_tc = data['fd_median_tc']
            fd_iqr_tc = data['fd_iqr_tc']
            fd_mean_icc = data['fd_mean_icc']
            fd_median_icc = data['fd_median_icc']
            fd_iqr_icc = data['fd_iqr_icc']
            # Generate image.
            image_file_mean = os.path.join(image_folder, os.path.basename(bam_file).replace('.bam', '.mean.png'))
            image_file_median = os.path.join(image_folder, os.path.basename(bam_file).replace('.bam', '.median.png'))
            image_file_iqr = os.path.join(image_folder, os.path.basename(bam_file).replace('.bam', '.iqr.png'))
            trs_image_generation(image_file_mean, fd_mean_tc, fd_mean_icc)
            trs_image_generation(image_file_median, fd_median_tc, fd_median_icc)
            trs_image_generation(image_file_iqr, fd_iqr_tc, fd_iqr_icc)


def main():
    # NEGATIVE_DATA_FOLDER_1 = '/data/tinhnh/NIPT/data/negatives'
    # NEGATIVE_DATA_FOLDER_2 = '/home/tinhnh/positives2'
    # POSITIVE_DATA_FOLDER = '/data/tinhnh/NIPT/data/positives'
    NEGATIVE_DATA_FOLDER_1 = 'data/negatives'
    NEGATIVE_DATA_FOLDER_2 = 'data/negatives'
    POSITIVE_DATA_FOLDER = 'data/positives'

    negative_data_folder = NEGATIVE_DATA_FOLDER_1
    negative_value_folder = 'data/negative_values'
    negative_image_folder = 'data/negative_images'
    handle(negative_data_folder, negative_value_folder, negative_image_folder, 'chr13', ['chr4', 'chr5', 'chr6'])

    negative_data_folder = NEGATIVE_DATA_FOLDER_2
    negative_value_folder = 'data/negative_values'
    negative_image_folder = 'data/negative_images'
    handle(negative_data_folder, negative_value_folder, negative_image_folder, 'chr13', ['chr4', 'chr5', 'chr6'])

    positive_data_folder = POSITIVE_DATA_FOLDER
    positive_value_folder = 'data/positive_values'
    positive_image_folder = 'data/positive_images'
    handle(positive_data_folder, positive_value_folder, positive_image_folder, 'chr13', ['chr4', 'chr5', 'chr6'])

if __name__ == '__main__':
    main()