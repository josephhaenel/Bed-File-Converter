
'''
    Joseph Haenel
    10/3/2024
    CS 490 Computing in Healthcare
    Project 1 - Creating Negative Bed File
    (Used SeqIO so I don't have to use a Linux VM)
'''

import os
# import sys
from Bio import SeqIO
from collections import Counter

# Function to determine the optimal sequence length X based on the input data file
def find_optimal_sequence_length(file_path):
    """
    Determines the most common sequence length (optimal length X) from a data file.

    Args:
        file_path (str): Path to the input data file containing nucleotide sequences.

    Returns:
        int: The most common sequence length found in the file.
    """
    lengths = []  # List to store the lengths of all sequences
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):  # Ignore header lines
                sequence = line.strip()  # Get the nucleotide sequence
                lengths.append(len(sequence))  # Store the length of the sequence

    # Find the most common sequence length
    length_counts = Counter(lengths)  # Count the occurrences of each length
    optimal_length = length_counts.most_common(1)[0][0]  # Get the most common length
    return optimal_length  # Return the optimal length X for cropping


# Function to filter sequences to length X and append them to the combined positive file
def append_to_positive_file(file_path, positive_file, optimal_length):
    """
    Appends sequences of length X (optimal length) from the input file to the combined positive sequences file.

    Args:
        file_path (str): Path to the input data file containing nucleotide sequences.
        positive_file (file object): File object to which the positive sequences will be appended.
        optimal_length (int): The optimal sequence length to which sequences should be trimmed.
    """
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):  # Ignore header lines
                sequence = line.strip()  # Get the nucleotide sequence
                if len(sequence) >= optimal_length:
                    # Write the sequence trimmed to optimal length
                    positive_file.write(f"{sequence[:optimal_length]}\n")


# Function to generate the negative file based on gaps between regions
def append_to_negative_file(data_file, genome_file, negative_file, optimal_length):
    """
    Identifies gaps between regions in the input data file and appends negative sequences
    (extracted from the genome file) of length X to the combined negative sequences file.

    Args:
        data_file (str): Path to the input data file containing region information.
        genome_file (str): Path to the genome FASTA file to extract negative sequences.
        negative_file (file object): File object to which the negative sequences will be appended.
        optimal_length (int): The optimal sequence length to which sequences should be cropped.
    """
    # Load the genome sequences using SeqIO
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    regions_by_chrom = {}  # Dictionary to store regions by chromosome

    # First pass: Collect all regions by chromosome
    with open(data_file, 'r') as data:
        for line in data:
            if line.startswith('>'):  # This line contains the chromosome and region information
                chrom, start_end = line.strip()[1:].split(':')  # Extract chromosome and position
                start, end = map(int, start_end.split('-'))  # Split start and end positions into integers

                # Add the region to the chromosome's list of regions
                if chrom not in regions_by_chrom:
                    regions_by_chrom[chrom] = []
                regions_by_chrom[chrom].append((start, end))

    # Second pass: Process each chromosome's regions in sorted order
    for chrom, regions in regions_by_chrom.items():
        # Sort regions by their start position
        regions.sort()

        previous_end = None  # Store the end position of the previous region

        # Loop through the sorted regions
        for start, end in regions:
            if previous_end is not None:
                gap_start = previous_end + 1  # Start of the gap (right after the previous region)
                gap_end = start - 1  # End of the gap (just before the current region)

                if gap_end - gap_start >= optimal_length:  # If the gap is long enough for a sequence
                    negative_start = gap_start
                    negative_end = gap_start + optimal_length - 1  # Crop the gap to length X
                    sequence = genome[chrom].seq[negative_start:negative_end]  # Get the sequence
                    negative_file.write(f"{sequence}\n")  # Write the negative sequence to the file

            # Update the previous region's end position
            previous_end = end


# Main function that processes all files in /data and saves combined results to /output
def process_all_files(data_dir, genome_file, output_dir):
    """
    Processes all .txt files in the /data directory, extracts positive and negative sequences,
    and saves the combined positive and negative sequences in /output.

    Args:
        data_dir (str): Path to the directory containing the data files (.txt format).
        genome_file (str): Path to the genome FASTA file.
        output_dir (str): Path to the directory where the combined positive and negative files will be saved.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the combined positive and negative output files
    combined_positive_file = os.path.join(output_dir, 'all_positive_sequences.txt')
    combined_negative_file = os.path.join(output_dir, 'all_negative_sequences.txt')

    # Open the output files for appending
    with open(combined_positive_file, 'w') as positive_file, open(combined_negative_file, 'w') as negative_file:
        # Iterate over all files in the /data directory
        for data_file in os.listdir(data_dir):
            # Only process files with a .txt extension
            if data_file.endswith('.txt'):
                data_file_path = os.path.join(data_dir, data_file)
                print(f"Processing file: {data_file_path}")

                # Step 1: Determine the optimal sequence length X for this file
                optimal_length = find_optimal_sequence_length(data_file_path)
                print(f"Optimal sequence length for {data_file}: {optimal_length}")

                # Step 2: Append positive sequences from this file to the combined positive file
                append_to_positive_file(data_file_path, positive_file, optimal_length)
                print(f"Positive sequences from {data_file} appended to {combined_positive_file}")

                # Step 3: Append negative sequences from this file to the combined negative file
                append_to_negative_file(data_file_path, genome_file, negative_file, optimal_length)
                print(f"Negative sequences from {data_file} appended to {combined_negative_file}")

    print(f"All files processed. Combined positive sequences saved to {combined_positive_file}")
    print(f"Combined negative sequences saved to {combined_negative_file}")


if __name__ == '__main__':
    # Define paths
    data_dir = 'data'  # Directory containing the data files
    genome_file = 'humanGenome/hg38.fa'  # Genome file
    output_dir = 'output'  # Directory to save combined positive and negative files

    # Process all files and generate combined outputs
    process_all_files(data_dir, genome_file, output_dir)
