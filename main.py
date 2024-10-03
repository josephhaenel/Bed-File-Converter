import os
import sys
from Bio import SeqIO
from collections import Counter


# Function to determine the optimal sequence length X based on the input data file
def find_optimal_sequence_length(file_path):
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


# Function to create the "positive" file by filtering sequences to the optimal length X
def generate_positive_file(file_path, output_path, optimal_length):
    with open(file_path, 'r') as file, open(output_path, 'w') as positive_file:
        for line in file:
            if not line.startswith('>'):  # Ignore header lines
                sequence = line.strip()  # Get the nucleotide sequence
                if len(sequence) >= optimal_length:
                    # Trim the sequence to the optimal length X and write to the positive file
                    positive_file.write(f"{sequence[:optimal_length]}\n")


# Function to create the "negative" file by identifying gaps between regions in the data file
def generate_negative_file(data_file, genome_file, output_path, optimal_length):
    # Load the entire genome from the genome file (FASTA format) using SeqIO
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    previous_end = None  # Store the end position of the previous region
    previous_chrom = None  # Store the chromosome of the previous region

    with open(data_file, 'r') as data, open(output_path, 'w') as negative_file:
        for line in data:
            if line.startswith('>'):  # This line contains the chromosome and region information
                chrom, start_end = line.strip()[1:].split(':')  # Extract chromosome and position
                start, end = map(int, start_end.split('-'))  # Split start and end positions

                # Only calculate gaps within the same chromosome
                if previous_end is not None and previous_chrom == chrom:
                    gap_start = previous_end + 1  # Start of the gap
                    gap_end = start - 1  # End of the gap

                    # If the gap size is larger than the optimal length X, extract a region
                    if gap_end - gap_start >= optimal_length:
                        negative_start = gap_start
                        negative_end = gap_start + optimal_length - 1  # Crop to length X
                        sequence = genome[chrom].seq[negative_start:negative_end]  # Get the sequence

                        # Write the negative region sequence to the file (no header)
                        negative_file.write(f"{sequence}\n")

                # Update the previous_end and previous_chrom for the next iteration
                previous_end = end
                previous_chrom = chrom


# Main function that coordinates generating the positive and negative sequence files
def process_files(data_file, genome_file):
    # Step 1: Determine the optimal sequence length X from the data file
    optimal_length = find_optimal_sequence_length(data_file)
    print(f"Optimal sequence length: {optimal_length}")

    # Step 2: Create the "positive" file by trimming sequences to the optimal length X
    positive_file = f"{os.path.splitext(data_file)[0]}_pos.txt"  # Positive file name
    generate_positive_file(data_file, positive_file, optimal_length)
    print(f"Positive sequences saved to {positive_file}")
    
    # Step 3: Create the "negative" file by identifying gaps between regions
    negative_file = f"{os.path.splitext(data_file)[0]}_neg.txt"  # Negative file name
    generate_negative_file(data_file, genome_file, negative_file, optimal_length)
    print(f"Negative sequences saved to {negative_file}")


if __name__ == '__main__':
    # Ensure that the correct number of command-line arguments are passed
    if len(sys.argv) != 3:
        print("Usage: python main.py <data_file> <genome_file>")
        sys.exit(1)  # Exit the program if arguments are incorrect

    # Get the data file and genome file from the command-line arguments
    data_file = sys.argv[1]
    genome_file = sys.argv[2]

    # Call the process_files function to generate positive and negative sequence files
    process_files(data_file, genome_file)
