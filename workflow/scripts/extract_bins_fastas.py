#!/usr/bin/env python3

# import sys
import os

def load_contig_bins(contig_bins_file):
    contig_to_bin = {}
    with open(contig_bins_file, 'r') as file:
        for line in file:
            contig, bin_name = line.strip().split()
            # Convert contig names and bin names to lowercase to handle case sensitivity
            contig = contig.strip().lower()
            bin_name = bin_name.strip().lower()
            if bin_name != "unbinned":
                contig_to_bin[contig] = bin_name
    return contig_to_bin

def extract_fasta_sequences(fasta_file, contig_to_bin):
    current_contig = None
    sequences = {}
    
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_contig = line[1:].strip().lower()  # Get contig name, convert to lowercase
                if current_contig in contig_to_bin:
                    bin_name = contig_to_bin[current_contig]
                    if bin_name not in sequences:
                        sequences[bin_name] = []
                    sequences[bin_name].append(line)  # Store header
                else:
                    current_contig = None
            elif current_contig:
                sequences[contig_to_bin[current_contig]].append(line)  # Store sequence line
    
    return sequences

def write_sequences_to_files(sequences, output_dir):
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)
    
    for bin_name, seq_lines in sequences.items():
        output_file_path = os.path.join(output_dir, f"{bin_name}.fa")
        with open(output_file_path, 'w') as out_file:
            out_file.write("\n".join(seq_lines) + "\n")

def main():
    # if len(sys.argv) != 4:
    #     print("Usage: python3 extract_fasta_by_bins.py <contig_bins.tsv> <input.fasta> <output_dir>")
    #     sys.exit(1)
    
    # contig_bins_file = sys.argv[2]
    # fasta_file = sys.argv[1]
    # output_dir = sys.argv[3]

    contig_bins_file = snakemake.input["contig_bin"]
    fasta_file = snakemake.input["fasta_file"]
    output_dir = snakemake.output["outdir"]

    os.makedirs(output_dir, exist_ok=True)
    
    contig_to_bin = load_contig_bins(contig_bins_file)
    sequences = extract_fasta_sequences(fasta_file, contig_to_bin)
    write_sequences_to_files(sequences, output_dir)

if __name__ == "__main__":
    main()
