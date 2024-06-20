import os
import sys

# Specify the directory containing your FASTA files
directory_path = sys.argv[1]

# Function to parse a FASTA file
def parse_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ""
        header = ""
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences.append((header, sequence))
                header = line[1:]
                sequence = ""
            else:
                sequence += line
        if header:
            sequences.append((header, sequence))
    return sequences

# Iterate through files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna"):
        file_path = os.path.join(directory_path, filename)

        bin=filename.rstrip("\..+")
        # Iterate through the sequences in the FASTA file
        sequences = parse_fasta(file_path)
        bin = filename.rsplit('.', 1)[0]

        for header, sequence in sequences:
            print(f"{header}\t{bin}")

