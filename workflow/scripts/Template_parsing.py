import re
import pandas as pd
"""
Functions to parse template names
"""

def extract_species_name(string):
    # general pattern for scientific species name
    pattern = re.compile(r'[A-Z][a-z]+ [a-z]+')
    # handle na
    species = string
    if pd.isna(species):
        species = "Unmapped"
    
    species = re.sub(r'\[|\]|\(|\)', '', species)
    species = re.sub('_', ' ', species)
    species = re.sub(r'Synthetic|scf', '', species)
    species = re.sub('E.coli', 'Escherichia coli', species)
    species = re.sub("veillonella", "Veillonella", species)
    species = re.sub("Clostridium difficil+e", "Clostridioides difficile", species)
    species = re.sub("Limosilactobacillus fermentum", "Lactobacillus fermentum", species)
    if not re.search(pattern, species):
        species = "Unmapped"
    # handle virus and phages: \W\s required as e.g. Treponema phagedenis exists
    elif re.search(r"([V|v]irus[\W\s])|([P|p]hage[\W\s])", species):
        if re.search(r"[V|v]irus", species):
            species = re.search(r"[A-z\s]+[V|v]irus", species).group()
        else:
            # phage with full name continuing behind "phage"
            if re.search(r"[A-z\s]+[P|p]hage[\W\s][A-z0-9]+", species):
                species = species.rstrip(", complete sequence")
                species = re.search(r"[A-z\s]+[P|p]hage[\W\s][A-z0-9]+", species).group()
            # phage with full name ending with "phage,"
            elif re.search(r"[A-z\s]+[P|p]hage,", species):
                species = re.search(r"[A-z\s]+[P|p]hage", species).group()
            else:
                # naming unclear --> keep complete name
                pass
                # species = re.search(pattern, species).group()
    elif re.search(r"[P|p]lasmid", species):
        if re.search(r"\|kraken:taxid\|\d+\s", species):
            species = species.rstrip(",\t")
            species = species.rstrip(", complete sequence")
            species = re.search(r"\|kraken:taxid\|\d+\s.+", species).group()
            species = re.sub(r"\|kraken:taxid\|\d+\s", "", species)
            
        else:
            # plasmid naming conventions unclear --> just keep the complete plasmid name
            pass
            # species = re.search(pattern, species).group()
    else:
        species = re.search(pattern, species).group()

    return species


def extract_arg_name(string):
    ARG = re.sub(r"_.+", "", string)
    ARG = re.sub(r"-\d+.+$", "", ARG)
    return ARG
