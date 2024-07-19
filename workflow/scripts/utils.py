import json
from pickle import FALSE
import pysam
import argparse, sys, csv, re, io
import pandas as pd
from collections import defaultdict, namedtuple
from contextlib import redirect_stdout

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


"""
Linking AMR genes to taxonomic classification results - read mapping based
"""
class AMRlinker:

    Query_taxamap = namedtuple('Query_taxamap', ['Template', 'Query_Length_taxa', "Temp_Align_Start", "Temp_Align_End", "Observed_Templ_Length", "ExactMatch_bp_taxa", "Template_Length"])
    Query_NDAROmap = namedtuple('Query_NDAROmap', ["Query_Length_NDARO", "ExactMatch_bp_NDARO", "ARG_Length"])
    Query_taxamap_NDAROmap = namedtuple("Query_taxamap_NDAROmap", Query_taxamap._fields + Query_NDAROmap._fields)


    def __init__(self, taxa_bam, NDARO_amr_bam, mode="assembly", NDARO_mappingfile='/db/gene_detection/NCBI_AMR/mapping_full.json', genomad_dir=None):
        self.taxa_bam = taxa_bam
        self.NDARO_amr_bam = NDARO_amr_bam
        self.mode = mode
        self.NDARO_mappingfile = NDARO_mappingfile
        self.genomad_dir = genomad_dir

    def _parse_taxa(self):
        align_taxa = pysam.AlignmentFile(self.taxa_bam, "rb")
        align_taxa_refseq = {ref["SN"]: ref["LN"] for ref in align_taxa.header.to_dict()["SQ"]}

        taxa_mapping = defaultdict(lambda: self.Query_taxamap)

        for query in align_taxa:
            n = query.query_name
            query_l = query.query_length
            Template1 = query.reference_name
            start = query.reference_start
            end = query.reference_end
            l = query.reference_length
            cigarEQ = query.get_cigar_stats()[0][7]
            Template1_l = align_taxa_refseq.get(Template1, None)

            taxa_mapping[n] = self.Query_taxamap(Template1, query_l, start, end, l, cigarEQ, Template1_l)

        return taxa_mapping

    def _parse_reads_NDARO_bam(self, tid_arg=90.0, tc_arg=60.0):
        
        NDARO_readmapping = defaultdict(lambda: defaultdict(lambda: self.Query_NDAROmap))

        # load mapping of NDARO ARG template names to ARGs and alleles in NDARO
        def load_mapping(filepath):
            with open(filepath, 'r') as file:
                return json.load(file)
        mapping = load_mapping(self.NDARO_mappingfile)

        pattern = r"seq_\d{6}$"
        align_amr = pysam.AlignmentFile(self.NDARO_amr_bam, "rb")
        align_amr_refseq = {ref["SN"]: ref["LN"] for ref in align_amr.header.to_dict()["SQ"]}

        for alignment in align_amr:
            n = alignment.query_name
            if alignment.reference_name is None:
                continue
            
            l= alignment.query_length
            Template2 = alignment.reference_name
            cigarEQ = alignment.get_cigar_stats()[0][7]
            Template2_l = align_amr_refseq.get(Template2, None)
            if 100 * alignment.reference_length < tc_arg or 100 * cigarEQ / alignment.reference_length < tid_arg:
                continue
            
            arg_match = re.search(pattern, Template2)
            if arg_match:
                arg = arg_match.group()
                arg_name = mapping[arg]["gene"]
                allele = mapping[arg]["allele"]
                # ARG mapping is not 1-to-1 e.g. contig can hold 2 or more ARGs -> nested dictionary
                NDARO_readmapping[n][(arg_name, allele)] = self.Query_NDAROmap(l, cigarEQ, Template2_l)

        return NDARO_readmapping  


    def match(self):

        self.AMRlinks = defaultdict(lambda: defaultdict(lambda: self.Query_taxamap_NDAROmap))

        taxa_mapping = self._parse_taxa()
        NDARO_mapping = self._parse_reads_NDARO_bam()

        for n in NDARO_mapping.keys():
            
            if n in taxa_mapping.keys():
                for allele in NDARO_mapping[n].keys():
                    self.AMRlinks[n][allele] = self.Query_taxamap_NDAROmap(*(taxa_mapping[n] + NDARO_mapping[n][allele]))
            else:
                for allele in NDARO_mapping[n].keys():
                    taxa_mapping[n] = self.Query_taxamap(None, NDARO_mapping[n][allele].Query_Length_NDARO, None, None, None, None, None)
                    self.AMRlinks[n][allele] = self.Query_taxamap_NDAROmap(*(taxa_mapping[n] + NDARO_mapping[n][allele]))
        
        flattened_data = []
        for query, arg_dict in self.AMRlinks.items():
            for (arg_name, allele), values in arg_dict.items():
                flattened_data.append((query, arg_name, allele) + values)

        columns = ['Query', 'ARG', 'Allele'] + list(self.Query_taxamap_NDAROmap._fields)
        self.AMRlinks = pd.DataFrame(flattened_data, columns=columns)
        
        if self.mode == "reads":
            grouped = self.AMRlinks.groupby(["Template", "Template_Length", "ARG", "Allele", "ARG_Length"], dropna=False)
                
            n_links = grouped.size()
            Total_Query_Length = grouped['Query_Length_NDARO'].sum()
            Observed_Templ_Depth = grouped['Observed_Templ_Length'].sum()
            ExactMatch_Taxa = grouped['ExactMatch_bp_taxa'].sum()
            ExactMatch_NDARO = grouped['ExactMatch_bp_NDARO'].sum()

            self.AMRlinks = pd.DataFrame({'n_links': n_links, 'Total_Query_Length': Total_Query_Length, 'Observed_Templ_Depth' : Observed_Templ_Depth, 'ExactMatch_Taxa' : ExactMatch_Taxa, 'ExactMatch_NDARO' : ExactMatch_NDARO})
            self.AMRlinks = self.AMRlinks.reset_index()


    def get_covinfo(self, coverage_bam):

        cov_output = pysam.coverage(coverage_bam)
        coverage_data = []
        for line in cov_output.split('\n'):
        # for line in cov_output:
            if line and not line.startswith('#'):
                fields = line.split('\t')
                rname = fields[0]
                converted_fields = [rname] + [int(x) for x in fields[1:5]] + [float(x) for x in fields[5:]]
                coverage_data.append(converted_fields)

        columns = ['rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq']
        coverage_data = pd.DataFrame(coverage_data, columns=columns)

        # if in assembly mode, we want coverage of contigs (=queries)
        if self.mode == "assembly":
            merged_df = self.AMRlinks.merge(coverage_data, left_on="Query", right_on="rname", how="outer")
            self.AMRlinks = merged_df
            self.AMRlinks["Query"] = self.AMRlinks["Query"].fillna(self.AMRlinks["rname"])
        # in read mode, match to the taxonomic template mapping
        else:
            merged_df = self.AMRlinks.merge(coverage_data, left_on="Template", right_on="rname", how="outer")
            self.AMRlinks = merged_df

    
    def get_genomad_info(self):
        if self.genomad_dir is None:
            raise ValueError("Genomad output directory is not available.")
        if self.mode == "reads":
            raise ValueError("Genomad not implemented for reads")
        
        plasmid_summary = self.genomad_dir.rstrip("/") + "/consensus_summary/consensus_plasmid_summary.tsv"
        virus_summary = self.genomad_dir.rstrip("/") + "/consensus_summary/consensus_virus_summary.tsv"


