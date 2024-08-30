import json
import pysam
from tqdm import tqdm
import argparse, sys, csv, re, io
import pandas as pd
import numpy as np
import networkx as nx
import pytaxonkit
import matplotlib.pyplot as plt
from collections import defaultdict, namedtuple


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
Linking AMR genes to taxonomic classification results
"""
class AMRlinker:

    Query_taxamap = namedtuple('Query_taxamap', ['Template', 'Query_Length_taxa', "Temp_Align_Start", "Temp_Align_End", "Observed_Templ_Length", "ExactMatch_bp_taxa", "Template_Length"])
    Query_ARGmap = namedtuple('Query_ARGmap', ["Query_Align_Start_ARG", "Query_Align_End_ARG", "Query_Length_ARG", "ExactMatch_bp_ARG", "ARG_Length"])
    Query_taxamap_ARGmap = namedtuple("Query_taxamap_ARGmap", Query_taxamap._fields + Query_ARGmap._fields)


    def __init__(self, taxa_bam, ARG_amr_bam, mode="assembly", 
    genomad_dir=None, taxonkit_db = "/data/brbloemen/tools/.taxonkit/"):
        self.taxa_bam = taxa_bam
        self.ARG_amr_bam = ARG_amr_bam
        self.mode = mode
        self.genomad_dir = genomad_dir
        self.taxonkit_db = taxonkit_db

    def _parse_taxa(self):
        align_taxa = pysam.AlignmentFile(self.taxa_bam, "rb")
        align_taxa_refseq = {ref["SN"]: ref["LN"] for ref in align_taxa.header.to_dict()["SQ"]}

        taxa_mapping = defaultdict(lambda: self.Query_taxamap)

        for query in tqdm(align_taxa, desc="Parsing read mapping to taxonomic references"):
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

    def _parse_reads_ARG_bam(self, tid_arg=90.0, tc_arg=80.0):
        
        ARG_readmapping = defaultdict(lambda: defaultdict(lambda: self.Query_ARGmap))
        align_amr = pysam.AlignmentFile(self.ARG_amr_bam, "rb")
        align_amr_refseq = {ref["SN"]: ref["LN"] for ref in align_amr.header.to_dict()["SQ"]}

        for alignment in tqdm(align_amr, desc="Parsing read mapping to ARG database"):
            n = alignment.query_name
            if alignment.reference_name is None:
                continue
            
            allele = alignment.reference_name
            arg_name = extract_arg_name(alignment.reference_name)
            l= alignment.query_length
            cigarEQ = alignment.get_cigar_stats()[0][7]
            Template2_l = align_amr_refseq.get(allele, None)
            start = alignment.query_alignment_start
            end = alignment.query_alignment_end
            if 100 * alignment.reference_length < tc_arg or 100 * cigarEQ / alignment.reference_length < tid_arg:
                continue

            ARG_readmapping[n][(arg_name, allele)] = self.Query_ARGmap(start, end, l, cigarEQ, Template2_l)

        return ARG_readmapping  


    def match(self):

        self.AMRlinks = defaultdict(lambda: defaultdict(lambda: self.Query_taxamap_ARGmap))

        taxa_mapping = self._parse_taxa()
        ARG_mapping = self._parse_reads_ARG_bam()

        for n in ARG_mapping.keys():
            
            if n in taxa_mapping.keys():
                for allele in ARG_mapping[n].keys():
                    self.AMRlinks[n][allele] = self.Query_taxamap_ARGmap(*(taxa_mapping[n] + ARG_mapping[n][allele]))
            else:
                for allele in ARG_mapping[n].keys():
                    taxa_mapping[n] = self.Query_taxamap(None, ARG_mapping[n][allele].Query_Length_ARG, None, None, None, None, None)
                    self.AMRlinks[n][allele] = self.Query_taxamap_ARGmap(*(taxa_mapping[n] + ARG_mapping[n][allele]))
        
        # add contigs/reads not mapped to any AMR gene
        for n in taxa_mapping.keys():

            if n in ARG_mapping.keys():
                continue
            else:
                arg = None
                allele = None
                ARG_mapping[n][(arg, allele)] = self.Query_ARGmap(None, None, None, None, None)
                self.AMRlinks[n][(arg, allele)] = self.Query_taxamap_ARGmap(*(taxa_mapping[n] + ARG_mapping[n][(arg, allele)]))

        flattened_data = []
        for query, arg_dict in self.AMRlinks.items():
            for (arg_name, allele), values in arg_dict.items():
                flattened_data.append((query, arg_name, allele) + values)

        columns = ['Query', 'ARG', 'Allele'] + list(self.Query_taxamap_ARGmap._fields)
        self.AMRlinks = pd.DataFrame(flattened_data, columns=columns)
        
        # when analyzing read-based AMR detection, group results per ARG and taxonomic tempalte
        if self.mode == "reads":
            grouped = self.AMRlinks.groupby(["Template", "Template_Length", "ARG", "Allele", "ARG_Length"], dropna=False)
                
            n_links = grouped.size()
            Total_Query_Length = grouped['Query_Length_ARG'].sum()
            Observed_Templ_Depth = grouped['Observed_Templ_Length'].sum()
            ExactMatch_Taxa = grouped['ExactMatch_bp_taxa'].sum()
            ExactMatch_ARG = grouped['ExactMatch_bp_ARG'].sum()
            # start and end pos of mapping are dropped for reads --> multiple covered intervals on templates = multiple start and end points

            self.AMRlinks = pd.DataFrame({'n_links': n_links, 'Total_Query_Length': Total_Query_Length, 'Observed_Templ_Depth' : Observed_Templ_Depth, 'ExactMatch_Taxa' : ExactMatch_Taxa, 'ExactMatch_ARG' : ExactMatch_ARG})
            self.AMRlinks = self.AMRlinks.reset_index()


    def get_covinfo(self, coverage_file):

        if coverage_file.endswith(".bam"):
            cov_output = pysam.coverage(coverage_file)
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

        elif coverage_file.endswith(".res"):
            coverage_data = pd.read_csv(coverage_file, sep="\t").drop(["q_value", "p_value", "Score", "Expected"], axis=1)
            coverage_data = coverage_data.add_suffix("_kma_res")
            coverage_data = coverage_data.rename(columns={"#Template_kma_res": "rname"})

        # if in assembly mode, we want coverage of contigs (=queries)
        if self.mode == "assembly":
            merged_df = self.AMRlinks.merge(coverage_data, left_on="Query", right_on="rname", how="outer")
            self.AMRlinks = merged_df
            self.AMRlinks["Query"] = self.AMRlinks["Query"].fillna(self.AMRlinks["rname"])
            self.AMRlinks = self.AMRlinks.drop("rname", axis=1)

        # in read mode, match to the taxonomic template mapping
        elif self.mode == "reads":
            merged_df = self.AMRlinks.merge(coverage_data, left_on="Template", right_on="rname", how="outer")
            self.AMRlinks = merged_df
            # in case res file is used instead of bam file (will be faster as res summarizes mapping results)
            if coverage_file.endswith(".res"):
                self.AMRlinks["Template"] = self.AMRlinks["Template"].fillna(self.AMRlinks["rname"])
                self.AMRlinks["Template_Length"] = self.AMRlinks["Template_Length"].fillna(self.AMRlinks["Template_length_kma_res"])
                self.AMRlinks = self.AMRlinks.drop(["rname", "Template_length_kma_res"], axis= 1)

    # TODO: split the class in two: create subclasses with assembly/read-specific methods which inherit from AMRlinker class
    def get_genomad_info(self):
        if self.genomad_dir is None:
            raise ValueError("Genomad output directory is not available.")
        if self.mode == "reads":
            raise ValueError("Genomad not implemented for reads")
        
        plasmid_sum_dir = self.genomad_dir.rstrip("/") + "/consensus_summary/consensus_plasmid_summary.tsv"
        virus_sum_dir = self.genomad_dir.rstrip("/") + "/consensus_summary/consensus_virus_summary.tsv"

        plasmid_sum = pd.read_csv(plasmid_sum_dir, sep="\t")
        plasmid_sum = plasmid_sum.drop("length", axis=1).add_suffix("_plasmid")
        self.genomad_plasmids = plasmid_sum

        virus_sum =pd.read_csv(virus_sum_dir, sep="\t").drop("length", axis=1).add_suffix("_virus")
        # split seq_name_virus from contig, take start and end position for provirus sequences
        virus_sum[["contig_name", "seq_name_virus"]] = virus_sum["seq_name_virus"].str.split("|", expand=True)
        virus_sum[["seq_name_virus", "virus_startp", "virus_endp"]] = virus_sum["seq_name_virus"].str.split("_", expand=True)
        self.genomad_virus = virus_sum


    def get_scapp_info(self, scapp_fp):
        if self.mode == "reads":
            raise ValueError("SCAPP not implemented for reads")
        scapp_paths = scapp_fp.rstrip("/") + "/intermediate_files/assembly_graph.cycs.paths.txt"
        scapp_plasmids = scapp_fp.rstrip("/") + "/assembly_graph.confident_cycs.fasta"
        
        lencov_pattern = re.compile(r"_length_\d+_cov_\d+.+")

        true_plasmids = set()
        with open(scapp_plasmids, 'r') as scapp_plasmidfile:
            for line in scapp_plasmidfile:
                if line.startswith(">"):
                    plasmid_name = line.strip(">\n")
                    plasmid_name = re.sub(lencov_pattern, "", plasmid_name)
                    true_plasmids.add(plasmid_name)

        len_pattern = re.compile(r"length_\d+")
        with open(scapp_paths, 'r') as scapp_pathfile:
            all_lines = scapp_pathfile.readlines()
            plasmid_names = [line.rstrip() for line in all_lines[0::3]]
            plasmid_lengths = [int(re.search(len_pattern, name).group().lstrip("length_")) for name in plasmid_names]
            plasmid_names = [re.sub(lencov_pattern, "", name) for name in plasmid_names]

        self.SCAPP_contigs = []
        for i in range(len(plasmid_names)):
            ctgs = {"contig_" + str(ctg.strip()) for ctg in all_lines[2::3][i].strip("[]\n").split(sep=",") if plasmid_lengths[i] < 1000000}

            if plasmid_names[i] in true_plasmids:
                for ctg in ctgs:
                    self.SCAPP_contigs.append([ctg, plasmid_names[i]])
            else:
                continue

        self.SCAPP_contigs = pd.DataFrame(self.SCAPP_contigs, columns=["SCAPP_contig", "SCAPP_plasmid"])


    def get_binning_info(self, binning_fp):
        if binning_fp is None:
            raise ValueError("Binning results tsv is not available.")
        if self.mode == "reads":
            raise ValueError("Binning not implemented for reads")

        binning_res = pd.read_csv(binning_fp, sep="\t", names=["Contig", "Bin"])

        merged_df_bin = self.AMRlinks.merge(binning_res, left_on="Query", right_on="Contig", how="outer")
        self.AMRlinks = merged_df_bin


    def get_bin_classification_kma(self):
        if self.mode == "reads":
            raise ValueError("Binning not implemented for reads")

        if "Bin" in self.AMRlinks.columns:         
            
            bin_taxa = self.AMRlinks[["Query", "Template", "Query_Length_taxa", "ExactMatch_bp_taxa", "Bin"]].dropna(subset="Bin")
            bin_taxa = bin_taxa[(bin_taxa["Bin"] != "") & (bin_taxa["Bin"] != None) & (bin_taxa["Bin"] != "unbinned")]

            taxid_pattern = re.compile(r"kraken:taxid\|(\d+)")
            def get_taxid(template):
                if isinstance(template, str):
                    taxid = re.search(taxid_pattern, template).group(1)
                    return taxid
                else:
                    return None
            
            bin_taxa["TaxID"] = bin_taxa["Template"].apply(get_taxid)
            Query_Length_taxa_sum = bin_taxa.groupby('Bin')["Query_Length_taxa"].sum().reset_index()
            Query_Length_taxa_sum = Query_Length_taxa_sum.rename(columns={'Query_Length_taxa': 'sum_Query_Length_taxa'})

            # FIXME: some species have multiple Taxids (e.g. Ecoli), this is not properly adressed here
            exactmatch_sum = bin_taxa.groupby(['Bin', 'TaxID'])['ExactMatch_bp_taxa'].sum().reset_index()
            exactmatch_sum = exactmatch_sum.rename(columns={'ExactMatch_bp_taxa': 'sum_ExactMatch_bp_taxa'})

            bin_taxa = pd.merge(bin_taxa, Query_Length_taxa_sum, on='Bin', how='left')
            bin_taxa = pd.merge(bin_taxa, exactmatch_sum, on=['Bin', 'TaxID'], how='left')
            bin_taxa["pid_taxa_query"] = bin_taxa["sum_ExactMatch_bp_taxa"]/bin_taxa["sum_Query_Length_taxa"]
            bin_taxa = bin_taxa.dropna(subset="pid_taxa_query")

            # bin_taxa = bin_taxa[(bin_taxa["pid_taxa_query"] >= 0.75) & (bin_taxa["sum_Query_Length_taxa"] >= 200000)]
            bin_taxa.to_csv(sys.stdout, sep="\t")
            bin_taxa = bin_taxa.loc[bin_taxa.groupby('Bin')['pid_taxa_query'].idxmax()]

            taxids = bin_taxa["TaxID"].dropna().to_list() 
            species_names = pytaxonkit.lineage(taxids, data_dir=self.taxonkit_db)[["TaxID", "Name"]]
            species_names.columns = ["TaxID", "Taxon"]
            species_names["TaxID"] = species_names["TaxID"].astype(str)

            bin_taxa = bin_taxa.merge(species_names, on="TaxID", how="left")
            # bin_taxa.to_csv(sys.stdout, sep="\t")

            self.AMRlinks = pd.merge(self.AMRlinks, bin_taxa[["Bin", "Taxon"]], on="Bin", how="left")

        else:
            # TODO: to first execute get_binning_information if binning info not yet present
            pass
    

    def get_bin_classification_gtdb(self,gtdb_bac_summary_fp):
        if self.mode == "reads":
            raise ValueError("Binning not implemented for reads")

        if "Bin" in self.AMRlinks.columns:         
            self.gtdb_classification = pd.read_csv(gtdb_bac_summary_fp, sep="\t")

            prefix_pattern = re.compile(r"[a-z]__")
            suffix_pattern = re.compile(r"_[A-Z]")
            def get_gtdb_taxon(classif):
                taxon = classif.split(";")[-1]
                taxon = re.sub(prefix_pattern, "", taxon)
                taxon = re.sub(suffix_pattern, "", taxon)
                if not taxon.startswith("Unclassified"):
                    return taxon
                else:
                    return None
            
            self.gtdb_classification["Taxon"] = self.gtdb_classification["classification"].apply(get_gtdb_taxon)
            self.gtdb_classification["Bin"] = self.gtdb_classification["user_genome"]
            self.AMRlinks = pd.merge(self.AMRlinks, self.gtdb_classification[["Bin", "Taxon"]], on="Bin", how="left")


    def get_nanomotif_include(self, include_contig_fp):
        """
        Required for retrieving methylation-based linking of contigs(plasmids) to given bin
        """
        if include_contig_fp is None:
            raise ValueError("nanomotif include_contig.tsv is not available.")
        if self.mode == "reads":
            raise ValueError("Binning not implemented for reads")

        self.nanomotif_links = pd.read_csv(include_contig_fp, sep="\t")
        self.nanomotif_links = self.nanomotif_links[self.nanomotif_links["contig_bin"] == "unbinned"]
        self.nanomotif_links = self.nanomotif_links[["bin", "contig"]]
        self.nanomotif_links.columns = ["Bin", "Contig"]


    def get_mobileOGdb_info(self, mobileOGdb_fp):
        if mobileOGdb_fp is None:
            raise ValueError("Binning results tsv is not available.")
        if self.mode == "reads":
            raise ValueError("Binning not implemented for reads")
        
        self.mobileOGdb = pd.read_csv(mobileOGdb_fp)
        # mobileOGdb = pd.read_csv(mobileOGdb_fp)
        # self.AMRlinks_MGEs = self.AMRlinks.merge(mobileOGdb, left_on="Query", right_on="Specific Contig", how="outer")
    

    def create_graph_ARGlinks(self):
        """
        Create graph of ARG gene sharing network:
            - if ARG on binned contig --> link to Bin, indicate ARG as plasmid, virus (if within pro(virus) sequence) or chromosomal node
            - if ARG not on binned contig, but on genomad/SCAPP plasmid contig --> link to either genomad contig or SCAPP plasmid
            - else: ARG present on unbinned chromosomal contig
        TODO: run helper methods first
        TODO: Split of plasmid or (fully!) viral contigs from the bins, create edge to bin instead
        """
        # data cleanup and formatting
        AMRlinks_temp = self.AMRlinks
        AMRlinks_temp["Bin"] = np.where((AMRlinks_temp["Bin"] == "unbinned") | (AMRlinks_temp["Bin"].isna()), AMRlinks_temp["Query"], AMRlinks_temp["Bin"])
        
        # Merge AMRlinks with plasmid information and determine if ARG is on an MGE
        # 1: Identify and classify plasmid-associated ARGs
        self.ARG_plasmids = AMRlinks_temp.merge(self.genomad_plasmids, left_on="Query", right_on="seq_name_plasmid", how="right")
        self.ARG_plasmids["ARG_on_MGE"] = np.where(self.ARG_plasmids["ARG"].isna(), 0, 1)
        self.ARG_plasmids["MGE"] = "Plasmid"
        self.ARG_plasmids["Node"] = self.ARG_plasmids["ARG"] + "_plasmid"
        self.ARG_plasmids = self.ARG_plasmids[["Node", "ARG", "Bin", "Query", "ARG_on_MGE", "MGE", "Taxon"]]

        # Remove plasmid-associated ARGs from the remaining set
        plasmid_ARGs = self.ARG_plasmids["Query"].to_list()
        AMRlinks_remaining = AMRlinks_temp[~AMRlinks_temp["Query"].isin(plasmid_ARGs)]

        # re-iterate process but now for plasmid contigs identified by scapp
        scapp_plasmid_ctgs = AMRlinks_remaining.merge(self.SCAPP_contigs, left_on="Query", right_on="SCAPP_contig", how="right")
        scapp_plasmid_ctgs["Bin"] = np.where(scapp_plasmid_ctgs["Bin"] == scapp_plasmid_ctgs["Query"], scapp_plasmid_ctgs["SCAPP_plasmid"], scapp_plasmid_ctgs["Bin"])
        scapp_plasmid_ctgs["ARG_on_MGE"] = np.where(scapp_plasmid_ctgs["ARG"].isna(), 0, 1)
        scapp_plasmid_ctgs["MGE"] = "Plasmid"
        scapp_plasmid_ctgs["Node"] = scapp_plasmid_ctgs["ARG"] + "_plasmid"
        scapp_plasmid_ctgs = scapp_plasmid_ctgs[["Node", "ARG", "Bin", "Query", "ARG_on_MGE", "MGE", "Taxon"]]
        self.ARG_plasmids = pd.concat([self.ARG_plasmids, scapp_plasmid_ctgs])

        plasmid_ARGs = self.ARG_plasmids["Query"].to_list()
        AMRlinks_remaining = AMRlinks_temp[~AMRlinks_temp["Query"].isin(plasmid_ARGs)]

        # 2: Identify and classify virus-associated ARGs
        self.ARG_virus = AMRlinks_remaining.merge(self.genomad_virus, left_on="Query", right_on="contig_name", how="right")
        self.ARG_virus[["Query_Align_Start_ARG", "virus_startp", "Query_Align_End_ARG", "virus_endp"]] = self.ARG_virus[["Query_Align_Start_ARG", "virus_startp", "Query_Align_End_ARG", "virus_endp"]].apply(pd.to_numeric, errors='coerce')
        
        # Determine if the ARG is located inside the viral sequence
        # TODO: modify here to label exclusively viral contigs (not proviruses)
        arg_in_virus = (
        (self.ARG_virus["Query_Align_Start_ARG"].fillna(float('-inf')) >= self.ARG_virus["virus_startp"].fillna(float('-inf'))) & 
        (self.ARG_virus["Query_Align_End_ARG"].fillna(float('inf')) <= self.ARG_virus["virus_endp"].fillna(float('inf')))
        )
        self.ARG_virus["ARG_on_MGE"] = np.where(arg_in_virus, 1, 0)
        self.ARG_virus["MGE"] = np.where(arg_in_virus, "Virus", "Chromosome")
        self.ARG_virus["Node"] = np.where(arg_in_virus, self.ARG_virus["ARG"] + "_virus", self.ARG_virus["ARG"] + "_chromosome") 
        self.ARG_virus = self.ARG_virus[["Node", "ARG", "Bin", "Query", "ARG_on_MGE", "MGE", "Taxon"]]

        # Remove virus-associated ARGs from the remaining set
        virus_ARGs = self.ARG_virus[self.ARG_virus["MGE"] == "Virus"]["Query"].to_list()
        AMRlinks_remaining = AMRlinks_remaining[~AMRlinks_remaining["Query"].isin(virus_ARGs)]

        # 3: Classify remaining ARGs as chromosomal
        self.ARG_remaining = AMRlinks_remaining.copy()
        self.ARG_remaining["ARG_on_MGE"] = 0
        self.ARG_remaining["MGE"] = "Chromosome"
        self.ARG_remaining["Node"] = self.ARG_remaining["ARG"] + "_chromosome"
        self.ARG_remaining = self.ARG_remaining[["Node", "ARG", "Bin", "Query", "ARG_on_MGE", "MGE", "Taxon"]]

        # 4: Concatenate all the processed DataFrames
        self.AMRlinks_simple = pd.concat([self.ARG_plasmids, self.ARG_virus, self.ARG_remaining])

        # 5: Create the graph from the simplified DataFrame
        AMRlinks_simple_tmp = self.AMRlinks_simple[self.AMRlinks_simple["ARG"].notna()]
        # Determine if ARGs occur on multiple MGEs and label them
        ARGs_with_multiple_mges = AMRlinks_simple_tmp.groupby("ARG")["MGE"].nunique()
        multi_mge_ARGs = ARGs_with_multiple_mges[ARGs_with_multiple_mges > 1].index.tolist()
        self.AMRlinks_graph = nx.Graph()

        # Add nodes for Bins and contigs
        for _, row in AMRlinks_simple_tmp.iterrows():
            # handle unbinned contigs --> no label
            if row['Bin'] == row['Query'] or row['Bin'].startswith("RNODE"):
                shape = 'o' if row['MGE'] == 'Plasmid' else ('D' if row['MGE'] == 'Virus' else 's')
                # label = row['Taxon'] if not pd.isna(row['Taxon']) else row['Bin']
                # self.AMRlinks_graph.add_node(row['Bin'], label=label, shape=shape)
                self.AMRlinks_graph.add_node(row['Bin'], label="", shape=shape)
            # add binned plasmid nodes --> no label
            elif row["MGE"] == 'Plasmid':
                self.AMRlinks_graph.add_node(row['Query'], label="", shape='o')
                label = row['Taxon'] if not pd.isna(row['Taxon']) else row['Bin']
                self.AMRlinks_graph.add_node(row['Bin'], label=label, shape='s')
            # add bins
            else:
                label = row['Taxon'] if not pd.isna(row['Taxon']) else row['Bin']
                self.AMRlinks_graph.add_node(row['Bin'], label=label, shape='s')

        # Add nodes for ARGs with their corresponding shapes
        for _, row in AMRlinks_simple_tmp.iterrows():
            if row['ARG'] in multi_mge_ARGs:
                shape = 'h'  # New shape for multiple MGEs (e.g., hexagon)
            else:
                shape = 'o' if row['MGE'] == 'Plasmid' else ('D' if row['MGE'] == 'Virus' else 's')
            self.AMRlinks_graph.add_node(row['ARG'], label=row['ARG'], shape=shape)

        # Add edges between Bins, plasmids and ARGs
        for _, row in AMRlinks_simple_tmp.iterrows():
            
            if row['MGE'] == 'Plasmid':
                self.AMRlinks_graph.add_edge(row['Query'], row['ARG'], MGE=row['MGE'])
                self.AMRlinks_graph.add_edge(row['Query'], row['Bin'], MGE=row['MGE']) if row['Query'] != row['Bin'] else None
            else:
                self.AMRlinks_graph.add_edge(row['Bin'], row['ARG'], MGE=row['MGE'])

        # # Add edges between ARGs on the same Query (contig)
        # for query, group in AMRlinks_simple_tmp.groupby('Query'):
        #     nodes = group['ARG'].tolist()
        #     for i in range(len(nodes)):
        #         for j in range(i + 1, len(nodes)):
        #             if nodes[i] != nodes[j] and nodes[i]:
        #                 self.AMRlinks_graph.add_edge(nodes[i], nodes[j])
        
    def plot_ARG_graph(self, output, seed=37, figsize=20):

        G = self.AMRlinks_graph
        # Calculate PageRank centrality for each node
        pagerank = nx.pagerank(G)

        # Normalize PageRank values to a suitable range for node sizes
        min_size = 300  # Minimum node size
        max_size = 3000  # Maximum node size
        node_sizes = np.array([pagerank[node] for node in G.nodes()])
        norm_node_sizes = min_size + (node_sizes - node_sizes.min()) * (max_size - min_size) / (node_sizes.max() - node_sizes.min())
        # Get positions for nodes
        # pos = nx.spring_layout(G, k=figsize/(4*np.sqrt(G.order())), iterations=120, scale=1.5, seed=seed)
        pos = nx.nx_agraph.graphviz_layout(G, prog='neato', args='-Goverlap=prism  -Gsplines=true')

        # Separate nodes by shape
        node_styles = {
            'plasmid': {
                'nodes': [node for node, attr in G.nodes(data=True) if attr['shape'] == 'o'],
                'color': 'red',
                'shape': 'o'
            },
            'virus': {
                'nodes': [node for node, attr in G.nodes(data=True) if attr['shape'] == 'D'],
                'color': 'green',
                'shape': 'D'
            },
            'chromosome': {
                'nodes': [node for node, attr in G.nodes(data=True) if attr['shape'] == 's'],
                'color': 'blue',
                'shape': 's'
            },
            'multiple_mge': {
                'nodes': [node for node, attr in G.nodes(data=True) if attr['shape'] == 'h'],
                'color': 'orange',
                'shape': 'h'
            }
        }

        # Draw nodes with the specified shapes
        plt.figure(figsize=(figsize, figsize))

        for style in node_styles.values():
            nx.draw_networkx_nodes(
                G, pos, 
                nodelist=style['nodes'], 
                node_shape=style['shape'], 
                label=style['shape'], 
                node_color=style['color'], 
                node_size=[norm_node_sizes[list(G.nodes).index(node)] for node in style['nodes']])

        # Draw edges
        nx.draw_networkx_edges(G, pos)

        # Draw labels using the 'label' attribute of nodes
        node_labels = {node: data['label'] for node, data in G.nodes(data=True)}
        nx.draw_networkx_labels(G, pos, labels=node_labels)

        # Add a legend
        plt.legend(scatterpoints=1)

        # Save the plot
        plt.savefig(output)