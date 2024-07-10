import sys, argparse, re, csv, time, logging, pysam, math
from datetime import datetime
from datetime import timedelta
from collections import defaultdict
import numpy as np
import pandas as pd
from Template_parsing import extract_species_name, extract_arg_name


logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')


def parse_arguments():
    parser = argparse.ArgumentParser(description='Retrieve alignments from alignment file 1 in alignment file 2; \nPrint results to standard output')
    parser.add_argument("-t", metavar="taxa.bam", dest="bam_taxa", type=str, help="Filepath of bam with alignment to taxa database")
#   parser.add_argument("-r", dest="taxa_res", type=str, help="Filepath of res file of KMA alignmnent to taxa database")
    parser.add_argument("-a", metavar="AMR.bam", dest="AMR_bam", type=str, help="Filepath of bam with alignments to resistance database")
    parser.add_argument("-o", metavar="output", dest="output", type=str, help="Name/filepath of output plots", default=None)
    parser.add_argument("--threads", type=int, help="Number of threads used to process SAM/BAM file", default=1)
    parser.add_argument("--rlen", type=float, help="Minimum read length", default=0.0)
    parser.add_argument("--tid_arg", type=float, help="Minimum Percentage Identity to AMR gene template", default=90.0)
    parser.add_argument("--tc_arg", type=float, help="Minimum template coverage for a arg template to be accepted", default=100.0)

    return parser



def parse_taxa(taxa_bam, threads):
    
    """
    Loop over alignment to taxa, filter out unwanted templates, determine start time and duration
    """
    
    logging.info(f"Parsing alignment to taxa: {taxa_bam}")
    align_taxa = pysam.AlignmentFile(taxa_bam, "rb", threads=threads)
    align_taxa_refseq = {ref["SN"]:ref["LN"] for ref in align_taxa.header.to_dict()["SQ"]}
    align_taxa_n_lines = 0

    # align_taxa dictionary: {read:[Template1, readlength, matching bases, Template1 total reference length]}
    align_taxa_reads = defaultdict(lambda: [str, 0, 0, 0])
    for read in align_taxa:
        align_taxa_n_lines += 1

        
        # store read identifiers
        n = read.query_name
        Template1 = read.reference_name
        l = read.query_length
        cigarEQ = read.get_cigar_stats()[0][7]
        Template1_l = align_taxa_refseq.get(Template1, None)

        align_taxa_reads[n] = [Template1, l, cigarEQ, Template1_l]

    logging.info(f"Done parsing alignment 1. Processed {align_taxa_n_lines} alignments")
    return align_taxa_reads
    

def parse_AMR(amr_bam, processed_taxa_reads, threads):
    """
    Loop over alignment to resistance database
    """


    # summary dictionary. structure: {hour: {AMR gene: {Species: [#reads, total readlength, #exactly matched bases, readlength mapped to AMR, exactly matched bases to AMR, AMR template length]}}}
    alignment_links = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0, 0, 0]))
    logging.info(f"Parsing alignment to AMR: {amr_bam}")
    align_amr = pysam.AlignmentFile(amr_bam, "rb", threads=threads)
    align_amr_refseq = {ref["SN"]:ref["LN"] for ref in align_amr.header.to_dict()["SQ"]}
    align_amr_n_lines = 0

    for read in align_amr:
        align_amr_n_lines += 1
        
        # read identifiers
        n = read.query_name
        # alignment stats
        # not aligned to AMR --> skip
        if read.reference_name == None: 
            continue
        

        Template2 = read.reference_name
        l = read.query_length
        cigarEQ = read.get_cigar_stats()[0][7] # EQ: exact sequence match (base in query matches base in template)
        Template2_l = align_amr_refseq.get(Template2, None)
        if 100*read.reference_length < args.tc_arg:
            continue
        if 100*(cigarEQ)/read.reference_length < args.tid_arg:
            continue


        if n in processed_taxa_reads:
            Template1, Template1_readl, Template1_cigarEQ, Template1_l = processed_taxa_reads.pop(n) # pop matched reads from amr_reads: unmapped reads will remain
            alignment_links[Template1][Template2][0] += 1
            alignment_links[Template1][Template2][1] += Template1_readl # will be total length of all reads aligned to both template 1 and 2
            alignment_links[Template1][Template2][2] += Template1_cigarEQ
            alignment_links[Template1][Template2][3] += cigarEQ 
            alignment_links[Template1][Template2][4] = Template1_l
            alignment_links[Template1][Template2][5] = Template2_l
        else:
            # read not in parsed reads => seq time not calculated


            alignment_links["Unmapped"][Template2][0] += 1
            alignment_links["Unmapped"][Template2][1] += l # will be total length of all reads aligned to template 2, with no match for temp
            alignment_links["Unmapped"][Template2][2] = None
            alignment_links["Unmapped"][Template2][3] += cigarEQ 
            alignment_links["Unmapped"][Template2][4] = None
            alignment_links["Unmapped"][Template2][5] = Template2_l

    logging.info(f"Done parsing alignment 2. Processed {align_amr_n_lines} alignments")

    return  alignment_links

def create_output(output, alignment_links):

    header = ["Species", "ARG", "n_reads", "Total_readlength", "n_match_bases1", "n_match_bases2", "Template1_length", "Template2_length"]
    # writer = csv.writer(sys.stdout, delimiter="\t")
    # writer.writerow(header)

    ARG_links = []

    for spec, arg_stats in alignment_links.items():
        for arg, stats in arg_stats.items():
            # writer.writerow([t, spec, arg] + stats)
            ARG_links.append([spec, arg] + stats)

    """
    summarize output
    """

    ARG_links = pd.DataFrame(ARG_links, columns=header)
    ARG_links["Species"] = ARG_links["Species"].apply(extract_species_name)
    ARG_links["ARG"] = ARG_links["ARG"].apply(extract_arg_name)
    
    grouped = ARG_links.groupby(['Species', 'ARG'])

    # Calculate other summary statistics
    n_reads = grouped['n_reads'].sum()
    Total_readlength = grouped['Total_readlength'].sum()
    n_match_bases1 = grouped['n_match_bases1'].sum()
    n_match_bases2 = grouped['n_match_bases2'].sum()
    Template1_length = grouped['Template1_length'].mean()
    Template2_length = grouped['Template2_length'].mean()
    Organism_QID = n_match_bases1/Total_readlength
    ARG_TID = n_match_bases2/(Template2_length * n_reads)

    ARG_links = pd.DataFrame({'n_reads': n_reads,
                        'Total_readlength': Total_readlength,
                        'n_match_bases1': n_match_bases1,
                        'n_match_bases2': n_match_bases2,
                        'Template1_length': Template1_length,
                        'Template2_length': Template2_length,
                        'Organism_QID': Organism_QID,
                        'ARG_TID': ARG_TID})

    ARG_links = ARG_links.reset_index()
    ARG_links.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == '__main__':
    logging.info(f"Process arguments, load and preprocess reads")
    args = parse_arguments().parse_args()

    processed_taxa_reads = parse_taxa(args.bam_taxa, args.threads)
        
    alignment_links = parse_AMR(args.AMR_bam, processed_taxa_reads, args.threads)

    create_output(args.output, alignment_links)
    logging.info(f"Done")
