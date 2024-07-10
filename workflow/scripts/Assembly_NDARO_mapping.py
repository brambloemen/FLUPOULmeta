import json
import pysam
import argparse, sys, csv, re
import pandas as pd
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Retrieve alignments from alignment file 1 in alignment file 2; \nPrint results to standard output')
    parser.add_argument("-t", metavar="taxa.bam", dest="bam_taxa", type=str, help="Filepath of bam with alignment to taxa database")
    parser.add_argument("-a", metavar="AMR.bam", dest="AMR_bam", type=str, help="Filepath of bam with alignments to resistance database")
    parser.add_argument("--tid_arg", type=float, help="Minimum Percentage Identity to AMR gene template", default=90.0)
    parser.add_argument("--tc_arg", type=float, help="Minimum template coverage for a arg template to be accepted", default=100.0)
    return parser

def load_mapping(filepath):
    with open(filepath, 'r') as file:
        return json.load(file)

def parse_taxa(taxa_bam):
    align_taxa = pysam.AlignmentFile(taxa_bam, "rb")
    align_taxa_refseq = {ref["SN"]: ref["LN"] for ref in align_taxa.header.to_dict()["SQ"]}
    align_taxa_contigs = {}

    for contig in align_taxa:
        n = contig.query_name
        contig_l = contig.query_length
        Template1 = contig.reference_name
        l = contig.reference_length
        cigarEQ = contig.get_cigar_stats()[0][7]
        Template1_l = align_taxa_refseq.get(Template1, None)

        align_taxa_contigs[n] = [Template1, contig_l, l, cigarEQ, Template1_l]

    return align_taxa_contigs

def parse_AMR(amr_bam, processed_taxa_reads, tid_arg, tc_arg, mapping):
    rows = []
    pattern = r"seq_\d{6}$"
    align_amr = pysam.AlignmentFile(amr_bam, "rb")
    align_amr_refseq = {ref["SN"]: ref["LN"] for ref in align_amr.header.to_dict()["SQ"]}

    for alignment in align_amr:
        n = alignment.query_name
        if alignment.reference_name is None:
            continue

        Template2 = alignment.reference_name
        cigarEQ = alignment.get_cigar_stats()[0][7]
        Template2_l = align_amr_refseq.get(Template2, None)
        if 100 * alignment.reference_length < tc_arg or 100 * cigarEQ / alignment.reference_length < tid_arg:
            continue

        if n in processed_taxa_reads:
            Template1, contig_l, Aligned_Contig_l, Template1_cigarEQ, Template1_l = processed_taxa_reads[n]
            arg_match = re.search(pattern, Template2)
            if arg_match:
                arg = arg_match.group()
                arg_name = mapping[arg]["gene"]
                allele = mapping[arg]["allele"]
                rows.append([n, Template1, arg_name, allele, contig_l, Aligned_Contig_l, Template1_cigarEQ, cigarEQ, Template1_l, Template2_l])
        else:
            arg_match = re.search(pattern, Template2)
            if arg_match:
                arg = arg_match.group()
                arg_name = mapping[arg]["gene"]
                allele = mapping[arg]["allele"]
                rows.append([n, "Unmapped", arg_name, allele, contig_l, None, None, cigarEQ, None, Template2_l])

    df = pd.DataFrame(rows, columns=["Contig", "Taxonomic_mapping", "ARG", "Allele", "contig_length",
                                     "l_alignment_taxonomic", "n_match_bases_taxonomic", 
                                     "n_match_bases_ARG", "Taxonomic_ref_length", "ARG_ref_length"])
    return df

def main():
    args = parse_arguments().parse_args()
    mapping = load_mapping('/db/gene_detection/NCBI_AMR/mapping_full.json')
    
    processed_taxa_reads = parse_taxa(args.bam_taxa)
    df = parse_AMR(args.AMR_bam, processed_taxa_reads, args.tid_arg, args.tc_arg, mapping)

    df.sort_values(by="Contig", inplace=True)
    
    df.to_csv(sys.stdout, sep="\t", index=False)

if __name__ == '__main__':
    main()