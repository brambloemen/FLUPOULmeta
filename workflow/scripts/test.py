from utils import AMRlinker
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
import pygraphviz

assemb_bam="/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/FlyeClassify/F4D4_assembly.bam"
ResF_bam="/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/FlyeClassify/F4D4_ResF.bam"
mobileOG_fp = "/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/mobileOGdb/DASTool_proteins.mobileOG.Alignment.Out.csv"
binning_fp = "/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/NanoMotif/bin/new_contig_bin.tsv"
scapp_fp = "/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/SCAPP"
gtdb_fp = "/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/GTDBtk/gtdbtk.bac120.summary.tsv"

test_parser = AMRlinker(assemb_bam, ResF_bam, genomad_dir = "/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/genomad_output")

test_parser.match()

# test_parser.get_covinfo("/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/MapToAssemb/smallsample.bam")
test_parser.get_binning_info(binning_fp)
# test_parser.AMRlinks.to_csv("/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/test_AMRlink_summary.tsv", sep="\t")

test_parser.get_genomad_info()
# print(test_parser.genomad_plasmids.head(5))
# print(test_parser.genomad_virus.head(5))


test_parser.get_mobileOGdb_info(mobileOG_fp)

test_parser.get_scapp_info(scapp_fp)
# test_parser.get_bin_classification_kma()
test_parser.get_bin_classification_gtdb(gtdb_fp)
test_parser.create_graph_ARGlinks()
# print(test_parser.ARG_virus[test_parser.ARG_virus["Allele"] != ""].head(5))

test_parser.AMRlinks.to_csv("/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/test_AMRlink.tsv", sep="\t")
test_parser.AMRlinks_simple.to_csv("/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/test_AMRlink_simple.tsv", sep="\t")

print(test_parser.AMRlinks_graph.order())
test_parser.plot_ARG_graph('test.png', figsize=30)

# # test read-based ARG linking
# readmap_bam_taxa="/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/F4D4.kmadb.bam"
# readmap_NDARO_bam="/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/NDARO/F4D4_NDARO.bam"
# test_parser = AMRlinker(readmap_bam_taxa, readmap_NDARO_bam, mode="reads", genomad_dir="/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/genomad_output")
# # test_parser = AMRlinker(assemb_bam, NDARO_bam, mode="reads")

# test_parser.match()

# print(test_parser.AMRlinks)

# test_parser.get_covinfo("/scratch/brbloemen/FLUPOUL/FLUPOULmeta/workflow/results/F4D4/F4D4_kma.res")
# print(test_parser.AMRlinks.head(10))



