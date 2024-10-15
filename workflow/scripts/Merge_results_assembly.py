# from snakemake.script import snakemake
from utils import AMRlinker
import networkx as nx

try:
    # Verify paths
    print(f"Assembled BAM: {snakemake.input['assemb_bam']}")
    print(f"ResF BAM: {snakemake.input['ResF_bam']}")
    print(f"Genomad output directory: {snakemake.input['genomad_outputdir']}")
    print(f"Covinfo: {snakemake.input['covinfo']}")
    print(f"Binning results: {snakemake.input['binning_results']}")
    print(f"SCAPP results: {snakemake.input['scapp']}")
    print(f"gtdb results: {snakemake.input['gtdb']}")
    print(f"Output file: {snakemake.output[0]}, {snakemake.output[1]}")

    # Initialize the AMRlinker object
    test_parser = AMRlinker(snakemake.input["assemb_bam"], snakemake.input["ResF_bam"], genomad_dir=snakemake.input["genomad_outputdir"])

    # Perform matching
    test_parser.match()

    # Get coverage information
    test_parser.get_covinfo(snakemake.input["covinfo"])
    
    # Get Genomad information
    test_parser.get_genomad_info()
    
    # Get binning information
    test_parser.get_binning_info(snakemake.input["binning_results"])

    # Scapp results
    scapp_fp = snakemake.input['scapp'].rstrip("/assembly_graph.confident_cycs.fasta")
    test_parser.get_scapp_info(scapp_fp)

    # classification
    test_parser.get_bin_classification_gtdb(f"{snakemake.input['gtdb']}/gtdbtk.bac120.summary.tsv")

    # Output to CSV
    test_parser.AMRlinks.to_csv(snakemake.output[0], sep="\t")

    test_parser.create_graph_ARGlinks()

    test_parser.plot_ARG_graph(snakemake.output[1], figsize=30)
    nx.write_graphml(test_parser.AMRlinks_graph, snakemake.output[2])

    read_parser = AMRlinker(snakemake.input["taxo_mapping_bam"], snakemake.input["reads_resf_bam"], mode="reads")
    read_parser.match()
    read_parser.get_covinfo(snakemake.input["reads_resf_res"])  
    read_parser.AMRlinks.to_csv(snakemake.output[3], sep="\t")  


except Exception as e:
    print(f"An error occurred: {e}")