# from snakemake.script import snakemake
from utils import AMRlinker

try:
    # Verify paths
    print(f"Assembled BAM: {snakemake.input['assemb_bam']}")
    print(f"NDARO BAM: {snakemake.input['NDARO_bam']}")
    print(f"Genomad output directory: {snakemake.input['genomad_outputdir']}")
    print(f"Covinfo: {snakemake.input['covinfo']}")
    print(f"Binning results: {snakemake.input['binning_results']}")
    print(f"Output file: {snakemake.output[0]}")

    # Initialize the AMRlinker object
    test_parser = AMRlinker(snakemake.input["assemb_bam"], snakemake.input["NDARO_bam"], genomad_dir=snakemake.input["genomad_outputdir"])

    # Perform matching
    test_parser.match()

    # Get coverage information
    test_parser.get_covinfo(snakemake.input["covinfo"])
    
    # Get Genomad information
    test_parser.get_genomad_info()
    
    # Get binning information
    test_parser.get_binning_info(snakemake.input["binning_results"])

    # Output to CSV
    test_parser.AMRlinks.to_csv(snakemake.output[0], sep="\t")

except Exception as e:
    print(f"An error occurred: {e}")