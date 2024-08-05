rule merge_assembly_results:
    input:
        assemb_bam="results/{sample}/FlyeClassify/{sample}_assembly.bam",
        NDARO_bam="results/{sample}/FlyeClassify/{sample}_NDARO.bam",
        covinfo="results/{sample}/MapToAssemb/{sample}_assembly.bam",
        genomad_outputdir="results/{sample}/genomad_output",
        binning_results="results/{sample}/NanoMotif/bin/new_contig_bin.tsv"
    output:
        "results/{sample}/Summary/Assembly_AMRlinks.tsv"
    log: "logs/merge_assemb_results/{sample}.log"
    resources:
        cpus_per_task=1,
        mem_mb=100000
    # shell:
    #     """
    #     echo test
    #     """
    script:
        "../scripts/Merge_results_assembly.py"
