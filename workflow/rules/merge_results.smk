rule merge_assembly_results:
    input:
        assemb_bam="results/{sample}/FlyeClassify/{sample}_assembly.bam",
        ResF_bam="results/{sample}/FlyeClassify/{sample}_ResF.bam",
        covinfo="results/{sample}/MapToAssemb/{sample}_assembly.bam",
        genomad_outputdir="results/{sample}/genomad_output",
        binning_results="results/{sample}/NanoMotif/bin/new_contig_bin.tsv",
        scapp="results/{sample}/SCAPP/assembly_graph.confident_cycs.fasta",
        gtdb="results/{sample}/GTDBtk_NM"
    output:
        "results/{sample}/Summary/Assembly_AMRlinks.tsv",
        "results/{sample}/Summary/Assembly_AMRlinks.png"
    log: "logs/merge_assemb_results/{sample}.log"
    resources:
        cpus_per_task=1,
        mem_mb=min(30000, config['memory'])
    # shell:
    #     """
    #     echo test
    #     """
    script:
        "../scripts/Merge_results_assembly.py"
