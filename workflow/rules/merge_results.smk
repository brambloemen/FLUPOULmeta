rule merge_assembly_results:
    input:
        assemb_bam="results/{sample}/FlyeClassify/{sample}_assembly.bam",
        ResF_bam="results/{sample}/FlyeClassify/{sample}_ResF.bam",
        covinfo="results/{sample}/MapToAssemb/{sample}_assembly.bam",
        genomad_outputdir="results/{sample}/genomad_output",
        binning_results="results/{sample}/NanoMotif/bin/new_contig_bin.tsv",
        scapp="results/{sample}/SCAPP/assembly_graph.confident_cycs.fasta",
        gtdb="results/{sample}/GTDBtk_NM",
        taxo_mapping_bam="results/{sample}/{sample}.kmadb.bam",
        reads_resf_bam="results/{sample}/{sample}_rep.resf.bam",
        reads_resf_res="results/{sample}/{sample}_resf.res"
    output:
        "results/{sample}/Summary/Assembly_AMRlinks.tsv",
        "results/{sample}/Summary/Assembly_AMRlinks.png",
        "results/{sample}/Summary/{sample}_Assembly_AMRlinks.xml",
        "results/{sample}/Summary/Reads_AMRlinks.tsv"
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
