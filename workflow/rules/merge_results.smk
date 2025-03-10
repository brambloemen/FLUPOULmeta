output_dir = config.get("output_dir", "results")

rule merge_assembly_results:
    input:
        assemb_bam=f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_assembly.bam",
        ResF_bam=f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_ResF.bam",
        covinfo=f"{output_dir}/{{sample}}/MapToAssemb/{{sample}}_assembly.bam",
        genomad_outputdir=f"{output_dir}/{{sample}}/genomad_output",
        binning_results=f"{output_dir}/{{sample}}/NanoMotif/bin/new_contig_bin.tsv",
        scapp=f"{output_dir}/{{sample}}/SCAPP/assembly_graph.confident_cycs.fasta",
        gtdb=f"{output_dir}/{{sample}}/GTDBtk_NM",
        taxo_mapping_bam=f"{output_dir}/{{sample}}/KMA_reads/{{sample}}.kmadb.bam",
        reads_resf_bam=f"{output_dir}/{{sample}}/KMA_reads/{{sample}}_rep.resf.bam",
        reads_taxo_res=f"{output_dir}/{{sample}}/KMA_reads/{{sample}}_kma"
    output:
        f"{output_dir}/{{sample}}/Summary/Assembly_AMRlinks.tsv",
        f"{output_dir}/{{sample}}/Summary/Assembly_AMRlinks_summary.tsv",
        f"{output_dir}/{{sample}}/Summary/Assembly_AMRlinks.png",
        f"{output_dir}/{{sample}}/Summary/{{sample}}_Assembly_AMRlinks.xml",
        f"{output_dir}/{{sample}}/Summary/Reads_AMRlinks.tsv"
    log: f"logs/merge_assemb_results/{{sample}}.log"
    resources:
        cpus_per_task=1,
        mem_mb=min(30000, config['memory'])
    conda:
        path.join(base_dir, "envs/FLUPOUL.yaml")
    script:
        path.join(base_dir, "scripts/Merge_results_assembly.py")
