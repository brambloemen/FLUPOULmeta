minimap2=TOOLS["minimap2"]["path"]
samtools=TOOLS["Samtools"]["path"]
output_dir=config.get("output_dir", "results")
skip_kraken2=config.get('skip_kraken2', False)

rule map_reads_to_graph:
    input:
        fastq_proka=fastq_proka,
        gfa=f"{output_dir}/{{sample}}/Flye/assembly_graph.gfa"
    output:
        bam=temp(f"{output_dir}/{{sample}}/SCAPP/{{sample}}_graphmap.bam"),
        graph=f"{output_dir}/{{sample}}/SCAPP/assembly_graph.fastg"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: f"logs/SCAPP/{{sample}}_mapping.log"
    params:
        script=path.join(base_dir, "scripts/metaflye_gfa2fastg.py")
    conda:
        path.join(base_dir, "envs/scapp.yaml")
    shell:
        """
        export PATH={minimap2}:{samtools}:$PATH
        python {params.script} {input.gfa} {output.graph} 2>{log}
        minimap2 -ax map-ont -t {resources.cpus_per_task} {output.graph} {input.fastq_proka} 2>>{log} | \
        samtools view -hb -F4 -@{resources.cpus_per_task} - | samtools sort -@{resources.cpus_per_task} > {output.bam} 2>>{log}
        """
#!!TODO: scapp uses unpolished gfa file (converted to fastg above); resulting plasmids are thus unpolished
#        Fix this by retrieving replacing unpolished contigs with the ones from Medaka/consensus.fasta where possible
rule scapp:
    input:
        graph=f"{output_dir}/{{sample}}/SCAPP/assembly_graph.fastg",
        bam=f"{output_dir}/{{sample}}/SCAPP/{{sample}}_graphmap.bam"
    output:
        plasmids=f"{output_dir}/{{sample}}/SCAPP/assembly_graph.confident_cycs.fasta",
        scapp=f"{output_dir}/{{sample}}/SCAPP/intermediate_files/assembly_graph.cycs.paths.txt"
    log: f"logs/SCAPP/{{sample}}_scapp.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    conda:
        path.join(base_dir, "envs/scapp.yaml")
    shell:
        """
        scapp -g {input.graph} -o {output_dir}/{wildcards.sample}/SCAPP -k 0 -b {input.bam} -p {resources.cpus_per_task} 2>{log}
        """