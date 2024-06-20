rule map_reads_to_graph:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz",
        gfa="results/{sample}/Flye/assembly_graph.gfa"
    output:
        bam=temp("results/{sample}/SCAPP/{sample}_graphmap.bam"),
        graph="results/{sample}/SCAPP/assembly_graph.fastg"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/SCAPP/{sample}_mapping.log"
    conda:
        "../envs/scapp.yaml"
    shell:
        """
        ml load minimap2/2.26 samtools/1.17
        python scripts/metaflye_gfa2fastg.py {input.gfa} {output.graph} 2>{log}
        minimap2 -ax map-ont -t {resources.cpus_per_task} {output.graph} {input.fastq_proka} 2>>{log} | \
        samtools view -hb -F4 -@{resources.cpus_per_task} > {output.bam} 2>>{log}
        """

rule scapp:
    input:
        graph="results/{sample}/SCAPP/assembly_graph.fastg",
        bam="results/{sample}/SCAPP/{sample}_graphmap.bam"
    output:
        plasmids="results/{sample}/SCAPP/assembly_graph.confident_cycs.fasta",
        scapp="results/{sample}/SCAPP/intermediate_files/assembly_graph.cycs.paths.txt"
    log: "logs/SCAPP/{sample}_scapp.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    conda:
        "../envs/scapp.yaml"
    shell:
        """
        scapp -g {input.graph} -o results/{wildcards.sample}/SCAPP -k 0 -b {input.bam} -p {resources.cpus_per_task} 2>{log}
        """