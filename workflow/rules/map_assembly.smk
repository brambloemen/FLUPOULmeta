rule minimap_assembly:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz",
        assembly="results/{sample}/Flye/assembly.fasta"
    output:
        sam=temp("results/{sample}/Flye/{sample}_assembly.sam")
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/map_assembly/{sample}_minimap.log"
    shell:
        """
        ml load minimap2/2.26
        minimap2 -ax map-ont -t {resources.cpus_per_task} {input.assembly} {input.fastq_proka} > {output.sam}
        """

rule samtools_assembly_mapping:
    input:
        sam="results/{sample}/Flye/{sample}_assembly.sam",
        assembly="results/{sample}/Flye/assembly.fasta"
    output:
        bam="results/{sample}/Flye/{sample}_assembly.bam",
        bai="results/{sample}/Flye/{sample}_assembly.bam.bai",
        fai="results/{sample}/Flye/assembly.fasta.fai"
    params:
        filt="-hb -F0x904"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/map_assembly/{sample}_samtools.log"
    shell:
        """
        ml load samtools/1.17
        samtools view {params.filt} -@{resources.cpus_per_task} {input.sam} | samtools sort -@{resources.cpus_per_task} > {output.bam}
        samtools index -@{resources.cpus_per_task} {output.bam} -o {output.bai}
        samtools faidx {input.assembly} -o {output.fai}
        """