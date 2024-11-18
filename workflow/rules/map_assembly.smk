minimap2=TOOLS["minimap2"]["path"]
samtools=TOOLS["Samtools"]["path"]

rule minimap_assembly:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz",
        assembly="results/{sample}/Medaka/consensus.fasta"
    output:
        sam=temp("results/{sample}/MapToAssemb/{sample}_assembly.sam")
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/map_assembly/{sample}_minimap.log"
    shell:
        """
        export PATH={minimap2}:{samtools}:$PATH
        minimap2 -ax map-ont -y -t {resources.cpus_per_task} {input.assembly} {input.fastq_proka} > {output.sam}
        """

rule samtools_assembly_mapping:
    input:
        sam="results/{sample}/MapToAssemb/{sample}_assembly.sam",
        assembly="results/{sample}/Medaka/consensus.fasta"
    output:
        bam="results/{sample}/MapToAssemb/{sample}_assembly.bam",
        bai="results/{sample}/MapToAssemb/{sample}_assembly.bam.bai",
        fai="results/{sample}/Medaka/consensus.fasta.fai"
    params:
        filt="-hb -F0x904"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/map_assembly/{sample}_samtools.log"
    shell:
        """
        export PATH={samtools}:$PATH
        samtools view {params.filt} -@{resources.cpus_per_task} {input.sam} | samtools sort -@{resources.cpus_per_task} > {output.bam}
        samtools index -@{resources.cpus_per_task} {output.bam} -o {output.bai}
        samtools faidx {input.assembly} -o {output.fai}
        """