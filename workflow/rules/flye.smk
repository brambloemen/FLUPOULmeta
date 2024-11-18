flye=TOOLS["Flye"]["path"]
minimap2=TOOLS["minimap2"]["path"]
samtools=TOOLS["Samtools"]["path"]
htslib=TOOLS["Medaka"]["htslib_path"]
medaka=TOOLS["Medaka"]["path"]

rule flye_assemble:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz"
    output:
        fasta="results/{sample}/Flye/assembly.fasta",
        graph="results/{sample}/Flye/assembly_graph.gfa",
        info="results/{sample}/Flye/assembly_info.txt"
    params:
        "--meta -i 2 --read-error 0.03"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/flye/{sample}.log"
    conda:
        "../envs/Flye.yaml"
    shell:
        """
        export PATH={flye}:$PATH
        flye --nano-hq {input.fastq_proka} {params} --out-dir results/{wildcards.sample}/Flye -t {resources.cpus_per_task}
        """

rule medaka_mini_align:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz",
        draft="results/{sample}/Flye/assembly.fasta"
    output:
        bam="results/{sample}/Medaka/calls_to_draft.bam",
        bai="results/{sample}/Medaka/calls_to_draft.bam.bai"
    params:
        bam=temp("results/{sample}/Medaka/calls_to_draft.tmp")
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
        """
        export PATH={minimap2}:{samtools}:{htslib}:{medaka}:$PATH
        source {medaka}/activate
        mini_align -i {input.fastq_proka} -r {input.draft} -m -p {params.bam} -t {resources.cpus_per_task}
        samtools sort {params.bam}.bam -@{resources.cpus_per_task} -o {output.bam}
        samtools index {output.bam}
        """

# use no more than 6 threads, as more doesn't improve performance of medaka consensus
cpus=min(config['threads'],12)
rule medaka_consensus:
    input:
        bam="results/{sample}/Medaka/calls_to_draft.bam",
        info="results/{sample}/Flye/assembly_info.txt"
    output:
        "results/{sample}/Medaka/consensus_probs.hdf"
    resources:
        cpus_per_task=cpus,
        mem_mb=config['memory']
    params:
        model=config['model']
    shell:
        """
        export PATH={minimap2}:{samtools}:{htslib}:{medaka}:$PATH
        source {medaka}/activate
        #only polish regions >9999bp, otherwise the consensus stage is slowed down significantly
        REGIONS=$(awk 'NR>1 && $2> 9999 {{print $1}}' {input.info} | tr '\n' ' ')
        medaka inference {input.bam} {output} --model {params.model} --batch 200 --threads {resources.cpus_per_task} --region $REGIONS
        """

rule medaka_stitch:
    input:
        probs="results/{sample}/Medaka/consensus_probs.hdf",
        draft="results/{sample}/Flye/assembly.fasta"
    output:
        "results/{sample}/Medaka/consensus.fasta"
    resources:
        cpus_per_task=1,
        mem_mb=24000
    shell:
        """
        export PATH={minimap2}:{samtools}:{htslib}:{medaka}:$PATH
        source {medaka}/activate
        medaka sequence {input.probs} {input.draft} {output}
        """
