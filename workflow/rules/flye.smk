flye=TOOLS["Flye"]["path"]
minimap2=TOOLS["minimap2"]["path"]
samtools=TOOLS["Samtools"]["path"]
htslib=TOOLS["Medaka"]["htslib_path"]
medaka=TOOLS["Medaka"]["path"]
output_dir=config.get("output_dir", "results")
skip_kraken2=config.get('skip_kraken2', False)

rule flye_assemble:
    input:
        fastq_proka=fastq_proka
    output:
        fasta=f"{output_dir}/{{sample}}/Flye/assembly.fasta",
        graph=f"{output_dir}/{{sample}}/Flye/assembly_graph.gfa",
        info=f"{output_dir}/{{sample}}/Flye/assembly_info.txt"
    params:
        "--meta -i 2 --read-error 0.03"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: f"logs/flye/{{sample}}.log"
    conda:
        path.join(base_dir, "envs/Flye.yaml")
    shell:
        """
        export PATH={flye}:$PATH
        flye --nano-hq {input.fastq_proka} {params} --out-dir {output_dir}/{wildcards.sample}/Flye -t {resources.cpus_per_task}
        """

rule medaka_mini_align:
    input:
        fastq_proka=fastq_proka,
        draft=f"{output_dir}/{{sample}}/Flye/assembly.fasta"
    output:
        bam=f"{output_dir}/{{sample}}/Medaka/calls_to_draft.bam",
        bai=f"{output_dir}/{{sample}}/Medaka/calls_to_draft.bam.bai"
    params:
        bam=temp(f"{output_dir}/{{sample}}/Medaka/calls_to_draft.tmp")
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
        bam=f"{output_dir}/{{sample}}/Medaka/calls_to_draft.bam",
        info=f"{output_dir}/{{sample}}/Flye/assembly_info.txt"
    output:
        f"{output_dir}/{{sample}}/Medaka/consensus_probs.hdf"
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
        probs=f"{output_dir}/{{sample}}/Medaka/consensus_probs.hdf",
        draft=f"{output_dir}/{{sample}}/Flye/assembly.fasta"
    output:
        f"{output_dir}/{{sample}}/Medaka/consensus.fasta"
    resources:
        cpus_per_task=1,
        mem_mb=24000
    shell:
        """
        export PATH={minimap2}:{samtools}:{htslib}:{medaka}:$PATH
        source {medaka}/activate
        medaka sequence {input.probs} {input.draft} {output}
        """
