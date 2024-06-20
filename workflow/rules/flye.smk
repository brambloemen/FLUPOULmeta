rule flye_assemble:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz"
    output:
        fasta="results/{sample}/Flye/assembly.fasta",
        graph="results/{sample}/Flye/assembly_graph.gfa"
    params:
        "--meta -i 2 --read-error 0.03"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/flye/{sample}.log"
    shell:
        """
        source /data/brbloemen/mambaforge/etc/profile.d/conda.sh
        conda deactivate
        ml load flye/2.9.3
        flye --nano-hq {input.fastq_proka} {params} --out-dir results/{wildcards.sample}/Flye -t {resources.cpus_per_task}
        """