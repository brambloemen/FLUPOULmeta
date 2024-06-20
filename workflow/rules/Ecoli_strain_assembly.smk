rule extract_ecoli_reads:
    input:
        # kraken="results/{sample}/{sample}.kraken",
        # report="results/{sample}/{sample}.kreport",
        mapped_assembly="results/{sample}/FlyeClassify/{sample}_assembly.bam",
        readmapping_assemb="results/{sample}/Flye/{sample}_assembly.bam"
    output:
        ecolcontiglist="results/{sample}/FlyeClassify/{sample}_EcolContigs.ls",
        filtered_reads="results/{sample}/Ecoli/{sample}_ecolreads.fastq.gz"
    params:
        ""
    resources:
        cpus_per_task=8,
        mem_mb=50000
    log:
        "logs/extract_ecoli_reads/{sample}.log"
    shell:
        """
        ml load samtools/1.17
        samtools view {input.mapped_assembly} | cut -f1-5 | grep "Escherichia coli" | cut -f1 > {output.ecolcontiglist}
        ecolcontigs=$(awk '{{$1=$1}}1' OFS=' ' {output.ecolcontiglist} | tr '\n' ' ')
        samtools view -h {input.readmapping_assemb} $ecolcontigs | samtools fastq | pigz -p{resources.cpus_per_task} -c > {output.filtered_reads}
        """

rule assemb_ecol_reads:
    input:
        reads="results/{sample}/Ecoli/{sample}_ecolreads.fastq.gz"
    output:
        fasta="results/{sample}/Ecoli/Flye_{sample}/assembly.fasta",
        graph="results/{sample}/Ecoli/Flye_{sample}/assembly_graph.gfa"
    params:
        "--meta --no-alt-contigs --keep-haplotypes -i 0 --read-error 0.03"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/Ecoli_flye/{sample}.log"
    shell:
        """
        source /data/brbloemen/mambaforge/etc/profile.d/conda.sh
        conda deactivate
        ml load flye/2.9.3
        flye --nano-hq {input.reads} {params} --out-dir results/{wildcards.sample}/Ecoli/Flye_{wildcards.sample} -t {resources.cpus_per_task}
        """