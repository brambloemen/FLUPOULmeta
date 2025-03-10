seqkit=config["tools"]["Seqkit"]["path"]
output_dir=config.get("output_dir", "results")

rule seqkit:
    input:
        lambda wildcards: config.get('fastq', '')
    output:
        fastq=f"{output_dir}/{{sample}}/fastq_hq/{{sample}}.HQ.fastq.gz"
    params:
        db="/db/kraken2_full/20240113/"
    resources:
        cpus_per_task=20,
        mem_mb=100000
    log:
        "logs/seqkit/{sample}.log"
    shell:
        """
        export PATH={seqkit}:$PATH
        seqkit seq -m 500 -Q 10 -j 20 {input} -o {output.fastq}
        """