rule kraken2:
    input:
        lambda wildcards: config.get('fastq', '')
    output:
        report="results/{sample}/{sample}.kreport",
        kraken="results/{sample}/{sample}.kraken",
        krona="results/{sample}/{sample}.krona.html"
    params:
        db="/db/kraken2_full/20240113/"
    resources:
        cpus_per_task=1,
        mem_mb=50000
    log:
        "logs/kraken2/{sample}.log"
    shell:
        """
        ml load kraken2_remote krona/2.7
        mkdir -p results/{wildcards.sample}
        kraken2 --db {params.db} --threads 20 --report {output.report} --output {output.kraken} {input}
        ktImportTaxonomy {output.report} -t 5 -m 3 -o {output.krona}
        ml -kraken2_remote 
        """

rule extract_kraken_reads:
    input:
        kraken="results/{sample}/{sample}.kraken",
        report="results/{sample}/{sample}.kreport",
        reads=lambda wildcards: config.get('fastq', '')
    output:
        filtered_reads="results/{sample}/{sample}_filtered.fastq.gz"
    params:
        "-t 2759 --exclude --include-children --fastq-output"
    resources:
        cpus_per_task=9,
        mem_mb=50000
    log:
        "logs/extract_kraken_reads/{sample}.log"
    shell:
        """
        python scripts/extract_kraken_reads.py -k {input.kraken} -s {input.reads} -o tmp.fastq {params} -r {input.report}
        pigz -p{resources.cpus_per_task} tmp.fastq 
        mv tmp.fastq.gz {output.filtered_reads}
        """