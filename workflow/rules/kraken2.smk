kraken2=config["tools"]["Kraken2"]["path"]
kraken2_env=config["tools"]["Kraken2"]["venv"]
krakendb=config["tools"]["Kraken2"]["db"]
krona=config["tools"]["Krona"]["path"]
output_dir=config.get("output_dir", "results")

rule kraken2:
    input:
        fastq=f"{output_dir}/{{sample}}/fastq_hq/{{sample}}.HQ.fastq.gz"
    output:
        report=f"{output_dir}/{{sample}}/{{sample}}.kreport",
        kraken=f"{output_dir}/{{sample}}/{{sample}}.kraken",
        krona=f"{output_dir}/{{sample}}/{{sample}}.krona.html"
    params:
        db=krakendb
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log:
        f"logs/kraken2/{{sample}}.log"
    shell:
        """
        export PATH={kraken2}:{kraken2_env}:{krona}:$PATH
        source {kraken2_env}/activate
        mkdir -p {output_dir}/{wildcards.sample}
        kraken2 --db {params.db} --threads 20 --report {output.report} --output {output.kraken} {input}
        ktImportTaxonomy {output.report} -t 5 -m 3 -o {output.krona}
        """

rule extract_kraken_reads:
    input:
        kraken=f"{output_dir}/{{sample}}/{{sample}}.kraken",
        report=f"{output_dir}/{{sample}}/{{sample}}.kreport",
        reads=lambda wildcards: config.get('fastq', '')
    output:
        filtered_reads=f"{output_dir}/{{sample}}/{{sample}}_filtered.fastq.gz"
    params:
        extract="-t 2759 --exclude --include-children --fastq-output",
        script=path.join(base_dir, "scripts/extract_kraken_reads.py")
    resources:
        cpus_per_task=9,
        mem_mb=50000
    log:
        f"logs/extract_kraken_reads/{{sample}}.log"
    conda:
        path.join(base_dir, "envs/FLUPOUL.yaml")
    shell:
        """
        python {params.script} -k {input.kraken} -s {input.reads} -o tmp.fastq {params.extract} -r {input.report}
        pigz -p{resources.cpus_per_task} tmp.fastq && mv tmp.fastq.gz {output.filtered_reads}
        """