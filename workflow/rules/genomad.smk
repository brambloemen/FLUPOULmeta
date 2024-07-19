rule genomad:
    input:
        assembly="results/{sample}/Medaka/consensus.fasta"
    output:
        "results/{sample}/genomad_output/genomad_done"
    log: "logs/genomad/{sample}_scapp.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    conda:
        "../envs/genomad.yaml"
    shell:
        """
        genomad end-to-end -t {resources.cpus_per_task} {input.assembly} results/{wildcards.sample}/genomad_output /data/brbloemen/db/genomad_db
        touch {output}
        """