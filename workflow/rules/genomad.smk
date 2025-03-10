genomad_db=config["tools"]["Genomad"]["db"]
output_dir=config.get("output_dir", "results")

rule genomad:
    input:
        assembly=f"{output_dir}/{{sample}}/Medaka/consensus.fasta"
    output:
        genomaddir=directory(f"{output_dir}/{{sample}}/genomad_output")
    log: f"logs/genomad/{{sample}}_scapp.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        db=genomad_db
    conda:
        path.join(base_dir, "envs/genomad.yaml")
    shell:
        """
        genomad end-to-end -t {resources.cpus_per_task} {input.assembly} {output.genomaddir} {params.db}
        """