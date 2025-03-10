output_dir = config.get("output_dir", "results")

rule mobtyper:
    input:
        assembly=f"{output_dir}/{{sample}}/Medaka/consensus.fasta"
    output:
        outdir=directory(f"{output_dir}/{{sample}}/MobTyper"),
        mobres=f"{output_dir}/{{sample}}/MobTyper/mobtyper_results.txt",
        mgeres=f"{output_dir}/{{sample}}/MobTyper/mobtyper_mge_results.txt"
    log: f"logs/mobtyper/{{sample}}_mobtyper.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    conda:
        path.join(base_dir, "envs/mob_suite.yaml")
    shell:
        """
        mob_typer --multi -i {input.assembly} -o {output.mobres} -g {output.mgeres} -n {resources.cpus_per_task}
        """