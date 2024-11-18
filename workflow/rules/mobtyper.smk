rule mobtyper:
    input:
        assembly="results/{sample}/Medaka/consensus.fasta"
    output:
        outdir=directory("results/{sample}/MobTyper"),
        mobres="results/{sample}/MobTyper/mobtyper_results.txt",
        mgeres="results/{sample}/MobTyper/mobtyper_mge_results.txt"
    log: "logs/mobtyper/{sample}_mobtyper.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    conda:
        "../envs/mob_suite.yaml"
    shell:
        """
        mob_typer --multi -i {input.assembly} -o {output.mobres} -g {output.mgeres} -n {resources.cpus_per_task}
        """