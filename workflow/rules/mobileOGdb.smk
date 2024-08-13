rule mobileOGdb:
    input:
        "results/{sample}/DAS_Tool/DASTool_DASTool_bins/"
    output:
        "results/{sample}/mobileOGdb/DASTool_proteins.mobileOG.Alignment.Out.csv"
    log:
      "results/{sample}/logs/mobileOGdb.log"
    conda:
      "../envs/mobileOGdb.yaml"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        params="-k 15 -e 1e-20 -p 90 -q 90",
        db="/data/brbloemen/db/mobileOG-db/mobileOG-db_beatrix-1.6.dmnd",
        metadata="/data/brbloemen/db/mobileOG-db/mobileOG-db-beatrix-1.6-All.csv"
    shell:
      """
      scripts/mobileOGs-pl-kyanite.sh -i results/{wildcards.sample}/DAS_Tool/DASTool_proteins.faa -d {params.db} -m {params.metadata} {params.params}
      """