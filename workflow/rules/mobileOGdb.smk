mobileOGdb=config["tools"]["MobileOGdb"]["db"]
metadata=config["tools"]["MobileOGdb"]["metadata"]
output_dir=config.get("output_dir", "results")

rule mobileOGdb:
  input:
    f"{output_dir}/{{sample}}/DAS_Tool/DASTool_DASTool_bins/"
  output:
    f"{output_dir}/{{sample}}/mobileOGdb/DASTool_proteins.mobileOG.Alignment.Out.csv"
  log:
    f"{output_dir}/{{sample}}/logs/mobileOGdb.log"
  conda:
    "../envs/mobileOGdb.yaml"
  resources:
    cpus_per_task=config['threads'],
    mem_mb=config['memory']
  params:
    params="-k 15 -e 1e-20 -p 90 -q 90",
    db=mobileOGdb,
    metadata=metadata,
    script=path.join(base_dir, "scripts/mobileOGs-pl-kyanite.sh")
  shell:
    """
    {params.script} -i {output_dir}/{wildcards.sample}/DAS_Tool/DASTool_proteins.faa -d {params.db} -m {params.metadata} {params.params}
    """