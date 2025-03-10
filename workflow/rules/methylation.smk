output_dir=config.get("output_dir", "results")

rule modkit_pileup:
  input:
    bam=f"{output_dir}/{{sample}}/MapToAssemb/{{sample}}_assembly.bam"
  output:
    bed=f"{output_dir}/{{sample}}/NanoMotif/{{sample}}_assembly_modpileup.bed"
  log:
    f"{output_dir}/{{sample}}/logs/modkit.log"
  conda:
    path.join(base_dir, "envs/modkit.yaml")
  resources:
    cpus_per_task=config['threads'],
    mem_mb=config['memory']
  params:
    params="--write_bins --write_unbinned --score_threshold=0",
    outdir=f"{output_dir}/{{sample}}/DAS_Tool/DASTool"
  shell:
    """
    modkit pileup --filter-threshold 0.7 {input.bam} {output.bed} -t {resources.cpus_per_task}
    """