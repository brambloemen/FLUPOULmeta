rule modkit_pileup:
    input:
        bam="results/{sample}/MapToAssemb/{sample}_assembly.bam"
    output:
        bed="results/{sample}/NanoMotif/{sample}_assembly_modpileup.bed"
    log:
      "results/{sample}/logs/modkit.log"
    conda:
      "../envs/modkit.yaml"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        params="--write_bins --write_unbinned --score_threshold=0",
        outdir="results/{sample}/DAS_Tool/DASTool"
    shell:
      """
       modkit pileup --filter-threshold 0.7 {input.bam} {output.bed} -t {resources.cpus_per_task}
      """