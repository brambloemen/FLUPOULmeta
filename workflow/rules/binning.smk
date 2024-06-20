rule metabat:
    input:
        assembly="results/{sample}/Flye/assembly.fasta",
        bam="results/{sample}/Flye/{sample}_assembly.bam"
    output:
        txt="results/{sample}/MetaBAT/depths.txt",
        binning_done="results/{sample}/MetaBAT/metabat_done"
    log:
        "results/{sample}/logs/metabat_done.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        outputdir="results/{sample}/MetaBAT/BIN/bin"
    shell:
        """
        ml load metabat2
        jgi_summarize_bam_contig_depths --outputDepth {output.txt} {input.bam} 2>{log}
        metabat2 -v -t {resources.cpus_per_task} -o {params.outputdir} -i {input.assembly} -a {output.txt} 2>>{log}
        touch {output.binning_done}
        """

rule contigbintsvs:
    input:
        metabat="results/{sample}/MetaBAT/metabat_done"
    output:
        metabatbins="results/{sample}/DAS_Tool/MetaBATBins.tsv"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        inputdir="results/{sample}/MetaBAT/BIN"
    shell:
        """
        python scripts/ContigBinsTSV.py {params.inputdir} > {output.metabatbins}
        """

rule dastool:
    input:
        fasta= "results/{sample}/Flye/assembly.fasta",
        metabatbins="results/{sample}/DAS_Tool/MetaBATBins.tsv"
    output:
        DASdone="results/{sample}/DAS_Tool/DASdone",
        ctg2bin="results/{sample}/DAS_Tool/DASTool_DASTool_contig2bin.tsv"
    log:
      "results/{sample}/logs/DAStool.log"
    conda:
      "../envs/DAS_Tool.yaml"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        params="--write_bins --write_unbinned --score_threshold=0",
        outdir="results/{sample}/DAS_Tool/DASTool"
    shell:
      """
      DAS_Tool -i {input.metabatbins} -c {input.fasta} \
      -o {params.outdir} --threads {resources.cpus_per_task} {params.params} 2>{log} && touch {output.DASdone}
      """