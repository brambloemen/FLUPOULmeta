rule metabat:
    input:
        assembly="results/{sample}/Medaka/consensus.fasta",
        bam="results/{sample}/MapToAssemb/{sample}_assembly.bam"
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

# rule graphmb:
#     input:
#         txt="results/{sample}/MetaBAT/depths.txt",
#         fasta="results/{sample}/Flye/assembly.fasta",
#         graph="results/{sample}/Flye/assembly.gfa",
#         bam="results/{sample}/Flye/{sample}_assembly.bam"
#     output:
#         binning_done="results/{sample}/GraphMB/GraphMB_done"
#     log:
#         "results/{sample}/logs/GraphMB_done.log"
#     conda:
#         "../envs/graphmb.yaml"
#     resources:
#         cpus_per_task=config['threads'],
#         mem_mb=config['memory']
#     params:
#         assembdir="results/{sample}/Flye/",
#         outputdir="results/{sample}/GraphMB/"
#     shell:
#         """
#         graphmb --assembly {params.assembdir} --outdir {params.outputdir} --assembly_name {input.fasta} --graph_file {input.graph} --depth {input.txt} --seed 37 --numcores {resources.cpus_per_task}
#         touch {output.binning_done}
#         """

rule medaka_contigbintsvs:
    input:
        metabat="results/{sample}/MetaBAT/metabat_done"
    output:
        metabatbins="results/{sample}/MetaBAT/MetaBATBins.tsv"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        inputdir="results/{sample}/MetaBAT/BIN"
    shell:
        """
        python scripts/ContigBinsTSV.py {params.inputdir} > {output.metabatbins}
        """

rule nanomotif_discovery:
    input:
        fasta="results/{sample}/Medaka/consensus.fasta",
        bed="results/{sample}/NanoMotif/{sample}_assembly_modpileup.bed",
        metabatbins="results/{sample}/MetaBAT/MetaBATBins.tsv"
    output:
        binmotifs="results/{sample}/NanoMotif/bin-motifs.tsv",
        motifsscored="results/{sample}/NanoMotif/motifs-scored.tsv",
        motifs="results/{sample}/NanoMotif/motifs.tsv"
    log:
      "results/{sample}/logs/nanomotif_disc.log"
    conda:
      "../envs/modkit.yaml"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        outdir="results/{sample}/NanoMotif/"
    shell:
      """
      nanomotif motif_discovery {input} --out {params.outdir} -t {resources.cpus_per_task}
      """

rule nanomotif_include:
    input:
        binmotifs="results/{sample}/NanoMotif/bin-motifs.tsv",
        motifsscored="results/{sample}/NanoMotif/motifs-scored.tsv",
        metabatbins="results/{sample}/MetaBAT/MetaBATBins.tsv"
    output:
        newbins="results/{sample}/NanoMotif/bin/new_contig_bin.tsv"
    log:
      "results/{sample}/logs/nanomotif_include.log"
    conda:
      "../envs/modkit.yaml"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        outdir="results/{sample}/NanoMotif/bin"
    shell:
      """
      nanomotif include_contigs --motifs_scored {input.motifsscored} \
      --bin_motifs {input.binmotifs} --contig_bins {input.metabatbins} \
      --out {params.outdir} -t {resources.cpus_per_task} --run_detect_contamination
      """

rule dastool:
    input:
        fasta="results/{sample}/Medaka/consensus.fasta",
        metabatbins="results/{sample}/MetaBAT/MetaBATBins.tsv",
        nanomotifbins="results/{sample}/NanoMotif/bin/new_contig_bin.tsv"
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
      DAS_Tool -i {input.metabatbins},{input.nanomotifbins} -c {input.fasta} \
      -o {params.outdir} --threads {resources.cpus_per_task} {params.params} 2>{log} && touch {output.DASdone}
      """