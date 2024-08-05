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

rule metabat_contigbintsvs:
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

rule semibin:
    input:
        assembly="results/{sample}/Medaka/consensus.fasta",
        bam="results/{sample}/MapToAssemb/{sample}_assembly.bam"
    output:
        contig_bins="results/{sample}/semibin/contig_bins.tsv",
        contig_bins_filt="results/{sample}/semibin/contig_bins_binned.tsv"
    log:
      "results/{sample}/logs/semibin.log"
    conda:
      "../envs/SemiBin.yaml"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
      """
      SemiBin2 single_easy_bin -i {input.assembly} -b {input.bam} -o results/{wildcards.sample}/semibin \
      --compression=none --sequencing-type long_read --environment global -t {resources.cpus_per_task}
      grep "\\-1" {output.contig_bins} -v | tail -n +2 > {output.contig_bins_filt}
      """

rule dastool:
    input:
        fasta="results/{sample}/Medaka/consensus.fasta",
        metabatbins="results/{sample}/MetaBAT/MetaBATBins.tsv",
        semibin_bins="results/{sample}/semibin/contig_bins_binned.tsv"
    output:
        DASdone="results/{sample}/DAS_Tool/DASdone",
        bins=directory("results/{sample}/DAS_Tool/DASTool_DASTool_bins/"),
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
      DAS_Tool -i {input.metabatbins},{input.semibin_bins} -c {input.fasta} \
      -o {params.outdir} --threads {resources.cpus_per_task} {params.params} 2>{log} && touch {output.DASdone}
      cd {output.bins}
      if test -f "unbinned.fa"; then
        awk '/^>/ {{if (seq) close(out); out=substr($0, 2) ".fa"; seq=1}} {{print > out}}' unbinned.fa
        mv unbinned.fa unbinned.fasta
      fi
      """

rule nanomotif_discovery:
    input:
        fasta="results/{sample}/Medaka/consensus.fasta",
        bed="results/{sample}/NanoMotif/{sample}_assembly_modpileup.bed",
        bins="results/{sample}/DAS_Tool/DASTool_DASTool_contig2bin.tsv"
    output:
        binmotifs="results/{sample}/NanoMotif/bin-motifs.tsv",
        motifsscored="results/{sample}/NanoMotif/motifs-scored.tsv",
        motifs="results/{sample}/NanoMotif/motifs.tsv"
    log:
      "results/{sample}/logs/nanomotif_disc.log"
    conda:
      "../envs/nanomotif.yaml"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        outdir="results/{sample}/NanoMotif/"
    shell:
      """
      nanomotif motif_discovery {input} --out {params.outdir} -t {resources.cpus_per_task} --threshold_methylation_general 0.4 --threshold_methylation_confident 0.5
      """
    #  defaults thresholds: general 0.7, confident 0.8
    # --threshold_methylation_general 0.4
    # --threshold_methylation_confident 0.5

rule nanomotif_include:
    input:
        binmotifs="results/{sample}/NanoMotif/bin-motifs.tsv",
        motifsscored="results/{sample}/NanoMotif/motifs-scored.tsv",
        bins="results/{sample}/DAS_Tool/DASTool_DASTool_contig2bin.tsv"
    output:
        newbins="results/{sample}/NanoMotif/bin/new_contig_bin.tsv"
    log:
      "results/{sample}/logs/nanomotif_include.log"
    conda:
      "../envs/nanomotif.yaml"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        outdir="results/{sample}/NanoMotif/bin"
    shell:
      """
      nanomotif include_contigs --motifs_scored {input.motifsscored} \
      --bin_motifs {input.binmotifs} --contig_bins {input.bins} \
      --out {params.outdir} -t {resources.cpus_per_task} --run_detect_contamination --save_scores
      mv {params.outdir}/new_contig_bin.tsv {params.outdir}/new_contig_bin.tmp.tsv
      tail -n +2 {params.outdir}/new_contig_bin.tmp.tsv > {output.newbins}
      """
#  --min_motif_comparisons 5 --mean_methylation_cutoff 0.1 --n_motif_contig_cutoff 10 --n_motif_bin_cutoff 500
# --ambiguous_motif_percentage_cutoff 0.4
# 
rule checkm:
    input:
        "results/{sample}/DAS_Tool/DASTool_DASTool_bins/"
    output:
        checkm="results/{sample}/CheckM/CheckM_summary.txt"
    log:
      "results/{sample}/logs/CheckM.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
      """
      ml load hmmer prodigal pplacer CheckM
      checkm lineage_wf -t{resources.cpus_per_task} -x fa {input} results/{wildcards.sample}/CheckM --pplacer_threads {resources.cpus_per_task} > {output}
      """