output_dir=config.get("output_dir", "results")

metabat=TOOLS["Metabat"]["path"]
metabat_bin=TOOLS["Metabat"]["bin_path"]
checkm=TOOLS["CheckM"]["path"]
checkm_dependencies=TOOLS["CheckM"]["dependencies"]
gtdbtk=TOOLS["GTDBtk"]["path"]
gtdbtk_dependencies=TOOLS["GTDBtk"]["dependencies"]

rule metabat:
    input:
        assembly=f"{output_dir}/{{sample}}/Medaka/consensus.fasta",
        bam=f"{output_dir}/{{sample}}/MapToAssemb/{{sample}}_assembly.bam"
    output:
        txt=f"{output_dir}/{{sample}}/MetaBAT/depths.txt",
        outputdir=directory(f"{output_dir}/{{sample}}/MetaBAT/BIN/")
    log:
        f"{output_dir}/{{sample}}/logs/metabat_done.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
        """
        export PATH={metabat}:{metabat_bin}:$PATH
        jgi_summarize_bam_contig_depths --outputDepth {output.txt} {input.bam} 2>{log}
        metabat2 -v -t {resources.cpus_per_task} -o {output.outputdir}/bin -i {input.assembly} -a {output.txt} 2>>{log}
        """

rule graphmb:
    input:
        txt=f"{output_dir}/{{sample}}/MetaBAT/depths.txt",
        fasta=f"{output_dir}/{{sample}}/Flye/assembly.fasta",
        bam=f"{output_dir}/{{sample}}/MapToAssemb/{{sample}}_assembly.bam"
    output:
        outputdir=directory(f"{output_dir}/{{sample}}/GraphMB/"),
        bins_dir=directory(f"{output_dir}/{{sample}}/GraphMB/_bins/"),
        contigbins=f"{output_dir}/{{sample}}/GraphMB/grapmb_contig2bin.tsv"
    log:
        f"{output_dir}/{{sample}}/logs/GraphMB.log"
    conda:
        path.join(base_dir, "envs/graphmb.yaml")
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        assembdir=f"{output_dir}/{{sample}}/Flye/",
        params="--assembly_type flye --vamb --minbin 250000 --mincontig 3000 --seed 37"
    shell:
        """
        cp {input.txt} {params.assembdir}/depths.txt
        # to update:
        graphmb --assembly {params.assembdir} --outdir {output.outputdir} --assembly_name assembly.fasta --depth depths.txt \
        --contignodes --numcores {resources.cpus_per_task} {params.params}
        grep -v '@' {output.outputdir}/_best_contig2bin.tsv > {output.contigbins}
        """

rule metabat_contigbintsvs:
    # TODO: process bins from other binning tools with the same script to generate contig_bins.tsv
    input:
        metabat=f"{output_dir}/{{sample}}/MetaBAT/BIN/"
    output:
        metabatbins=f"{output_dir}/{{sample}}/MetaBAT/MetaBATBins.tsv"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        inputdir=f"{output_dir}/{{sample}}/MetaBAT/BIN",
        script=path.join(base_dir, "scripts/ContigBinsTSV.py")
    conda:
        path.join(base_dir, "envs/FLUPOUL.yaml")
    shell:
        """
        python {params.script} {params.inputdir} > {output.metabatbins}
        """

rule semibin:
    input:
        assembly=f"{output_dir}/{{sample}}/Medaka/consensus.fasta",
        bam=f"{output_dir}/{{sample}}/MapToAssemb/{{sample}}_assembly.bam"
    output:
        contig_bins=f"{output_dir}/{{sample}}/semibin/contig_bins.tsv",
        contig_bins_filt=f"{output_dir}/{{sample}}/semibin/contig_bins_binned.tsv",
        outdir=directory(f"{output_dir}/{{sample}}/semibin")
    log:
      f"{output_dir}/{{sample}}/logs/semibin.log"
    conda:
        path.join(base_dir, "envs/SemiBin.yaml")
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
      """
      SemiBin2 single_easy_bin -i {input.assembly} -b {input.bam} -o {output.outdir} \
      --compression=none --sequencing-type long_read --environment global -t {resources.cpus_per_task}
      grep "\\-1" {output.contig_bins} -v | tail -n +2 > {output.contig_bins_filt}
      """

rule dastool:
    input:
        fasta=f"{output_dir}/{{sample}}/Medaka/consensus.fasta",
        metabatbins=f"{output_dir}/{{sample}}/MetaBAT/MetaBATBins.tsv",
        semibin_bins=f"{output_dir}/{{sample}}/semibin/contig_bins_binned.tsv",
        graphmb_bins=f"{output_dir}/{{sample}}/GraphMB/grapmb_contig2bin.tsv"
    output:
        DASdone=f"{output_dir}/{{sample}}/DAS_Tool/DASdone",
        bins=directory(f"{output_dir}/{{sample}}/DAS_Tool/DASTool_DASTool_bins/"),
        ctg2bin=f"{output_dir}/{{sample}}/DAS_Tool/DASTool_DASTool_contig2bin.tsv"
    log:
      f"{output_dir}/{{sample}}/logs/DAStool.log"
    conda:
       path.join(base_dir, "envs/DAS_Tool.yaml")
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        params="--write_bins --write_unbinned --score_threshold=0",
        outdir=f"{output_dir}/{{sample}}/DAS_Tool/DASTool"
    shell:
      """
      DAS_Tool -i {input.metabatbins},{input.semibin_bins},{input.graphmb_bins} -c {input.fasta} \
      -o {params.outdir} --threads {resources.cpus_per_task} {params.params} 2>{log} && touch {output.DASdone}
      cd {output.bins}
      if test -f "unbinned.fa"; then
        awk '/^>/ {{if (seq) close(out); out=substr($0, 2) ".fa"; seq=1}} {{print > out}}' unbinned.fa
        mv unbinned.fa unbinned.fasta
      fi
      """

rule nanomotif_discovery:
    input:
        fasta=f"{output_dir}/{{sample}}/Medaka/consensus.fasta",
        bed=f"{output_dir}/{{sample}}/NanoMotif/{{sample}}_assembly_modpileup.bed",
        bins=f"{output_dir}/{{sample}}/DAS_Tool/DASTool_DASTool_contig2bin.tsv"
    output:
        binmotifs=f"{output_dir}/{{sample}}/NanoMotif/bin-motifs.tsv",
        motifsscored=f"{output_dir}/{{sample}}/NanoMotif/motifs-scored.tsv",
        motifs=f"{output_dir}/{{sample}}/NanoMotif/motifs.tsv"
    log:
      f"{output_dir}/{{sample}}/logs/nanomotif_disc.log"
    conda:
        path.join(base_dir, "envs/nanomotif.yaml")
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        outdir=f"{output_dir}/{{sample}}/NanoMotif/"
    shell:
      """
      nanomotif motif_discovery {input} --out {params.outdir} -t {resources.cpus_per_task} --min_motif_score 0.2 --read_level_methylation --threshold_valid_coverage 1
      """
    #  defaults thresholds: general 0.7, confident 0.8
    # --threshold_methylation_general 0.4
    # --threshold_methylation_confident 0.5

rule nanomotif_include:
    input:
        assembly=f"{output_dir}/{{sample}}/Medaka/consensus.fasta",
        binmotifs=f"{output_dir}/{{sample}}/NanoMotif/bin-motifs.tsv",
        pileup=f"{output_dir}/{{sample}}/NanoMotif/{{sample}}_assembly_modpileup.bed",
        bins=f"{output_dir}/{{sample}}/DAS_Tool/DASTool_DASTool_contig2bin.tsv"
    output:
        newbins=f"{output_dir}/{{sample}}/NanoMotif/bin/new_contig_bin.tsv"
    log:
      f"{output_dir}/{{sample}}/logs/nanomotif_include.log"
    conda:
        path.join(base_dir, "envs/nanomotif.yaml")
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    params:
        outdir=f"{output_dir}/{{sample}}/NanoMotif/bin"
    shell:
      """
      nanomotif include_contigs --pileup {input.pileup} --assembly {input.assembly} \
      --bin_motifs {input.binmotifs} --contig_bins {input.bins} \
      --out {params.outdir} -t {resources.cpus_per_task} --run_detect_contamination
      mv {params.outdir}/new_contig_bin.tsv {params.outdir}/new_contig_bin.tmp.tsv
      tail -n +2 {params.outdir}/new_contig_bin.tmp.tsv > {output.newbins}
      """
#  --min_motif_comparisons 5 --mean_methylation_cutoff 0.1 --n_motif_contig_cutoff 10 --n_motif_bin_cutoff 500
# --ambiguous_motif_percentage_cutoff 0.4
# 
# TODO: update bin fastas with nanomotif binning; perform checkM and gtdb on these updated bins
rule checkm:
    input:
        f"{output_dir}/{{sample}}/DAS_Tool/DASTool_DASTool_bins/"
    output:
        outdir=directory(f"{output_dir}/{{sample}}/CheckM"),
        checkm=f"{output_dir}/{{sample}}/CheckM/CheckM_summary.txt"
    log:
      f"{output_dir}/{{sample}}/logs/CheckM.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
      """
      export PATH={checkm}:{checkm_dependencies}:$PATH
      source {checkm}/activate
      checkm lineage_wf -t{resources.cpus_per_task} -x fa {input} {output.outdir} --pplacer_threads {resources.cpus_per_task} > {output.checkm}
      """
    
rule gtdbtk:
    input:
        f"{output_dir}/{{sample}}/DAS_Tool/DASTool_DASTool_bins/"
    output:
        directory(f"{output_dir}/{{sample}}/GTDBtk")
    log:
      f"{output_dir}/{{sample}}/logs/gtdbtk.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
      """
      export PATH={gtdbtk}:{gtdbtk_dependencies}:$PATH
      export GTDBTK_DATA_PATH=/db/gtdb/release214
      source {gtdbtk}/activate
      gtdbtk classify_wf --genome_dir {input} --out_dir {output} --cpus {resources.cpus_per_task} --pplacer_cpus 1 \
      --extension fa --skip_ani_screen
      """

rule nanomotif_bins:
    input:
        fasta_file = f"{output_dir}/{{sample}}/Medaka/consensus.fasta",
        contig_bin = f"{output_dir}/{{sample}}/NanoMotif/bin/new_contig_bin.tsv"
    output:
        outdir = directory(f"{output_dir}/{{sample}}/NanoMotif/bin_fastas/")
    log:
        f"{output_dir}/{{sample}}/logs/nanomotif_bins.log"
    resources:
        cpus_per_task=1,
        mem_mb=10000
    script:
        path.join(base_dir, "scripts/extract_bins_fastas.py")

rule checkm_NM:
    input:
        f"{output_dir}/{{sample}}/NanoMotif/bin_fastas/"
    output:
        outdir=directory(f"{output_dir}/{{sample}}/CheckM_NM"),
        checkm=f"{output_dir}/{{sample}}/CheckM_NM/CheckM_summary.txt"
    log:
      f"{output_dir}/{{sample}}/logs/CheckM_NM.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
      """
      export PATH={checkm}:{checkm_dependencies}:$PATH
      source {checkm}/activate
      checkm lineage_wf -t{resources.cpus_per_task} -x fa {input} {output.outdir} --pplacer_threads {resources.cpus_per_task} > {output.checkm}
      """
    
rule gtdbtk_NM:
    input:
        f"{output_dir}/{{sample}}/NanoMotif/bin_fastas/"
    output:
        directory(f"{output_dir}/{{sample}}/GTDBtk_NM")
    log:
      f"{output_dir}/{{sample}}/logs/gtdbtk_NM.log"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    shell:
      """
      export PATH={gtdbtk}:{gtdbtk_dependencies}:$PATH
      export GTDBTK_DATA_PATH=/db/gtdb/release214
      source {gtdbtk}/activate
      gtdbtk classify_wf --genome_dir {input} --out_dir {output} --cpus {resources.cpus_per_task} --pplacer_cpus 1 \
      --extension fa --skip_ani_screen
      """