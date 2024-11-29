KMA=TOOLS["KMA"]["path"]
samtools=TOOLS["Samtools"]["path"]
taxa_db=TOOLS["KMA"]["resfinder_db"]
resf_db=TOOLS["KMA"]["taxonomic_db"]

rule Flye_kma_classify:
    input:
        assembly="results/{sample}/Medaka/consensus.fasta"
    output:
        sam=temp("results/{sample}/FlyeClassify/{sample}_assembly.sam"),
        res="results/{sample}/FlyeClassify/{sample}"
    params:
        database="/scratch/alvanuffelen/kma_v2/kma_db",
        map="-tmp $(pwd)/ -mem_mode -ID 0.0 -ef -proxi 0.9 -na -nc -nf -1t1 -verbose 2"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/Flye_kma_{sample}.log"
    shell:
        """
        export PATH={KMA}:$PATH
        kma -i {input} -o {output.res} -t {resources.cpus_per_task} -t_db {params.database} {params.map} -sam > {output.sam} 2>{log}
        touch {output.res}
        """

rule Flye_kma_classify_samtools_filter:
    input:
        sam="results/{sample}/FlyeClassify/{sample}_assembly.sam"
    output:
        bam="results/{sample}/FlyeClassify/{sample}_assembly.bam"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=20000
    log: "logs/samtools_{sample}.log"
    shell:
        """
        export PATH={samtools}:$PATH
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 {input.sam} | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        samtools index -@{resources.cpus_per_task} {output.bam}
        """

rule Flye_kma_align_resfinder:
    input:
        assembly="results/{sample}/Medaka/consensus.fasta"
    output:
        sam=temp("results/{sample}/FlyeClassify/{sample}_ResF.sam"),
        res="results/{sample}/FlyeClassify/{sample}_ResF"
    params:
        database="/db/resfinder4/latest/all",
        map="-tmp $(pwd)/ -mem_mode -ID 0.0 -cge -ef -na -nc -nf -verbose 2"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/kma_{sample}.log"
    shell:
        """
        export PATH={KMA}:$PATH
        kma -i {input} -o {output.res} -t {resources.cpus_per_task} -t_db {params.database} {params.map} -sam > {output.sam} 2>{log}
        touch {output.res}
        """

rule Flye_kma_ResF_repair_header:
    input:
        sam="results/{sample}/FlyeClassify/{sample}_ResF.sam"
    output:
        bam="results/{sample}/FlyeClassify/{sample}_ResF.bam"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=20000
    log: "logs/repair_ResF_{sample}.log"
    shell:
        """
        export PATH={samtools}:$PATH
        samtools view -@{resources.cpus_per_task} -H --no-PG {input.sam} | sort | uniq > header.txt
        grep -v '^@' {input.sam} > reads.txt
        cat header.txt reads.txt > temp.sam
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 temp.sam | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        rm reads.txt header.txt temp.sam
        samtools index -@{resources.cpus_per_task} {output.bam}
        """

rule AMRFinder:
    input:
        assembly="results/{sample}/Medaka/consensus.fasta"
    output:
        "results/{sample}/AMRFinder/{sample}_amrfinder.txt"
    resources:
        cpus_per_task=1,
        mem_mb=50000
    log: "logs/AMRFinder_{sample}.log"
    conda:
      "../envs/amrfinder.yaml"
    shell:
        """
        amrfinder -n {input.assembly} > {output}
        """