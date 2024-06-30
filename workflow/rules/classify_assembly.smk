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
        ml load kma/1.4.12a
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
        ml load samtools/1.17
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
        ml load kma/1.4.12a
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
        ml load samtools/1.17
        samtools view -@{resources.cpus_per_task} -H --no-PG {input.sam} | sort | uniq > header.txt
        grep -v '^@' {input.sam} > reads.txt
        cat header.txt reads.txt > temp.sam
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 temp.sam | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        rm reads.txt header.txt temp.sam
        samtools index -@{resources.cpus_per_task} {output.bam}
        """

rule Flye_kma_align_pointfinder:
    input:
        assembly="results/{sample}/Medaka/consensus.fasta"
    output:
        sam=temp("results/{sample}/FlyeClassify/{sample}_PointF.sam"),
        res="results/{sample}/FlyeClassify/{sample}_PointF"
    params:
        database="/db/resfinder4/pointfinder/escherichia_coli/escherichia_coli",
        map="-tmp $(pwd)/ -mem_mode -ID 0.0 -cge -ef -na -nc -nf -verbose 2"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/kma_{sample}.log"
    shell:
        """
        ml load kma/1.4.12a
        kma -i {input} -o {output.res} -t {resources.cpus_per_task} -t_db {params.database} {params.map} -sam > {output.sam} 2>{log}
        touch {output.res}
        """

rule Flye_kma_repair_header_pointf:
    input:
        sam="results/{sample}/FlyeClassify/{sample}_PointF.sam"
    output:
        bam="results/{sample}/FlyeClassify/{sample}_PointF.bam"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=20000
    log: "logs/repair_PointF_{sample}.log"
    shell:
        """
        ml load samtools/1.17
        samtools view -@{resources.cpus_per_task} -H --no-PG {input.sam} | sort | uniq > header.txt
        grep -v '^@' {input.sam} > reads.txt
        cat header.txt reads.txt > temp_pointf.sam
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 temp_pointf.sam | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        rm reads.txt header.txt temp_pointf.sam
        samtools index -@{resources.cpus_per_task} {output.bam}
        """
