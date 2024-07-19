rule kma_align:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz"
    output:
        sam=temp("results/{sample}/{sample}_kma.sam"),
        res="results/{sample}/{sample}_kma"
    params:
        database="/scratch/alvanuffelen/kma_v2/kma_db",
        map="-tmp $(pwd)/ -mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nc -nf -1t1 -ca -verbose 2"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=config['memory']
    log: "logs/kma/{sample}.log"
    shell:
        """
        ml load kma/1.4.15 samtools/1.17
        kma -i {input} -o {output.res} -t {resources.cpus_per_task} -t_db {params.database} {params.map} -sam > {output.sam} 2>{log}
        touch {output.res}
        """

rule samtools_filter:
    input:
        sam="results/{sample}/{sample}_kma.sam"
    output:
        bam="results/{sample}/{sample}.kmadb.bam",
        bai="results/{sample}/{sample}.kmadb.bam.bai"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=50000
    log: "logs/kma_filter/{sample}.log"
    shell:
        """
        ml load kma/1.4.15 samtools/1.17
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 {input.sam} | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        samtools index -@{resources.cpus_per_task} {output.bam}
        """

rule kma_align_resfinder:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz"
    output:
        sam=temp("results/{sample}/{sample}_resf.sam"),
        res="results/{sample}/{sample}_resf"
    params:
        database="/db/resfinder4/latest/all",
        map="-tmp $(pwd)/ -mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nc -nf -ca -verbose 2"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=50000
    log: "logs/kma/{sample}.log"
    shell:
        """
        ml load kma/1.4.15 samtools/1.17
        kma -i {input} -o {output.res} -t {resources.cpus_per_task} -t_db {params.database} {params.map} -sam > {output.sam} 2>{log}
        touch {output.res}
        """

rule repair_header:
    input:
        sam="results/{sample}/{sample}_resf.sam"
    output:
        bam="results/{sample}/{sample}_rep.resf.bam"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=50000
    log: "logs/repair_ResF/{sample}.log"
    shell:
    # use grep -v below to work around samtools error caused by duplicate header entries (created by kma mapper)
        """
        ml load kma/1.4.15 samtools/1.17
        samtools view -@{resources.cpus_per_task} -H --no-PG {input.sam} | sort | uniq > header.txt
        grep -v '^@' {input.sam} > reads.txt
        cat header.txt reads.txt > temp.sam
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 temp.sam | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        rm reads.txt header.txt temp.sam
        """

rule kma_align_NDARO:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz"
    output:
        sam=temp("results/{sample}/NDARO/{sample}_NDARO.sam"),
        res="results/{sample}/NDARO/{sample}_NDARO"
    params:
        database="/db/gene_detection/NCBI_AMR/kma/ncbi_amr-clustered_80",
        map="-tmp $(pwd)/ -mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nc -nf -ca -verbose 2"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=50000
    log: "logs/kma/{sample}.log"
    shell:
        """
        ml load kma/1.4.15 samtools/1.17
        kma -i {input} -o {output.res} -t {resources.cpus_per_task} -t_db {params.database} {params.map} -sam > {output.sam} 2>{log}
        touch {output.res}
        # copy database metadata to see how recent the used db was
        cp /db/gene_detection/NCBI_AMR/db_metadata.txt results/{wildcards.sample}/NDARO/NDARO_db_metadata.txt
        """

rule samtools_NDARO:
    input:
        sam="results/{sample}/NDARO/{sample}_NDARO.sam"
    output:
        bam="results/{sample}/NDARO/{sample}_NDARO.bam"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=50000
    log: "logs/kma_filter/{sample}.log"
    shell:
        """
        ml load kma/1.4.15 samtools/1.17
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 {input.sam} | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        """

rule kma_align_pointfinder:
    input:
        fastq_proka="results/{sample}/{sample}_filtered.fastq.gz"
    output:
        sam=temp("results/{sample}/{sample}_pointf.sam"),
        res="results/{sample}/{sample}_pointf"
    params:
        database="/db/resfinder4/pointfinder/escherichia_coli/escherichia_coli",
        map="-tmp $(pwd)/ -mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nc -nf -ca -verbose 2"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=50000
    log: "logs/kma/{sample}.log"
    shell:
        """
        ml load kma/1.4.15 samtools/1.17
        kma -i {input} -o {output.res} -t {resources.cpus_per_task} -t_db {params.database} {params.map} -sam > {output.sam} 2>{log}
        touch {output.res}
        """

rule repair_header_pointf:
    input:
        sam="results/{sample}/{sample}_pointf.sam"
    output:
        bam="results/{sample}/{sample}_rep.pointf.bam"
    resources:
        cpus_per_task=config['threads'],
        mem_mb=50000
    log: "logs/repair_PointF/{sample}.log"
    shell:
        """
        ml load kma/1.4.15 samtools/1.17
        samtools view -@{resources.cpus_per_task} -H --no-PG {input.sam} | sort | uniq > header.txt
        grep -v '^@' {input.sam} > reads.txt
        cat header.txt reads.txt > temp_pointf.sam
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 temp_pointf.sam | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        rm reads.txt header.txt temp_pointf.sam
        """
