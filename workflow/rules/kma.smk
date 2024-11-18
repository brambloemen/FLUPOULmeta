KMA=TOOLS["KMA"]["path"]
samtools=TOOLS["Samtools"]["path"]
taxa_db=TOOLS["KMA"]["resfinder_db"]
resf_db=TOOLS["KMA"]["taxonomic_db"]

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
        export PATH={KMA}:$PATH
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
        export PATH={samtools}:$PATH
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
        export PATH={KMA}:$PATH
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
        export PATH={samtools}:$PATH
        samtools view -@{resources.cpus_per_task} -H --no-PG {input.sam} | sort | uniq > header.txt
        grep -v '^@' {input.sam} > reads.txt
        cat header.txt reads.txt > temp.sam
        samtools view -hb -@{resources.cpus_per_task} -F 0x904 temp.sam | \
        samtools sort -@{resources.cpus_per_task} > {output.bam} 2>{log}
        rm reads.txt header.txt temp.sam
        """