KMA=TOOLS["KMA"]["path"]
samtools=TOOLS["Samtools"]["path"]
taxa_db=TOOLS["KMA"]["resfinder_db"]
resf_db=TOOLS["KMA"]["taxonomic_db"]
output_dir=config.get("output_dir", "results")

rule Flye_kma_classify:
    input:
        assembly=f"{output_dir}/{{sample}}/Medaka/consensus.fasta"
    output:
        sam=temp(f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_assembly.sam"),
        res=f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}"
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
        sam=f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_assembly.sam"
    output:
        bam=f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_assembly.bam"
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
        assembly=f"{output_dir}/{{sample}}/Medaka/consensus.fasta"
    output:
        sam=temp(f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_ResF.sam"),
        res=f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_ResF"
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
        sam=f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_ResF.sam"
    output:
        bam=f"{output_dir}/{{sample}}/FlyeClassify/{{sample}}_ResF.bam"
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
        assembly=f"{output_dir}/{{sample}}/Medaka/consensus.fasta"
    output:
        report=f"{output_dir}/{{sample}}/AMRFinder/{{sample}}_amrfinder.txt"
    resources:
        cpus_per_task=1,
        mem_mb=50000
    log: "logs/AMRFinder_{sample}.log"
    params:
        db_flag = path.join(base_dir, "resources/amrfinder_db/.setup_complete")
    conda:
        path.join(base_dir, "envs/amrfinder.yaml")
    shell:
        """
        # Setup DB if not already done
        if [ ! -f {params.db_flag} ]; then
            amrfinder -u
        fi
        amrfinder -n {input.assembly} > {output.report}
        """