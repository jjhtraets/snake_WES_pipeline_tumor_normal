rule bwa_mem:
    input:
        reads = get_fastq_files,
        ref = config["params"]["bwa"]["ref"]
    output:
        temp(config["output_folder"]+"/mapped/{sample}_hg38.sam")
    singularity:
        config["SIF"]["bwa"]
    threads:
        config["params"]["bwa"]["threads"]
    log:
        config["output_folder"]+"/logs/GATK/{sample}_bwa_mem.log"
    shell:
        """
        bwa mem -M -t {threads} {input.ref} {input.reads} > {output} 2> {log}
        """

rule samtools_sort: # coordinate
    input:
        config["output_folder"]+"/mapped/{sample}_hg38.sam"
    output:
        temp(config["output_folder"]+"/mapped/{sample}_sorted_hg38.bam")
    singularity:
        config["SIF"]["samtools"]
    threads:
        1
    shell:
        """
        samtools view -Sb {input} | samtools sort > {output}
        """

rule picard_groups:
    input:
        config["output_folder"]+"/mapped/{sample}_sorted_hg38.bam"
    output:
        temp(config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG.bam")
    threads:
        config["params"]["picard"]["threads"]
    params:
        groups = config["params"]["picard"]["groups"]
    singularity:
        config["SIF"]["picard"]
    log:
        config["output_folder"]+"/logs/GATK/{sample}_picard_ARRG.log"
    shell:
        """
        java -XX:ParallelGCThreads={threads} -jar /usr/picard/picard.jar AddOrReplaceReadGroups I= {input} O= {output} {params.groups} RGSM={wildcards.sample} VALIDATION_STRINGENCY=LENIENT 2> {log}
        """

rule picard_dedup:
    input:
        config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG.bam"
    output:
        dup = temp(config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup.bam"),
        met = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup.bam_metrix.txt"
    threads:
        config["params"]["picard"]["threads"]
    singularity:
        config["SIF"]["picard"]
    log:
        config["output_folder"]+"/logs/GATK/{sample}_picard_dedup.log"
    shell:
        """
        java -XX:ParallelGCThreads={threads} -jar /usr/picard/picard.jar MarkDuplicates I= {input} O= {output.dup} METRICS_FILE= {output.met} 2> {log}
        """

rule samtools_index:
    input:
        config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup.bam"
    output:
        temp(config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup.bam.bai")
    singularity:
        config["SIF"]["samtools"]
    threads:
        1
    shell:
        """
        samtools index {input}
        """
        
rule recal_base_qualities:
    input:
        bam = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup.bam",
        genome = config["params"]["gatk"]["ref"],
        known_sites = config["params"]["gatk"]["known_sites"],
        bai = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup.bam.bai"
    output:
        config["output_folder"]+"/mapped/{sample}_recal_data.table"
    singularity:
        config["SIF"]["gatk"]
    log:
        config["output_folder"]+"/logs/GATK/{sample}_gatk_bq.log"
    threads:
        config["params"]["gatk"]["threads"]
    shell:
        """
        gatk BaseRecalibrator -R {input.genome} -I {input.bam} -O {output} --known-sites {input.known_sites} 2> {log}
        """

rule apply_bqscore:
    input:
        bam = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup.bam",
        genome = config["params"]["gatk"]["ref"],
        table_recal = config["output_folder"]+"/mapped/{sample}_recal_data.table"
    output:
        config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam"
    singularity:
        config["SIF"]["gatk"]
    log:
        config["output_folder"]+"/logs/GATK/{sample}_gatk_recal.log"
    threads:
        config["params"]["gatk"]["threads"]
    shell:
        """
        gatk ApplyBQSR -R {input.genome} --bqsr-recal-file {input.table_recal} -I {input.bam} -O {output} && samtools index {output} 2> {log}
        """
