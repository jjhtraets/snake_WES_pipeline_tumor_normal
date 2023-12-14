rule index_strelka:
    input:
        bam = config["output_folder"]+ "/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam",
    output:
        config["output_folder"]+ "/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam.bai",
    singularity:
        config["SIF"]["samtools"]
    threads:
        1
    shell:
        """
        samtools index {input.bam}
        """

rule config_strelka:
    input:
        normal_bam = config["output_folder"]+ "/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",
        normal_bai = config["output_folder"]+ "/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam.bai",
        tumor_bam = config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam",
        tumor_bai = config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam.bai",
        ref = config["params"]["strelka"]["ref"]
    output:
        config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/runWorkflow.py"
    params:
        run_dir = config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/"
    threads:
        config["params"]["strelka"]["threads"]
    singularity:
        "docker://biocontainers/strelka:v2.9.7_cv2"
    log:
        config["output_folder"]+"/logs/strelka/{tumor}-vs-{normal}_strelka_conf.log"
    shell:
        """
        configureStrelkaSomaticWorkflow.py --normalBam {input.normal_bam} --tumorBam {input.tumor_bam} --referenceFasta {input.ref} --runDir {params.run_dir} --exome 2> {log}
        """

rule run_strelka:
    input:
        normal_bam = config["output_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",
        tumor_bam = config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam",
        ref = config["params"]["strelka"]["ref"],
        run = config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/runWorkflow.py"
    output:
        snvs = config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/results/variants/somatic.snvs.vcf.gz",
        indels = config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/results/variants/somatic.indels.vcf.gz"
    singularity:
        "docker://biocontainers/strelka:v2.9.7_cv2"
    threads:
        config["params"]["strelka"]["threads"]
    log:
        config["output_folder"]+"/logs/strelka/{tumor}-vs-{normal}_strelka_run.log"
    shell:
        """
        {input.run} -m local -j {threads} 2> {log}
        """
