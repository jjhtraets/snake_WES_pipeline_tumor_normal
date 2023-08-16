rule manta_workflow:
    input:
        normal_bam = config["output_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",
        tumor_bam = config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam"
    params:
        ref = config["params"]["manta"]["ref"],
        run_dir = config["output_folder"]+"/Manta/{tumor}-vs-{normal}/"
    output:
        config["output_folder"]+"/Manta/{tumor}-vs-{normal}/runWorkflow.py"
    conda:
        "../envs/manta.yaml"
    threads:
        5
    shell:
        "configManta.py --exome --normalBam {input.normal_bam} --tumorBam {input.tumor_bam} --referenceFasta {params.ref} --runDir {params.run_dir}" 

rule run_manta_workflow:
    input:
        normal_bam = config["output_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",
        tumor_bam = config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam",
        workflow = config["output_folder"]+"/Manta/{tumor}-vs-{normal}/runWorkflow.py"
    output:
        diploid = config["output_folder"]+"/Manta/{tumor}-vs-{normal}/results/variants/diploidSV.vcf.gz",
        somatic = config["output_folder"]+"/Manta/{tumor}-vs-{normal}/results/variants/somaticSV.vcf.gz",
        candidateSV = config["output_folder"]+"/Manta/{tumor}-vs-{normal}/results/variants/candidateSV.vcf.gz",
        candidateIndel = config["output_folder"]+"/Manta/{tumor}-vs-{normal}/results/variants/candidateSmallIndels.vcf.gz",
    conda:
       "../envs/manta.yaml"
    threads: 4
    shell:
        "python2 {input.workflow}"
