rule fastqc:
    input:
        reads = [config["input_folder"]+"/{sample}"+ config["read_name"]["R1"] +".fastq.gz",config["input_folder"]+"/{sample}"+ config["read_name"]["R2"] +".fastq.gz"]
    output:
        f1=config["output_folder"]+ "/fastqc/{sample}"+ config["read_name"]["R1"] +"_fastqc.html",
        f2=config["output_folder"]+ "/fastqc/{sample}"+ config["read_name"]["R2"] +"_fastqc.html"
    singularity:
        config["SIF"]["fastqc"]
    threads:
        1
    params:
        input_f = config["output_folder"]+"/fastqc/"
    shell:
        """
        fastqc -o {params.input_f} -t {threads} {input.reads}
        """

rule multiqc:
    input:
        f1=expand(config["output_folder"]+"/fastqc/{sample}"+ config["read_name"]["R2"] +"_fastqc.html",sample=samples["sample_ID"]),
        f2=expand(config["output_folder"]+"/fastqc/{sample}"+ config["read_name"]["R2"] +"_fastqc.html",sample=samples["sample_ID"])
    output:
        config["output_folder"]+"/multiqc/multiqc_report.html"
    singularity:
        config["SIF"]["multiqc"]
    params:
        input_f = config["output_folder"]+"/fastqc/",
        output_f = config["output_folder"]+"/multiqc/"
    shell:
        """
        multiqc {params.input_f} -o {params.output_f} -n multiqc_report
        """

rule samtools_stat:
    input:
        bam = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam",
        genome = config["params"]["gatk"]["ref"]
    output:
        config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.text"
    singularity:
        config["SIF"]["samtools"]
    shell:
        """
        samtools stats {input.bam} --ref-seq {input.genome} > {output}
        """

rule bcftools_stat:
    input:
        vcf = config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered-pass.vcf",
        genome = config["params"]["gatk"]["ref"]
    output:
        config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered-pass.text"
    singularity:
        config["SIF"]["bcftools"]
    shell:
        """
        bcftools stats {input.vcf} > {output}
        """
        
rule mosdepth:
    input:
        bam = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam",
        genome = config["params"]["gatk"]["ref"]
    output:
        config["output_folder"]+"/mapped/{sample}.mosdepth.global.dist.txt"
    params:
        bed = config["params"]["CNVkit"]["bait_bed"],
        output = config["output_folder"]+"/mapped/{sample}"
    singularity:
        config["SIF"]["mosdepth"]
    threads:
        config["params"]["gatk"]["threads"]
    shell:
        """
        mosdepth --threads {threads} --by {params.bed} {params.output} {input.bam}
        """

rule multiqc_gatk:
    input:
        f1=expand(config["output_folder"]+"/fastqc/{sample}"+ config["read_name"]["R2"] +"_fastqc.html",sample=samples["sample_ID"]),
        f2=expand(config["output_folder"]+"/fastqc/{sample}"+ config["read_name"]["R2"] +"_fastqc.html",sample=samples["sample_ID"]),
        bam = expand(config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup.bam_metrix.txt",sample=samples["sample_ID"]),
        bam_stat =  expand(config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.text",sample=samples["sample_ID"]),
        bam_depth = expand(config["output_folder"]+"/mapped/{sample}.mosdepth.global.dist.txt",sample=samples["sample_ID"]), 
        vcf_stat =  expand(config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered-pass.text",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"]),
        gatk = expand(config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered-pass.vcf",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])
    output:
        config["output_folder"]+"/multiqc/multiqc_report_gatk.html"
    singularity:
        config["SIF"]["multiqc"]
    params:
        input_f = config["output_folder"]+"/fastqc/",
        bam_f = config["output_folder"] + "/mapped/",
        gatk_f = config["output_folder"] + "/GATK_out/",
        output_f = config["output_folder"]+"/multiqc/"
    shell:
        """
        multiqc {params.input_f} {params.bam_f} {params.gatk_f} -o {params.output_f} -n multiqc_report_gatk
        """

rule bamqc:
    input:
        config["output_folder"]+"/mapped/{sample}_sorted_hg38.bam"
    output:
        config["output_folder"]+"/bamqc/{sample}/qualimapReport.html"
    singularity:
        config["SIF"]["qualimap"]
    params:
        output_f = config["output_folder"]+"/bamqc/{sample}/",
    threads:
        config["params"]["bamqc"]["threads"]
    shell:
        "qualimap bamqc -bam {input} -nt {threads} -outdir {params.output_f} --java-mem-size=120000M"

rule multiqc_bam:
    input:
        expand(config["output_folder"]+"/bamqc/{sample}/qualimapReport.html",sample=samples["sample_ID"])
    output:
        config["output_folder"]+"/multiqc_bam/multiqc_report.html"
    singularity:
        config["SIF"]["multiqc"]
    params:
        input_f = config["output_folder"]+"/bamqc/",
        output_f = config["output_folder"]+"/multiqc_bam/"
    shell:
        "multiqc {params.input_f} -o {params.output_f} -n multiqc_report"
