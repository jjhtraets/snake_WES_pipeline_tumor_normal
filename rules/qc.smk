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
        "fastqc -o {params.input_f} -t {threads} {input.reads}"

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
        "multiqc {params.input_f} -o {params.output_f} -n multiqc_report"

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
