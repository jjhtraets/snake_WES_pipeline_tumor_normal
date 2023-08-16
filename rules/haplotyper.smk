rule bam_index:
    input:
        bam = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam",
        genome = config["params"]["gatk"]["ref"]
    output:
        config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bai"
    conda:
        "../envs/mapping.yaml"
    params:
        par1 = config["params"]["picard"]["java_set"],
        par2 = config["params"]["picard"]["jar_loc"],
    shell:
        "{params.par1} {params.par2} BuildBamIndex INPUT={input.bam} OUTPUT={output}"

rule gatk_haplocal:
    input:
        bam = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam",
        genome = config["params"]["gatk"]["ref"],
        known_sites = config["params"]["gatk"]["known_sites"]
    output:
        config["output_folder"]+"/haplo/{sample}_sorted_hg38_ARRG_dedup_recal.vcf"
    conda:
        "../envs/calling.yaml"
    log:
        config["output_folder"]+"/logs/GATK/{sample}_gatk_haplo.log"
    shell:
        "gatk HaplotypeCaller -R {input.genome} -I {input.bam} --dbsnp {input.known_sites} --output {output} 2> {log}"

rule gatk_haplo_filter:
    input:
        vcf = config["output_folder"]+"/haplo/{sample}_sorted_hg38_ARRG_dedup_recal.vcf"
    output:
        config["output_folder"]+"/haplo/{sample}_sorted_hg38_ARRG_dedup_recal_snp.vcf"
    conda:
        "../envs/calling.yaml"
    shell:
        "gatk SelectVariants -R {input.genome} -V {input.vcf} --output {output} --select-type-to-include SNP"

rule gatk_haplo_hard_filter:
    input:
        vcf = config["output_folder"]+"/haplo/{sample}_sorted_hg38_ARRG_dedup_recal_snp.vcf"
    output:
        config["output_folder"]+"/haplo/{sample}_sorted_hg38_ARRG_dedup_recal_snp_hard_filter.vcf"
    conda:
        "../envs/calling.yaml"
    params:
        filter = config["params"]["gatk"]["filtering"]["snvs"]
    shell:
        "gatk VariantFiltration -R {input.genome} -V {input.vcf} --output {output} --filterExpression {params.filter}"

rule recalibrate_calls:
    input:
        vcf = config["output_folder"]+"/haplo/{sample}_sorted_hg38_ARRG_dedup_recal_snp.vcf"
    output:
        config["output_folder"]+"/haplo/{sample}_sorted_hg38_ARRG_dedup_recal_snp_recal.vcf"
    conda:
        "../envs/calling.yaml"
    shell:
        "gatk VariantRecalibrator -R {input.genome} -V {input.vcf} --output {output} -mode SNP"


