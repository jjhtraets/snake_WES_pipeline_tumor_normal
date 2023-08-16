rule gatk_mutect2:
    input:
        bam_tumor = config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam",
        bam_normal = config["output_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",
        genome = config["params"]["gatk"]["ref"],
        ref_gatk = config["params"]["gatk"]["gnomad"],
    output:
        config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}.vcf"
    singularity:
        "docker://broadinstitute/gatk:4.4.0.0"
    threads:
        config["params"]["gatk"]["threads"]
    log:
        config["output_folder"]+"/logs/GATK/{tumor}-vs-{normal}_mutect2.log"
    shell:
        """
        gatk Mutect2 -R {input.genome} -I {input.bam_tumor} -tumor {wildcards.tumor} -I {input.bam_normal} -normal {wildcards.normal} -O {output} --germline-resource {input.ref_gatk}  2> {log}
        """

rule gatk_filter:
    input:
        VCFs = config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}.vcf",
        genome = config["params"]["gatk"]["ref"] 
    output:
        config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered.vcf"
    singularity:
        "docker://broadinstitute/gatk:4.4.0.0"
    threads:
        config["params"]["gatk"]["threads"]
    log:
        config["output_folder"]+"/logs/GATK/{tumor}-vs-{normal}_mutect2_filter.log"
    shell:
        """
        gatk FilterMutectCalls -V {input.VCFs} -O {output} -R {input.genome} 2> {log}
        """

rule keep_pass_SNPs:
    input:
        VCFs = config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered.vcf",
        genome = config["params"]["gatk"]["ref"]  
    output:
        config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered-pass.vcf"
    singularity:
        "docker://broadinstitute/gatk:4.4.0.0"
    threads:
        config["params"]["gatk"]["threads"]
    log:
        config["output_folder"]+"/logs/GATK/{tumor}-vs-{normal}_mutect2_filter_pass.log"
    shell:
        """
        gatk SelectVariants -R {input.genome} --variant {input.VCFs} -O {output} --exclude-filtered 2> {log}
        """

