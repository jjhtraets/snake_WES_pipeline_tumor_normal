
  
rule mpileup_ngs:
    input:
        reads = config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam"
    output:
        config["output_folder"]+"/NGScheckmate/{sample}_output_matched.vcf"
    singularity:
        config["SIF"]["checkmate"]
    params:
        SNP_ngs = config["params"]["ngscheck"]["SNP"],
        ref = config["params"]["gatk"]["ref"]
    shell:
        "samtools mpileup -I -uf {params.ref} -l {params.SNP_ngs} {input.reads} | bcftools call -c - > {output}"


rule test_mates:
    input:
        bam_files = [config["output_folder"]+"/NGScheckmate/{normal}_output_matched.vcf",config["output_folder"]+"/NGScheckmate/{tumor}_output_matched.vcf"],
        ref = config["params"]["gatk"]["ref"]
    output:
        config["output_folder"]+"/NGScheckmate/{tumor}-vs-{normal}/output_matched.txt"
    singularity:
        config["SIF"]["checkmate"]
    params:
        file_l = write_file,
        output_d = config["output_folder"]+"/NGScheckmate/{tumor}-vs-{normal}/",
        SNP_ngs = config["params"]["ngscheck"]["SNP"],
    shell:
        "python ./NGS_resources/ncm.py -V -l {params.file_l} -bed {params.SNP_ngs} -O {params.output_d}"

