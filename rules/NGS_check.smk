
rule test_ngs:
    input:
        reads = [config["output_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam"]
    output:
        config["output_folder"]+"/NGScheckmate/{tumor}-vs-{normal}/output_matched.txt"
    conda:
        "../envs/NGS_check.yaml"
    params:
        file_l = write_file,
        output_d = config["output_folder"]+"/NGScheckmate/{tumor}-vs-{normal}/"
    shell:
        "python /DATA/j.traets/Tools/NGSCheckMate/NGSCheckMate/ncm.py -B -l {params.file_l} -bed /DATA/j.traets/Tools/NGSCheckMate/NGSCheckMate/SNP/SNP_GRCh38_hg38_woChr.bed -O {params.output_d}"
