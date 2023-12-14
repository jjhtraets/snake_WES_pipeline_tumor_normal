rule CNVkit_access:
    input:
        genome = config["params"]["CNVkit"]["ref"]
    output:
        access = config["output_folder"]+"/CNVkit/access.hg38.bed"
    singularity:
        config["SIF"]["cnvkit"]
    log:
        config["output_folder"] + "/logs/CNVkit/access_CNVkit_run.log"
    shell:
        """
        cnvkit.py access {input.genome} -o {output.access} 2> {log}
        """

rule CNVkit_autobin:
    input:
        bam_files = expand(config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam",sample=samples["sample_ID"]),
        access = config["output_folder"]+"/CNVkit/access.hg38.bed"
    params:
        bait = config["params"]["CNVkit"]["bait_bed"],
        refflat = config["params"]["CNVkit"]["refflat"],
    output:
        target = config["output_folder"]+"/CNVkit/target.bed",
        antitarget = config["output_folder"]+"/CNVkit/antitarget.bed",
    singularity:
        config["SIF"]["cnvkit"]
    log:
        config["output_folder"] + "/logs/CNVkit/autobin_CNVkit_run.log"
    shell:
        """
        cnvkit.py autobin {input.bam_files} -t {params.bait} -g {input.access} --annotate {params.refflat} --target-output-bed {output.target} --antitarget-output-bed {output.antitarget} 2> {log}
        """

rule CNVkit_reference:
    input:
        bam_normal = normals_list,
        genome = config["params"]["CNVkit"]["ref"],
        target = config["output_folder"]+"/CNVkit/target.bed",
        antitarget = config["output_folder"]+"/CNVkit/antitarget.bed",
        access = config["output_folder"]+"/CNVkit/access.hg38.bed"
    output:
        config["output_folder"]+"/CNVkit/results_normal/normals_reference.cnn"
    singularity: 
        config["SIF"]["cnvkit"]
    params:
        output_f = config["output_folder"]+"/CNVkit/results_normal/"
    threads:
        config["params"]["CNVkit"]["threads"]
    log:
        config["output_folder"]+"/logs/CNVkit/CNVkit_normal.log"
    shell:
        """
        cnvkit.py batch --normal {input.bam_normal} --output-reference {output} --targets {input.target} -g {input.access} -f {input.genome} -a {input.antitarget} --output-dir {params.output_f} -p {threads} 2> {log}
        """       

rule CNVkit_run:
    input:
        bam_tumor = config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam",
        ref = config["output_folder"]+"/CNVkit/results_normal/normals_reference.cnn"
    output:
        out = config["output_folder"]+"/CNVkit/results/{tumor}/{tumor}_sorted_hg38_ARRG_dedup_recal.cnr"
    singularity: 
        config["SIF"]["cnvkit"]
    log:
        config["output_folder"]+"/logs/CNVkit/{tumor}_CNVkit_run.log"
    params:
        output_f = config["output_folder"]+"/CNVkit/results/{tumor}/"
    shell:
        """
        cnvkit.py batch {input.bam_tumor} -r {input.ref} --output-dir {params.output_f} --diagram --scatter 2> {log}
        """

rule CNVkit_purity:
    input:
        cns = config["output_folder"]+"/CNVkit/results/{tumor}/{tumor}_sorted_hg38_ARRG_dedup_recal.cns",
        facets = config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz_fitted.csv"
    output:
        config["output_folder"]+"/CNVkit/results/{tumor}/{tumor}-vs-{normal}_sorted_hg38_ARRG_dedup_recal_purity_thr.call.cns"
    params:
        purity = get_tumor_purity,
        ploidy = get_tumor_ploidy,
        thresholds = get_threshold
    singularity: 
        config["SIF"]["cnvkit"]
    shell:
        """
        cnvkit.py call {input.cns} -y -m threshold -t={params.thresholds} --purity {params.purity} -o {output}
        """

rule export_CNVkit:
    input:
        config["output_folder"]+"/CNVkit/results/{tumor}-vs-{normal}/{tumor}_sorted_hg38_ARRG_dedup_recal.cnr"
    output:
        config["output_folder"]+"/CNVkit/results/{tumor}-vs-{normal}/{tumor}_sorted_hg38_ARRG_dedup_recal.seg"
    singularity: 
        config["SIF"]["cnvkit"]
    shell:
        """
        cnvkit.py export seg {input} --enumerate-chroms -o {output}
        """
