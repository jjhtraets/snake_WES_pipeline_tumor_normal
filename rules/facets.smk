rule run_snppile:
        input:
                tumor_bam = config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam",
                normal_bam = config["output_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam"
        params:
                common_snps = config["params"]["facets"]["common_snps"],
                snppile = config["params"]["facets"]["snppile"],
        singularity:
                config["SIF"]["snppile"]
        output:
                out_file = config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz"
        shell:
                """
		snp-pileup {params.snppile} {params.common_snps} {output.out_file} {input.normal_bam} {input.tumor_bam}
                """

rule run_facets_r:
        input:
                tumor_snp = config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz"
        params:
                config["output_folder"] + "/FACETS/fitted/",
                snppile = config["params"]["facets"]["snppile"],
                targets = config["params"]["facets"]["targets"],
                build = "hg38"
        singularity:
                config["SIF"]["facets"]
        log: 
                config["output_folder"]+"/logs/FACETS/{tumor}-vs-{normal}_run.log"
        output:
                out_file = config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz_fitted.csv"
        shell:
                """
                Rscript scripts/FACETS_run.R -i {input.tumor_snp} -o {output.out_file} -t {params.targets} -b {params.build}
                """

rule facets_pyclone:
        input:
                inp = config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz_fitted.csv"
        params:
        #        run_dir = config["params"]["pyclone"]["path_to_run"]
                output_dir = config["output_folder"] + "/FACETS/",
        #        indel_maf = config["output_folder"]+"/maf_files/{tumor}-vs-{normal}-strelka_indels.vep.maf",
                maf = config["output_folder"]+"/maf_files/filtered_{tumor}-vs-{normal}-GATKFiltered-pass.vep.maf",
        #        inp = config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz_fitted.csv"
        output:
                config["output_folder"] + "/FACETS/{tumor}-vs-{normal}_pyclone_maf_filtered_ready.tsv"
        shell:
                "python scripts/FACETS_to_Pyclone_input.py {params.maf} --facets_output {input.inp} {params.output_dir} --sample_id {wildcards.tumor}-vs-{wildcards.normal}"


