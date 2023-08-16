rule run_pyclone_facets_input:
        input:
                tsv = config["output_folder"] + "/FACETS/{tumor}-vs-{normal}_pyclone_maf_filtered_ready.tsv"
        params:
                tumor_purity = get_tumor_purity,
                work_dir = config["output_folder"] + "/pyclone/{tumor}-vs-{normal}/"
        # fix environment, graph plotting issue
        #conda:
        #        "../envs/pyclone.yaml"
        output:
                config["output_folder"] + "/pyclone/{tumor}-vs-{normal}/tables/cluster.tsv"
        shell:
                "PyClone run_analysis_pipeline --in_files {input.tsv} --working_dir {params.work_dir} --tumour_contents {params.tumor_purity} --burnin 1000 --seed 100 --samples {wildcards.tumor}"

