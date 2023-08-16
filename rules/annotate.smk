
#rule get_vep_cache:
#    output:
#        directory("../vep/cache"),
#    params:
#        species="Homo_sapiens",
#        build="GRCh38",
#        release="105"
#    shell:
#        "vep_install -a cf -s {params.species} -y {params.build} -c {output} --CONVERT"

rule download_vep_plugins:
    output:
        directory("../vep/plugins")
    params:
        release=100
    wrapper:
        "v1.1.0/bio/vep/plugins"

rule annotate_variants:
    input:
        calls=config["output_folder"]+"/GATK_out/{tumor}_{normal}_GATKFiltered_pass.vcf",
        plugins="../vep/plugins",
        # optionally add reference genome fasta
        fasta=config["params"]["gatk"]["ref"] 
    output:
        calls=config["output_folder"]+"/annotated/{tumor}_{normal}_variants.annotated.vcf",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["Frameshift","Wildtype"],
    conda:
        "../envs/test_anno.yaml"
    threads: 4
    shell:
        "vep --format vcf --input_file {input.calls} --dir_cache /DATA/j.traets/Reference_files/homo_sapiens/GRCh38_VEP_direct/ --output_file {output.calls} --fasta {input.fasta} --offline --plugin Frameshift --plugin Wildtype --vcf --symbol --terms SO --tsl --hgvs"
