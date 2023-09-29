rule vcf2maf: # TODO change into singularity
    input:
        calls=config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered-pass.vcf",
        fasta=config["params"]["gatk"]["ref"]
    output:
        mafs=config["output_folder"]+"/maf_files/{tumor}-vs-{normal}-GATKFiltered-pass.vep.maf",
    threads: config["params"]["vcf2maf"]["threads"]
    params:
        ref=config["params"]["vcf2maf"]["ref"]
    #singularity:
    #    config["SIF"]["vcf2maf"]
    conda:
        "../envs/anno.yaml"
    shell:
        "perl scripts/vcf2maf.pl --input-vcf {input.calls} --output-maf {output.mafs} --normal-id {wildcards.normal} --tumor-id {wildcards.tumor} --ncbi-build GRCh38 --ref-fasta {input.fasta} --vep-path /opt/vep/ --vep-path $CONDA_PREFIX/bin/ --vep-data {params.ref}"
        #"""
        #perl scripts/vcf2maf.pl --input-vcf {input.calls} --output-maf {output.mafs} --normal-id {wildcards.normal} --tumor-id {wildcards.tumor} --ncbi-build GRCh38 --ref-fasta {input.fasta} --vep-path /opt/vep/.vep  --vep-data {params.ref}
        #"""

rule unzip_indels:
    input:
        in_t=config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/results/variants/somatic.indels.vcf.gz",
    output:
        config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/results/variants/somatic.indels.vcf",
    shell:
        "gzip -kd {input.in_t}"

rule vcf2maf_indels:
    input:
        in_t=config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/results/variants/somatic.indels.vcf",
        fasta=config["params"]["strelka"]["ref"]
    output:
        mafs=config["output_folder"]+"/maf_files/{tumor}-vs-{normal}-strelka_indels.vep.maf",
    threads: 4
    conda:
        "../envs/anno.yaml"
    shell:
        "perl /DATA/j.traets/Tools/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf {input.in_t} --output-maf {output.mafs} --normal-id {wildcards.normal} --tumor-id {wildcards.tumor} --ncbi-build GRCh38 --ref-fasta {input.fasta} --vep-path $CONDA_PREFIX/bin/ --vep-data /DATA/j.traets/Reference_files/homo_sapiens/GRCh38_VEP_direct/"


