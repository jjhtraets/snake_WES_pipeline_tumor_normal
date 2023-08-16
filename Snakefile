###################################
#### Snakemake pipeline DNAseq ####
###################################

container: "docker://continuumio/miniconda3:latest"

include: "rules/common.smk"
include: "rules/qc.smk"
include: "rules/mapping.smk"
include: "rules/variant.smk"
include: "rules/strelka.smk"
include: "rules/CNVkit.smk"
include: "rules/vcf2maf.smk"
include: "rules/NGS_check.smk"
include: "rules/facets.smk"
include: "rules/manta.smk"
#include: "rules/haplotyper.smk"

rule all:
    input:
        output_rules_all(),
