# common script
import glob
import itertools
import pandas as pd
import numpy as np

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"],sep=";")
#print(set(samples["sample_ID"]))

# function to create new data frame, tumor vs normal
def creat_matched(samples):
    norm_list = []
    tumor_list = []
    for x_p in set(samples["Patient"]):
        sel_pat = samples[samples["Patient"]==x_p]
        norm_list.append(list(list(sel_pat[sel_pat["Tumor_yes"]=="NO"]["sample_ID"])*len(sel_pat[sel_pat["Tumor_yes"]=="YES"]["sample_ID"])))
        tumor_list.append(list(sel_pat[sel_pat["Tumor_yes"]=="YES"]["sample_ID"]))
    norm_list = list(itertools.chain(*norm_list))
    tumor_list = list(itertools.chain(*tumor_list))
    matched_col = pd.DataFrame(list(zip(norm_list,tumor_list)),columns =['Normal', 'Tumor'])
    return(matched_col)

samples_matched = creat_matched(samples)
#print("Check data frame!")
#print(samples_matched["Tumor"][0])

# collect output
def output_rules_all():
    samples_matched = creat_matched(samples)

    bam_output = expand(config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam",sample=set(samples["sample_ID"])),

    qc_output = config["output_folder"]+"/multiqc/multiqc_report.html",
    qc_output_gatk = config["output_folder"]+"/multiqc/multiqc_report_gatk.html",
    qc_output_b = config["output_folder"]+"/multiqc_bam/multiqc_report.html",

    variant_GATK_output = expand(config["output_folder"]+"/GATK_out/{tumor}-vs-{normal}-GATKFiltered-pass.vcf",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])    
    variant_strelka_output = expand(config["output_folder"]+"/strelka_out/{tumor}-vs-{normal}/results/variants/somatic.snvs.vcf.gz",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])
   
    cnv_normal =config["output_folder"]+"/CNVkit/results_normal/normals_reference.cnn"
    cnv_output = expand(config["output_folder"]+"/CNVkit/results/{tumor}/{tumor}_sorted_hg38_ARRG_dedup_recal.cnr",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])
    cnv_purity = expand(config["output_folder"]+"/CNVkit/results/{tumor}/{tumor}-vs-{normal}_sorted_hg38_ARRG_dedup_recal_purity_thr.call.cns",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])

    maf_files = expand(config["output_folder"]+"/maf_files/{tumor}-vs-{normal}-GATKFiltered-pass.vep.maf",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])
    maf_files_indels = expand(config["output_folder"]+"/maf_files/{tumor}-vs-{normal}-strelka_indels.vep.maf",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])

    NGS_check_files = expand(config["output_folder"]+"/NGScheckmate/{tumor}-vs-{normal}/output_matched.txt",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])

    facets = expand(config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz_fitted.csv",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])

    # requires filtering of maf, TODO add r script
    pyclone_input = expand(config["output_folder"] + "/FACETS/{tumor}-vs-{normal}_pyclone_maf_filtered_ready.tsv",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])
    # TODO fix env pyclone
    pyclone_output = expand(config["output_folder"] + "/pyclone/{tumor}-vs-{normal}/tables/cluster.tsv",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])

    haplo = expand(config["output_folder"]+"/haplo/{sample}_sorted_hg38_ARRG_dedup_recal_snp_recal.vcf",sample=set(samples["sample_ID"]))

    manta = expand(config["output_folder"]+"/Manta/{tumor}-vs-{normal}/results/variants/somaticSV.vcf.gz",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])
    
    facets = expand(config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz_fitted.csv",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])

    modes = list()
    if config["run_modes"]["gatk"] == True:
      modes.append([variant_GATK_output,maf_files])
    if config["run_modes"]["strelka"] == True:
      modes.append(variant_strelka_output)
    if config["run_modes"]["QCs"] == True:
        if config["run_modes"]["QCs"] == True and config["run_modes"]["gatk"] == False:
            modes.append(qc_output)
        if config["run_modes"]["QCs"] == True and config["run_modes"]["gatk"] == True:
            modes.append(qc_output_gatk)
    if config["run_modes"]["CNV"] == True:
      modes.append([cnv_normal,cnv_output,facets])

    return modes
    
# possible to add more rules to output
# examples:
#    return bam_output, qc_output, qc_output_b, variant_GATK_output, variant_strelka_output, variant_strelka_output_idx,cnv_normal, cnv_output, maf_files,maf_files_indels,NGS_check_files,facets,pyclone_input, variant_varscan_mpile,variant_varscan_output
#    return bam_output, qc_output, variant_GATK_output, maf_files, variant_strelka_output, variant_strelka_output_idx, variant_GATK_output, qc_output,qc_output_b, variant_strelka_output, cnv_output,variant_varscan_mpile,variant_varscan_output, facets, pyclone_output


### common functions

def get_fastq_files(wildcards):
    R1 = config["read_name"]["R1"]
    R2 = config["read_name"]["R2"]
    fastq_1 = config["input_folder"]+"/"+wildcards[0] + R1 +".fastq.gz"
    fastq_2 = config["input_folder"]+"/"+wildcards[0] + R2 +".fastq.gz"      
    return([fastq_1,fastq_2])

def get_output_qc_fastq_R1(wildcards):
    R1 = config["read_name"]["R1"]
    fastq_1 = config["output_folder"]+"/fastqc/"+wildcards[0] + R1 +"_fastqc.html"
    return(fastq_1)

def get_output_qc_fastq_R2(wildcards):
    R2 = config["read_name"]["R2"]
    fastq_2 = config["output_folder"]+"/fastqc/"+wildcards[0] + R2 +"_fastq.html"
    return(fastq_2)

def write_file(wildcards):
    temp_file = open("/DATA/j.traets/Tools/NGSCheckMate/bam_list.txt","w")
    temp_file.write(config["output_folder"]+"/mapped/"+wildcards[0]+"_sorted_hg38_ARRG_dedup_recal.bam"+"\n")
    temp_file.write(config["output_folder"]+"/mapped/"+wildcards[1]+"_sorted_hg38_ARRG_dedup_recal.bam")
    temp_file.close()
    return("/DATA/j.traets/Tools/NGSCheckMate/bam_list.txt")

def get_tumor_purity(wildcards):
    temp_file = pd.read_csv(config["output_folder"]+"/FACETS/fitted/"+wildcards[0]+"-vs-"+wildcards[1]+".snppile.csv.gz_fitted.csv",sep=',')
    purity = temp_file["Purity"][0]
    return(purity)

def get_tumor_ploidy(wildcards):
    temp_file = pd.read_csv(config["output_folder"]+"/FACETS/fitted/"+wildcards[0]+"-vs-"+wildcards[1]+".snppile.csv.gz_fitted.csv",sep=',')
    ploidy = temp_file["Ploidy"][0]
    return(ploidy)

def get_target(wildcards):
    loc_bed_file = config["params"]["CNVkit"]["bait_bed"]
    name_bed = loc_bed_file.split("/")[len(loc_bed_file.split("/"))-1]
    #loc_bed = loc_bed_file.split("/")[1:(len(loc_bed_file.split("/"))-1)]
    #loc_bed = "/".join(loc_bed)
    return( str(name_bed[:-4]) + ".target.bed")
    
def get_antitarget(wildcards):
    loc_bed_file = config["params"]["CNVkit"]["bait_bed"]
    name_bed = loc_bed_file.split("/")[len(loc_bed_file.split("/"))-1]
    #loc_bed = loc_bed_file.split("/")[1:(len(loc_bed_file.split("/"))-1)]
    #loc_bed = "/".join(loc_bed)
    return( str(name_bed[:-4]) + ".antitarget.bed")

def get_threshold(wildcards):
    temp_file = pd.read_csv(config["output_folder"]+"/FACETS/fitted/"+wildcards[0]+"-vs-"+wildcards[1]+".snppile.csv.gz_fitted.csv",sep=',')
    ploidy = temp_file["Ploidy"][0]
    purity = temp_file["Purity"][0]
    thresholds = np.log2((1 - float(purity)) + float(purity) * (np.arange(4) + 0.5) / float(ploidy))
    return(str(thresholds[0])+","+str(thresholds[1])+","+str(thresholds[2])+","+str(thresholds[3]))

def normals_list(wildcards):
    bam_normals = expand(config["output_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",normal=samples_matched["Normal"])
    return(list(set(bam_normals)))

