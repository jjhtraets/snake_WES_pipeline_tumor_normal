# DNAseq pipeline #
April 2023
version v2, alpha
Conda > Docker 

### TODO
Use Docker/Singularity instead of conda
- GATK/bwa
- Sterlka
- CNVkit/FACETS/Manta
- Pyclone
- loFreq
- varDict
- QCs

Add two more variant callers
- loFreq
- varDict

Make use of wrappers to clean up scripts

FIX:
- NGS match check on all samples
- FACETS R script clean up


### Old
Input: WES data
Output: QC's (multiQC on fastq's + BAMs & match check), CNVs (CNVkit & FACETS), SNVs (GATK & Strelka & Varscan), Indels (Strelka), clones (Pyclone), germline SNPs(+short indels), SVs (Manta)


### renew how to run
### TODO use shell/python script combine with singularity volume assign

#### How to run
1. create a new conda env (or use one in which only snakemake is installed)
conda create -n snakemake_run snakemake

2. activate new environment
conda activate snakemake_run

3. create annotation file, see example "sample_list.txt". Important: keep column names!

4. configure config file in config/config.yaml

5. configure output files in rules/common.smk, output_rules_all()

6. run snakemake in (copied) snakemake DNAseq pipeline folder, adjust amount of cores
snakemake --use-conda --cores 40

(remove .git folder from copied snakemake folder, rm -r if .git)

#### Data manangement
- to save space on the server, only keeps sorted.bam, dedup.bam and dedup_recal.bam.

#### Known issues
- new mamba integration in snakemake can be buggy from time to time, use --conda-frontend conda if it is causing issues
- if continuously stuck on "solving env", try to update conda

### Annotation file format
Annotation file should contain "sample_ID" (name of the fastq files without [READ]fastq.gz), "Tumor_yes" (YES or NO), "Patient" (same for matched tumor and normal) column names. The format of the file should be ";" separated.


## Latest fixes/additions:
- February 2023, adding pyclone, WIP
- November 2022, fixed bugs
- September 2022, added haplotypecaller GATK, Manta SVs
- August 2022, debugging pipeline, stable parts listed on top
- June 2022, CNVkit purity, ploidy, threshold for call correction, added logs, somewhat more userfriendly
- May 2022, added FACETS
- March 2022, fixed sample name issues with "_", replaced with "-vs-"


## TODO
- Fix FACETS environment conda
- Fix pyclone environment conda, remove graph plotting
- Create different modes for running pipeline

## Gitlab
https://gitlab.rhpc.nki.nl/j.traets/dna_snakemake_pipeline

