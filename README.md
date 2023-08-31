# WES pipeline (snakemake), for matched tumor & normal
version v2.0 (singularity/conda)
last update, August 2023

Input: paired end  reads, fastq or bam files, WES data
Output: SNVs/indels (Mutect2 GATK, Strelka), CNVs (CNVkit, FACETS), annotation (VEP, maf), QCs (checkmate, bamqc, fastqc, multiqc)

Optional, not integrated(!): Haplotyper, Pyclone, Vardict, Manta

## ⚠️⚠️ work in process! ⚠️⚠️ ##

## Prepare sif containers
- run `Buid_sif_containers.py` (default package versions used in WES pipeline 23/06)
- or build your own sif containers with selfbuild or pulled docker images


## How to run with the Run_DNApipeline.py shell
1. create a new conda env (or use one in which only snakemake is installed)
```bash
conda create -n snakemake_run snakemake
```
2. activate new environment
```bash
conda activate snakemake_run
```
3. configure config file:
`config/config.yaml`
4. run "Run_RNApipeline.py"
```bash
  python Run_RNApipeline.py
```

```bash

  options:
  -h, --help            show this help message and exit
  -c CORES, --cores CORES
                        Number of cores to be used (int) (default: 1)
```
Example:
```bash
python Run_RNApipeline.py -c 5
```

## Alternative: How to run without the Run_RNApipeline.py shell
1. activate snakemake environment
2. configure config files:
`config/config.yaml`
`var/config.yaml`
3. run snakemake in snakemake RNAseq pipeline folder, adjust amount of cores
etc. 
```bash
snakemake --use-conda --use-singularity --singularity-args "-B $PATH_INPUT -B $PATH_OUTPUT -B $PATH_REF" --cores 1 -k
```

### Sample file format 
Annotation file should contain "sample_ID" (name of the fastq files without [READ]fastq.gz), "Tumor_yes" (YES or NO), "Patient" (same for matched tumor and normal) column names. The format of the file should be ";" separated. See example `samples_test.txt`

### Information on config/config.yaml
"input_folder" in the `config/config/yaml` file locates the folder containing the fastq files. And "output_folder" locates the folder for the output files. Setting the modes to True/False determines the output.

### Updates
- February 2023, adding pyclone
- November 2022, fixed bugs
- September 2022, added haplotypecaller GATK, Manta SVs
- August 2022, debugging pipeline
- June 2022, CNVkit purity, ploidy, threshold for call correction, added logs
- May 2022, added FACETS
- March 2022, fixed sample name issues

### TODO
- make sample file more userfriendly
- add slurm options
- replace all conda environments with singularity containers
- add all quality checks
- add modes (haplotype caller etc)
