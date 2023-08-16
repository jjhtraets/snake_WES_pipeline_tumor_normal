# build containers from docker images when not present
# default containers listed for the DNAseq pipeline v2 (2307)

import os

if __name__ == "__main__":
  print("build containers from docker images")
  print("WARING: requires sudo rights")

  if (os.path.isfile("sif/fastqc_v0.11.9_cv8.sif")) == (not True):
    os.system("singularity build sif/fastqc_v0.11.9_cv8.sif docker://biocontainers/fastqc:v0.11.9_cv8")
  else:
    print("fastqc_v0.11.9_cv8.sif already present")

  if (os.path.isfile("sif/multiqc_v1.8_staphb.sif")) == (not True):
    os.system("singularity build sif/multiqc_v1.8_staphb.sif docker://staphb/multiqc:1.8")
  else:
    print("multiqc_v1.8_staphb.sif already present")

  if (os.path.isfile("sif/bwa_v0.7.17_biocontainers.sif")) == (not True):
    os.system("singularity build sif/bwa_v0.7.17_biocontainers.sif docker://biocontainers/bwa:v0.7.17_cv1")
  else:
    print("bwa_v0.7.17_biocontainers.sif already present")

  if (os.path.isfile("sif/samtools_v1.9.4deb_biocontainers.sif")) == (not True):
    os.system("singularity build sif/samtools_v1.9.4deb_biocontainers.sif docker://biocontainers/samtools:v1.9-4-deb_cv1")
  else:
    print("samtools_v1.9.4deb_biocontainers.sif already present")
    
  if (os.path.isfile("sif/picard_v2.27.5_broadinstitute.sif")) == (not True):
    os.system("singularity build sif/picard_v2.27.5_broadinstitute.sif docker://broadinstitute/picard:2.27.5")
  else:
    print("picard_v2.27.5_broadinstitute.sif already present")

  if (os.path.isfile("sif/gatk_v4.4.0_broadinstitute.sif")) == (not True):
    os.system("singularity build sif/gatk_v4.4.0_broadinstitute.sif docker://broadinstitute/gatk:4.4.0.0")
  else:
    print("gatk_v4.4.0_broadinstitute.sif already present")
    
  if (os.path.isfile("sif/cnvkit_0.9.10_etal.sif")) == (not True):
    os.system("singularity build sif/cnvkit_0.9.10_etal.sif docker://etal/cnvkit:0.9.10")
  else:
    print("cnvkit_0.9.10_etal.sif already present")

  if (os.path.isfile("sif/qualimap_v2.2.1_pegi3s.sif")) == (not True):
    os.system("singularity build sif/qualimap_v2.2.1_pegi3s.sif docker://pegi3s/qualimap:2.2.1")
  else:
    print("qualimap_v2.2.1_pegi3s.sif already present")
    
  if (os.path.isfile("sif/snppile_v0.6.2_biocontainers.sif")) == (not True):
    os.system("singularity build sif/snppile_v0.6.2_biocontainers.sif docker://quay.io/biocontainers/snp-pileup:0.6.2--h6b7c446_4")
  else:
    print("snppile_v0.6.2_biocontainers.sif already present")

  if (os.path.isfile("sif/facets_v0.6.2_biocontainers.sif")) == (not True):
    os.system("singularity build sif/facets_v0.6.2_biocontainers.sif docker://quay.io/biocontainers/r-facets:0.6.2--r42h1c9e865_3")
  else:
    print("facets_v0.6.2_biocontainers.sif already present")


