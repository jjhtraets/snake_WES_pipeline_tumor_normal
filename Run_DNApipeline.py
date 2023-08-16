#                              WES snakemake pipeline tumor/normal                             #
#                              ### script to run DNA pipeline ###                              #
# WES data
# Calling SNVs, CNVs, indels, SVs
# Tumor - normal samples
# Singularity

import sys
import subprocess
import argparse
import os
import yaml
import pandas as pd
import logging
from datetime import datetime

# config file location fixed
def load_yaml():
    with open("./config/config.yaml", "r") as stream:
        try:
            safe_config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return safe_config

def check_path_format():
    if yaml_config["input_folder"][-1] != "/":
        yaml_config["input_folder"] = yaml_config["input_folder"]+"/"
    if yaml_config["output_folder"][-1] != "/":
        yaml_config["output_folder"] = yaml_config["output_folder"]+"/"

def check_sample_file():
    if os.path.exists(yaml_config["samples"]):
        with open(yaml_config["samples"], "r") as file:
            header = file.readline().strip()
            if "sample_ID" in header:
                return "OK"
            else:
                return "sample_ID column name is missing!"
            if "Tumor_yes" in header:
                return "OK"
            else:
                return "Tumor_yes column name is missing!"
            if "Patient" in header:
                return "OK"
            else:
                return "Patient column name is missing!"
    else:
        return "Sample file missing!"
      
def check_mode_config_yaml():
    mode_config = ", ".join(key for key, value in yaml_config["run_modes"].items() if value == True)
    return mode_config
  
# run snakemakepipeline with OS command based on user input
def run_snakemake(yaml_config,config):
    # prepare bind volumes
    bind_d1 = yaml_config["params"]["bwa"]["ref"]
    bind_d1 = bind_d1.split("/")
    bind_d1.pop()
    bind_d1 = "/".join(bind_d1)
    
    bind_d2 = yaml_config["input_folder"]
    bind_d2 = bind_d2.split("/")
    bind_d2.pop()
    bind_d2 = "/".join(bind_d2)
    
    bind_d3 = yaml_config["output_folder"]
    bind_d3 = bind_d3.split("/")
    bind_d3.pop()
    bind_d3 = "/".join(bind_d3)
    
    bind_d4 = yaml_config["params"]["CNVkit"]["bait_bed"]
    bind_d4 = bind_d4.split("/")
    bind_d4.pop()
    bind_d4 = "/".join(bind_d4)

    log.info("bind singularity:")
    log.info(bind_d1)
    log.info(bind_d2)
    log.info(bind_d3)
    log.info(bind_d4)
    
    logging.shutdown()

    # run snakemake
    os.system('snakemake --use-singularity --singularity-args "-B ' + bind_d1 + ' -B ' + bind_d2  +  ' -B ' +  bind_d3 + ' -B ' + bind_d4 + ' " --cores ' + str(config["cores"]) + ' ' + str(yaml_config["snakemake"]))
#--use-conda --conda-frontend conda

if __name__ == "__main__":
  os.system("cat banner_cmd.txt")
  
    # prepare logging, stdout and log file
  now = datetime.now()
  file_path = now.strftime("%Y-%m-%d_%H:%M:%S.log")
  print(file_path)

  log = logging.getLogger()
  logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
  log.setLevel(logging.DEBUG)

  fh = logging.FileHandler(file_path)
  fh.setFormatter(logFormatter)
  log.addHandler(fh)

  ch = logging.StreamHandler()
  ch.setFormatter(logFormatter)
  log.addHandler(ch)
  
  # TODO: add slurm in the future
  parser = argparse.ArgumentParser(description="Snakemake WES TN - Pipeline",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-c", "--cores",default=1, help="Number of cores to be used (int)")
  args = parser.parse_args()
  config = vars(args)
    
  # print input user
  log.info("\n \n ### User input ### ")
  log.info("Nr of cores: " + str(config["cores"]))
  
  yaml_config = load_yaml()
  check_path_format()
    
  # print most important parameters in config file
  log.info("\n \n ### Parameters in config.yaml file ### ")
  log.info("Sample file: " + yaml_config["samples"])
  log.info("Read name 1: *" + yaml_config["read_name"]["R1"] + ".fastq.gz")
  log.info("Read name 2: *" + yaml_config["read_name"]["R2"] + ".fastq.gz")
  log.info("Reference: " + yaml_config["params"]["bwa"]["ref"])
  log.info("Target file: " + yaml_config["params"]["CNVkit"]["bait_bed"]) 
  log.info("Access file: " + yaml_config["params"]["CNVkit"]["access_bed"])
  log.info("Known sites (GATK): " + yaml_config["params"]["gatk"]["known_sites"])
  log.info("Gnomad (GATK): " + yaml_config["params"]["gatk"]["gnomad"] )
  log.info("Additional arguments Snakemake: " + yaml_config["snakemake"] )
  run_modes_yaml = check_mode_config_yaml()
  log.info("Run mode: " + str(run_modes_yaml))
  
  # general checks on yaml file
  check_samples_return = check_sample_file()
  log.info("\n \n ### Checks ### ")
  log.info("Sample file: " + str(check_samples_return))
  
  # ask user to continue
  startRun = True
  while startRun:
      var = input("Do you want to continue with the current settings? (y/n): ")
      if var != 'y' and var != 'n':
          print("Sorry, I didn't catch that!")
      else:
          startRun = False
  
  # TODO bind more volumes when bed files etc are outside the current bind location
  if var == "y":
      run_snakemake(yaml_config,config)
  elif var == "n":
      sys.exit("Change config/config.yaml file first before rerunning")


