########
# Import required packages
########

import pandas as pd
from pathlib import Path
import multiprocessing
import hashlib
import sys
import tempfile
import os
import subprocess
import re
import glob

######
# Check genopred conda env is activated
######

conda_env_name = os.getenv('CONDA_DEFAULT_ENV')
if not conda_env_name == 'genopred':
  print("Error: The genopred conda environment must be active when running the pipeline.\nFor more information: https://opain.github.io/GenoPred/pipeline_readme.html#Step_2:_Create_conda_environment_for_pipeline")
  sys.exit(1)

########
# Create required functions
########

# Create function to check whether path of gwas or score exist
def check_list_paths(df):
  for index, row in df.iterrows():
    file_path = row['path']
    if pd.isna(file_path):
        continue
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

# Create function to return the range of chromosomes requested
def get_chr_range(testing):
  if testing != 'NA':
    val = int(testing[-2:])
    return range(val, val + 1)
  else:
    return range(1, 23)

# Create function to check whether path of gwas or score exist
def check_target_paths(df, chr):
  for index, row in df.iterrows():
    if row['type'] == '23andMe':
      file_path = row['path']
    if row['type'] == 'plink1':
      file_path =  row['path'] + ".chr" + chr + ".bed"
    if row['type'] == 'plink2':
      file_path =  row['path'] + ".chr" + chr + ".pgen"
    if row['type'] == 'bgen':
      file_path =  row['path'] + ".chr" + chr + ".bgen"
    if row['type'] == 'vcf':
      file_path =  row['path'] + ".chr" + chr + ".vcf.gz"
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    else :
      return []

# Create function that checks the type column given in the target_list
def check_target_type(df, column='type'):
    valid_formats = {'plink1', 'plink2', 'vcf', 'bgen', '23andMe'}
    invalid_formats = df[~df[column].isin(valid_formats)]

    if not invalid_formats.empty:
        raise ValueError(f"Invalid format entries found in column '{column}': {invalid_formats[column].unique()}. Must be either 'plink1', 'plink2', 'vcf', 'bgen' or '23andMe'")

def check_target_type(df, column='type'):
    valid_formats = {'plink1', 'plink2', 'vcf', 'bgen', '23andMe'}

    # Check if the dataframe is empty
    if df.empty:
        return

    # Check if the column exists
    if column not in df.columns:
        return

    # Check for invalid formats
    invalid_formats = df[~df[column].isin(valid_formats)]

    if not invalid_formats.empty:
        raise ValueError(f"Invalid format entries found in column '{column}': {invalid_formats[column].unique()}. Must be either 'plink1', 'plink2', 'vcf', 'bgen' or '23andMe'")

######
# Check config file
######

def check_config_parameters(config):
    required_params = [
        'outdir', 'resdir', 'refdir', 'config_file', 'gwas_list', 'target_list',
        'score_list', 'pgs_methods', 'ptclump_pts', 'dbslmm_h2f', 'prscs_phi',
        'prscs_ldref', 'ldpred2_model', 'ldpred2_inference', 'ancestry_prob_thresh',
        'testing'
    ]

    missing_params = []
    for param in required_params:
        if config.get(param) is None:
            missing_params.append(param)

    if missing_params:
        print("Error: Missing parameters in user-specified and default config files:", missing_params)
        sys.exit(1)

# Check the config
check_config_parameters(config)

# Set outdir parameter
outdir=config['outdir']

#######
# Check config files
#######

# Function to check for duplicate names in a dataframe
def check_for_duplicates(df, name_col, list_name):
    duplicate_names = df[df[name_col].duplicated(keep=False)]
    if not duplicate_names.empty:
        raise ValueError(f"Duplicate values found in '{name_col}' column of {list_name}: {', '.join(duplicate_names[name_col].unique())}")


###
# target_list
###

if 'target_list' in config and config["target_list"] != 'NA':
  target_list_df = pd.read_table(config["target_list"], sep=r'\s+')
  if 'unrel' not in target_list_df.columns:
    target_list_df['unrel'] = 'NA'  # Adding a column with string 'NA' values
  target_list_df_23andMe = target_list_df.loc[target_list_df['type'] == '23andMe']
  samp_types = ['plink1', 'plink2', 'bgen', 'vcf']
  target_list_df_samp = target_list_df[target_list_df['type'].isin(samp_types)]
  target_list_df_indiv_report = target_list_df.loc[(target_list_df['indiv_report'].isin(['T', 'TRUE', True]))]
else:
  target_list_df = pd.DataFrame(columns = ["name", "path" "type", "indiv_report","unrel"])
  target_list_df_23andMe = pd.DataFrame(columns = ["name", "path" "type", "indiv_report","unrel"])
  target_list_df_samp = pd.DataFrame(columns = ["name", "path" "type", "indiv_report","unrel"])
  target_list_df_indiv_report = pd.DataFrame(columns = ["name", "path" "type", "indiv_report","unrel"])

# Check for duplicate values in the 'name' column
check_for_duplicates(target_list_df, 'name', 'target_list')

# Check specific target paths exist
check_target_type(df = target_list_df)

# Check specific target paths exist
check_target_paths(df = target_list_df, chr = str(get_chr_range(config['testing'])[0]))

###
# gwas_list
###

# Read in the gwas_list or make an empty version
if 'gwas_list' in config and config["gwas_list"] != 'NA':
  gwas_list_df = pd.read_table(config["gwas_list"], sep=r'\s+')
else:
  gwas_list_df = pd.DataFrame(columns = ["name", "path", "population", "n", "sampling", "prevalence", "mean", "sd", "label"])

# Remove commas in the 'n' column and convert to numeric
gwas_list_df['n'] = gwas_list_df['n'].replace({',': ''}, regex=True)

# Check for duplicate values in the 'name' column
check_for_duplicates(gwas_list_df, 'name', 'gwas_list')

# Check whether gwas_list paths exist
check_list_paths(gwas_list_df)

# Identify gwas_list with population == 'EUR'
gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

###
# score_list
###

# Read in score_list or create empty score_list
if 'score_list' in config and config["score_list"] != 'NA':
  score_list_df = pd.read_table(config["score_list"], sep=r'\s+')
  pgs_methods = config['pgs_methods']
  pgs_methods_all = list(config['pgs_methods'])
  pgs_methods_all.append('external')

  # Check whether score_list paths exist
  check_list_paths(score_list_df)
else:
  score_list_df = pd.DataFrame(columns = ["name", "path", "label"])
  pgs_methods = config['pgs_methods']
  pgs_methods_all = config['pgs_methods']

# Check for duplicate values in the 'name' column
check_for_duplicates(score_list_df, 'name', 'score_list')

# Check whether score_list paths exist
check_list_paths(score_list_df)

###
# gwas_groups
###

# Read in the gwas_groups or make an empty version
if 'gwas_groups' in config and config["gwas_groups"] != 'NA':
  gwas_groups_df = pd.read_table(config["gwas_groups"], sep=r'\s+')
else:
  gwas_groups_df = pd.DataFrame(columns = ["name", "gwas", "label"])

# Check for duplicate values in the 'name' column
check_for_duplicates(gwas_groups_df, 'name', 'gwas_groups')

# Function to get the list of GWAS names for a given group
def get_gwas_names(gwas_group):
    gwas_names_str = gwas_groups_df[gwas_groups_df['name'] == gwas_group]['gwas'].iloc[0]
    return gwas_names_str.split(',')

# Function to generate comma-separated list of populations for each name
def get_populations(gwas_group):
    gwas_names = get_gwas_names(gwas_group)
    sumstats_populations = []
    for gwas in gwas_names:
        gwas_info = gwas_list_df[gwas_list_df['name'] == gwas].iloc[0]
        sumstats_populations.append(gwas_info['population'])
    return sumstats_populations

# Check whether gwas_groups contains gwas that are not in the gwas_list
gwas_groups_gwas = gwas_groups_df['gwas'].str.split(',', expand=True).stack().unique()
gwas_list_names = gwas_list_df['name'].unique()
missing_gwas = set(gwas_groups_gwas) - set(gwas_list_names)
if missing_gwas:
    raise ValueError(f"The following GWAS are in gwas_groups but missing in gwas_list: {', '.join(missing_gwas)}")

# Subset gwas_groups to those with 2 GWAS specified
gwas_groups_df_two = gwas_groups_df[gwas_groups_df['gwas'].str.count(',') == 1]

###
# Check there are no duplicate values in name columns of gwas_list, score_list, gwas_groups
###

def check_for_duplicates_across_lists(df_list, name_col, list_names):
    combined_names = pd.concat([df[name_col] for df in df_list])

    # Find duplicates across all lists
    duplicate_names = combined_names[combined_names.duplicated(keep=False)]

    if not duplicate_names.empty:
        raise ValueError(f"Duplicate values found across {', '.join(list_names)}: {', '.join(duplicate_names.unique())}")

check_for_duplicates_across_lists(
    df_list=[gwas_list_df, score_list_df, gwas_groups_df],
    name_col='name',
    list_names=['gwas_list', 'score_list', 'gwas_groups']
)

#########
# Set PRS-CS ld reference path
if config['prscs_ldref'] == 'ukb':
    prscs_ldref='ukbb'
elif config['prscs_ldref'] == '1kg':
    prscs_ldref='1kg'

# Set resdir parameter
# If resdir is NA, set resdir to 'resources'
if config['resdir'] == 'NA':
  resdir='resources'
else:
  resdir=config['resdir']

# Set ldpred2 reference path
if config['ldpred2_ldref'] == 'NA':
  ldpred2_ldref=f"{resdir}/data/ldpred2_ref"
else:
  ldpred2_ldref=config['ldpred2_ldref']

# Check the ldpred2 ldref data is present for the required populations in the pgwas_list
if 'ldpred2' in config['pgs_methods']:
  for pop in gwas_list_df['population'].unique():
    path = f"{ldpred2_ldref}/{pop}"
    # Check if map.rds file exists
    map_file = os.path.join(path, "map.rds")
    if not os.path.exists(map_file):
      print(f"File not found: {map_file}")
      raise FileNotFoundError(f"Required file not found: {map_file}. LDpred2 reference data must include map.rds for all populations.")

    # Check if LD_with_blocks_chr${chr}.rds files exist for chr 1 to 22
    for chr in range(1, 23):
      ld_file = os.path.join(path, f"LD_with_blocks_chr{chr}.rds")
      if not os.path.exists(ld_file):
        print(f"File not found: {ld_file}")
        raise FileNotFoundError(f"Required file not found: {ld_file}. LDpred2 reference data must include files for all chromosomes.")

# Set sbayesr reference path
if config['sbayesr_ldref'] == 'NA':
  sbayesr_ldref=f"{resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr"
else:
  sbayesr_ldref=config['sbayesr_ldref']

# Check the sbayesr ldref data is present for the required populations in the gwas_list
if 'sbayesr' in config['pgs_methods']:
  for pop in gwas_list_df['population'].unique():
    path = f"{sbayesr_ldref}/{pop}"
    # Check if map.rds file exists
    map_file = os.path.join(path, "map.rds")
    if not os.path.exists(map_file):
      print(f"File not found: {map_file}")
      raise FileNotFoundError(f"Required file not found: {map_file}. SBayesR reference data must include map.rds for all populations.")

    # Check if LD_with_blocks_chr${chr}.rds files exist for chr 1 to 22
    for chr in range(1, 23):
      ld_file = os.path.join(path, f"LD_with_blocks_chr{chr}.rds")
      if not os.path.exists(ld_file):
        print(f"File not found: {ld_file}")
        raise FileNotFoundError(f"Required file not found: {ld_file}. SBayesR reference data must include files for all chromosomes.")

# Set quickprs reference path
if (config["leopard_methods"] and config["leopard_methods"] != "NA") or "quickprs" in config["pgs_methods"]:
  if config['quickprs_ldref'] == 'NA':
    quickprs_ldref=f"{resdir}/data/quickprs"
    
    # Check if gwas_list contains invalid populations
    valid_pops = {'EUR', 'EAS', 'AFR', 'SAS'}
    invalid_pops = set(gwas_list_df['population'].unique()) - valid_pops
  
    if invalid_pops:
      raise ValueError(
        f"Default quickprs reference data is only available for EUR, EAS, AFR and SAS populations. For other populations, please provide your own quickprs reference data using the quickprs_ldref parameter."
      )
  else:
    quickprs_ldref=config['quickprs_ldref']
  
    # Check the quickprs ldref data is present for the required populations in the gwas_list
    for pop in gwas_list_df['population'].unique():
      path = f"{quickprs_ldref}/{pop}"
      # Check if required files exists
      cors_file = os.path.join(path, f"{pop}.cors.bin")
      if not os.path.exists(cors_file):
        print(f"File not found: {cors_file}")
        raise FileNotFoundError(f"Required file not found: {cors_file}. quickprs reference data must include .cors.bin for all populations when quickprs_ldref is specified.")

# Set quickprs_multi reference path
quickprs_multi_ldref=config['quickprs_multi_ldref']
if config["leopard_methods"] and config["leopard_methods"] != "NA":
  missing_files = []
  for pop in gwas_list_df['population'].unique():
    path = f"{quickprs_multi_ldref}/{pop}"
    # Check if required files exists
    if not os.path.exists(f"{path}/{pop}.subset_1.bed"):
      missing_files.append(f"{path}/{pop}.subset_1.bed")
    if not os.path.exists(f"{path}/{pop}.subset_2.bed"):
      missing_files.append(f"{path}/{pop}.subset_2.bed")
    if not os.path.exists(f"{path}/{pop}.subset_3.bed"):
      missing_files.append(f"{path}/{pop}.subset_3.bed")
    if missing_files:
      raise FileNotFoundError(f"The following quickprs_multi reference data are missing: {', '.join(missing_files)}")

# Set sbayesrc reference path
if config['sbayesrc_ldref'] == 'NA':
  sbayesrc_ldref=f"{resdir}/data/sbayesrc_ref"
else:
  sbayesrc_ldref=config['sbayesrc_ldref']

# Check the sbayesrc ldref data is present for the required populations in the gwas_list
if 'sbayesrc' in config['pgs_methods']:
  for pop in gwas_list_df['population'].unique():
    path = f"{sbayesrc_ldref}/{pop}"
    # Check if required files exists
    cors_file = os.path.join(path, f"ldm.info")
    if not os.path.exists(cors_file):
      print(f"File not found: {cors_file}")
      raise FileNotFoundError(f"Required file not found: {cors_file}. sbayesrc reference data must include ldm.info for all populations.")

####
# Check reference data
####
if config['refdir'] == 'NA':
    refdir = f"{resdir}/data/ref"
    ref_input=f"{refdir}/ref.pop.txt"
else:
    refdir = config['refdir']

    ref_input = [os.path.join(refdir, f"ref.chr{i}.{ext}") for i in get_chr_range(testing=config['testing']) for ext in ['pgen', 'pvar', 'psam', 'rds']]
    ref_input.append(os.path.join(refdir, 'ref.pop.txt'))
    
    # Read populations from ref.pop.txt
    populations = set()
    ref_pop_file = os.path.join(refdir, 'ref.pop.txt')
    if os.path.exists(ref_pop_file):
        with open(ref_pop_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) == 2:
                    populations.add(parts[1])
    
    # Check keep files for populations in ref.pop.txt
    keep_dir = os.path.join(refdir, "keep_files")
    for pop in populations:
        keep_file = os.path.join(keep_dir, f"{pop}.keep")
        ref_input.append(keep_file)
    
    # Check frequency files for populations in ref.pop.txt and TRANS
    freq_dir = os.path.join(refdir, "freq_files")
    for pop in list(populations) + ['TRANS']:
        for i in range(1, 23):  # Chromosomes 1-22
            freq_file = os.path.join(freq_dir, pop, f"ref.{pop}.chr{i}.afreq")
            ref_input.append(freq_file)
    
    # Verify that all required files exist
    for full_path in ref_input:
        if not os.path.exists(full_path):
            raise FileNotFoundError(f"File not found: {full_path}. Check reference data format.")

#####

# Check valid pgs_methods are specified
def check_pgs_methods(x):
    # If pgs_methods is NA (None) or an empty list, return early without error
    if x is None or x == "NA" or not x:
        return

    valid_pgs_methods = {
        "ptclump", "dbslmm", "prscs", "sbayesr","sbayesrc", "lassosum", "ldpred2", "megaprs", "quickprs", "xwing", "prscsx"
    }

    invalid_methods = [method for method in x if method not in valid_pgs_methods]

    if invalid_methods:
        raise ValueError(f"Invalid pgs_methods specified: {', '.join(invalid_methods)}. "
                         f"Valid methods are: {', '.join(valid_pgs_methods)}.")

check_pgs_methods(config['pgs_methods'])

# Check valid tlprs_methods are specified
def check_tlprs_methods(config):
    valid_tlprs_methods = {
        "ptclump", "dbslmm", "prscs", "sbayesrc", "lassosum", "ldpred2", "megaprs", "quickprs"
    }

    # Check if 'tlprs_methods' is empty
    if config["tlprs_methods"] and config["tlprs_methods"] != "NA":
        # Check for invalid methods
        invalid_methods = [method for method in config["tlprs_methods"] if method not in valid_tlprs_methods]

        if invalid_methods:
            raise ValueError(f"Invalid tlprs_methods specified: {', '.join(invalid_methods)}. "
                             f"Valid methods are: {', '.join(valid_tlprs_methods)}.")

check_tlprs_methods(config)

# Check valid leopard_methods are specified
def check_leopard_methods(config):
    valid_leopard_methods = {
        "ptclump", "dbslmm", "prscs", "sbayesrc", "lassosum", "ldpred2", "megaprs","quickprs"
    }

    # Check if 'leopard_methods' is empty
    if config["leopard_methods"] and config["leopard_methods"] != "NA":
        # Check for invalid methods
        invalid_methods = [method for method in config["leopard_methods"] if method not in valid_leopard_methods]

        if invalid_methods:
            raise ValueError(f"Invalid leopard_methods specified: {', '.join(invalid_methods)}. "
                             f"Valid methods are: {', '.join(valid_leopard_methods)}.")

check_leopard_methods(config)

########
# Check for repo version updates
########

# If there has been a change to the major or minor version numbers, we will rerun the entire pipeline

# Define the path for storing the last known version
os.makedirs(resdir, exist_ok=True)
last_version_file = f"{resdir}/last_version.txt"

def get_current_version():
    cmd = "git describe --tags"
    tag = subprocess.check_output(cmd, shell=True).decode().strip()
    match = re.match(r"v?(\d+)\.(\d+)", tag)
    if match:
        return int(match.group(1)), int(match.group(2))  # Major, Minor
    else:
        raise ValueError("Git tag does not contain a valid version format.")

def read_last_version():
    if os.path.exists(last_version_file):
        with open(last_version_file, "r") as file:
            major, minor = file.read().strip().split('.')
            return int(major), int(minor)
    return 0, 0  # Default to 0.0 if file does not exist

def write_last_version(major, minor):
    with open(last_version_file, "w") as file:
        file.write(f"{major}.{minor}")

# Access overwrite flag from config
overwrite = config.get("overwrite", "false").lower() == "true"

# Main logic to check version and decide on execution
current_major, current_minor = get_current_version()
last_major, last_minor = read_last_version()

# Check if the last version is 0.0, proceed without requiring overwrite
if last_major == 0 and last_minor == 0:
    print(f"Initial version setup detected. Updating to v{current_major}.{current_minor}.")
    write_last_version(current_major, current_minor)
else:
    # Check for both major and minor version changes
    if current_major != last_major or current_minor != last_minor:
        if not overwrite:
            print(f"Change in version of GenoPred detected from v{last_major}.{last_minor} to v{current_major}.{current_minor}. Use --config overwrite=true to proceed.")
            sys.exit(1)
        else:
            print("Proceeding with version update due to overwrite=true config.")
            write_last_version(current_major, current_minor)  # Update the stored version

########
# Download dependencies
########

# Download impute2_data
rule download_impute2_data:
  output:
    directory(f"{resdir}/data/impute2/1000GP_Phase3")
  benchmark:
    f"{resdir}/data/benchmarks/download_impute2_data.txt"
  log:
    f"{resdir}/data/logs/download_impute2_data.log"
  shell:
    """
    {{
      rm -r -f {resdir}/data/impute2; \
      mkdir -p {resdir}/data/impute2/; \
      wget --no-check-certificate -O {resdir}/data/impute2/1000GP_Phase3.tgz https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz; \
      tar -zxvf {resdir}/data/impute2/1000GP_Phase3.tgz -C {resdir}/data/impute2/; \
      rm {resdir}/data/impute2/1000GP_Phase3.tgz; \
      wget --no-check-certificate -O {resdir}/data/impute2/1000GP_Phase3_chrX.tgz https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz; \
      tar -zxvf {resdir}/data/impute2/1000GP_Phase3_chrX.tgz -C {resdir}/data/impute2/1000GP_Phase3/; \
      rm {resdir}/data/impute2/1000GP_Phase3_chrX.tgz
    }} > {log} 2>&1
    """

# Download PLINK. DBSLMM requires the binary to be specified, which is challenging with conda environments. I have tried to avoid this again but no joy. The conda environment may not exist when the snakemake is executed which will cause problems if trying to access the conda environment manually.
rule download_plink:
  output:
    f"{resdir}/software/plink/plink"
  benchmark:
    f"{resdir}/data/benchmarks/download_plink.txt"
  log:
    f"{resdir}/data/logs/download_plink.log"
  shell:
    """
    {{
      rm -r -f {resdir}/software/plink; \
      mkdir -p {resdir}/software/plink; \
      wget --no-check-certificate -O {resdir}/software/plink/plink_linux_x86_64_20210606.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip; \
      unzip {resdir}/software/plink/plink_linux_x86_64_20210606.zip -d {resdir}/software/plink; \
      rm {resdir}/software/plink/plink_linux_x86_64_20210606.zip
    }} > {log} 2>&1
    """

# Download LDSC
rule download_ldsc:
  output:
    f"{resdir}/software/ldsc/ldsc.py"
  benchmark:
    f"{resdir}/data/benchmarks/download_ldsc.txt"
  log:
    f"{resdir}/data/logs/download_ldsc.log"
  shell:
    """
    {{
      rm -r -f {resdir}/software/ldsc; \
      git clone https://github.com/bulik/ldsc.git {resdir}/software/ldsc/; \
      cd {resdir}/software/ldsc/; \
      git reset --hard aa33296abac9569a6422ee6ba7eb4b902422cc74
    }} > {log} 2>&1
    """

# Download ld scores from PanUKB
rule download_ldscores_panukb:
  output:
    directory(f"{resdir}/data/ld_scores")
  benchmark:
    f"{resdir}/data/benchmarks/download_ldscores_panukb.txt"
  log:
    f"{resdir}/data/logs/download_ldscores_panukb.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data/ld_scores; \
      wget --no-check-certificate -O {resdir}/data/ld_scores.tar.gz https://zenodo.org/records/10666891/files/ld_scores.tar.gz?download=1; \
      tar -xvf {resdir}/data/ld_scores.tar.gz -C {resdir}/data/; \
      rm {resdir}/data/ld_scores.tar.gz
    }} > {log} 2>&1
    """

# Download hapmap3 snplist
rule download_hm3_snplist:
  output:
    f"{resdir}/data/hm3_snplist/w_hm3.snplist"
  benchmark:
    f"{resdir}/data/benchmarks/download_hm3_snplist.txt"
  log:
    f"{resdir}/data/logs/download_hm3_snplist.log"
  shell:
    """
    {{
      rm -r -f {resdir}/data/hm3_snplist; \
      mkdir -p {resdir}/data/hm3_snplist; \
      wget --no-check-certificate -O {resdir}/data/hm3_snplist/w_hm3.snplist.gz https://zenodo.org/record/7773502/files/w_hm3.snplist.gz?download=1; \
      gunzip {resdir}/data/hm3_snplist/w_hm3.snplist.gz
    }} > {log} 2>&1
    """

# Download DBSLMM
rule download_dbslmm:
  output:
    directory(f"{resdir}/software/dbslmm/")
  benchmark:
    f"{resdir}/data/benchmarks/download_dbslmm.txt"
  conda:
    "../envs/analysis.yaml"
  log:
    f"{resdir}/data/logs/download_dbslmm.log"
  shell:
    """
    {{
      git clone https://github.com/biostat0903/DBSLMM.git {output}; \
      cd {output}; \
      git reset --hard d158a144dd2f2dc883ad93d0ea71e8fc48e80bd3; \
      wget --no-check-certificate -O software/dbslmm "https://drive.usercontent.google.com/download?id=1eAbEyhF8rO_faOFL3jqRo9LmfgJNRH6K&export=download&authuser=0"; \
      chmod a+x software/dbslmm
    }} > {log} 2>&1
    """

# Download LD block data
rule download_ld_blocks:
  output:
    directory(f"{resdir}/data/ld_blocks/")
  benchmark:
    f"{resdir}/data/benchmarks/download_ld_blocks.txt"
  log:
    f"{resdir}/data/logs/download_ld_blocks.log"
  shell:
    """
    {{
      git clone https://bitbucket.org/nygcresearch/ldetect-data.git {output}; \
      mv {resdir}/data/ld_blocks/ASN {resdir}/data/ld_blocks/EAS
    }} > {log} 2>&1
    """

# Download UKB-based PRScs reference
# Note. Using Yengo Height GWAS with OpenSNP target, using the UKB reference performed significantly worse than the 1KG reference.
prscs_ref_ukb_dropbox = {
    'eur': 'https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=0',
    'eas': 'https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=0',
    'afr': 'https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=0',
    'amr': 'https://www.dropbox.com/s/y7ruj364buprkl6/ldblk_ukbb_amr.tar.gz?dl=0',
    'sas': 'https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz?dl=0',
}

rule download_prscs_ref_ukb:
  output:
    f"{resdir}/data/prscs_ref/ukbb/ldblk_ukbb_{{population}}/ldblk_ukbb_chr1.hdf5"
  benchmark:
    f"{resdir}/data/benchmarks/download_prscs_ref_ukb-{{population}}.txt"
  log:
    f"{resdir}/data/logs/download_prscs_ref_ukb-{{population}}.log"
  params:
    url=lambda w: prscs_ref_ukb_dropbox.get(w.population)
  shell:
    """
    {{
      mkdir -p {resdir}/data/prscs_ref/ukbb; \
      rm -r -f {resdir}/data/prscs_ref/ukbb/ldblk_ukbb_{wildcards.population}; \
      wget --no-check-certificate -O {resdir}/data/prscs_ref/ukbb/ldblk_ukbb_{wildcards.population}.tar.gz {params.url}; \
      tar -zxvf {resdir}/data/prscs_ref/ukbb/ldblk_ukbb_{wildcards.population}.tar.gz -C {resdir}/data/prscs_ref/ukbb/; \
      rm {resdir}/data/prscs_ref/ukbb/ldblk_ukbb_{wildcards.population}.tar.gz
    }} > {log} 2>&1
    """

rule download_prscs_ref_ukb_all:
  input:
    lambda w: expand(f"{resdir}/data/prscs_ref/ukbb/ldblk_ukbb_{{population}}/ldblk_ukbb_chr1.hdf5", population=['eur','eas','afr','amr','sas'])

# Download 1KG-based PRScs reference
prscs_ref_1kg_dropbox = {
    'eur': 'https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0',
    'eas': 'https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=0',
    'afr': 'https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0',
    'amr': 'https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz?dl=0',
    'sas': 'https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz?dl=0',
}

rule download_prscs_ref_1kg:
  output:
    f"{resdir}/data/prscs_ref/1kg/ldblk_1kg_{{population}}/ldblk_1kg_chr1.hdf5"
  benchmark:
    f"{resdir}/data/benchmarks/download_prscs_ref_1kg-{{population}}.txt"
  log:
    f"{resdir}/data/logs/download_prscs_ref_1kg-{{population}}.log"
  params:
    url=lambda w: prscs_ref_1kg_dropbox.get(w.population)
  shell:
    """
    {{
      mkdir -p {resdir}/data/prscs_ref/1kg; \
      rm -r -f {resdir}/data/prscs_ref/1kg/ldblk_1kg_{wildcards.population}; \
      wget --no-check-certificate -O {resdir}/data/prscs_ref/1kg/ldblk_1kg_{wildcards.population}.tar.gz {params.url}; \
      tar -zxvf {resdir}/data/prscs_ref/1kg/ldblk_1kg_{wildcards.population}.tar.gz -C {resdir}/data/prscs_ref/1kg/; \
      rm {resdir}/data/prscs_ref/1kg/ldblk_1kg_{wildcards.population}.tar.gz
    }} > {log} 2>&1
    """

rule download_prscs_ref_1kg_all:
  input:
    lambda w: expand(f"{resdir}/data/prscs_ref/1kg/ldblk_1kg_{{population}}/ldblk_1kg_chr1.hdf5", population=['eur','eas','afr','amr','sas'])

# Download PRScs software
rule download_prscs_software:
  output:
    directory(f"{resdir}/software/prscs/")
  benchmark:
    f"{resdir}/data/benchmarks/download_prscs_software.txt"
  log:
    f"{resdir}/data/logs/download_prscs_software.log"
  shell:
    """
    {{
      git clone https://github.com/getian107/PRScs.git {output}; \
      cd {output}; \
      git reset --hard 621fdc80daac56c93d9528eb1a1187f7b1fc9afb
    }} > {log} 2>&1
    """

# Download PRS-CSx software
rule download_prscsx_software:
  output:
    directory(f"{resdir}/software/prscsx/")
  benchmark:
    f"{resdir}/data/benchmarks/download_prscsx_software.txt"
  log:
    f"{resdir}/data/logs/download_prscsx_software.log"
  shell:
    """
    {{
      git clone https://github.com/getian107/PRScsx.git {output}; \
      cd {output}; \
      git reset --hard 29a1148875f6ae3f2594b25579f40d4b587c5691
    }} > {log} 2>&1
    """

####
# Download PRS-CSx SNP data for reference
####

rule download_prscs_snp_data_ukb:
  output:
    f"{resdir}/data/prscs_ref/ukbb/snpinfo_mult_ukbb_hm3"
  benchmark:
    f"{resdir}/data/benchmarks/download_prscs_snp_data_ukb.txt"
  log:
    f"{resdir}/data/logs/download_prscs_snp_data_ukb.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data/prscs_ref/ukbb; \
      rm -r -f {resdir}/data/prscs_ref/ukbb/snpinfo_mult_ukbb_hm3; \
      wget --no-check-certificate -O {resdir}/data/prscs_ref/ukbb/snpinfo_mult_ukbb_hm3 https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3?dl=0; \
    }} > {log} 2>&1
    """

rule download_prscs_snp_data_1kg:
  output:
    f"{resdir}/data/prscs_ref/1kg/snpinfo_mult_1kg_hm3"
  benchmark:
    f"{resdir}/data/benchmarks/download_prscs_snp_data_ukb.txt"
  log:
    f"{resdir}/data/logs/download_prscs_snp_data_ukb.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data/prscs_ref/1kg; \
      rm -r -f {resdir}/data/prscs_ref/1kg/snpinfo_mult_1kg_hm3; \
      wget --no-check-certificate -O {resdir}/data/prscs_ref/1kg/snpinfo_mult_1kg_hm3 https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3?dl=0; \
    }} > {log} 2>&1
    """

# Download gctb reference
rule download_gctb_ref:
  output:
    directory(f"{resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse")
  benchmark:
    f"{resdir}/data/benchmarks/download_gctb_ref.txt"
  log:
    f"{resdir}/data/logs/download_gctb_ref.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data/gctb_ref; \
      wget --no-check-certificate -O {resdir}/data/gctb_ref/ukbEURu_hm3_sparse.zip https://zenodo.org/record/3350914/files/ukbEURu_hm3_sparse.zip?download=1; \
      unzip {resdir}/data/gctb_ref/ukbEURu_hm3_sparse.zip -d {resdir}/data/gctb_ref; \
      for chr in $(seq 1 22);do \
      mv {resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k.ldm.sparse.bin {resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr${{chr}}.ldm.sparse.bin; \
      mv {resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k.ldm.sparse.info {resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr${{chr}}.ldm.sparse.info; \
      mv {resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k_sparse.log {resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_sparse_chr${{chr}}.log; \
      done; \
      rm {resdir}/data/gctb_ref/ukbEURu_hm3_sparse.zip
    }} > {log} 2>&1
    """
# Download GCTB
rule download_gctb_software:
  output:
    f"{resdir}/software/gctb/gctb_2.03beta_Linux/gctb"
  benchmark:
    f"{resdir}/data/benchmarks/download_gctb_software.txt"
  log:
    f"{resdir}/data/logs/download_gctb_software.log"
  shell:
    """
    {{
      rm -r -f {resdir}/software/gctb; \
      mkdir -p {resdir}/software/gctb; \
      wget --no-check-certificate -O {resdir}/software/gctb/gctb_2.03beta_Linux.zip https://cnsgenomics.com/software/gctb/download/gctb_2.03beta_Linux.zip; \
      unzip {resdir}/software/gctb/gctb_2.03beta_Linux.zip -d {resdir}/software/gctb; \
      rm {resdir}/software/gctb/gctb_2.03beta_Linux.zip
    }} > {log} 2>&1
    """

# Download GCTB v2.5.2 for SBayesRC
rule download_gctb252_software:
  output:
    f"{resdir}/software/gctb_2.5.2/gctb_2.5.2_Linux/gctb"
  benchmark:
    f"{resdir}/data/benchmarks/download_gctb252_software.txt"
  log:
    f"{resdir}/data/logs/download_gctb252_software.log"
  shell:
    """
    {{
      rm -r -f {resdir}/software/gctb_2.5.2; \
      mkdir -p {resdir}/software/gctb_2.5.2; \
      wget --no-check-certificate -O {resdir}/software/gctb_2.5.2/gctb_2.5.2_Linux.zip https://cnsgenomics.com/software/gctb/download/gctb_2.5.2_Linux.zip; \
      unzip {resdir}/software/gctb_2.5.2/gctb_2.5.2_Linux.zip -d {resdir}/software/gctb_2.5.2; \
      rm {resdir}/software/gctb_2.5.2/gctb_2.5.2_Linux.zip
    }} > {log} 2>&1
    """

# Download annotations for SBayesRC
rule download_sbayesrc_annot:
  output:
    f"{resdir}/data/sbayesrc_annot/annot_baseline2.2.txt"
  benchmark:
    f"{resdir}/data/benchmarks/download_sbayesrc_annot.txt"
  log:
    f"{resdir}/data/logs/download_sbayesrc_annot.log"
  shell:
    """
    {{
      rm -r -f {resdir}/data/sbayesrc_annot; \
      mkdir -p {resdir}/data/sbayesrc_annot; \
      wget --no-check-certificate -O {resdir}/data/sbayesrc_annot/annot_baseline2.2.zip https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/annot_baseline2.2.zip; \
      unzip {resdir}/data/sbayesrc_annot/annot_baseline2.2.zip -d {resdir}/data/sbayesrc_annot; \
      rm {resdir}/data/sbayesrc_annot/annot_baseline2.2.zip
    }} > {log} 2>&1
    """

# Download SBayesRC R package
rule install_sbayesrc:
  input:
    "envs/sbayesrc.yaml"
  output:
    touch("resources/software/install_sbayesrc.done")
  conda:
    "../envs/sbayesrc.yaml"
  benchmark:
    "resources/data/benchmarks/install_sbayesrc.txt"
  log:
    "resources/data/logs/install_sbayesrc.log"
  shell:
    """
    {{
      Rscript -e 'install.packages(\"https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.6/SBayesRC_0.2.6.tar.gz\", repos=NULL, type=\"source\")'
    }} > {log} 2>&1
    """

# Install GenoUtils in SBayesRC environment
rule install_genoutils_sbayesrc:
  input:
    rules.install_sbayesrc.output
  output:
    touch("resources/software/install_genoutils_sbayesrc.done")
  conda:
    "../envs/sbayesrc.yaml"
  benchmark:
    "resources/data/benchmarks/install_genoutils_sbayesrc.txt"
  log:
    "resources/data/logs/install_genoutils_sbayesrc.log"
  shell:
    """
    {{
      Rscript -e 'devtools::install_github(\"opain/GenoUtils@6334159ab5d95ce936896e6938a1031c38ed4f30\")'
    }} > {log} 2>&1
    """

# Download LDpred2 reference
rule download_ldpred2_ref:
  output:
    f"{resdir}/data/ldpred2_ref/EUR/map.rds"
  benchmark:
    f"{resdir}/data/benchmarks/download_ldpred2_ref.txt"
  log:
    f"{resdir}/data/logs/download_ldpred2_ref.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data/ldpred2_ref/EUR; \
      wget --no-check-certificate -O {resdir}/data/ldpred2_ref/EUR/download.zip https://figshare.com/ndownloader/articles/19213299/versions/2; \
      unzip {resdir}/data/ldpred2_ref/EUR/download.zip -d {resdir}/data/ldpred2_ref/EUR/; \
      rm {resdir}/data/ldpred2_ref/EUR/download.zip; \
      unzip {resdir}/data/ldpred2_ref/EUR/ldref_with_blocks.zip -d {resdir}/data/ldpred2_ref/EUR/; \
      mv {resdir}/data/ldpred2_ref/EUR/ldref/* {resdir}/data/ldpred2_ref/EUR/; \
      rm {resdir}/data/ldpred2_ref/EUR/ldref_with_blocks.zip; \
      rm -r {resdir}/data/ldpred2_ref/EUR/ldref
    }} > {log} 2>&1
    """
# Download LDAK
rule download_ldak:
  output:
    f"{resdir}/software/ldak/ldak5.1.linux"
  benchmark:
    f"{resdir}/data/benchmarks/download_ldak.txt"
  log:
    f"{resdir}/data/logs/download_ldak.log"
  shell:
    """
    {{
      rm -r -f {resdir}/software/ldak; \
      mkdir -p {resdir}/software/ldak; \
      wget --no-check-certificate -O {resdir}/software/ldak/ldak5.1.linux_.zip https://dougspeed.com/wp-content/uploads/ldak5.1.linux_.zip; \
      unzip {resdir}/software/ldak/ldak5.1.linux_.zip -d {resdir}/software/ldak/; \
      chmod a+x {resdir}/software/ldak/ldak5.1.linux; \
      rm {resdir}/software/ldak/ldak5.1.linux_.zip
    }} > {log} 2>&1
    """
# Download LDAK map data
rule download_ldak_map:
  output:
    f"{resdir}/data/ldak_map/genetic_map_b37/genetic_map_chr22_combined_b37.txt"
  benchmark:
    f"{resdir}/data/benchmarks/download_ldak_map.txt"
  log:
    f"{resdir}/data/logs/download_ldak_map.log"
  shell:
    """
    {{
      rm -r {resdir}/data/ldak_map; \
      mkdir -p {resdir}/data/ldak_map; \
      wget --no-check-certificate -O {resdir}/data/ldak_map/genetic_map_b37.zip https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip; \
      unzip {resdir}/data/ldak_map/genetic_map_b37.zip -d {resdir}/data/ldak_map/; \
      rm {resdir}/data/ldak_map/genetic_map_b37.zip
    }} > {log} 2>&1
    """
# Download LDAK bld snp annotations
rule download_ldak_bld:
  output:
    f"{resdir}/data/ldak_bld/README.txt"
  benchmark:
    f"{resdir}/data/benchmarks/download_ldak_bld.txt"
  log:
    f"{resdir}/data/logs/download_ldak_bld.log"
  shell:
    """
    {{
      rm -r {resdir}/data/ldak_bld; \
      mkdir -p {resdir}/data/ldak_bld; \
      wget --no-check-certificate -O {resdir}/data/ldak_bld/bld.zip https://genetics.ghpc.au.dk/doug/bld.zip; \
      unzip {resdir}/data/ldak_bld/bld.zip -d {resdir}/data/ldak_bld/; \
      rm {resdir}/data/ldak_bld/bld.zip
    }} > {log} 2>&1
    """
# Download LDAK high ld regions file
rule download_ldak_highld:
  output:
    f"{resdir}/data/ldak_highld/highld.txt"
  benchmark:
    f"{resdir}/data/benchmarks/download_ldak_highld.txt"
  log:
    f"{resdir}/data/logs/download_ldak_highld.log"
  shell:
    """
    {{
      rm -r {resdir}/data/ldak_highld; \
      mkdir -p {resdir}/data/ldak_highld; \
      wget --no-check-certificate -O {resdir}/data/ldak_highld/highld.txt https://dougspeed.com/wp-content/uploads/highld.txt
    }} > {log} 2>&1
    """

# Download LDAK V5.2 for QuickPRS
# Only this version works for QuickPRS
rule download_ldak5_2:
  output:
    f"{resdir}/software/ldak5.2/ldak5.2.linux"
  benchmark:
    f"{resdir}/data/benchmarks/download_ldak5_2.txt"
  log:
    f"{resdir}/data/logs/download_ldak5_2.log"
  shell:
    """
    {{
      rm -r -f {resdir}/software/ldak5.2; \
      mkdir -p {resdir}/software/ldak5.2; \
      wget --no-check-certificate -O {resdir}/software/ldak5.2/ldak5.2.linux "https://drive.google.com/uc?export=download&id=19knXZnbPNDz3J5dBKeVyZZoe6iZnLPEk"; \
      chmod a+x {resdir}/software/ldak5.2/ldak5.2.linux
    }} > {log} 2>&1
    """

# NOTE. This doesn't currently work as the reference data on LDAK website isn't in the right format for LDAK 5.1, 5.2, or 6
# In due course we will update this to download from our own repo. For the timebeing. The user must specify the reference data themselves
rule download_quickprs_ref:
  output:
    f"{resdir}/data/quickprs/{{population}}/{{population}}.cors.bin"
  benchmark:
    f"{resdir}/data/benchmarks/download_quickprs_ref-{{population}}.txt"
  log:
    f"{resdir}/data/logs/download_quickprs_ref-{{population}}.log"
  params:
    pop_code=lambda wildcards: {'EUR': 'gbr', 'SAS': 'sas', 'EAS': 'eas', 'AFR': 'afr'}[wildcards.population]
  shell:
    """
    {{
      mkdir -p {resdir}/data/quickprs; \
      rm -r -f {resdir}/data/quickprs/{wildcards.population}; \
      wget --no-check-certificate -O {resdir}/data/quickprs/{wildcards.population}.hapmap.tar.gz https://genetics.ghpc.au.dk/doug/{params.pop_code}.hapmap.tar.gz; \
      tar -zxvf {resdir}/data/quickprs/{wildcards.population}.hapmap.tar.gz -C {resdir}/data/quickprs/; \
      mv {resdir}/data/quickprs/{params.pop_code}.hapmap {resdir}/data/quickprs/{wildcards.population}; \
      find {resdir}/data/quickprs/{wildcards.population} -type f -name '*{params.pop_code}*' -exec bash -c 'mv \"$0\" \"${{0//{params.pop_code}/{wildcards.population}}}\"' {{}} \; \
      find {resdir}/data/quickprs/{wildcards.population} -type f -name '*hapmap*' -exec bash -c 'mv \"$0\" \"${{0//.hapmap./.}}\"' {{}} \; \
      rm {resdir}/data/quickprs/{wildcards.population}.hapmap.tar.gz
    }} > {log} 2>&1
    """

rule download_quickprs_ref_all:
  input:
    lambda w: expand(f"{resdir}/data/quickprs/{{population}}/{{population}}.cors.bin", population=['EUR', 'SAS', 'EAS', 'AFR'])

# Download preprocessed reference data (1KG+HGDP HapMap3)
rule download_default_ref:
  output:
    f"{resdir}/data/ref/ref.pop.txt"
  benchmark:
    f"{resdir}/data/benchmarks/download_default_ref.txt"
  log:
    f"{resdir}/data/logs/download_default_ref.log"
  shell:
    """
    {{
      rm -r {resdir}/data/ref; \
      mkdir -p {resdir}/data/ref; \
      gdown --id 1vYH6V-7F68Ji1vy9TaH0ysjmdYJFef-f -O resources/data/ref/genopred_1kg_hgdp.tar.gz; \
      tar -xzvf {resdir}/data/ref/genopred_1kg_hgdp.tar.gz -C {resdir}/data/ref/; \
      mv {resdir}/data/ref/ref/* {resdir}/data/ref/; \
      rm -r {resdir}/data/ref/ref; \
      rm {resdir}/data/ref/genopred_1kg_hgdp.tar.gz
    }} > {log} 2>&1
    """

# install ggchicklet
rule install_ggchicklet:
  input:
    "envs/analysis.yaml"
  output:
    touch("resources/software/install_ggchicklet.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    "resources/data/benchmarks/install_ggchicklet.txt"
  log:
    "resources/data/logs/install_ggchicklet.log"
  shell:
    """
    {{
      Rscript -e 'remotes::install_github(\"hrbrmstr/ggchicklet@64c468dd0900153be1690dbfc5cfb35710da8183\")'
    }} > {log} 2>&1
    """

# install lassosum
rule install_lassosum:
  input:
    "envs/analysis.yaml"
  output:
    touch("resources/software/install_lassosum.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    "resources/data/benchmarks/install_lassosum.txt"
  log:
    "resources/data/logs/install_lassosum.log"
  shell:
    """
    {{
      Rscript -e 'remotes::install_github(\"tshmak/lassosum@v0.4.5\")'
    }} > {log} 2>&1
    """

# Install GenoUtils
rule install_genoutils:
  input:
    "envs/analysis.yaml"
  output:
    touch("resources/software/install_genoutils.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    "resources/data/benchmarks/install_genoutils.txt"
  log:
    "resources/data/logs/install_genoutils.log"
  shell:
    """
    {{
      Rscript -e 'devtools::install_github(\"opain/GenoUtils@6334159ab5d95ce936896e6938a1031c38ed4f30\")'
    }} > {log} 2>&1
    """

# Install R packages (handy function for when conda env updates erroneously)
rule install_r_packages:
  input:
    rules.install_ggchicklet.output,
    rules.install_lassosum.output,
    rules.install_genoutils.output

# Download pgscatalog_utils
rule download_pgscatalog_utils:
  output:
    "resources/software/pgscatalog_utils/download_pgscatalog_utils.done"
  conda:
    "../envs/pgscatalog_utils.yaml"
  benchmark:
    "resources/data/benchmarks/download_pgscatalog_utils.txt"
  log:
    "resources/data/logs/download_pgscatalog_utils.log"
  shell:
    """
    # Create the log file and mark the rule as done.
    echo "Environment setup triggered" > {log}
    touch {output}
    """

# Download XPASS for X-wing dependencies
rule install_xpass:
  input:
    "envs/xwing.yaml"
  output:
    touch("resources/software/install_xpass.done")
  conda:
    "../envs/xwing.yaml"
  benchmark:
    "resources/data/benchmarks/install_xpass.txt"
  log:
    "resources/data/logs/install_xpass.log"
  shell:
    """
    {{
      Rscript -e 'devtools::install_github(\"YangLabHKUST/XPASS@65877ffba60dce69e0a6aa31c2e61045bf36dc40\")'
    }} > {log} 2>&1
    """

# Install GenoUtils in X-wing environment
rule install_genoutils_xwing:
  input:
    rules.install_xpass.output
  output:
    touch("resources/software/install_genoutils_xwing.done")
  conda:
    "../envs/xwing.yaml"
  benchmark:
    "resources/data/benchmarks/install_genoutils_xwing.txt"
  log:
    "resources/data/logs/install_genoutils_xwing.log"
  shell:
    """
    {{
      Rscript -e 'devtools::install_github(\"opain/GenoUtils@6334159ab5d95ce936896e6938a1031c38ed4f30\")'
    }} > {log} 2>&1
    """


# Download X-wing repo
rule download_xwing_software:
  input:
    rules.install_xpass.output,
    rules.install_genoutils_xwing.output
  output:
    "resources/software/xwing/block_partition.txt"
  conda:
    "../envs/xwing.yaml"
  benchmark:
    "resources/data/benchmarks/download_xwing_software.txt"
  log:
    "resources/data/logs/download_xwing_software.log"
  shell:
    """
    {{
      rm -r -f resources/software/xwing; \
      git clone https://github.com/opain/X-Wing resources/software/xwing; \
      cd resources/software/xwing; \
      git reset --hard e9fcc264266e0e884323311816bfe20053fd3f7a
    }} > {log} 2>&1
    """

# Download LOGODetect (X-wing) reference data
rule download_logodetect_ref:
  output:
    f"{resdir}/data/LOGODetect_1kg_ref/{{population}}/1000G_{{population}}_QC.bim"
  benchmark:
    f"{resdir}/data/benchmarks/logodetect_ref-{{population}}.txt"
  log:
    f"{resdir}/data/logs/logodetect_ref-{{population}}.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data; \
      wget --no-check-certificate -O {resdir}/data/LOGODetect_1kg_{wildcards.population}.tar.gz ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_{wildcards.population}.tar.gz; \
      tar -zxvf {resdir}/data/LOGODetect_1kg_{wildcards.population}.tar.gz -C {resdir}/data/; \
      rm {resdir}/data/LOGODetect_1kg_{wildcards.population}.tar.gz
    }} > {log} 2>&1
    """

rule download_logodetect_ref_all:
  input:
    lambda w: expand(f"{resdir}/data/LOGODetect_1kg_ref/{{population}}/1000G_{{population}}_QC.bim", population=['EUR','EAS','AFR','SAS','AMR'])

# Download PANTHER (X-wing) reference data
# The reference data is the same as the PRS-CS reference data
# The PRS-CS ref parameter will also affect the X-WING/PANTHER analysis

# Download LEOPARD (X-wing) reference data
rule download_leopard_ref:
  output:
    f"{resdir}/data/LEOPARD_1kg_ref/{{population}}/{{population}}_part1.bed"
  benchmark:
    f"{resdir}/data/benchmarks/download_leopard_ref-{{population}}.txt"
  log:
    f"{resdir}/data/logs/download_leopard_ref-{{population}}.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data; \
      wget --no-check-certificate -O {resdir}/data/LEOPARD_1kg_hm3_{wildcards.population}.tar.gz ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LEOPARD/LEOPARD_1kg_hm3_{wildcards.population}.tar.gz; \
      tar -zxvf {resdir}/data/LEOPARD_1kg_hm3_{wildcards.population}.tar.gz -C {resdir}/data/; \
      rm {resdir}/data/LEOPARD_1kg_hm3_{wildcards.population}.tar.gz
    }} > {log} 2>&1
    """

rule download_leopard_ref_all:
  input:
    lambda w: expand(f"{resdir}/data/LEOPARD_1kg_ref/{{population}}/{{population}}_part1.bed", population=['EUR','EAS','AFR','SAS','AMR'])

# Download LEOPARD and subsampled PANTHER (X-wing) reference data
rule download_leopard_panther_ref:
  output:
    f"{resdir}/data/PANTHER_LEOPARD_1kg_ref/ldblk_1kg_{{population}}/ldblk_1kg_chr13.hdf5"
  benchmark:
    f"{resdir}/data/benchmarks/download_leopard_panther_ref-{{population}}.txt"
  log:
    f"{resdir}/data/logs/download_leopard_panther_ref-{{population}}.log"
  params:
    pop_upper=lambda w: w.population.upper()
  shell:
    """
    {{
      mkdir -p {resdir}/data; \
      wget --no-check-certificate -O {resdir}/data/PANTHER_LEOPARD_1kg_{wildcards.population}.tar.gz ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LEOPARD/PANTHER_LEOPARD_1kg_{params.pop_upper}.tar.gz; \
      tar -zxvf {resdir}/data/PANTHER_LEOPARD_1kg_{wildcards.population}.tar.gz -C {resdir}/data/; \
      rm {resdir}/data/PANTHER_LEOPARD_1kg_{wildcards.population}.tar.gz
    }} > {log} 2>&1
    """

rule download_leopard_panther_ref_all:
  input:
    lambda w: expand(f"{resdir}/data/PANTHER_LEOPARD_1kg_ref/ldblk_1kg_{{population}}/ldblk_1kg_chr13.hdf5", population=['eur','eas','afr','sas','amr'])

rule download_leopard_panther_snp_data:
  output:
    f"{resdir}/data/PANTHER_LEOPARD_1kg_ref/snpinfo_mult_1kg_hm3"
  benchmark:
    f"{resdir}/data/benchmarks/download_leopard_panther_snp_data.txt"
  log:
    f"{resdir}/data/logs/download_leopard_panther_snp_data.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data; \
      wget --no-check-certificate -O {resdir}/data/snpinfo_mult_1kg_hm3_PANTHER_LEOPARD.tar.gz ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LEOPARD/snpinfo_mult_1kg_hm3_PANTHER_LEOPARD.tar.gz; \
      tar -zxvf {resdir}/data/snpinfo_mult_1kg_hm3_PANTHER_LEOPARD.tar.gz -C {resdir}/data/; \
      rm {resdir}/data/snpinfo_mult_1kg_hm3_PANTHER_LEOPARD.tar.gz
    }} > {log} 2>&1
    """

############

# Install TL-PRS
rule install_tlprs:
  output:
    touch("resources/software/install_tlprs.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{resdir}/data/benchmarks/install_tlprs.txt"
  log:
    f"{resdir}/data/logs/install_tlprs.log"
  shell:
    """
    {{
      Rscript -e 'devtools::install_github(\"opain/TLPRS@5a5528a3f709ca7d627381a3f09ccdcb923b50f4\")'
    }} > {log} 2>&1
    """

############

# Install GenoUtils in BridgePRS environment
rule install_genoutils_bridgeprs:
  output:
    touch("resources/software/install_genoutils_bridgeprs.done")
  conda:
    "../envs/bridgeprs.yaml"
  benchmark:
    f"{resdir}/data/benchmarks/install_genoutils_bridgeprs.txt"
  log:
    f"{resdir}/data/logs/install_genoutils_bridgeprs.log"
  shell:
    """
    {{
      Rscript -e 'devtools::install_github(\"opain/GenoUtils@6334159ab5d95ce936896e6938a1031c38ed4f30\")'
    }} > {log} 2>&1
    """

# Download BridgePRS
rule download_bridgeprs_software:
  input:
    rules.install_genoutils_bridgeprs.output
  output:
    f"{resdir}/software/bridgeprs/bridgePRS"
  benchmark:
    f"{resdir}/data/benchmarks/download_bridgeprs_software.txt"
  log:
    f"{resdir}/data/logs/download_bridgeprs_software.log"
  shell:
    """
    {{
      rm -r -f {resdir}/software/bridgeprs; \
      git clone https://github.com/clivehoggart/BridgePRS.git {resdir}/software/bridgeprs; \
      cd {resdir}/software/bridgeprs; \
      git reset --hard 1e1c9477d42d44ed312759751adbc25d146aeb9f
    }} > {log} 2>&1
    """

############
# Check all dependencies are available
############

# Rule to install and download commonly required software and data dependencies
rule get_dependencies:
  input:
    rules.download_plink.output,
    rules.download_ldscores_panukb.output,
    rules.download_ldsc.output,
    rules.download_hm3_snplist.output,
    rules.download_default_ref.output,
    rules.download_dbslmm.output,
    rules.download_ld_blocks.output,
    rules.install_ggchicklet.output,
    rules.install_lassosum.output,
    rules.install_genoutils.output,
    rules.download_pgscatalog_utils.output
  output:
    touch(f"{resdir}/software/get_dependencies.done")

############
# Rules for preparing offline resources
############

rule get_all_resources:
  input:
    rules.download_plink.output,
    rules.download_ldsc.output,
    rules.download_dbslmm.output,
    rules.download_prscs_software.output,
    rules.download_gctb_software.output,
    rules.download_ldak.output,
    rules.download_ldscores_panukb.output,
    rules.download_hm3_snplist.output,
    rules.download_ld_blocks.output,
    rules.download_prscs_ref_ukb_all.input,
    rules.download_prscs_ref_1kg_all.input,
    rules.download_gctb_ref.output,
    rules.download_ldpred2_ref.output,
    rules.download_ldak_map.output,
    rules.download_ldak_bld.output,
    rules.download_ldak_highld.output,
    rules.download_default_ref.output,
    rules.download_quickprs_ref_all.output
  output:
    touch(f"{resdir}/software/get_all_resources.done")

rule get_key_resources:
  input:
    rules.download_plink.output,
    rules.download_ldsc.output,
    rules.download_dbslmm.output,
    rules.download_ldak.output,
    rules.download_ldscores_panukb.output,
    rules.download_hm3_snplist.output,
    rules.download_ld_blocks.output,
    rules.download_ldak_map.output,
    rules.download_ldak_bld.output,
    rules.download_ldak_highld.output,
    rules.download_default_ref.output
  output:
    touch(f"{resdir}/software/get_key_resources.done")

rule get_prscs_resources:
  input:
    rules.get_key_resources.output,
    rules.download_prscs_software.output,
    rules.download_prscs_ref_ukb_all.input,
    rules.download_prscs_ref_1kg_all.input
  output:
    touch(f"{resdir}/software/get_prscs_resources.done")

rule get_ldpred2_resources:
  input:
    rules.get_key_resources.output,
    rules.download_ldpred2_ref.output
  output:
    touch(f"{resdir}/software/get_ldpred2_resources.done")

rule get_sbayesr_resources:
  input:
    rules.get_key_resources.output,
    rules.download_gctb_software.output,
    rules.download_gctb_ref.output
  output:
    touch(f"{resdir}/software/get_sbayesr_resources.done")
