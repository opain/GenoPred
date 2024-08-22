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

# Set refdir parameter
# If refdir is NA, set refdir to '${resdir}/data/ref'
if config['refdir'] == 'NA':
  refdir=f"{resdir}/data/ref"
  ref_input=f"{refdir}/ref.pop.txt"
else:
  refdir=config['refdir']
  ref_input = [os.path.join(refdir, f"ref.chr{i}.{ext}") for i in get_chr_range(testing = config['testing']) for ext in ['pgen', 'pvar', 'psam', 'rds']] + \
                 [os.path.join(refdir, file_name) for file_name in ['ref.pop.txt', 'ref.keep.list']]

  for full_path in ref_input:
      if not os.path.exists(full_path):
          raise FileNotFoundError(f"File not found: {full_path}. Check reference data format.")

# Check valid pgs_methods are specified
def check_pgs_methods(x):
    valid_pgs_methods = {
        "ptclump", "dbslmm", "prscs", "sbayesr", "lassosum", "ldpred2", "megaprs", "xwing", "prscsx", "tlprs"
    }

    invalid_methods = [method for method in x if method not in valid_pgs_methods]

    if invalid_methods:
        raise ValueError(f"Invalid pgs_methods specified: {', '.join(invalid_methods)}. "
                         f"Valid methods are: {', '.join(valid_pgs_methods)}.")

check_pgs_methods(config['pgs_methods'])

# Check valid tlprs_methods are specified
def check_tlprs_methods(config):
    valid_tlprs_methods = {
        "ptclump", "dbslmm", "prscs", "sbayesr", "lassosum", "ldpred2", "megaprs"
    }

    # Check if 'tlprs' is in the pgs_methods list
    if 'tlprs' in config.get('pgs_methods', []):
        # Check if tlprs_methods is defined and not None/NA
        tlprs_methods = config.get('tlprs_methods')

        if tlprs_methods is None or tlprs_methods == 'NA':
            raise ValueError("tlprs_methods must be specified when 'tlprs' is included in pgs_methods.")

        # Check for invalid methods
        invalid_methods = [method for method in tlprs_methods if method not in valid_tlprs_methods]

        if invalid_methods:
            raise ValueError(f"Invalid tlprs_methods specified: {', '.join(invalid_methods)}. "
                             f"Valid methods are: {', '.join(valid_tlprs_methods)}.")

check_tlprs_methods(config)

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
# Download LDpred2 reference
rule download_ldpred2_ref:
  output:
    directory(f"{resdir}/data/ldpred2_ref")
  benchmark:
    f"{resdir}/data/benchmarks/download_ldpred2_ref.txt"
  log:
    f"{resdir}/data/logs/download_ldpred2_ref.log"
  shell:
    """
    {{
      mkdir -p {resdir}/data/ldpred2_ref; \
      wget --no-check-certificate -O {resdir}/data/ldpred2_ref/download.zip https://figshare.com/ndownloader/articles/19213299/versions/2; \
      unzip {resdir}/data/ldpred2_ref/download.zip -d {resdir}/data/ldpred2_ref/; \
      rm {resdir}/data/ldpred2_ref/download.zip; \
      unzip {resdir}/data/ldpred2_ref/ldref_with_blocks.zip -d {resdir}/data/ldpred2_ref/; \
      mv {resdir}/data/ldpred2_ref/ldref/* {resdir}/data/ldpred2_ref/; \
      rm {resdir}/data/ldpred2_ref/ldref_with_blocks.zip; \
      rm -r {resdir}/data/ldpred2_ref/ldref
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
      rm -r {resdir}/software/ldak; \
      mkdir -p {resdir}/software/ldak; \
      wget --no-check-certificate -O {resdir}/software/ldak/ldak5.1.linux_.zip https://dougspeed.com/wp-content/uploads/ldak5.1.linux_.zip; \
      unzip {resdir}/software/ldak/ldak5.1.linux_.zip -d {resdir}/software/ldak/; \
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
      wget --no-check-certificate -O {resdir}/data/ref/genopred_1kg_hgdp.tar.gz https://zenodo.org/records/10666983/files/genopred_1kg_hgdp.tar.gz?download=1; \
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
      Rscript -e 'devtools::install_github(\"opain/GenoUtils@50ac8a2078226c8c2349064f904031576fbfe606\")'
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
    {{
      rm -r -f resources/software/pgscatalog_utils; \
      git clone https://github.com/PGScatalog/pgscatalog_utils.git resources/software/pgscatalog_utils; \
      cd resources/software/pgscatalog_utils; \
      git reset --hard 6da7eb0e157ba4e73f941233ee8d8ae4fb5e3926; \
      poetry install; \
      poetry build; \
      pip3 install --user dist/*.whl; \
      download_scorefiles -h > download_pgscatalog_utils.done
    }} > {log} 2>&1
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
      Rscript -e 'devtools::install_github(\"opain/GenoUtils@50ac8a2078226c8c2349064f904031576fbfe606\")'
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
    "resources/data/benchmarks/install_tlprs.txt"
  log:
    "resources/data/logs/install_tlprs.log"
  shell:
    """
    {{
      Rscript -e 'devtools::install_github(\"opain/TLPRS@fb71076267d405d5f7df97e445ab0de73d76bc0f\")'
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
    rules.download_default_ref.output
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
