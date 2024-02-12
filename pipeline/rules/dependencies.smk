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

######
# Check config file
######

# Check for the presence of at least one of gwas_list, target_list, or score_list
lists_to_check = ['gwas_list', 'target_list', 'score_list']
if not any(lst in config for lst in lists_to_check):
  print("Error: At least one of 'gwas_list', 'score_list', or 'target_list' must be specified in the config file.")
  sys.exit(1)

# Check for missing required configuration parameters
required_config_params = ['outdir', 'pgs_methods', 'config_file']
missing_config_params = [param for param in required_config_params if param not in config]
if missing_config_params:
  # Print an informative message
  print(f"Missing required configuration parameters: {', '.join(missing_config_params)}. Please specify these in the configuration file.")

  # Exit Snakemake gracefully
  sys.exit(1)

# Set outdir parameter
outdir=config['outdir']

# Set refdir parameter
# If refdir is NA, set refdir to resources/data/ref
if config['refdir'] == 'NA':
  refdir='resources/data/ref'
  ref_input="resources/data/ref/ref.pop.txt"
else:
  refdir=config['refdir']
  ref_input = [os.path.join(refdir, f"ref.chr{i}.{ext}") for i in range(1, 23) for ext in ['pgen', 'pvar', 'psam', 'rds']] + \
                 [os.path.join(refdir, file_name) for file_name in ['ref.pop.txt', 'ref.keep.list']]

  for full_path in ref_input:
      if not os.path.exists(full_path):
          raise FileNotFoundError(f"File not found: {full_path}. Check reference data format.")

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

########
# Check for repo version updates
########

# If there has been a change to the major or minor version numbers, we will rerun the entire pipeline

# Define the path for storing the last known version
os.makedirs("resources", exist_ok=True)
last_version_file = "resources/last_version.txt"

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

# Check for both major and minor version changes
if current_major != last_major or current_minor != last_minor:
    if not overwrite:
        print(f"Change in version of GenoPred detected from v{last_major}.{last_minor} to v{current_major}.{current_minor}. Use --params overwrite=true to proceed.")
        sys.exit(1)
    else:
        print("Proceeding with version update due to --overwrite flag.")
        write_last_version(current_major, current_minor)  # Update the stored version


########
# Download dependencies
########

# Download impute2_data
rule download_impute2_data:
  output:
    directory("resources/data/impute2/1000GP_Phase3")
  benchmark:
    "resources/data/benchmarks/download_impute2_data.txt"
  log:
    "resources/data/logs/download_impute2_data.log"
  shell:
    """
    {{
      rm -r -f resources/data/impute2; \
      mkdir -p resources/data/impute2/; \
      wget --no-check-certificate -O resources/data/impute2/1000GP_Phase3.tgz https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz; \
      tar -zxvf resources/data/impute2/1000GP_Phase3.tgz -C resources/data/impute2/; \
      rm resources/data/impute2/1000GP_Phase3.tgz; \
      wget --no-check-certificate -O resources/data/impute2/1000GP_Phase3_chrX.tgz https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz; \
      tar -zxvf resources/data/impute2/1000GP_Phase3_chrX.tgz -C resources/data/impute2/1000GP_Phase3/; \
      rm resources/data/impute2/1000GP_Phase3_chrX.tgz
    }} > {log} 2>&1
    """

# Download PLINK. DBSLMM requires the binary to be specified, which is challenging with conda environments. I have tried to avoid this again but no joy. The conda environment may not exist when the snakemake is executed which will cause problems if trying to access the conda environment manually.
rule download_plink:
  output:
    "resources/software/plink/plink"
  benchmark:
    "resources/data/benchmarks/download_plink.txt"
  log:
    "resources/data/logs/download_plink.log"
  shell:
    """
    {{
      rm -r -f resources/software/plink; \
      mkdir -p resources/software/plink; \
      wget --no-check-certificate -O resources/software/plink/plink_linux_x86_64_20210606.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip; \
      unzip resources/software/plink/plink_linux_x86_64_20210606.zip -d resources/software/plink; \
      rm resources/software/plink/plink_linux_x86_64_20210606.zip
    }} > {log} 2>&1
    """

# Download LDSC
rule download_ldsc:
  output:
    "resources/software/ldsc/ldsc.py"
  benchmark:
    "resources/data/benchmarks/download_ldsc.txt"
  log:
    "resources/data/logs/download_ldsc.log"
  shell:
    """
    {{
      rm -r -f resources/software/ldsc; \
      git clone https://github.com/bulik/ldsc.git resources/software/ldsc/; \
      cd resources/software/ldsc/; \
      git reset --hard aa33296abac9569a6422ee6ba7eb4b902422cc74
    }} > {log} 2>&1
    """
# Download LDSC reference data
rule download_ldsc_ref:
  output:
    directory("resources/data/ldsc_ref/eur_w_ld_chr")
  benchmark:
    "resources/data/benchmarks/download_ldsc_ref.txt"
  log:
    "resources/data/logs/download_ldsc_ref.log"
  shell:
    """
    {{
      rm -r -f resources/data/ldsc_ref; \
      mkdir -p resources/data/ldsc_ref; \
      wget --no-check-certificate -O resources/data/ldsc_ref/eur_w_ld_chr.tar.gz https://zenodo.org/record/8182036/files/eur_w_ld_chr.tar.gz?download=1; \
      tar -xvf resources/data/ldsc_ref/eur_w_ld_chr.tar.gz -C resources/data/ldsc_ref/; \
      rm resources/data/ldsc_ref/eur_w_ld_chr.tar.gz
    }} > {log} 2>&1
    """

# Download hapmap3 snplist
rule download_hm3_snplist:
  output:
    "resources/data/hm3_snplist/w_hm3.snplist"
  benchmark:
    "resources/data/benchmarks/download_hm3_snplist.txt"
  log:
    "resources/data/logs/download_hm3_snplist.log"
  shell:
    """
    {{
      rm -r -f resources/data/hm3_snplist; \
      mkdir -p resources/data/hm3_snplist; \
      wget --no-check-certificate -O resources/data/hm3_snplist/w_hm3.snplist.gz https://zenodo.org/record/7773502/files/w_hm3.snplist.gz?download=1; \
      gunzip resources/data/hm3_snplist/w_hm3.snplist.gz
    }} > {log} 2>&1
    """
# Download DBSLMM
# Specifying old commit as developer has deleted dbslmm binary (accidentally?)
rule download_dbslmm:
  output:
    directory("resources/software/dbslmm/")
  benchmark:
    "resources/data/benchmarks/download_dbslmm.txt"
  log:
    "resources/data/logs/download_dbslmm.log"
  shell:
    """
    {{
      git clone https://github.com/biostat0903/DBSLMM.git {output}; \
      cd {output}; \
      git reset --hard aa6e7ad5b8a7d3b6905556a4007c4a1fa2925b7d; \
      chmod a+x software/dbslmm
    }} > {log} 2>&1
    """
# Download LD block data
rule download_ld_blocks:
  output:
    directory("resources/data/ld_blocks/")
  benchmark:
    "resources/data/benchmarks/download_ld_blocks.txt"
  log:
    "resources/data/logs/download_ld_blocks.log"
  shell:
    """
    {{
      git clone https://bitbucket.org/nygcresearch/ldetect-data.git {output}; \
      mv resources/data/ld_blocks/ASN resources/data/ld_blocks/EAS
    }} > {log} 2>&1
    """

# Download PRScs reference
prscs_ref_dropbox = {
    'eur': 'https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=0',
    'eas': 'https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=0',
    'afr': 'https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=0',
    'amr': 'https://www.dropbox.com/s/y7ruj364buprkl6/ldblk_ukbb_amr.tar.gz?dl=0',
    'sas': 'https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz?dl=0',
}

rule download_prscs_ref_ukb:
  output:
    "resources/data/prscs_ref/ldblk_ukbb_{population}/ldblk_ukbb_chr1.hdf5"
  benchmark:
    "resources/data/benchmarks/download_prscs_ref_ukb-{population}.txt"
  log:
    "resources/data/logs/download_prscs_ref_ukb-{population}.log"
  params:
    url=lambda w: prscs_ref_dropbox.get(w.population)
  shell:
    """
    {{
      mkdir -p resources/data/prscs_ref; \
      rm -r -f resources/data/prscs_ref/ldblk_ukbb_{wildcards.population}; \
      wget --no-check-certificate -O resources/data/prscs_ref/ldblk_ukbb_{wildcards.population}.tar.gz {params.url}; \
      tar -zxvf resources/data/prscs_ref/ldblk_ukbb_{wildcards.population}.tar.gz -C resources/data/prscs_ref/; \
      rm resources/data/prscs_ref/ldblk_ukbb_{wildcards.population}.tar.gz
    }} > {log} 2>&1
    """

# Download PRScs software
rule download_prscs_software:
  output:
    directory("resources/software/prscs/")
  benchmark:
    "resources/data/benchmarks/download_prscs_software.txt"
  log:
    "resources/data/logs/download_prscs_software.log"
  shell:
    """
    {{
      git clone https://github.com/getian107/PRScs.git {output}; \
      cd {output}; \
      git reset --hard 621fdc80daac56c93d9528eb1a1187f7b1fc9afb
    }} > {log} 2>&1
    """
# Download gctb reference
rule download_gctb_ref:
  output:
    directory("resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse")
  benchmark:
    "resources/data/benchmarks/download_gctb_ref.txt"
  log:
    "resources/data/logs/download_gctb_ref.log"
  shell:
    """
    {{
      mkdir -p resources/data/gctb_ref; \
      wget --no-check-certificate -O resources/data/gctb_ref/ukbEURu_hm3_sparse.zip https://zenodo.org/record/3350914/files/ukbEURu_hm3_sparse.zip?download=1; \
      unzip resources/data/gctb_ref/ukbEURu_hm3_sparse.zip -d resources/data/gctb_ref; \
      for chr in $(seq 1 22);do \
      mv resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k.ldm.sparse.bin resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr${{chr}}.ldm.sparse.bin; \
      mv resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k.ldm.sparse.info resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr${{chr}}.ldm.sparse.info; \
      mv resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k_sparse.log resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_sparse_chr${{chr}}.log; \
      done; \
      rm resources/data/gctb_ref/ukbEURu_hm3_sparse.zip
    }} > {log} 2>&1
    """
# Download GCTB
rule download_gctb_software:
  output:
    "resources/software/gctb/gctb_2.03beta_Linux/gctb"
  benchmark:
    "resources/data/benchmarks/download_gctb_software.txt"
  log:
    "resources/data/logs/download_gctb_software.log"
  shell:
    """
    {{
      rm -r -f resources/software/gctb; \
      mkdir -p resources/software/gctb; \
      wget --no-check-certificate -O resources/software/gctb/gctb_2.03beta_Linux.zip https://cnsgenomics.com/software/gctb/download/gctb_2.03beta_Linux.zip; \
      unzip resources/software/gctb/gctb_2.03beta_Linux.zip -d resources/software/gctb; \
      rm resources/software/gctb/gctb_2.03beta_Linux.zip
    }} > {log} 2>&1
    """
# Download LDpred2 reference
rule download_ldpred2_ref:
  output:
    directory("resources/data/ldpred2_ref")
  benchmark:
    "resources/data/benchmarks/download_ldpred2_ref.txt"
  log:
    "resources/data/logs/download_ldpred2_ref.log"
  shell:
    """
    {{
      mkdir -p resources/data/ldpred2_ref; \
      wget --no-check-certificate -O resources/data/ldpred2_ref/download.zip https://figshare.com/ndownloader/articles/19213299/versions/2; \
      unzip resources/data/ldpred2_ref/download.zip -d resources/data/ldpred2_ref/; \
      rm resources/data/ldpred2_ref/download.zip; \
      unzip resources/data/ldpred2_ref/ldref_with_blocks.zip -d resources/data/ldpred2_ref/; \
      mv resources/data/ldpred2_ref/ldref/* resources/data/ldpred2_ref/; \
      rm resources/data/ldpred2_ref/ldref_with_blocks.zip; \
      rm -r resources/data/ldpred2_ref/ldref
    }} > {log} 2>&1
    """
# Download LDAK
rule download_ldak:
  output:
    "resources/software/ldak/ldak5.1.linux"
  benchmark:
    "resources/data/benchmarks/download_ldak.txt"
  log:
    "resources/data/logs/download_ldak.log"
  shell:
    """
    {{
      rm -r resources/software/ldak; \
      mkdir -p resources/software/ldak; \
      wget --no-check-certificate -O resources/software/ldak/ldak5.1.linux_.zip https://dougspeed.com/wp-content/uploads/ldak5.1.linux_.zip; \
      unzip resources/software/ldak/ldak5.1.linux_.zip -d resources/software/ldak/; \
      rm resources/software/ldak/ldak5.1.linux_.zip
    }} > {log} 2>&1
    """
# Download LDAK map data
rule download_ldak_map:
  output:
    "resources/data/ldak_map/genetic_map_b37/genetic_map_chr22_combined_b37.txt"
  benchmark:
    "resources/data/benchmarks/download_ldak_map.txt"
  log:
    "resources/data/logs/download_ldak_map.log"
  shell:
    """
    {{
      rm -r resources/data/ldak_map; \
      mkdir -p resources/data/ldak_map; \
      wget --no-check-certificate -O resources/data/ldak_map/genetic_map_b37.zip https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip; \
      unzip resources/data/ldak_map/genetic_map_b37.zip -d resources/data/ldak_map/; \
      rm resources/data/ldak_map/genetic_map_b37.zip
    }} > {log} 2>&1
    """
# Download LDAK bld snp annotations
rule download_ldak_bld:
  output:
    "resources/data/ldak_bld/README.txt"
  benchmark:
    "resources/data/benchmarks/download_ldak_bld.txt"
  log:
    "resources/data/logs/download_ldak_bld.log"
  shell:
    """
    {{
      rm -r resources/data/ldak_bld; \
      mkdir -p resources/data/ldak_bld; \
      wget --no-check-certificate -O resources/data/ldak_bld/bld.zip https://genetics.ghpc.au.dk/doug/bld.zip; \
      unzip resources/data/ldak_bld/bld.zip -d resources/data/ldak_bld/; \
      rm resources/data/ldak_bld/bld.zip
    }} > {log} 2>&1
    """
# Download LDAK high ld regions file
rule download_ldak_highld:
  output:
    "resources/data/ldak_highld/highld.txt"
  benchmark:
    "resources/data/benchmarks/download_ldak_highld.txt"
  log:
    "resources/data/logs/download_ldak_highld.log"
  shell:
    """
    {{
      rm -r resources/data/ldak_highld; \
      mkdir -p resources/data/ldak_highld; \
      wget --no-check-certificate -O resources/data/ldak_highld/highld.txt https://dougspeed.com/wp-content/uploads/highld.txt
    }} > {log} 2>&1
    """

# Download preprocessed reference data (1KG+HGDP HapMap3)
rule download_default_ref:
  output:
    "resources/data/ref/ref.pop.txt"
  benchmark:
    "resources/data/benchmarks/download_default_ref.txt"
  log:
    "resources/data/logs/download_default_ref.log"
  shell:
    """
    {{
      rm -r resources/data/ref; \
      mkdir -p resources/data/ref; \
      wget --no-check-certificate -O resources/data/ref/genopred_1kg_hgdp.tar.gz https://zenodo.org/records/10640523/files/genopred_1kg_hgdp.tar.gz?download=1; \
      tar -xzvf resources/data/ref/genopred_1kg_hgdp.tar.gz -C resources/data/ref/; \
      mv resources/data/ref/ref/* resources/data/ref/; \
      rm -r resources/data/ref/ref; \
      rm resources/data/ref/genopred_1kg_hgdp.tar.gz
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
      Rscript -e 'devtools::install_github(\"opain/GenoUtils@edf5bec1be0e396d953eb8974488dc4e3d57c134\")'
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
############
# Check all dependencies are available
############

# Rule to install and download all dependencies
rule get_dependencies:
  input:
    rules.download_plink.output,
    rules.download_ldsc.output,
    rules.download_ldsc_ref.output,
    rules.download_hm3_snplist.output,
    rules.download_default_ref.output,
    rules.download_dbslmm.output,
    rules.download_ld_blocks.output,
    rules.install_ggchicklet.output,
    rules.install_lassosum.output,
    rules.install_genoutils.output,
    rules.download_pgscatalog_utils.output
  output:
    touch("resources/software/get_dependencies.done")

