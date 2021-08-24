##########
# Pipeline Preparation
##########

####
# Download dependencies
####

# Download qctool v2
rule download_qctool2:
  output:
    "resources/software/qctool2/qctool"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/software/qctool2/; \
     wget -O resources/software/qctool2/qctool2.tgz https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz; \
     tar -zxvf resources/software/qctool2/qctool2.tgz -C resources/software/qctool2/ --strip-components=1; \
     rm resources/software/qctool2/qctool2.tgz"

# Download impute2_data
rule download_impute2_data:
  output:
    directory("resources/data/impute2/1000GP_Phase3")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/data/impute2/; \
     wget -O resources/data/impute2/1000GP_Phase3.tgz https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz; \
     tar -zxvf resources/data/impute2/1000GP_Phase3.tgz -C resources/data/impute2/; \
     rm resources/data/impute2/1000GP_Phase3.tgz; \
     wget -O resources/data/impute2/1000GP_Phase3_chrX.tgz https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz; \
     tar -zxvf resources/data/impute2/1000GP_Phase3_chrX.tgz -C resources/data/impute2/1000GP_Phase3/; \
     rm resources/data/impute2/1000GP_Phase3_chrX.tgz"

# Download PLINK. DBSLMM requires the binary to be specified, which is challenging with conda environments.
rule download_plink:
  output:
    "resources/software/plink/plink"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/software/plink; \
     wget -O resources/software/plink/plink_linux_x86_64_20210606.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip; \
     unzip resources/software/plink/plink_linux_x86_64_20210606.zip -d resources/software/plink; \
     rm resources/software/plink/plink_linux_x86_64_20210606.zip"

# Download LDSC
rule download_ldsc:
  output:
    "resources/software/ldsc/ldsc.py"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "git clone git@github.com:bulik/ldsc.git resources/software/ldsc/"

# Download LDSC reference data
rule dowload_ldsc_ref:
  output:
    directory("resources/data/ldsc_ref/eur_w_ld_chr")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/data/ldsc_ref; \
     wget -O resources/data/ldsc_ref/eur_w_ld_chr.tar.bz2 https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2; \
     tar -jxvf resources/data/ldsc_ref/eur_w_ld_chr.tar.bz2 -C resources/data/ldsc_ref/; \
     rm resources/data/ldsc_ref/eur_w_ld_chr.tar.bz2"

# Download hapmap3 snplist
rule download_hm3_snplist:
  output:
    "resources/data/hm3_snplist/w_hm3.snplist"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/data/hm3_snplist; \
     wget -O resources/data/hm3_snplist/w_hm3.snplist.bz2 https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2; \
     bunzip2 resources/data/hm3_snplist/w_hm3.snplist.bz2"

# Download DBSLMM
# Specifying old commit as developer has deleted dbslmm binary (accidentally?)
rule download_dbslmm:
  output:
    directory("resources/software/dbslmm/")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "git clone git@github.com:biostat0903/DBSLMM.git {output}; \
     cd {output}; \
     git reset --hard aa6e7ad5b8a7d3b6905556a4007c4a1fa2925b7d; \
     chmod a+x software/dbslmm"

# Download LD block data
rule download_ld_blocks:
  output:
    directory("resources/data/ld_blocks/")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "git clone https://bitbucket.org/nygcresearch/ldetect-data.git {output}"

# install lassosum
rule install_lassosum:
  output:
    touch("resources/software/install_lassosum.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e 'remotes::install_github(\"tshmak/lassosum\")'"

# install bigsnpr
rule install_bigsnpr:
  output:
    touch("resources/software/install_bigsnpr.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e 'remotes::install_github(\"privefl/bigsnpr\")'"

# install ggchicklet
rule install_ggchicklet:
  output:
    touch("resources/software/install_ggchicklet.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e 'remotes::install_github(\"hrbrmstr/ggchicklet\")'"

# Download PRScs reference
rule download_prscs_ref_1kg_eur:
  output:
    "resources/data/prscs_ref/ldblk_1kg_eur/ldblk_1kg_chr1.hdf5"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/data/prscs_ref; \
     wget -O resources/data/prscs_ref/ldblk_1kg_eur.tar.gz https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0; \
     tar -zxvf resources/data/prscs_ref/ldblk_1kg_eur.tar.gz -C resources/data/prscs_ref/; \
     rm resources/data/prscs_ref/ldblk_1kg_eur.tar.gz"

# Download PRScs software
rule download_prscs_software:
  output:
    directory("resources/software/prscs/")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "git clone https://github.com/getian107/PRScs.git {output}"

# Download gctb reference
rule download_gctb_ref:
  output:
    directory("resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/data/gctb_ref; wget -O resources/data/gctb_ref/ukbEURu_hm3_sparse.zip https://zenodo.org/record/3350914/files/ukbEURu_hm3_sparse.zip?download=1; unzip resources/data/gctb_ref/ukbEURu_hm3_sparse.zip -d resources/data/gctb_ref; for chr in $(seq 1 22);do mv resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k.ldm.sparse.bin resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr${{chr}}.ldm.sparse.bin; mv resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k.ldm.sparse.info resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr${{chr}}.ldm.sparse.info; mv resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${{chr}}_v3_50k_sparse.log resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_sparse_chr${{chr}}.log; done; rm resources/data/gctb_ref/ukbEURu_hm3_sparse.zip"

# Download GCTB
rule download_gctb_software:
  output:
    "resources/software/gctb/gctb_2.03beta_Linux/gctb"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/software/gctb; \
     wget -O resources/software/gctb/gctb_2.03beta_Linux.zip https://cnsgenomics.com/software/gctb/download/gctb_2.03beta_Linux.zip; \
     unzip resources/software/gctb/gctb_2.03beta_Linux.zip -d resources/software/gctb; \
     rm resources/software/gctb/gctb_2.03beta_Linux.zip"

# Download 7z
rule download_7z:
  output:
    "resources/software/7z/7zz"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/software/7z; \
     wget -O resources/software/7z/7z2103-linux-x64.tar.xz https://www.7-zip.org/a/7z2103-linux-x64.tar.xz; \
     tar -xf resources/software/7z/7z2103-linux-x64.tar.xz -C resources/software/7z/; \
     rm resources/software/7z/7z2103-linux-x64.tar.xz"

# Download LDpred2 reference
rule download_ldpred2_ref:
  input:
    rules.download_7z.output
  output:
    directory("resources/data/ldpred2_ref/ukbEURu_hm3_shrunk_sparse")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/data/ldpred2_ref; wget -O resources/data/ldpred2_ref/download.zip https://ndownloader.figshare.com/articles/13034123/versions/3; resources/software/7z/7zz e resources/data/ldpred2_ref/download.zip -oresources/data/ldpred2_ref/; rm resources/data/ldpred2_ref/download.zip"

# Download and format 1kg reference data
rule prep_1kg:
  resources: 
    mem_mb=20000
  input:
    rules.download_hm3_snplist.output
  output:
    directory("resources/data/1kg/freq_files"),directory("resources/data/1kg/keep_files")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript scripts/prep_1kg.R"

# Create GW merged version of 1kg refrence data
rule merge_1kg_GW:
  input:
    rules.prep_1kg.output
  output:
    "resources/data/1kg/1KGPhase3.w_hm3.GW.bed"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "ls resources/data/1kg/1KGPhase3.w_hm3.chr*.bed | sed -e 's/\.bed//g' > resources/data/1kg/merge_list.txt; \
     plink \
       --merge-list resources/data/1kg/merge_list.txt \
       --make-bed \
       --out resources/data/1kg/1KGPhase3.w_hm3.GW; \
     rm resources/data/1kg/merge_list.txt"

# Download 1kg population code definitions
rule download_1kg_pop_codes:
  output:
    "resources/data/1kg/1kg_pop_codes.tsv",
    "resources/data/1kg/1kg_super_pop_codes.tsv"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "wget -O resources/data/1kg/1kg_pop_codes.tsv http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv; \
     wget -O resources/data/1kg/1kg_super_pop_codes.tsv http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.superpopulations.tsv"

####
# Principal Component (Ancestry) Scoring
####

# Create PC score files specific to each super population
rule super_pop_pc_scoring:
  input:
    ref=rules.prep_1kg.output
  output:
    "resources/data/1kg/pc_score_files/{population}/1KGPhase3.w_hm3.{population}.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_keep resources/data/1kg/keep_files/{wildcards.population}_samples.keep \
      --plink plink \
      --plink2 plink2 \
      --n_pcs 100 \
      --output resources/data/1kg/pc_score_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}"

populations=["AFR","AMR","EAS","EUR","SAS"]

rule run_super_pop_pc_scoring:
  input: expand("resources/data/1kg/pc_score_files/{population}/1KGPhase3.w_hm3.{population}.scale", population=populations)

####
# Polygenic Scoring
####

##
# QC and format GWAS summary statistics
##

import pandas as pd
gwas_list_df = pd.read_table(config["gwas_list"], sep=' ')
gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

rule sumstat_prep:
  input:
    rules.prep_1kg.output
  output:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    path= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
  shell:
    "Rscript ../Scripts/sumstat_cleaner/sumstat_cleaner.R \
      --sumstats {params.path} \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_freq_chr resources/data/1kg/freq_files/{params.population}/1KGPhase3.w_hm3.{params.population}.chr \
      --output resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned"
    
rule run_sumstat_prep:
  input: expand("resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz", gwas=gwas_list_df['name'])

##
# pT+clump (sparse, nested)
##

rule prs_scoring_ptclump:
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz"
  output:
    "resources/data/1kg/prs_score_files/pt_clump/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator/polygenic_score_file_creator.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink plink \
      --output resources/data/1kg/prs_score_files/pt_clump/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"
  
rule run_prs_scoring_ptclump:
  input: expand("resources/data/1kg/prs_score_files/pt_clump/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# DBSLMM
##

rule prs_scoring_dbslmm:
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_ldsc.output,
    rules.dowload_ldsc_ref.output,
    rules.download_dbslmm.output,
    rules.download_ld_blocks.output,
    rules.download_plink.output
  output:
    "resources/data/1kg/prs_score_files/dbslmm/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    sampling= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'sampling'].iloc[0],
    prevalence= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'prevalence'].iloc[0],
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator_DBSLMM/polygenic_score_file_creator_DBSLMM.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink resources/software/plink/plink \
      --ld_blocks resources/data/ld_blocks/{params.population} \
      --rscript Rscript \
      --dbslmm resources/software/dbslmm/software \
      --munge_sumstats resources/software/ldsc/munge_sumstats.py \
      --ldsc resources/software/ldsc/ldsc.py \
      --ldsc_ref resources/data/ldsc_ref/eur_w_ld_chr \
      --hm3_snplist resources/data/hm3_snplist/w_hm3.snplist \
      --sample_prev {params.sampling} \
      --pop_prev {params.prevalence} \
      --output resources/data/1kg/prs_score_files/dbslmm/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"

rule run_prs_scoring_dbslmm:
  input: expand("resources/data/1kg/prs_score_files/dbslmm/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# PRScs
##

rule prs_scoring_prscs:
  resources: 
    mem_mb=110000,
    cpus=6
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_plink.output,
    rules.download_prscs_ref_1kg_eur.output,
    rules.download_prscs_software.output
  output:
    "resources/data/1kg/prs_score_files/prscs/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "export MKL_NUM_THREADS=1; \
     export NUMEXPR_NUM_THREADS=1; \
     export OMP_NUM_THREADS=1; \
     Rscript ../Scripts/polygenic_score_file_creator_PRScs/polygenic_score_file_creator_PRScs.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink plink \
      --memory 5000 \
      --output resources/data/1kg/prs_score_files/prscs/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list \
      --PRScs_path resources/software/prscs/PRScs.py \
      --PRScs_ref_path resources/data/prscs_ref/ldblk_1kg_eur \
      --n_cores {resources.cpus} \
      --phi_param auto"
 
# --phi_param 1e-6,1e-4,1e-2,1,auto"
   
rule run_prs_scoring_prscs:
  input: expand("resources/data/1kg/prs_score_files/prscs/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# SBayesR
##

rule prs_scoring_sbayesr:
  resources: 
    mem_mb=50000,
    cpus=6
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_plink.output,
    rules.download_gctb_ref.output,
    rules.download_gctb_software.output
  output:
    "resources/data/1kg/prs_score_files/sbayesr/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator_SBayesR/polygenic_score_file_creator_SBayesR.R \
      --ref_plink resources/data/1kg/1KGPhase3.w_hm3.GW \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink plink \
      --gctb resources/software/gctb/gctb_2.03beta_Linux/gctb \
      --ld_matrix_chr resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr \
      --memory {resources.mem_mb} \
      --robust T \
      --n_cores {resources.cpus} \
      --output resources/data/1kg/prs_score_files/sbayesr/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"

rule run_prs_scoring_sbayesr:
  input: expand("resources/data/1kg/prs_score_files/sbayesr/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# lassosum
##

rule prs_scoring_lassosum:
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_plink.output,
    rules.install_lassosum.output
  output:
    "resources/data/1kg/prs_score_files/lassosum/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator_lassosum/polygenic_score_file_creator_lassosum.R \
     --ref_plink_gw resources/data/1kg/1KGPhase3.w_hm3.GW \
     --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
     --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
     --output resources/data/1kg/prs_score_files/lassosum/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
     --plink plink \
     --ref_pop_scale resources/data/1kg/super_pop_keep.list"
    
rule run_prs_scoring_lassosum:
  input: expand("resources/data/1kg/prs_score_files/lassosum/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# LDpred2
##

rule prs_scoring_ldpred2:
  resources: 
    mem_mb=70000,
    cpus=10
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_plink.output,
    rules.install_bigsnpr.output,
    rules.download_ldpred2_ref.output
  output:
    "resources/data/1kg/prs_score_files/ldpred2/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator_LDPred2/polygenic_score_file_creator_LDPred2_LDPredRef.R \
      --ref_plink resources/data/1kg/1KGPhase3.w_hm3.GW \
      --ref_keep resources/data/1kg/keep_files/EUR_samples.keep \
      --ldpred2_ref_dir resources/data/ldpred2_ref \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink plink \
      --memory {resources.mem_mb} \
      --n_cores {resources.cpus} \
      --output resources/data/1kg/prs_score_files/ldpred2/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"
    
rule run_prs_scoring_ldpred2:
  input: expand("resources/data/1kg/prs_score_files/ldpred2/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# Estimate R2/AUC of PRS using lassosum pseudovalidate
##

rule pseudovalidate_prs:
  input:
    rules.merge_1kg_GW.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.install_lassosum.output
  output:
    "resources/data/1kg/prs_pseudoval/{gwas}/lassosum_pseudo_{gwas}.pseudovalidate.png"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/lassosum_pseudovalidate/lassosum_pseudovalidate.R \
      --ref_plink_gw resources/data/1kg/1KGPhase3.w_hm3.GW \
      --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --prune_mhc T \
      --output resources/data/1kg/prs_pseudoval/{wildcards.gwas}/lassosum_pseudo_{wildcards.gwas} \
      --plink plink \
      --n_cores 1"

rule run_pseudovalidate_prs:
  input: expand("resources/data/1kg/prs_pseudoval/{gwas}/lassosum_pseudo_{gwas}.pseudovalidate.png", gwas=gwas_list_df['name'])

rule pipeline_prep:
  input:
    rules.run_super_pop_pc_scoring.input,
    rules.run_prs_scoring_ptclump.input,
    rules.run_prs_scoring_dbslmm.input,
    rules.run_pseudovalidate_prs.input

##########
# Target sample processing
##########

import pandas as pd
target_list_df = pd.read_table(config["target_list"], sep=' ')

target_list_df_23andMe = target_list_df.loc[target_list_df['type'] == '23andMe']
target_list_df_samp_imp_plink1 = target_list_df.loc[target_list_df['type'] == 'samp_imp_plink1']
target_list_df_samp_imp_bgen = target_list_df.loc[target_list_df['type'] == 'samp_imp_bgen']
target_list_df_samp_imp = target_list_df.loc[(target_list_df['type'] == 'samp_imp_plink1') | (target_list_df['type'] == 'samp_imp_bgen')]

target_list_df['pre_harm_path'] = target_list_df['path']
target_list_df.loc[target_list_df['type'] == '23andMe', 'pre_harm_path'] = "resources/data/target/" + target_list_df.loc[target_list_df['type'] == '23andMe', 'name'] + "/" + target_list_df.loc[target_list_df['type'] == '23andMe', 'name'] + ".1KGphase3"
target_list_df.loc[target_list_df['type'] == 'samp_imp_bgen', 'pre_harm_path'] = "resources/data/target/" + target_list_df.loc[target_list_df['type'] == 'samp_imp_bgen', 'name'] + "/" + target_list_df.loc[target_list_df['type'] == 'samp_imp_bgen', 'name'] + ".w_hm3"

####
# Format target data
####

##
# 23andMe
##
# Largely based on Impute.Me by Lasse Folkersen

rule format_impute_23andme_target:
  resources: 
    mem_mb=90000,
    cpus=6
  input:
    rules.download_impute2_data.output,
    rules.download_qctool2.output
  output:
    "resources/data/target/{name}/{name}.sample"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    name= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'name'].iloc[0],
    path= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0],
  shell:
    "Rscript ../Scripts/23andMe_imputer/23andMe_imputer.R \
      --FID {params.name} \
      --IID {params.name} \
      --geno {params.path} \
      --plink plink \
      --out resources/data/target/{params.name}/{params.name} \
      --ref resources/data/impute2/1000GP_Phase3 \
      --shapeit shapeit \
      --impute2 impute2 \
      --n_core {resources.cpus} \
      --qctool resources/software/qctool2/qctool"

rule run_format_impute_23andme_target:
  input: expand("resources/data/target/{name}/{name}.sample", name=target_list_df_23andMe['name'])

#opt$FID<-'Oliver_Pain'
#opt$IID<-'Oliver_Pain'
#opt$geno<-'/users/k1806347/brc_scratch/Data/23andMe_personal/Oliver_Pain/genome_Oliver_Pain_v5_Full_20210804054305.zip'
#opt$plink<-'plink'
#opt$out<-'resources/data/target/Oliver_Pain/Oliver_Pain'
#opt$ref<-'resources/data/impute2/1000GP_Phase3'
#opt$shapeit<-'shapeit'
#opt$impute2<-'impute2' 
#opt$n_core<-6
#opt$qctool<-'resources/software/qctool2/qctool'

##
# Note. one chunk failed for me due to insufficient memory
# Insert code to check for stings in summary files hat indicate completion
# i.e. grep -IRiL "There are no SNPs\|Imputation accuracy assessment" *_summary
##

##
# samp_imp_plink1
##

# Touch non-23andMe data
rule touch_imp:
  output:
    touch("resources/data/target/{name}/touch_imp.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    path=lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0]
  shell:
    "ls {params.path}*"

rule run_touch_imp:
  input: expand("resources/data/target/{name}/touch_imp.done", name=target_list_df_samp_imp_plink1['name'])

##
# samp_imp_bgen
##

# Estimate MAF and INFO
rule compute_snp_stats_target:
  input:
    rules.download_qctool2.output
  output:
    "resources/data/target/{name}/{name}_snp_stats_chr{chr}.txt"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    path=lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0]
  shell:
    "resources/software/qctool2/qctool \
      -g {params.path}.chr{wildcards.chr}.bgen \
      -s {params.path}.sample \
      -snp-stats \
      -osnp resources/data/target/{wildcards.name}/{wildcards.name}_snp_stats_chr{wildcards.chr}.txt"

rule run_compute_snp_stats_target:
  input: expand("resources/data/target/{name}/{name}_snp_stats_chr{chr}.txt", name=target_list_df_samp_imp_bgen['name'], chr=range(1, 23))

# NOTE. I do not implement the run_compute_snp_stats_target rule as target sample size may be too small to accurately estimate INFO and MAF. I suggest people perform this QC in advance, if at all. The use of HapMap3 variants in subsequent analyses should make these filters far less important.

# Convert to plink1 binary
rule convert_bgen_target:
  input:
    rules.download_qctool2.output
  output:
    "resources/data/target/{name}/{name}.w_hm3.chr{chr}.bed"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    path=lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0]
  shell:
    "resources/software/qctool2/qctool \
      -g {params.path}.chr{wildcards.chr}.bgen \
      -s {params.path}.sample \
      -incl-rsids resources/data/hm3_snplist/w_hm3.snplist \
      -ofiletype binary_ped \
      -og resources/data/target/{wildcards.name}/{wildcards.name}.w_hm3.chr{wildcards.chr} \
      -threshold 0.9; \
     awk \'$1=$2\' resources/data/target/{wildcards.name}/{wildcards.name}.w_hm3.chr{wildcards.chr}.fam > resources/data/target/{wildcards.name}/{wildcards.name}.w_hm3.chr{wildcards.chr}.fam.tmp; \
     mv resources/data/target/{wildcards.name}/{wildcards.name}.w_hm3.chr{wildcards.chr}.fam.tmp resources/data/target/{wildcards.name}/{wildcards.name}.w_hm3.chr{wildcards.chr}.fam"

rule run_convert_bgen_target:
  input: 
    lambda w: expand("resources/data/target/{name}/{name}.w_hm3.chr{chr}.bed", name=w.name, chr=range(1, 23))
  output:
    touch("resources/data/target/{name}/converted_to_plink.done")

rule run_convert_bgen_target_2:
  input: 
    expand("resources/data/target/{name}/converted_to_plink.done", name=target_list_df_samp_imp_bgen['name'])

####
# Harmonise with reference
####

rule harmonise_target_with_ref:
  input:
    rules.prep_1kg.output,
    lambda w: "resources/data/target/" + w.name + "/" + w.name + ".sample" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == '23andMe' else ("resources/data/target/{name}/touch_imp.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == 'samp_imp_plink1' else "resources/data/target/{name}/converted_to_plink.done")
  output:
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    path=lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'pre_harm_path'].iloc[0]
  shell:
    "Rscript ../Scripts/hm3_harmoniser/hm3_harmoniser.R \
      --target {params.path}.chr \
      --ref resources/data/1kg/1KGPhase3.w_hm3.chr \
      --plink plink \
      --out resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr"

rule run_harmonise_target_with_ref:
  input: expand("resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed", name=target_list_df['name'])

####
# Identify super_population
####

rule target_super_population:
  input:
    rules.prep_1kg.output,
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed"
  output:
    "resources/data/target/{name}/ancestry/ancestry_all/{name}.Ancestry.SAS.eigenvec"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --n_pcs 6 \
      --plink plink \
      --plink2 plink2 \
      --output resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list \
      --pop_data resources/data/1kg/ref_pop_dat.txt \
      --prob_thresh 0.5"

rule run_target_super_population:
  input: expand("resources/data/target/{name}/ancestry/ancestry_all/{name}.Ancestry.SAS.eigenvec", name=target_list_df['name'])

# Create a file listing target samples and super population assignments
checkpoint ancestry_reporter:
  input:
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry/ancestry_all/{name}.Ancestry.SAS.eigenvec"
  output:
    "resources/data/target/{name}/ancestry/ancestry_report.txt"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript scripts/ancestry_reporter.R {wildcards.name}"

rule run_ancestry_reporter:
  input: expand("resources/data/target/{name}/ancestry/ancestry_report.txt", name=target_list_df['name'])

####
# Identify population within super population
####

from pathlib import Path

def ancestry_munge(x):
    checkpoint_output = checkpoints.ancestry_reporter.get(name=x).output[0]
    ancestry_report_df = pd.read_table(checkpoint_output, sep=' ')
    return ancestry_report_df['population'].tolist()

rule target_population:
  input:
    rules.prep_1kg.output,
    "resources/data/target/{name}/ancestry/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed"
  output:
    touch("resources/data/target/{name}/ancestry/ancestry_{population}/ancestry_{population}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --n_pcs 6 \
      --plink plink \
      --plink2 plink2 \
      --output resources/data/target/{wildcards.name}/ancestry/ancestry_{wildcards.population}/{wildcards.name}.Ancestry \
      --ref_pop_scale resources/data/1kg/pop_keep_for_{wildcards.population}.list \
      --pop_data resources/data/1kg/ref_pop_dat_for_{wildcards.population}.txt \
      --prob_thresh 0.5"

rule run_target_population_all_pop:
  input: 
    lambda w: expand("resources/data/target/{name}/ancestry/ancestry_{population}/ancestry_{population}.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target/{name}/ancestry/ancestry_all_pop.done")

rule run_target_population_all:
  input: 
    expand("resources/data/target/{name}/ancestry/ancestry_all_pop.done", name=target_list_df['name'])

##########
# Create ancestry-only reports
##########

## 
# Individual ancestry reports
##

rule create_individual_ancestry_report:
  input:
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry/ancestry_all_pop.done",
    rules.download_1kg_pop_codes.output
  output:
    touch('resources/data/target/{name}/{name}_indiv_ancestry_report.done') 
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e \"rmarkdown::render(\'scripts/indiv_ancestry_report_creator.Rmd\', \
     output_file = \'../resources/data/target/{wildcards.name}/{wildcards.name}_ancestry_report.html\', \
     params = list(name = \'{wildcards.name}\'))\""

rule run_create_individual_ancestry_report:
  input: 
    expand('resources/data/target/{name}/{name}_indiv_ancestry_report.done', name=target_list_df_23andMe['name'])

##
# Sample ancestry reports
##

rule create_sample_ancestry_report:
  input:
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry/ancestry_all_pop.done",
    rules.download_1kg_pop_codes.output
  output:
    touch('resources/data/target/{name}/{name}_samp_ancestry_report.done')
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e \"rmarkdown::render(\'scripts/samp_ancestry_report_creator.Rmd\', \
    output_file = \'../resources/data/target/{wildcards.name}/{wildcards.name}_ancestry_report.html\', \
    params = list(name = \'{wildcards.name}\'))\""

rule run_create_sample_ancestry_report:
  input: 
    expand('resources/data/target/{name}/{name}_samp_ancestry_report.done', name=target_list_df_samp_imp['name'])

rule run_create_ancestry_reports:
  input: 
    rules.run_create_individual_ancestry_report.input,
    rules.run_create_sample_ancestry_report.input

##########
# Target sample scoring
##########

####
# PC scoring
####

rule target_pc:
  input:
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry/ancestry_report.txt",
    rules.run_super_pop_pc_scoring.input
  output:
    touch("resources/data/target/{name}/projected_pc/{population}/{name}_{population}_target_pc.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/scaled_ancestry_scorer/scaled_ancestry_scorer.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_freq_chr resources/data/1kg/freq_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.chr \
      --ref_score resources/data/1kg/pc_score_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.eigenvec.var \
      --ref_eigenvec resources/data/1kg/pc_score_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.eigenvec \
      --ref_scale resources/data/1kg/pc_score_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.scale \
      --plink2 plink2 \
      --output resources/data/target/{wildcards.name}/projected_pc/{wildcards.population}/{wildcards.name} \
      --pop_data resources/data/1kg/ref_pop_dat.txt"

rule run_target_pc_all_pop:
  input: 
    lambda w: expand("resources/data/target/{name}/projected_pc/{population}/{name}_{population}_target_pc.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target/{name}/projected_pc/{name}_target_pc.done")

rule run_target_pc_all:
  input: 
    expand("resources/data/target/{name}/projected_pc/{name}_target_pc.done", name=target_list_df['name'])

####
# Polygenic scoring
####

##
# pt+clump
##

rule target_prs_ptclump:
  input:
    "resources/data/target/{name}/ancestry/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/1kg/prs_score_files/pt_clump/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  output:
    touch("resources/data/target/{name}/prs/{population}/ptclump/{gwas}/{name}_prs_{population}_ptclump_{gwas}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/Scaled_polygenic_scorer/Scaled_polygenic_scorer.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_score resources/data/1kg/prs_score_files/pt_clump/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_scale resources/data/1kg/prs_score_files/pt_clump/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.{wildcards.population}.scale \
      --ref_freq_chr resources/data/1kg/freq_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.chr \
      --plink plink \
      --pheno_name {wildcards.gwas} \
      --output resources/data/target/{wildcards.name}/prs/{wildcards.population}/pt_clump/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.{wildcards.population}"

rule run_target_prs_ptclump_all_gwas:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/ptclump/{gwas}/{name}_prs_{population}_ptclump_{gwas}.done", name=w.name, gwas=gwas_list_df['name'], population=w.population)
  output:
    touch("resources/data/target/{name}/prs/{population}/ptclump/{name}_prs_{population}_ptclump.done")

rule run_target_prs_ptclump_all_pop:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/ptclump/{name}_prs_{population}_ptclump.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target/{name}/prs/{name}_prs_ptclump.done")

rule run_target_prs_ptclump_all_name:
  input: 
    expand("resources/data/target/{name}/prs/{name}_prs_ptclump.done", name=target_list_df['name'])
  output:
    touch("resources/data/target/prs_ptclump.done")

##
# DBSLMM
##

rule target_prs_dbslmm:
  input:
    "resources/data/target/{name}/ancestry/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/1kg/prs_score_files/dbslmm/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  output:
    touch("resources/data/target/{name}/prs/{population}/dbslmm/{gwas}/{name}_prs_{population}_dbslmm_{gwas}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/Scaled_polygenic_scorer_DBSLMM/Scaled_polygenic_scorer_DBSLMM.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_score resources/data/1kg/prs_score_files/dbslmm/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.dbslmm.GW.txt \
      --ref_scale resources/data/1kg/prs_score_files/dbslmm/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.{wildcards.population}.scale \
      --ref_freq_chr resources/data/1kg/freq_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.chr \
      --plink plink \
      --pheno_name {wildcards.gwas} \
      --output resources/data/target/{wildcards.name}/prs/{wildcards.population}/dbslmm/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.{wildcards.population}"

rule run_target_prs_dbslmm_all_gwas:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/dbslmm/{gwas}/{name}_prs_{population}_dbslmm_{gwas}.done", name=w.name, gwas=gwas_list_df['name'], population=w.population)
  output:
    touch("resources/data/target/{name}/prs/{population}/dbslmm/{name}_prs_{population}_dbslmm.done")

rule run_target_prs_dbslmm_all_pop:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/dbslmm/{name}_prs_{population}_dbslmm.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target/{name}/prs/{name}_prs_dbslmm.done")

rule run_target_prs_dbslmm_all_name:
  input: 
    expand("resources/data/target/{name}/prs/{name}_prs_dbslmm.done", name=target_list_df['name'])
  output:
    touch("resources/data/target/prs_dbslmm.done")

##
# PRScs
##

rule target_prs_prscs:
  input:
    "resources/data/target/{name}/ancestry/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/1kg/prs_score_files/prscs/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  output:
    touch("resources/data/target/{name}/prs/{population}/prscs/{gwas}/{name}_prs_{population}_prscs_{gwas}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/Scaled_polygenic_scorer_PRScs/Scaled_polygenic_scorer_PRScs.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_score resources/data/1kg/prs_score_files/prscs/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_scale resources/data/1kg/prs_score_files/prscs/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.{wildcards.population}.scale \
      --ref_freq_chr resources/data/1kg/freq_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.chr \
      --plink plink \
      --pheno_name {wildcards.gwas} \
      --output resources/data/target/{wildcards.name}/prs/{wildcards.population}/prscs/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.{wildcards.population}"

rule run_target_prs_prscs_all_gwas:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/prscs/{gwas}/{name}_prs_{population}_prscs_{gwas}.done", name=w.name, gwas=gwas_list_df['name'], population=w.population)
  output:
    touch("resources/data/target/{name}/prs/{population}/prscs/{name}_prs_{population}_prscs.done")

rule run_target_prs_prscs_all_pop:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/prscs/{name}_prs_{population}_prscs.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target/{name}/prs/{name}_prs_prscs.done")

rule run_target_prs_prscs_all_name:
  input: 
    expand("resources/data/target/{name}/prs/{name}_prs_prscs.done", name=target_list_df['name'])
  output:
    touch("resources/data/target/prs_prscs.done")

##
# lassosum
##

rule target_prs_lassosum:
  resources:
    cpus=5,
    mem_mb=10000
  input:
    "resources/data/target/{name}/ancestry/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/1kg/prs_score_files/lassosum/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  output:
    touch("resources/data/target/{name}/prs/{population}/lassosum/{gwas}/{name}_prs_{population}_lassosum_{gwas}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "for population in $(echo {params.population}); do \
     Rscript ../Scripts/Scaled_polygenic_scorer_lassosum/Scaled_polygenic_scorer_lassosum.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_score resources/data/1kg/prs_score_files/lassosum/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_scale resources/data/1kg/prs_score_files/lassosum/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.{wildcards.population}.scale \
      --ref_freq_chr resources/data/1kg/freq_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.chr \
      --plink plink \
      --pheno_name {wildcards.gwas} \
      --n_cores {resources.cpus} \
      --output resources/data/target/{wildcards.name}/prs/{wildcards.population}/lassosum/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.{wildcards.population}; \
     done"

rule run_target_prs_lassosum_all_gwas:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/lassosum/{gwas}/{name}_prs_{population}_lassosum_{gwas}.done", name=w.name, gwas=gwas_list_df['name'], population=w.population)
  output:
    touch("resources/data/target/{name}/prs/{population}/lassosum/{name}_prs_{population}_lassosum.done")

rule run_target_prs_lassosum_all_pop:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/lassosum/{name}_prs_{population}_lassosum.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target/{name}/prs/{name}_prs_lassosum.done")

rule run_target_prs_lassosum_all_name:
  input: 
    expand("resources/data/target/{name}/prs/{name}_prs_lassosum.done", name=target_list_df['name'])
  output:
    touch("resources/data/target/prs_lassosum.done")

##
# SBayesR
##

rule target_prs_sbayesr:
  input:
    "resources/data/target/{name}/ancestry/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/1kg/prs_score_files/sbayesr/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  output:
    touch("resources/data/target/{name}/prs/{population}/sbayesr/{gwas}/{name}_prs_{population}_sbayesr_{gwas}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "for population in $(echo {params.population}); do \
     Rscript ../Scripts/Scaled_polygenic_scorer_SBayesR/Scaled_polygenic_scorer_SBayesR.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_score resources/data/1kg/prs_score_files/sbayesr/{wildcards.gwas}/GWAS_sumstats_SBayesR.GW.snpRes \
      --ref_scale resources/data/1kg/prs_score_files/sbayesr/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.{wildcards.population}.scale \
      --ref_freq_chr resources/data/1kg/freq_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.chr \
      --plink plink \
      --pheno_name {wildcards.gwas} \
      --output resources/data/target/{wildcards.name}/prs/{wildcards.population}/sbayesr/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.{wildcards.population}; \
     done"

rule run_target_prs_sbayesr_all_gwas:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/sbayesr/{gwas}/{name}_prs_{population}_sbayesr_{gwas}.done", name=w.name, gwas=gwas_list_df['name'], population=w.population)
  output:
    touch("resources/data/target/{name}/prs/{population}/sbayesr/{name}_prs_{population}_sbayesr.done")

rule run_target_prs_sbayesr_all_pop:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/sbayesr/{name}_prs_{population}_sbayesr.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target/{name}/prs/{name}_prs_sbayesr.done")

rule run_target_prs_sbayesr_all_name:
  input: 
    expand("resources/data/target/{name}/prs/{name}_prs_sbayesr.done", name=target_list_df['name'])
  output:
    touch("resources/data/target/prs_sbayesr.done")

##
# LDpred2
##

rule target_prs_ldpred2:
  resources:
    cpus=5,
    mem_mb=10000
  input:
    "resources/data/target/{name}/ancestry/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/1kg/prs_score_files/ldpred2/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  output:
    touch("resources/data/target/{name}/prs/{population}/ldpred2/{gwas}/{name}_prs_{population}_ldpred2_{gwas}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "for population in $(echo {params.population}); do \
     Rscript ../Scripts/Scaled_polygenic_scorer_LDPred2/Scaled_polygenic_scorer_LDPred2.R \
      --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep resources/data/target/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_score resources/data/1kg/prs_score_files/ldpred2/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_scale resources/data/1kg/prs_score_files/ldpred2/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.{wildcards.population}.scale \
      --ref_freq_chr resources/data/1kg/freq_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}.chr \
      --plink plink \
      --pheno_name {wildcards.gwas} \
      --n_cores {resources.cpus} \
      --output resources/data/target/{wildcards.name}/prs/{wildcards.population}/ldpred2/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.{wildcards.population}; \
     done"
     
rule run_target_prs_ldpred2_all_gwas:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/ldpred2/{gwas}/{name}_prs_{population}_ldpred2_{gwas}.done", name=w.name, gwas=gwas_list_df['name'], population=w.population)
  output:
    touch("resources/data/target/{name}/prs/{population}/ldpred2/{name}_prs_{population}_ldpred2.done")

rule run_target_prs_ldpred2_all_pop:
  input: 
    lambda w: expand("resources/data/target/{name}/prs/{population}/ldpred2/{name}_prs_{population}_ldpred2.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target/{name}/prs/{name}_prs_ldpred2.done")

rule run_target_prs_ldpred2_all_name:
  input: 
    expand("resources/data/target/{name}/prs/{name}_prs_ldpred2.done", name=target_list_df['name'])
  output:
    touch("resources/data/target/prs_ldpred2.done")

##
# Calculate PRS using all methods
##

rule target_prs_all:
  input:
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_ptclump.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_dbslmm.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_prscs.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_sbayesr.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_lassosum.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_ldpred2.done", name=w.name)
  output:
    touch('resources/data/target/{name}/prs/target_prs_all.done')

rule run_target_prs_all:
  input: expand('resources/data/target/{name}/prs/target_prs_all.done', name=target_list_df['name'])

##########
# Create full reports
##########

## 
# Individual reports
##

rule create_individual_report:
  input:
    rules.install_ggchicklet.output,
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry/ancestry_all_pop.done",
    lambda w: expand("resources/data/target/{name}/projected_pc/{name}_target_pc.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_ptclump.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_dbslmm.done.done", name=w.name),
    rules.run_pseudovalidate_prs.input,
    rules.download_1kg_pop_codes.output
  output:
    touch('resources/data/target/{name}/{name}_indiv_report.done') 
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e \"rmarkdown::render(\'scripts/indiv_report_creator.Rmd\', \
     output_file = \'../resources/data/target/{wildcards.name}/{wildcards.name}_report.html\', \
     params = list(name = \'{wildcards.name}\'))\""

rule run_create_individual_report:
  input: expand('resources/data/target/{name}/{name}_indiv_report.done', name=target_list_df_23andMe['name'])

##
# Sample reports
##

rule create_sample_report:
  input:
    rules.install_ggchicklet.output,
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry/ancestry_all_pop.done",
    lambda w: expand("resources/data/target/{name}/projected_pc/{name}_target_pc.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_ptclump.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/{name}_prs_dbslmm.done", name=w.name),
    rules.download_1kg_pop_codes.output
  output:
    touch('resources/data/target/{name}/{name}_samp_report.done')
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e \"rmarkdown::render(\'scripts/samp_report_creator.Rmd\', \
     output_file = \'../resources/data/target/{wildcards.name}/{wildcards.name}_report.html\', \
     params = list(name = \'{wildcards.name}\'))\""

rule run_create_sample_report:
  input: expand('resources/data/target/{name}/{name}_samp_report.done', name=target_list_df_samp_imp['name'])

rule run_create_reports:
  input: 
    rules.run_create_individual_report.input,
    rules.run_create_sample_report.input
