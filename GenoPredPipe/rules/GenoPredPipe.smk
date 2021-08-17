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
    "mkdir -p resources/software/qctool2/; wget -O resources/software/qctool2/qctool2.tgz https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz; tar -zxvf resources/software/qctool2/qctool2.tgz -C resources/software/qctool2/ --strip-components=1; rm resources/software/qctool2/qctool2.tgz"

# Download impute2_data
rule download_impute2_data:
  output:
    directory("resources/data/impute2/1000GP_Phase3")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/data/impute2/; wget -O resources/data/impute2/1000GP_Phase3.tgz https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz; tar -zxvf resources/data/impute2/1000GP_Phase3.tgz -C resources/data/impute2/; rm resources/data/impute2/1000GP_Phase3.tgz; wget -O resources/data/impute2/1000GP_Phase3_chrX.tgz https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz; tar -zxvf resources/data/impute2/1000GP_Phase3_chrX.tgz -C resources/data/impute2/1000GP_Phase3/; rm resources/data/impute2/1000GP_Phase3_chrX.tgz"

# Download PLINK. DBSLMM requires the binary to be specified, which is challenging with conda environments.
rule download_plink:
  output:
    "resources/software/plink/plink"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/software/plink; wget -O resources/software/plink/plink_linux_x86_64_20210606.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip; unzip resources/software/plink/plink_linux_x86_64_20210606.zip -d resources/software/plink; rm resources/software/plink/plink_linux_x86_64_20210606.zip"

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
    "mkdir -p resources/data/ldsc_ref; wget -O resources/data/ldsc_ref/eur_w_ld_chr.tar.bz2 https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2; tar -jxvf resources/data/ldsc_ref/eur_w_ld_chr.tar.bz2 -C resources/data/ldsc_ref/; rm resources/data/ldsc_ref/eur_w_ld_chr.tar.bz2"

# Download hapmap3 snplist
rule download_hm3_snplist:
  output:
    "resources/data/hm3_snplist/w_hm3.snplist"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "mkdir -p resources/data/hm3_snplist; wget -O resources/data/hm3_snplist/w_hm3.snplist.bz2 https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2; bunzip2 resources/data/hm3_snplist/w_hm3.snplist.bz2"

# Download DBSLMM
# Specifying old commit as developer has deleted dbslmm binary (accidentally?)
rule download_dbslmm:
  output:
    directory("resources/software/dbslmm/")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "git clone git@github.com:biostat0903/DBSLMM.git {output}; cd {output}; git reset --hard aa6e7ad5b8a7d3b6905556a4007c4a1fa2925b7d; chmod a+x software/dbslmm"

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

# install ggchicklet
rule install_ggchicklet:
  output:
    touch("resources/software/install_ggchicklet.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e 'remotes::install_github(\"hrbrmstr/ggchicklet\")'"

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
    "ls resources/data/1kg/1KGPhase3.w_hm3.chr*.bed | sed -e 's/\.bed//g' > resources/data/1kg/merge_list.txt; plink --merge-list resources/data/1kg/merge_list.txt --make-bed --out resources/data/1kg/1KGPhase3.w_hm3.GW; rm resources/data/1kg/merge_list.txt"

# Download 1kg population code definitions
rule download_1kg_pop_codes:
  output:
    "resources/data/1kg/1kg_pop_codes.tsv",
    "resources/data/1kg/1kg_super_pop_codes.tsv"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "wget -O resources/data/1kg/1kg_pop_codes.tsv http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv; wget -O resources/data/1kg/1kg_super_pop_codes.tsv http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.superpopulations.tsv"

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
    "Rscript ../Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr --ref_keep resources/data/1kg/keep_files/{wildcards.population}_samples.keep --plink plink --plink2 plink2 --n_pcs 100 --output resources/data/1kg/pc_score_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}"

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
    "Rscript ../Scripts/sumstat_cleaner/sumstat_cleaner.R --sumstats {params.path} --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr --ref_freq_chr resources/data/1kg/freq_files/{params.population}/1KGPhase3.w_hm3.{params.population}.chr --output resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned"
    
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
    "Rscript ../Scripts/polygenic_score_file_creator/polygenic_score_file_creator.R --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz --plink plink --output resources/data/1kg/prs_score_files/pt_clump/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} --ref_pop_scale resources/data/1kg/super_pop_keep.list"

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
    "Rscript ../Scripts/polygenic_score_file_creator_DBSLMM/polygenic_score_file_creator_DBSLMM.R --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz --plink resources/software/plink/plink --ld_blocks resources/data/ld_blocks/{params.population} --rscript Rscript --dbslmm resources/software/dbslmm/software --munge_sumstats resources/software/ldsc/munge_sumstats.py --ldsc resources/software/ldsc/ldsc.py --ldsc_ref resources/data/ldsc_ref/eur_w_ld_chr --hm3_snplist resources/data/hm3_snplist/w_hm3.snplist --sample_prev {params.sampling} --pop_prev {params.prevalence} --output resources/data/1kg/prs_score_files/dbslmm/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} --ref_pop_scale resources/data/1kg/super_pop_keep.list"

rule run_prs_scoring_dbslmm:
  input: expand("resources/data/1kg/prs_score_files/dbslmm/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

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
    "Rscript ../Scripts/lassosum_pseudovalidate/lassosum_pseudovalidate.R --ref_plink_gw resources/data/1kg/1KGPhase3.w_hm3.GW --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz --prune_mhc T --output resources/data/1kg/prs_pseudoval/{wildcards.gwas}/lassosum_pseudo_{wildcards.gwas} --plink plink --n_cores 1"

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
target_list_df_samp_imp = target_list_df.loc[target_list_df['type'] == 'samp_imp']
target_list_df['pre_harm_path'] = target_list_df['path']
target_list_df.loc[target_list_df['type'] == '23andMe', 'pre_harm_path'] = "resources/data/target/" + target_list_df.loc[target_list_df['type'] == '23andMe', 'name'] + "/" + target_list_df.loc[target_list_df['type'] == '23andMe', 'name'] + ".1KGphase3.chr"

####
# Format and impute input data
####
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
    "Rscript ../Scripts/23andMe_imputer/23andMe_imputer.R --FID {params.name} --IID {params.name} --geno {params.path} --plink plink --out resources/data/target/{params.name}/{params.name} --ref resources/data/impute2/1000GP_Phase3 --shapeit shapeit --impute2 impute2 --n_core {resources.cpus} --qctool resources/software/qctool2/qctool"

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
  input: expand("resources/data/target/{name}/touch_imp.done", name=target_list_df_samp_imp['name'])

####
# Harmonise with reference
####

rule harmonise_target_with_ref:
  input:
    rules.prep_1kg.output,
    lambda w: "resources/data/target/" + w.name + "/" + w.name + ".sample" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == '23andMe' else "resources/data/target/{name}/touch_imp.done"
  output:
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    path=lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'pre_harm_path'].iloc[0]
  shell:
    "Rscript ../Scripts/hm3_harmoniser/hm3_harmoniser.R --target {params.path} --ref resources/data/1kg/1KGPhase3.w_hm3.chr --plink plink --out resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr"

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
    "resources/data/target/{name}/ancestry/{name}.Ancestry.SAS.eigenvec"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr --n_pcs 6 --plink plink --plink2 plink2 --output resources/data/target/{wildcards.name}/ancestry/{wildcards.name}.Ancestry --ref_pop_scale resources/data/1kg/super_pop_keep.list --pop_data resources/data/1kg/ref_pop_dat.txt --prob_thresh 0.5"

rule run_target_super_population:
  input: expand("resources/data/target/{name}/ancestry/{name}.Ancestry.SAS.eigenvec", name=target_list_df['name'])

# Create a file listing target samples and super population assignments
rule ancestry_reporter:
  input:
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry/{name}.Ancestry.SAS.eigenvec"
  output:
    "resources/data/target/{name}/ancestry_report.txt"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript scripts/ancestry_reporter.R {wildcards.name}"

rule run_ancestry_reporter:
  input: expand("resources/data/target/{name}/ancestry_report.txt", name=target_list_df['name'])

####
# Identify population within super population
####

from pathlib import Path
  
def ancestry_munge(x):
  my_file = Path("".join(["resources/data/target/",x,"/ancestry_report.txt"]))
  if my_file.is_file():
    ancestry_report_df = pd.read_table("".join(["resources/data/target/",x,"/ancestry_report.txt"]), sep=' ')
    return ancestry_report_df['population'].tolist()

rule target_population:
  input:
    rules.prep_1kg.output,
    "resources/data/target/{name}/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed"
  output:
    touch("resources/data/target/{name}/ancestry_pop.done")
  params:
    population=lambda w: ancestry_munge("{}".format(w.name))
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "for population in $(echo {params.population}); do Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr --target_keep resources/data/target/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.${{population}}.keep --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr --n_pcs 6 --plink plink --plink2 plink2 --output resources/data/target/{wildcards.name}/ancestry_${{population}}/{wildcards.name}.Ancestry --ref_pop_scale resources/data/1kg/pop_keep_for_${{population}}.list --pop_data resources/data/1kg/ref_pop_dat_for_${{population}}.txt --prob_thresh 0.5; done"

rule run_target_population:
  input: expand("resources/data/target/{name}/ancestry_pop.done", name=target_list_df['name'])

##########
# Target sample scoring
##########

####
# PC scoring
####

rule target_pc:
  input:
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry_report.txt",
    rules.run_super_pop_pc_scoring.input
  output:
    touch("resources/data/target/{name}/projected_pc/{name}_target_pc.done")
  params:
    population=lambda w: ancestry_munge("{}".format(w.name))
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "for population in $(echo {params.population}); do Rscript ../Scripts/scaled_ancestry_scorer/scaled_ancestry_scorer.R --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr --target_keep resources/data/target/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.${{population}}.keep --ref_freq_chr resources/data/1kg/freq_files/${{population}}/1KGPhase3.w_hm3.${{population}}.chr --ref_score resources/data/1kg/pc_score_files/${{population}}/1KGPhase3.w_hm3.${{population}}.eigenvec.var --ref_eigenvec resources/data/1kg/pc_score_files/${{population}}/1KGPhase3.w_hm3.${{population}}.eigenvec --ref_scale resources/data/1kg/pc_score_files/${{population}}/1KGPhase3.w_hm3.${{population}}.scale --plink2 plink2 --output resources/data/target/{wildcards.name}/projected_pc/${{population}}/{wildcards.name} --pop_data resources/data/1kg/ref_pop_dat.txt; done"

rule run_target_pc:
  input: expand("resources/data/target/{name}/projected_pc/{name}_target_pc.done", name=target_list_df['name'])

####
# Polygenic scoring
####

##
# pt+clump
##

rule target_prs_ptclump:
  input:
    "resources/data/target/{name}/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/1kg/prs_score_files/pt_clump/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  output:
    touch("resources/data/target/{name}/prs/target_prs_ptclump_{gwas}.done")
  params:
    population=lambda w: ancestry_munge("{}".format(w.name))
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "for population in $(echo {params.population}); do Rscript ../Scripts/Scaled_polygenic_scorer/Scaled_polygenic_scorer.R --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr --target_keep resources/data/target/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.${{population}}.keep --ref_score resources/data/1kg/prs_score_files/pt_clump/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} --ref_scale resources/data/1kg/prs_score_files/pt_clump/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.${{population}}.scale --ref_freq_chr resources/data/1kg/freq_files/${{population}}/1KGPhase3.w_hm3.${{population}}.chr --plink plink --pheno_name {wildcards.gwas} --output resources/data/target/{wildcards.name}/prs/${{population}}/pt_clump/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.${{population}}; done"

rule run_target_prs_ptclump:
  input: expand("resources/data/target/{name}/prs/target_prs_ptclump_{gwas}.done", name=target_list_df['name'], gwas=gwas_list_df['name'])

##
# DBSLMM
##

rule target_prs_dbslmm:
  input:
    "resources/data/target/{name}/ancestry_report.txt",
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/1kg/prs_score_files/dbslmm/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  output:
    touch("resources/data/target/{name}/prs/target_prs_dbslmm_{gwas}.done")
  params:
    population=lambda w: ancestry_munge("{}".format(w.name))
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "for population in $(echo {params.population}); do Rscript ../Scripts/Scaled_polygenic_scorer_DBSLMM/Scaled_polygenic_scorer_DBSLMM.R --target_plink_chr resources/data/target/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr --target_keep resources/data/target/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.${{population}}.keep --ref_score resources/data/1kg/prs_score_files/dbslmm/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.dbslmm.GW.txt --ref_scale resources/data/1kg/prs_score_files/dbslmm/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas}.${{population}}.scale --ref_freq_chr resources/data/1kg/freq_files/${{population}}/1KGPhase3.w_hm3.${{population}}.chr --plink plink --pheno_name {wildcards.gwas} --output resources/data/target/{wildcards.name}/prs/${{population}}/dbslmm/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.${{population}}; done"

rule run_target_prs_dbslmm:
  input: expand("resources/data/target/{name}/prs/target_prs_dbslmm_{gwas}.done", name=target_list_df_23andMe['name'], gwas=gwas_list_df['name'])

##########
# Create reports
##########

## 
# Individual reports
##

rule create_individual_report:
  input:
    rules.install_ggchicklet.output,
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry_pop.done",
    lambda w: expand("resources/data/target/{name}/projected_pc/{name}_target_pc.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/target_prs_ptclump_{gwas}.done", name=w.name, gwas=gwas_list_df['name']),
    lambda w: expand("resources/data/target/{name}/prs/target_prs_dbslmm_{gwas}.done", name=w.name, gwas=gwas_list_df['name']),
    rules.run_pseudovalidate_prs.input,
    rules.download_1kg_pop_codes.output
  output:
    touch('resources/data/target/{name}/{name}_indiv_report.done') 
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e \"rmarkdown::render(\'scripts/indiv_report_creator.Rmd\', output_file = \'../resources/data/target/{wildcards.name}/{wildcards.name}_report.html\', params = list(name = \'{wildcards.name}\'))\""

rule run_create_individual_report:
  input: expand('resources/data/target/{name}/{name}_indiv_report.done', name=target_list_df_23andMe['name'])

##
# Sample reports
##

rule create_sample_report:
  input:
    rules.install_ggchicklet.output,
    "resources/data/target/{name}/{name}.1KGphase3.hm3.chr22.bed",
    "resources/data/target/{name}/ancestry_pop.done",
    lambda w: expand("resources/data/target/{name}/projected_pc/{name}_target_pc.done", name=w.name),
    lambda w: expand("resources/data/target/{name}/prs/target_prs_ptclump_{gwas}.done", name=w.name, gwas=gwas_list_df['name']),
    lambda w: expand("resources/data/target/{name}/prs/target_prs_dbslmm_{gwas}.done", name=w.name, gwas=gwas_list_df['name']),
    rules.run_pseudovalidate_prs.input,
    rules.download_1kg_pop_codes.output
  output:
    touch('resources/data/target/{name}/{name}_samp_report.done')
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript -e \"rmarkdown::render(\'scripts/samp_report_creator.Rmd\', output_file = \'../resources/data/target/{wildcards.name}/{wildcards.name}_report.html\', params = list(name = \'{wildcards.name}\'))\""

rule run_create_sample_report:
  input: expand('resources/data/target/{name}/{name}_samp_report.done', name=target_list_df_samp_imp['name'])

rule run_create_reports:
  input: 
    rules.run_create_individual_report.input,
    rules.run_create_sample_report.input
