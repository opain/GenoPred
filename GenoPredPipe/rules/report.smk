##########
# Create full reports
##########

## 
# Individual reports
##

# Create individual level report for an individuals genotype dataset
rule create_individual_report:
  input:
    lambda w: expand("resources/data/target_checks/{name}/run_target_pc_all_pop.done", name=w.name),
    lambda w: expand("resources/data/target_checks/{name}/run_target_prs_ptclump_all_pop.done", name=w.name),
    lambda w: expand("resources/data/target_checks/{name}/run_target_prs_dbslmm_all_pop.done", name=w.name),
    rules.run_pseudovalidate_prs.input,
    rules.install_ggchicklet.output,
    "scripts/indiv_report_creator.Rmd"
  output:
    touch('resources/data/target_checks/{name}/create_individual_report.done') 
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0],
    report_output= lambda w: report_output_munge(w.name)
  shell:
    "mkdir -p {params.output}/{wildcards.name}/reports; Rscript -e \"rmarkdown::render(\'scripts/indiv_report_creator.Rmd\', \
     output_file = \'{params.report_output}/{wildcards.name}/reports/{wildcards.name}_report.html\', \
     params = list(name = \'{wildcards.name}\', output = \'{params.output}\'))\""

rule run_create_individual_report:
  input: expand('resources/data/target_checks/{name}/create_individual_report.done', name=target_list_df_23andMe['name'])

# Create individual level reports for all individuals in a sample
rule create_individual_report_for_sample:
  input:
    lambda w: expand("resources/data/target_checks/{name}/run_target_pc_all_pop.done", name=w.name),
    lambda w: expand("resources/data/target_checks/{name}/run_target_prs_ptclump_all_pop.done", name=w.name),
    lambda w: expand("resources/data/target_checks/{name}/run_target_prs_dbslmm_all_pop.done", name=w.name),
    rules.run_pseudovalidate_prs.input,
    rules.install_ggchicklet.output,
    "scripts/indiv_report_creator_for_sample.Rmd"
  output:
    touch('resources/data/target_checks/{name}/create_individual_report_for_sample_{id}.done') 
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0],
    report_output= lambda w: report_output_munge(w.name)
  shell:
    "mkdir -p {params.output}/{wildcards.name}/reports; Rscript -e \"rmarkdown::render(\'scripts/indiv_report_creator_for_sample.Rmd\', \
     output_file = \'{params.report_output}/{wildcards.name}/reports/{wildcards.name}_{wildcards.id}_report.html\', \
     params = list(name = \'{wildcards.name}\',id = \'{wildcards.id}\', output = \'{params.output}\'))\""

rule run_create_individual_report_for_sample_all_id:
  input:
    lambda w: expand('resources/data/target_checks/{name}/create_individual_report_for_sample_{id}.done', name=w.name, id=id_munge("{}".format(w.name)))
  output:
    touch('resources/data/target_checks/{name}/run_create_individual_report_for_sample_all_id.done')

rule run_create_individual_report_for_sample_all_indiv:
  input: expand('resources/data/target_checks/{name}/run_create_individual_report_for_sample_all_id.done', name=target_list_df_samp_imp_indiv_report['name'])

##
# Sample reports
##

rule create_sample_report:
  input:
    lambda w: expand("resources/data/target_checks/{name}/run_target_pc_all_pop.done", name=w.name),
    lambda w: expand("resources/data/target_checks/{name}/run_target_prs_ptclump_all_pop.done", name=w.name),
    lambda w: expand("resources/data/target_checks/{name}/run_target_prs_all_method.done", name=w.name),
    "scripts/samp_report_creator.Rmd"
  output:
    touch('resources/data/target_checks/{name}/create_sample_report.done')
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0],
    report_output= lambda w: report_output_munge(w.name)
  shell:
    "mkdir -p {params.output}/{wildcards.name}/reports; Rscript -e \"rmarkdown::render(\'scripts/samp_report_creator.Rmd\', \
     output_file = \'{params.report_output}/{wildcards.name}/reports/{wildcards.name}_report.html\', \
     params = list(name = \'{wildcards.name}\', output = \'{params.output}\'))\""

rule run_create_sample_report:
  input: expand('resources/data/target_checks/{name}/create_sample_report.done', name=target_list_df_samp_imp['name'])

rule run_create_reports:
  input: 
    rules.run_create_individual_report.input,
    rules.run_create_individual_report_for_sample_all_indiv.input,
    rules.run_create_sample_report.input

##########
# Create ancestry-only reports
##########

def report_output_munge(x):
    output = target_list_df.loc[target_list_df['name'] == x, 'output'].iloc[0]
    report_output=output if output[0] == "/" else "../" + output
    return report_output

## 
# Individual ancestry reports
##

rule create_individual_ancestry_report:
  input:
    "resources/data/target_checks/{name}/ancestry_reporter.done",
    "scripts/indiv_ancestry_report_creator.Rmd"
  output:
    touch('resources/data/target_checks/{name}/indiv_ancestry_report.done') 
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0],
    report_output= lambda w: report_output_munge(w.name)
  shell:
    "mkdir -p {params.output}/{wildcards.name}/reports; Rscript -e \"rmarkdown::render(\'scripts/indiv_ancestry_report_creator.Rmd\', \
     output_file = \'{params.report_output}/{wildcards.name}/reports/{wildcards.name}_ancestry_report.html\', \
     params = list(name = \'{wildcards.name}\', output = \'{params.output}\'))\""

rule run_create_individual_ancestry_report:
  input: 
    expand('resources/data/target_checks/{name}/indiv_ancestry_report.done', name=target_list_df_23andMe['name'])

# Create individual level reports for all individuals in a sample
def id_munge(x):
    checkpoint_output = checkpoints.ancestry_reporter.get(name=x).output[0]
    checkpoint_output = target_list_df.loc[target_list_df['name'] == x, 'output'].iloc[0] + "/" + x + "/" + x + ".1KGphase3.hm3.chr22.fam"
    fam_df = pd.read_table(checkpoint_output, delim_whitespace=True, usecols=[0,1], names=['FID', 'IID'], header=None)
    fam_df['id'] = fam_df.FID.apply(str) + '.' + fam_df.IID.apply(str)
    return fam_df['id'].tolist()

rule create_individual_ancestry_report_for_sample:
  input:
    "resources/data/target_checks/{name}/ancestry_reporter.done",
    "scripts/indiv_ancestry_report_creator_for_sample.Rmd"
  output:
    touch('resources/data/target_checks/{name}/create_individual_ancestry_report_for_sample_{id}.done') 
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0],
    report_output= lambda w: report_output_munge(w.name)
  shell:
    "mkdir -p {params.output}/{wildcards.name}/reports; Rscript -e \"rmarkdown::render(\'scripts/indiv_ancestry_report_creator_for_sample.Rmd\', \
     output_file = \'{params.report_output}/{wildcards.name}/reports/{wildcards.name}_{wildcards.id}_ancestry_report.html\', \
     params = list(name = \'{wildcards.name}\',id = \'{wildcards.id}\', output = \'{params.output}\'))\""

rule run_create_individual_ancestry_report_for_sample_all_id:
  input:
    lambda w: expand('resources/data/target_checks/{name}/create_individual_ancestry_report_for_sample_{id}.done', name=w.name, id=id_munge("{}".format(w.name)))
  output:
    touch('resources/data/target_checks/{name}/run_create_individual_ancestry_report_for_sample_all_id.done')

rule run_create_individual_ancestry_report_for_sample_all_indiv:
  input: expand('resources/data/target_checks/{name}/run_create_individual_ancestry_report_for_sample_all_id.done', name=target_list_df_samp_imp_indiv_report['name'])

##
# Sample ancestry reports
##

rule create_sample_ancestry_report:
  input:
    "resources/data/target_checks/{name}/ancestry_reporter.done",
    "scripts/samp_ancestry_report_creator.Rmd"
  output:
    touch('resources/data/target_checks/{name}/samp_ancestry_report.done')
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0],
    report_output= lambda w: report_output_munge(w.name)
  shell:
    "mkdir -p {params.output}/{wildcards.name}/reports; Rscript -e \"rmarkdown::render(\'scripts/samp_ancestry_report_creator.Rmd\', \
    output_file = \'{params.report_output}/{wildcards.name}/reports/{wildcards.name}_ancestry_report.html\', \
     params = list(name = \'{wildcards.name}\', output = \'{params.output}\'))\""

rule run_create_sample_ancestry_report:
  input: 
    expand('resources/data/target_checks/{name}/samp_ancestry_report.done', name=target_list_df_samp_imp['name'])

rule run_create_ancestry_reports:
  input: 
    rules.run_create_individual_ancestry_report.input,
    rules.run_create_individual_ancestry_report_for_sample_all_indiv.input,
    rules.run_create_sample_ancestry_report.input
