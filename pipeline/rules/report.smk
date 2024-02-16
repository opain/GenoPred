output_all_input = list()
label_list = pd.Series(dtype=object)

if 'target_list' in config and config["target_list"] != 'NA':
  output_all_input.append(rules.target_pgs.input)

if 'gwas_list' in config and config["gwas_list"] != 'NA':
  output_all_input.append(rules.prep_pgs.input)
  output_all_input.append(rules.sumstat_prep.input)
  label_list = pd.concat([label_list, gwas_list_df['label']])

if 'score_list' in config and config["score_list"] != 'NA':
  output_all_input.append(rules.prep_pgs.input)
  label_list = pd.concat([label_list, score_list_df['label']])

# Identify temp directory
tmpdir = 'resources/tmp'

# Make a list of all values in sampling prevalence mean sd columns in gwas_list
gwas_dist = [gwas_list_df['sampling'], gwas_list_df['prevalence'], gwas_list_df['mean'], gwas_list_df['sd']]
gwas_dist = [item for sublist in gwas_dist for item in sublist]

#####
# Create a report for each target sample
#####

rule sample_report_i:
  input:
    output_all_input
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/sample_report.done")
  benchmark:
    f"{outdir}/reference/benchmarks/sample_report_i-{{name}}.txt"
  log:
    f"{outdir}/reference/logs/sample_report_i-{{name}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    outdir=config["outdir"],
    labels=' '.join(label_list),
    gwas_dist=gwas_dist,
    cwd=os.getcwd(),
    tempdir= lambda w: tmpdir + "/" + w.name + "/",
    config_file = config["config_file"],
    report_out= lambda w: outdir if outdir[0] == "/" else os.getcwd() + "/" + outdir
  shell:
    """
    {{
      mkdir -p {params.tempdir} ; \
      cp ../Scripts/pipeline_reports/samp_report_creator.Rmd {params.tempdir}/samp_report_creator.Rmd ; \
      mkdir -p {outdir}/{wildcards.name}/reports; \
      Rscript -e \"rmarkdown::render(\'{params.tempdir}/samp_report_creator.Rmd\', \
      output_file = \'{params.report_out}/{wildcards.name}/reports/{wildcards.name}-report.html\', \
      params = list(name = \'{wildcards.name}\', config = \'{params.config_file}\', cwd = \'{params.cwd}\'))\" ; \
      rm -r {params.tempdir}
    }} > {log} 2>&1
    """

rule sample_report:
  input: expand(f"{outdir}/reference/target_checks/{{name}}/sample_report.done", name=target_list_df_samp['name'])

#####
# Create individual-level reports for each target sample
#####

def id_munge(name):
  if config['testing'] != 'NA':
    val = config['testing'][-2:]
    val = str(val)
  else:
    val = str(22)

  checkpoint_output = checkpoints.ancestry_reporter.get(name=name).output[0]
  checkpoint_output = checkpoints.score_reporter.get(name=name).output[0]
  checkpoint_output = f"{outdir}/{name}/geno/{name}.ref.chr{val}.psam"
  fam_df = pd.read_table(checkpoint_output, delim_whitespace=True, usecols=[0,1], names=['FID', 'IID'], header=0)
  fam_df['id'] = fam_df.FID.apply(str) + '.' + fam_df.IID.apply(str)

  return fam_df['id'].tolist()

rule indiv_report_i:
  input:
    rules.install_ggchicklet.output,
    rules.prep_pgs_lassosum.input,
    output_all_input
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/indiv_report-{{id}}-report.done")
  benchmark:
    f"{outdir}/reference/benchmarks/indiv_report_i-{{name}}-{{id}}.txt"
  log:
    f"{outdir}/reference/logs/indiv_report_i-{{name}}-{{id}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    outdir=config["outdir"],
    labels=' '.join(label_list),
    gwas_dist=gwas_dist,
    cwd=os.getcwd(),
    tempdir= lambda w: tmpdir + "/" + w.name + "." + w.id + "/",
    config_file = config["config_file"],
    report_out= lambda w: outdir if outdir[0] == "/" else os.getcwd() + "/" + outdir
  shell:
    """
    {{
      mkdir -p {params.tempdir} ; \
      cp ../Scripts/pipeline_reports/indiv_report_creator.Rmd {params.tempdir}/indiv_report_creator.Rmd ; \
      mkdir -p {outdir}/{wildcards.name}/reports/individual; \
      Rscript -e \"rmarkdown::render(\'{params.tempdir}/indiv_report_creator.Rmd\', \
      output_file = \'{params.report_out}/{wildcards.name}/reports/individual/{wildcards.name}-{wildcards.id}-report.html\', \
      params = list(name = \'{wildcards.name}\', id = \'{wildcards.id}\', config = \'{params.config_file}\', cwd = \'{params.cwd}\'))\" ; \
      rm -r {params.tempdir}
    }} > {log} 2>&1
    """

rule indiv_report_all:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/indiv_report-{{id}}-report.done", name=w.name, id=id_munge(name="{}".format(w.name)))
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/indiv_report.done")

rule indiv_report:
  input: expand(f"{outdir}/reference/target_checks/{{name}}/indiv_report.done", name= target_list_df_indiv_report['name'])

#####
# Create a rule that checks all defaults outputs given certain outputs are present
#####

rule output_all:
  input:
    output_all_input,
    rules.sample_report.input,
    rules.indiv_report.input
