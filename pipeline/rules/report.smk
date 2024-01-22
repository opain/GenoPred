output_all_input = list()

if 'target_list' in config:
  output_all_input.append(rules.target_pgs.input)
else:
  if 'gwas_list' in config or 'score_list' in config:
    output_all_input.append(rules.prep_pgs.input)

#####
# Create a report for each target sample
#####

rule sample_report_i:
  input:
    output_all_input
  output:
    touch('{outdir}/reference/target_checks/{name}/sample_report.done')
  conda:
    "../envs/analysis.yaml"
  params:
    config_file = config["config_file"],
    report_out= lambda w: outdir if outdir[0] == "/" else "../../pipeline/" + outdir
  shell:
    "mkdir -p {outdir}/{wildcards.name}/reports; \
     Rscript -e \"rmarkdown::render(\'../Scripts/pipeline_reports/samp_report_creator.Rmd\', \
     output_file = \'{params.report_out}/{wildcards.name}/reports/{wildcards.name}-report.html\', \
     params = list(name = \'{wildcards.name}\', config = \'{params.config_file}\'))\""

rule sample_report:
  input: expand('{outdir}/reference/target_checks/{name}/sample_report.done', name=target_list_df_samp['name'], outdir=outdir)

#####
# Create individual-level reports for each target sample
#####

def id_munge(name, outdir):
  if config['testing'] != 'NA':
    val = config['testing'][-2:]
    val = str(val)
  else:
    val = str(22)

  checkpoint_output = checkpoints.ancestry_reporter.get(name=name, outdir=outdir).output[0]
  checkpoint_output = outdir + "/" + name + "/geno/" + name + ".ref.chr" + val + ".fam"
  fam_df = pd.read_table(checkpoint_output, delim_whitespace=True, usecols=[0,1], names=['FID', 'IID'], header=None)
  fam_df['id'] = fam_df.FID.apply(str) + '.' + fam_df.IID.apply(str)

  return fam_df['id'].tolist()

rule indiv_report_i:
  input:
    rules.install_ggchicklet.output,
    rules.prep_pgs_lassosum.input,
    output_all_input
  output:
    touch('{outdir}/reference/target_checks/{name}/indiv_report-{id}-report.done')
  conda:
    "../envs/analysis.yaml"
  params:
    config_file = config["config_file"],
    report_out= lambda w: outdir if outdir[0] == "/" else "../../pipeline/" + outdir
  shell:
    "mkdir -p {outdir}/{wildcards.name}/reports; \
     Rscript -e \"rmarkdown::render(\'../Scripts/pipeline_reports/indiv_report_creator.Rmd\', \
     output_file = \'{params.report_out}/{wildcards.name}/reports/{wildcards.name}-{wildcards.id}-report.html\', \
     params = list(name = \'{wildcards.name}\', id = \'{wildcards.id}\', config = \'{params.config_file}\'))\""

rule indiv_report_all_id:
  input:
    lambda w: expand('{outdir}/reference/target_checks/{name}/indiv_report-{id}-report.done', name=w.name, id=id_munge(name="{}".format(w.name), outdir=w.outdir), outdir=w.outdir)
  output:
    touch('{outdir}/reference/target_checks/{name}/indiv_report_all_id.done')

rule indiv_report:
  input: expand('{outdir}/reference/target_checks/{name}/indiv_report_all_id.done', name= target_list_df_indiv_report['name'], outdir=outdir)

#####
# Create a rule that checks all defaults outputs given certain outputs are present
#####

rule output_all:
  input:
    output_all_input,
    rules.sample_report.input,
    rules.indiv_report.input
