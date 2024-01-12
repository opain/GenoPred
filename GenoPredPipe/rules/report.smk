# Create a rule that checks all defaults outputs given certain outputs are present
output_all_input = list()

if 'target_list' in config:
  output_all_input.append(expand("{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done", outdir=outdir, name=target_list_df['name']))
  if 'gwas_list' in config or 'score_list' in config:
    output_all_input.append(expand("{outdir}/resources/data/target_checks/{name}/run_target_pgs_all_method.done", outdir=outdir, name=target_list_df['name']))
    # For pseudovalidation if not already requested
    #output_all_input.append(rules.run_prep_pgs_lassosum.input)
else:
  if 'gwas_list' in config or 'score_list' in config:
    output_all_input.append(rules.pgs_methods_complete.input)

rule output_all:
  input: output_all_input
