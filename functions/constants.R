# Set range of chromosome numbers to use by default
CHROMS <- 1:22

# Coordinates of high ld regions in build GRCh37
long_ld_coord<-do.call(rbind, list(
  data.frame(CHR = 1, P0 = 48e6, P1 = 52e6),
  data.frame(CHR = 2, P0 = 86e6, P1 = 100.5e6),
  data.frame(CHR = 2, P0 = 134.5e6, P1 = 138e6),
  data.frame(CHR = 2, P0 = 183e6, P1 = 190e6),
  data.frame(CHR = 3, P0 = 47.5e6, P1 = 50e6),
  data.frame(CHR = 3, P0 = 83.5e6, P1 = 87e6),
  data.frame(CHR = 3, P0 = 89e6, P1 = 97.5e6),
  data.frame(CHR = 5, P0 = 44.5e6, P1 = 50.5e6),
  data.frame(CHR = 5, P0 = 98e6, P1 = 100.5e6),
  data.frame(CHR = 5, P0 = 129e6, P1 = 132e6),
  data.frame(CHR = 5, P0 = 135.5e6, P1 = 138.5e6),
  data.frame(CHR = 6, P0 = 25.5e6, P1 = 33.5e6),
  data.frame(CHR = 6, P0 = 57e6, P1 = 64e6),
  data.frame(CHR = 6, P0 = 140e6, P1 = 142.5e6),
  data.frame(CHR = 7, P0 = 55e6, P1 = 66e6),
  data.frame(CHR = 8, P0 = 8e6, P1 = 12e6),
  data.frame(CHR = 8, P0 = 43e6, P1 = 50e6),
  data.frame(CHR = 8, P0 = 112e6, P1 = 115e6),
  data.frame(CHR = 10, P0 = 37e6, P1 = 43e6),
  data.frame(CHR = 11, P0 = 46e6, P1 = 57e6),
  data.frame(CHR = 11, P0 = 87.5e6, P1 = 90.5e6),
  data.frame(CHR = 12, P0 = 33e6, P1 = 40e6),
  data.frame(CHR = 12, P0 = 109.5e6, P1 = 112e6),
  data.frame(CHR = 20, P0 = 32e6, P1 = 34.5e6)
))

# Make a data.frame giving labels to the 1KG reference populations
ref_pop <- data.frame(
  pop = c('AFR','AMR','EAS','EUR','CSA','MID'),
  label = c('African','American','East Asian','European','Central and South Asian','Middle Eastern')
)

# Make a data.frame giving labels to the 1KG reference populations
pgs_method_labels <- data.frame(
  method = c('ptclump','dbslmm','ldpred2','sbayesr','lassosum','prscs','megaprs','external'),
  label = c('pT+clump','DBSLMM','LDpred2','SBayesR','lassosum','PRS-CS','MegaPRS','External')
)
pgs_method_labels[order(pgs_method_labels$method),]

# Make vector indicating pgs_methods that can be applied to non-european GWAS
pgs_methods_noneur <- c('ptclump','lassosum','megaprs','prscs','dbslmm')

# Make vector indicating pgs_methods that are to be applied to gwas_groups
pgs_group_methods <- c('prscsx')

