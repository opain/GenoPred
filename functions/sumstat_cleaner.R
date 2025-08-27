#' Interpret the header of GWAS summary statistics
#'
#' This function interprets the header of GWAS (Genome-Wide Association Studies) summary statistics. It reads the header, matches column names to a predefined dictionary, and provides a data frame indicating the original and interpreted column names, whether each column is kept or ignored, the reason for ignoring any columns, and a description of each interpreted column.
#'
#' @param sub_ss A data frame or matrix containing GWAS summary statistics, where the column names represent different data fields.
#'
#' @return A data frame with columns 'Original' (original column names from `sub_ss`), 'Interpreted' (standardized column names based on a predefined dictionary), 'Keep' (logical indicating whether the column is kept), 'Reason' (reason for ignoring a column, if applicable), and 'Description' (description of the interpreted column).
#' @export
#'
#' @import data.table
#'
#' @examples
#' header_info <- head_interp(raw_sumstats_1)
#' print(header_info)
head_interp<-function(sub_ss){
    message("head_interp() function started")  
  # Convert data.frame to data.table if necessary
  if (is.data.frame(sub_ss)) {
    sub_ss <- data.table::as.data.table(sub_ss)
  }

  # Read in the header and interpret column names
  sub_header<-names(sub_ss)

  # Remove columns that are all NA
  sub_ss_comp<-sub_ss[,apply(sub_ss, 2, function(x) !(all(is.na(x)))), with=F]
  sub_header_comp<-names(sub_ss_comp)

  int_header <- sub_header_comp
  for(i in names(ss_head_dict)){
    int_header[toupper(int_header) %in% ss_head_dict[[i]]] <- i
  }
  int_header[!(toupper(int_header) %in% unlist(ss_head_dict))]<-NA

  # Show original and interpreted header
  header_interp <- data.frame(Original = sub_header_comp,
                              Interpreted = int_header)

  # Show columns that are ignored due to be irrelevant, duplicated, or all NA
  header_interp$Keep<-!(
    !(int_header %in% names(ss_head_dict)) |
      duplicated(int_header)
  )

  # Insert reason it was ignored
  header_interp$Reason<-NA
  header_interp$Reason[!(int_header %in% names(ss_head_dict))]<-'Not recognised'
  header_interp$Reason[duplicated(int_header)]<-'Duplicated'

  # Show columns ignored due to missingness
  for(i in sub_header[!(sub_header %in% header_interp$Original)]){
    header_interp<-rbind(header_interp, data.frame(Original = i,
                                                   Interpreted = NA,
                                                   Keep=F,
                                                   Reason = 'First 1000 rows NA'))
  }

  header_interp[is.na(header_interp)]<-'NA'
  header_interp$Keep<-as.character(header_interp$Keep)

  # Insert description of each column after interpretation
  header_labels<-data.frame(Interpreted=c('SNP','CHR','BP','A1','A2','P','OR','BETA','Z','SE','N','N_CAS','N_CON','NEF','FREQ','FRQ_A','FREQ_U','INFO'),
                            Description=c("RSID for variant",
                                          "Chromosome number",
                                          "Base pair position",
                                          "Allele 1 (effect allele)",
                                          "Allele 2",
                                          "P-value of association",
                                          "Odds ratio effect size",
                                          "BETA effect size",
                                          "Z-score",
                                          "Standard error of log(OR) or BETA",
                                          "Total sample size",
                                          "Number of cases",
                                          "Number of controls",
                                          "Effective sample size",
                                          "Allele frequency",
                                          "Allele frequency in cases",
                                          "Allele frequency in controls",
                                          "Imputation quality"))

  header_interp<-merge(header_interp, header_labels, by='Interpreted', all.x=T)
  header_interp$Keep<-factor(header_interp$Keep, levels=c('TRUE','FALSE'))
  header_interp<-header_interp[order(header_interp$Keep),]
  header_interp<-header_interp[,c('Original','Interpreted','Keep','Reason','Description')]

  return(header_interp)
}

#' Check GWAS sumstat header is valid and update based on interpreted meaning
#'
#' This function checks if the header of a GWAS (Genome-Wide Association Studies) summary statistics dataset is valid and updates it based on interpreted meanings. It first interprets the header using `head_interp`, then logs the interpretation, checks for the presence of required columns, logs any errors or warnings, and finally updates the summary statistics dataset with the interpreted header.
#'
#' @param sumstats A data frame containing GWAS summary statistics where the column names are expected to be in the header.
#' @param log_file An optional string specifying the path and name of the log file where messages about the header interpretation and validation will be recorded. If NULL, messages are printed to the console.
#'
#' @return A data frame of GWAS summary statistics with updated and validated column names based on the interpreted header.
#' @export
#'
#' @import data.table
#'
#' @examples
#' formatted_ss_data <- format_header(raw_sumstats_1)
#' print(formatted_ss_data)
format_header<-function(sumstats, log_file = NULL){

  # Convert data.frame to data.table if necessary
  if (is.data.frame(sumstats)) {
    sumstats <- data.table::as.data.table(sumstats)
  }

  header_interpretation<-head_interp(sumstats[1:1000,])
  header_interpretation<-header_interpretation[match(names(sumstats), header_interpretation$Original),]

  log_add(log_file = log_file, message = '---------------')

  if(is.null(log_file)){
    print.data.frame(header_interpretation, row.names=F, quote=F, right=F)
  } else {
    sink(file = log_file, append = T)
    print.data.frame(header_interpretation, row.names=F, quote=F, right=F)
    sink()
  }

  log_add(log_file = log_file, message = '---------------')

  # Throw an error if required columns are not available
  error<-F
  messages <- NULL
  if(!all(c('A1','A2') %in% header_interpretation$Interpreted)){
    messages <- c(messages, "Error: Missing A1 and A2 information.")
    error<-T
  }
  if(!any(c('OR','BETA','Z') %in% header_interpretation$Interpreted)){
    messages <- c(messages, "Error: Effect size (BETA, OR, Z) must be present.")
    error<-T
  }
  if(!('SNP' %in% header_interpretation$Interpreted) & !(all(c('CHR','BP') %in% header_interpretation$Interpreted))){
    messages <- c(messages, "Error: Either SNP, or CHR and BP, must be present.")
    error<-T
  }
  if(any(c('N_CAS','N_CON') %in% header_interpretation$Interpreted) & sum(c('N_CAS','N_CON') %in% header_interpretation$Interpreted) == 1 & !('N' %in% names(sumstats))){
    messages <- c(messages, "Error: Either N, or both N_CAS and N_CON must be present.")
    error<-T
  }
  if(any(c('FRQ_A','FRQ_U') %in% header_interpretation$Interpreted) & sum(c('FRQ_A','FRQ_U') %in% header_interpretation$Interpreted) == 1){
    messages <- c(messages, "Warning: Both FRQ_A and FRQ_U must be present for either to be considered.")
  }

  if(!is.null(messages)){
    log_add(log_file = log_file, message = messages)

    if(error){
      stop(paste(messages, sep='\n'))
    }
  }

  # Remove columns that are not interpreted
  sumstats <- sumstats[, as.logical(header_interpretation$Keep), with=F]

  # Update sumstat header
  names(sumstats)<-header_interpretation$Interpreted[as.logical(header_interpretation$Keep)]

  # Check whether SNP contains RSIDs or is CHR:BP:A1:A2
  if('SNP' %in% names(sumstats)){
    if(all(grepl("^\\d+:[0-9]+:.+:.+$", sumstats$SNP))){
      log_add(log_file = log_file, message = 'SNP column contains CHR:BP:A1:A2 information.')
      snp_dat<-data.table(do.call(rbind, strsplit(sumstats$SNP, ':')))
      if(!('CHR' %in% names(sumstats))){
        log_add(log_file = log_file, message = 'CHR information extracted from SNP column.')
        sumstats$CHR<-snp_dat$V1
      }
      if(!('BP' %in% names(sumstats))){
        log_add(log_file = log_file, message = 'BP information extracted from SNP column.')
        sumstats$BP<-snp_dat$V2
      }
      log_add(log_file = log_file, message = 'SNP column has been removed.')
      sumstats$SNP <- NULL
    }
  }

  return(sumstats)
}

#' Calculate average allele frequency across cases and controls
#'
#' This function calculates the average allele frequency across cases and controls in GWAS (Genome-Wide Association Study) summary statistics. It accounts for different scenarios based on the available data, such as when case and control counts are available, or when a sampling fraction is provided. The function also logs the method used for the calculation.
#'
#' @param sumstats A data frame containing GWAS summary statistics, which should include columns for allele frequencies in cases (`FRQ_A`), controls (`FRQ_U`), and optionally the number of cases (`N_CAS`) and controls (`N_CON`).
#' @param sampling Optional; a numeric value representing the sampling fraction of cases in the study population. If not provided or NA, and `N_CAS` and `N_CON` are missing, the function assumes equal case and control numbers.
#' @param log_file Optional; a string specifying the path and name of the log file to record the process of the function. If NULL or not provided, no logging is performed.
#'
#' @return Returns a numeric vector representing the average allele frequency across cases and controls. If allele frequency information is not available in `sumstats`, it returns NULL.
#' @export
#'
#' @examples
#' mean_freq <- calc_mean_freq(clean_sumstats_1, sampling = 0.5)
#' print(mean_freq)
calc_mean_freq<-function(sumstats, sampling = NA, log_file = NULL){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(sumstats)) {
    sumstats <- data.table::as.data.table(sumstats)
  }

  # Convert data.frame to data.table if necessary
  if (is.data.frame(sumstats)) {
    sumstats <- data.table::as.data.table(sumstats)
  }

  if(all(c('FRQ_A','FRQ_U') %in% names(sumstats))){
    if(all(c('FRQ_A','FRQ_U','N_CAS','N_CON') %in% names(sumstats))){
      FREQ<- ((sumstats$FRQ_A*sumstats$N_CAS) + (sumstats$FRQ_U*sumstats$N_CON))/(sumstats$N_CAS + sumstats$N_CON)
      log_add(log_file = log_file, message = 'Average allele frequency calculated weighted by N_CAS and N_CON.')
    }
    if(all(c('FRQ_A','FRQ_U') %in% names(sumstats)) & any(!(c('FRQ_A','FRQ_U') %in% names(sumstats))) & !(is.na(sampling))){
      FREQ<- (sumstats$FRQ_A*sampling) + (sumstats$FRQ_U*(1-sampling))
      log_add(log_file = log_file, message = 'Average allele frequency calculated weighted by sampling fraction.')
    }
    if(all(c('FRQ_A','FRQ_U') %in% names(sumstats)) & any(!(c('FRQ_A','FRQ_U') %in% names(sumstats))) & is.na(sampling)){
      FREQ<- mean(sumstats$FRQ_A, sumstats$FRQ_U)
      log_add(log_file = log_file, message = 'Average allele frequency calculated assuming N_CAS == N_CON.')
    }
  } else {
    FREQ<-NULL
  }

  return(FREQ)
}

#' Create IUPAC Codes for pairs of SNP alleles
#'
#' This function takes two vectors representing pairs of alleles for SNP (Single Nucleotide Polymorphism) pairs and converts them into IUPAC (International Union of Pure and Applied Chemistry) codes. It returns an error message if the input vectors are of different lengths.
#'
#' @param x A character vector representing the first nucleotide of the allele pairs. Possible values are 'A', 'T', 'C', 'G'.
#' @param y A character vector of the same length as `x`, representing the second nucleotide of the allele pairs. Possible values are 'A', 'T', 'C', 'G'.
#'
#' @return A character vector of IUPAC codes corresponding to the input SNP allele pairs. If `x` and `y` are different lengths, the function prints an error message and does not return a value.
#' @export
#'
#' @examples
#' iupac_codes <- snp_iupac(clean_sumstats_1$A1, clean_sumstats_1$A2)
snp_iupac<-function(x=NA, y=NA){
  if(length(x) != length(y)){
    print('x and y are different lengths')
  } else {
    iupac<-rep(NA, length(x))
    iupac[x == 'A' & y =='T' | x == 'T' & y =='A']<-'W'
    iupac[x == 'C' & y =='G' | x == 'G' & y =='C']<-'S'
    iupac[x == 'A' & y =='G' | x == 'G' & y =='A']<-'R'
    iupac[x == 'C' & y =='T' | x == 'T' & y =='C']<-'Y'
    iupac[x == 'G' & y =='T' | x == 'T' & y =='G']<-'K'
    iupac[x == 'A' & y =='C' | x == 'C' & y =='A']<-'M'
    return(iupac)
  }
}

#' Detect the Genome Build of Target Using Reference Chromosome and Base Pair Information
#'
#' This function compares target SNP data with reference data to determine the genome build of the target. It checks chromosome and base pair concordance across different builds.
#'
#' @param ref A data frame containing reference SNP data with chromosome, base pair, and other relevant information.
#' @param targ A data frame containing target SNP data to be compared against the reference.
#' @param log_file An optional path to a log file where messages will be recorded.
#' @param overlap_thresh A numeric value specifying the minimum overlap proportion required to consider a genome build as a good match. The overlap is typically expressed as the fraction of shared SNPs between the target and reference datasets. A higher threshold ensures stricter matching, while a lower threshold allows for more flexibility when data quality or completeness is lower.
#'
#' @return Returns the detected genome build as a character string, or NA if the build cannot be determined.
#' @export
#'
#' @examples
#' reference_data_path <- gsub('22.rds', '',
#'   system.file("extdata", "ref.chr22.rds", package = "GenoUtils"))
#' ref <- readRDS(paste0(reference_data_path,'22.rds'))
#' targ <- clean_sumstats_1[clean_sumstats_1$CHR == 22,]
#' detected_build <- detect_build(ref, targ)
detect_build<-function(ref, targ, log_file = NULL, overlap_thresh=0.2){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(ref)) {
    ref <- data.table::as.data.table(ref)
  }
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  if(!(all(c('CHR','BP') %in% names(targ)))){
    stop('CHR and BP columns must be present in targ.\n')
  }
  if(!(all(c("CHR","A1","A2") %in% names(targ)) & any(grepl('BP_GRCh', names(ref))))){
    stop('CHR, A1, A2 and BP coordinates must be present in ref.\n')
  }

  # Insert IUPAC codes if missing
  if(!('IUPAC' %in% names(targ))){
    targ$IUPAC <- snp_iupac(targ$A1, targ$A2)
  }
  if(!('IUPAC' %in% names(ref))){
    ref$IUPAC <- snp_iupac(ref$A1, ref$A2)
  }

  builds<-as.numeric(gsub('BP_GRCh','',names(ref)[grepl('BP_GRCh', names(ref))]))

  # Check target-ref condordance of BP across builds
  ref$CHR<-as.character(ref$CHR)
  targ$CHR<-as.character(targ$CHR)

  matched<-list()
  build_overlap<-NULL
  target_build<-NA
  for(build_i in builds){
    matched<-merge(targ, ref, by.x=c('CHR','BP','IUPAC'), by.y=c('CHR',paste0('BP_GRCh', build_i),'IUPAC'))
    overlap_target <- nrow(matched) / nrow(targ)  # Proportion of matched SNPs in target
    overlap_ref <- nrow(matched) / nrow(ref)      # Proportion of matched SNPs in reference

    overlap<-nrow(matched)/nrow(targ)
    build_overlap<-rbind(build_overlap, data.frame(build = build_i,
                                                   nsnp = nrow(matched),
                                                   overlap_target = overlap_target,
                                                   overlap_ref = overlap_ref))
    log_add(log_file = log_file, message = paste0('GRCh', build_i, ' match: ', round(overlap_target * 100, 2), '% (Target), ', round(overlap_ref * 100, 2), '% (Ref)'))
  }

  best_build <- build_overlap[which.max(build_overlap$overlap_target), ]

  # Set a threshold where either the target or reference overlap is greater than 10%
  if (best_build$overlap_target > overlap_thresh | best_build$overlap_ref > overlap_thresh) {
    target_build <- paste0('GRCh', best_build$build)
  }

  return(target_build)
}

#' Convert SNP Alleles to Their Complementary Alleles
#'
#' This function converts each allele in a given vector to its complementary allele.
#'
#' @param x A character vector containing alleles ('A', 'T', 'C', 'G').
#'
#' @return A character vector with the complementary alleles.
#' @export
#'
#' @examples
#' complementary_alleles <- snp_allele_comp(c('A', 'C', 'G', 'T'))
snp_allele_comp<-function(x=NA){
  x_new<-x
  x_new[x == 'A']<-'T'
  x_new[x == 'T']<-'A'
  x_new[x == 'G']<-'C'
  x_new[x == 'C']<-'G'
  x_new[!(x %in% c('A','T','G','C'))]<-NA
  return(x_new)
}

#' Calculate Effective Sample Size
#'
#' This function calculates the effective sample size for a study based on the number of cases and controls.
#'
#' @param ncas The number of cases in the study.
#' @param ncon The number of controls in the study.
#'
#' @return Returns the effective sample size as a numeric value.
#' @export
#'
#' @examples
#' effective_size <- neff(1000, 1500)
neff<-function(ncas, ncon){
  return(4/((1/ncas)+(1/ncon)))
}

#' Remove Ambiguous SNPs from Data
#'
#' This function filters out SNPs with ambiguous IUPAC codes (R, Y, K, M) from the given data.
#'
#' @param dat A data frame containing SNP data with IUPAC codes or alleles (A1, A2).
#'
#' @return A data frame with ambiguous SNPs removed.
#' @export
#'
#' @examples
#' cleaned_data <- remove_ambig(clean_sumstats_1)
remove_ambig<-function(dat){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(dat)) {
    dat <- data.table::as.data.table(dat)
  }

  if(!('IUPAC' %in% names(dat))){
    if(!(all(c('A1','A2') %in% names(dat)))){
      stop('Either IUPAC, or A1 and A2 must be present')
    } else {
      iupac_codes <- snp_iupac(dat$A1, dat$A2)
      subset_condition <- iupac_codes %in% c('R', 'Y', 'K', 'M')
    }
  } else {
    subset_condition <- dat$IUPAC %in% c('R', 'Y', 'K', 'M')
  }

  return(dat[subset_condition, ])
}

#' Detect Strand Flips in SNP Data
#'
#' This function identifies SNPs in the target that are on the opposite strand compared to the reference.
#'
#' @param targ A vector of IUPAC codes for the target SNPs.
#' @param ref A vector of IUPAC codes for the reference SNPs.
#'
#' @return A logical vector indicating which SNPs are flipped.
#' @export
#'
#' @examples
#' flipped_snps <- detect_strand_flip(iupac_codes$target, iupac_codes$reference)
detect_strand_flip<-function(targ, ref){
  flipped<-(
      (!is.na(targ) & !is.na(ref)) &
      (
        (targ == 'R' & ref == 'Y') |
        (targ == 'Y' & ref == 'R') |
        (targ == 'K' & ref == 'M') |
        (targ == 'M' & ref == 'K')
      )
    )

  return(flipped)
}

#' Harmonize Target Data with Reference Data
#'
#' This function harmonizes target SNP data with reference SNP data by alleles, and either chromosome and base pair, or SNP ID. It accounts for strand flips between the reference and target.
#'
#' @param targ A data frame of target SNP data to be harmonized.
#' @param ref_rds The path to reference data files in RDS format.
#' @param population The reference population matching the GWAS sample. This option is used to determine ancestry-matched reference allele frequencies.
#' @param log_file An optional path to a log file where messages will be recorded.
#' @param chr An optional numeric vector indicating chromosomes to be processed.
#'
#' @return A harmonized data frame of target SNP data.
#' @export
#'
#' @import data.table
#'
#' @examples
#' # Get path and prefix to example ref_rds data
#' reference_data_path <- gsub( '22.rds','',
#'   system.file("extdata", "ref.chr22.rds", package = "GenoUtils"))
#' harmonised_data <- ref_harmonise(clean_sumstats_1, reference_data_path, 'EUR')
#' print(head(harmonised_data))
ref_harmonise<-function(targ, ref_rds, population, log_file = NULL, chr = 1:23){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  # Insert IUPAC in targ if not already present
  if(!any(names(targ) == 'IUPAC')){
    targ$IUPAC <- snp_iupac(targ$A1, targ$A2)
  }

  ref_tmp<-readRDS(file = paste0(ref_rds, max(chr), '.rds'))
  if(!(all(c('CHR','A1','A2','IUPAC') %in% names(ref_tmp)) & any(grepl('BP_GRCh', names(ref_tmp))))){
    stop('CHR, A1, A2, IUPAC and BP coordinates must be present in ref_rds files.\n')
  }

  if(!(any(grepl('REF.FRQ', names(ref_tmp))))){
    stop('There are no REF.FRQ.<POP> columns.\n')
  }

  ref_pops <- gsub('REF.FRQ.', '', names(ref_tmp)[grepl('REF.FRQ', names(ref_tmp))])
  if(!(population %in% ref_pops)){
    stop(paste0('Specified reference population must be one of the following: ', paste0(ref_pops, collapse=', '),'\n'))
  }

  # Check whether CHR and BP information are present
  chr_bp_avail<-sum(c('CHR','BP') %in% names(targ)) == 2

  # Check whether RSIDs are available for majority of SNPs in GWAS
  rsid_avail<-(sum(grepl('rs', targ$SNP)) > 0.9*length(targ$SNP))

  # Remove REF.FREQ column from targ if present
  if(any(names(targ) == 'REF.FREQ')){
    targ$REF.FREQ<-NULL
  }

  targ_matched<-NULL
  flip_logical_all<-NULL
  target_build<-NA

  if(chr_bp_avail){
    log_add(log_file = log_file, message = 'Merging sumstats with reference using CHR, BP, A1, and A2')

    ###
    # Determine build
    ###

    ref_tmp<-readRDS(file = paste0(ref_rds, max(targ$CHR), '.rds'))
    target_build <- detect_build( ref = ref_tmp,
                                  targ = targ[targ$CHR == max(targ$CHR),],
                                  log_file = log_file)

    if(!is.na(target_build)){
      for(i in chr){

        # Read reference data
        ref_i<-readRDS(file = paste0(ref_rds,i,'.rds'))

        # Retain only non-ambiguous SNPs
        ref_i<-remove_ambig(ref_i)

        # Rename columns prior to merging with target
        names(ref_i)<-paste0('REF.',names(ref_i))
        names(ref_i)[names(ref_i) == paste0('REF.REF.FRQ.',population)]<-'REF.FREQ'
        ref_i$BP<-ref_i[[paste0('REF.BP_',target_build)]]
        ref_i<-ref_i[, c('REF.CHR','REF.SNP','BP','REF.BP_GRCh37','REF.A1','REF.A2','REF.IUPAC','REF.FREQ'), with=F]

        # Subset chromosome i from target
        targ_i<-targ[targ$CHR == i,]

        # Merge target and reference by BP
        ref_target<-merge(targ_i, ref_i, by='BP')

        # Identify targ-ref strand flips, and flip target
        flip_logical<-detect_strand_flip(targ = ref_target$IUPAC, ref = ref_target$REF.IUPAC)
        flip_logical_all<-c(flip_logical_all, flip_logical)

        flipped<-ref_target[flip_logical,]
        flipped$A1<-snp_allele_comp(flipped$A1)
        flipped$A2<-snp_allele_comp(flipped$A2)
        flipped$IUPAC<-snp_iupac(flipped$A1, flipped$A2)

        # Identify SNPs that have matched IUPAC
        matched<-ref_target[ref_target$IUPAC == ref_target$REF.IUPAC,]
        matched<-rbind(matched, flipped)

        # Flip REF.FREQ if alleles are swapped
        matched$REF.FREQ[matched$A1 != matched$REF.A1]<-1-matched$REF.FREQ[matched$A1 != matched$REF.A1]

        # Retain reference SNP, REF.FREQ, and REF.BP_GRCh37 data
        matched<-matched[, names(matched) %in% c('CHR','REF.BP_GRCh37','REF.SNP','A1','A2','BETA','SE','OR','Z','FREQ','REF.FREQ','N','INFO','P'), with=F]
        names(matched)[names(matched) == 'REF.SNP']<-'SNP'
        names(matched)[names(matched) == 'REF.BP_GRCh37']<-'BP'

        targ_matched<-rbind(targ_matched, matched)
      }
    }

    if(is.na(target_build) & !(rsid_avail)){
      log_add(log_file = log_file, message = 'Error: Target build could not be determined and SNP IDs unavailable.')
      stop('Target build could not be determined and SNP IDs unavailable.\n')
    }
    if(is.na(target_build) & rsid_avail){
      log_add(log_file = log_file, message = 'Target build could not be determined from CHR and BP data.')
    }
  }

  if(is.na(target_build) & rsid_avail){
    log_add(log_file = log_file, message = 'Using SNP, A1 and A2 to merge with the reference.')

    for(i in chr){

      # Read reference data
      ref_i<-readRDS(file = paste0(ref_rds,i,'.rds'))

      # Retain only non-ambiguous SNPs
      ref_i<-remove_ambig(ref_i)

      # Rename columns prior to merging with target
      names(ref_i)<-paste0('REF.',names(ref_i))
      names(ref_i)[names(ref_i) == paste0('REF.REF.FRQ.',population)]<-'REF.FREQ'
      ref_i<-ref_i[, c('REF.CHR','REF.SNP','REF.BP_GRCh37','REF.A1','REF.A2','REF.IUPAC','REF.FREQ'), with=F]

      # Merge target and reference by SNP ID
      ref_target<-merge(targ, ref_i, by.x='SNP', by.y='REF.SNP')

      # Identify targ-ref strand flips, and flip target
      flip_logical<-detect_strand_flip(targ = ref_target$IUPAC, ref = ref_target$REF.IUPAC)
      flip_logical_all<-c(flip_logical_all, flip_logical)

      flipped<-ref_target[flip_logical,]
      flipped$A1<-snp_allele_comp(flipped$A1)
      flipped$A2<-snp_allele_comp(flipped$A2)
      flipped$IUPAC<-snp_iupac(flipped$A1, flipped$A2)

      # Identify SNPs that have matched IUPAC
      matched<-ref_target[ref_target$IUPAC == ref_target$REF.IUPAC,]
      matched<-rbind(matched, flipped)

      # Flip REF.FREQ if alleles are swapped
      matched$REF.FREQ[matched$A1 != matched$REF.A1]<-1-matched$REF.FREQ[matched$A1 != matched$REF.A1]

      # Retain reference CHR and BP_GRCh37 data
      matched<-matched[, names(matched) %in% c('REF.CHR','REF.BP_GRCh37','SNP','A1','A2','BETA','SE','OR','Z','FREQ','REF.FREQ','N','INFO','P'), with=F]
      names(matched)[names(matched) == 'REF.CHR']<-'CHR'
      names(matched)[names(matched) == 'REF.BP_GRCh37']<-'BP'

      targ_matched<-rbind(targ_matched, matched)
    }
  }

  log_add(log_file = log_file, message = paste0('After matching variants to the reference, ',nrow(targ_matched),' variants remain.'))
  log_add(log_file = log_file, message = paste0(sum(flip_logical_all), ' variants were flipped to match reference.'))

  return(targ_matched)
}

#' Filter Variants Based on INFO Score Threshold
#'
#' This function removes variants from the target data that have an INFO score below the specified threshold.
#'
#' @param targ A data frame of SNP data with an INFO column.
#' @param thresh A numeric threshold for the INFO score.
#' @param log_file An optional path to a log file where messages will be recorded.
#'
#' @return A filtered data frame of SNP data.
#' @export
#'
#' @import data.table
#'
#' @examples
#' filtered_data <- filter_info(clean_sumstats_1, 0.8)
filter_info<-function(targ, thresh, log_file = NULL){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  if(sum(names(targ) == 'INFO') == 1){
    targ<-targ[targ$INFO >= thresh,]
    log_add(log_file = log_file, message = paste0('After removal of SNPs with INFO < ',thresh,', ',nrow(targ),' variants remain.'))
  } else {
    log_add(log_file = log_file, message = 'INFO column is not present.')
  }
  return(targ)
}

#' Filter Variants Based on Minor Allele Frequency (MAF) Threshold
#'
#' This function removes variants from the target dataset that have a minor allele frequency (MAF) below a specified threshold. It uses either the reported MAF (`FREQ` column) or the reference MAF (`REF.FREQ` column), depending on the `ref` argument.
#'
#' @param targ A data frame containing SNP data.
#' @param thresh A numeric threshold for MAF.
#' @param ref Logical, if TRUE, uses the `REF.FREQ` column for filtering, otherwise uses `FREQ`.
#' @param log_file Optional; path to a log file where messages will be recorded.
#'
#' @return A filtered data frame with variants meeting the MAF threshold.
#' @export
#'
#' @import data.table
#'
#' @examples
#' filtered_data <- filter_maf(targ = clean_sumstats_1, thresh = 0.05, ref = FALSE)
filter_maf<-function(targ, thresh, ref = F, log_file = NULL){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  if(ref == F){
    if(any(names(targ) == 'FREQ')){
      targ<-targ[targ$FREQ >= thresh & targ$FREQ <= (1-thresh),]
      log_add(log_file = log_file, message = paste0('After removal of SNPs with reported MAF < ',thresh,', ',nrow(targ),' variants remain.'))
    } else {
      log_add(log_file = log_file, message = 'FREQ  column is not present.')
    }
  } else {
    if(any(names(targ) == 'REF.FREQ')){
      targ<-targ[targ$REF.FREQ >= thresh & targ$REF.FREQ <= (1-thresh),]
      log_add(log_file = log_file, message = paste0('After removal of SNPs with reference MAF < ',thresh,', ',nrow(targ),' variants remain.'))
    } else {
      log_add(log_file = log_file, message = 'REF.FREQ  column is not present.')
    }
  }

  return(targ)
}

#' Remove Variants with Discrepancy in Reported and Reference Allele Frequencies
#'
#' This function filters out variants from the target dataset where the absolute difference between reported and reference allele frequencies exceeds a specified threshold. Optionally, it can generate a plot of these discrepancies.
#'
#' @param targ A data frame containing SNP data with both reported and reference allele frequencies.
#' @param thresh A numeric threshold for the absolute frequency difference.
#' @param log_file Optional; path to a log file where messages will be recorded.
#' @param plot_file Optional; path to save a plot file visualizing the frequency discrepancies.
#'
#' @return A filtered data frame with variants within the specified frequency discrepancy threshold.
#' @export
#'
#' @import data.table
#' @import grDevices
#' @import graphics

#' @examples
#' filtered_data <- discord_maf( targ = clean_sumstats_1,
#'                               thresh = 0.01,
#'                               plot_file = paste0(tempdir(),"/plot.png"))
discord_maf<-function(targ, thresh, log_file = NULL, plot_file = NA){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  if(!(any(names(targ) %in% 'REF.FREQ'))) {
    stop('REF.FREQ column must be present in targ')
  }

  if(sum(names(targ) == 'FREQ') == 1){
    targ$diff<-abs(targ$FREQ-targ$REF.FREQ)

    if(!is.na(plot_file)){
      png(plot_file, units='px', res=300, width=1200, height=1200)
      plot(targ$REF.FREQ[targ$diff > thresh],targ$FREQ[targ$diff > thresh], xlim=c(0,1), ylim=c(0,1), xlab='Reference Allele Frequency', ylab='Sumstat Allele Frequency')
      abline(coef = c(0,1))
      dev.off()
    }

    targ<-targ[targ$diff < thresh,]
    targ$diff<-NULL

    log_add(log_file = log_file, message = paste0('After removal of SNPs with absolute MAF difference of < ',thresh,', ',nrow(targ),' variants remain.'))
  } else {
    log_add(log_file = log_file, message = 'Reported MAF column is not present, so discordance with reference cannot be determined.')
  }
  return(targ)
}

#' Filter Variants Based on Sample Size Deviation
#'
#' This function removes variants from the target dataset where the sample size (`N`) is more than a specified number of standard deviations away from the median sample size.
#'
#' @param targ A data frame containing SNP data with a sample size column (`N`).
#' @param n_sd The number of standard deviations from the median to use as the threshold.
#' @param log_file Optional; path to a log file where messages will be recorded.
#'
#' @return A filtered data frame with variants within the specified sample size range.
#' @export
#'
#' @import data.table
#' @import stats
#'
#' @examples
#' filtered_data <- filter_n(targ = clean_sumstats_1, n_sd = 3)
filter_n <- function(targ, n_sd = 3, log_file = NULL){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  if(length(unique(targ$N)) > 1){
    thresh <- n_sd*sd(targ$N)
    targ <- targ[targ$N < median(targ$N) + thresh & targ$N > median(targ$N) - thresh,]

    log_add(log_file = log_file, message = paste0('After removal of SNPs with N > ',median(targ$N)+thresh,' or < ',median(targ$N)-thresh,', ',nrow(targ),' variants remain.'))
  } else {
    log_add(log_file = log_file, message = 'N column is not present or invariant.')
  }
  return(targ)
}

#' Convert Odds Ratios or Z-scores to Beta Coefficients
#'
#' This function calculates Beta coefficients for SNP data, using either available odds ratios (OR) or Z-scores and allele frequencies, if Beta coefficients are not already present.
#'
#' @param targ A data frame containing SNP data with OR or Z-scores.
#' @param log_file Optional; path to a log file where messages will be recorded.
#'
#' @return The input data frame with added `BETA` column.
#' @export
#'
#' @import data.table
#'
#' @examples
#' updated_data <- insert_beta(targ = clean_sumstats_1)
insert_beta <- function(targ, log_file = NULL){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  # If OR is present but BETA is not, convert OR to logOR and name it BETA
  if('OR' %in% names(targ) & !('BETA' %in% names(targ))){
    targ$BETA<-log(targ$OR)
    log_add(log_file = log_file, message = 'BETA column inserted based on log(OR).')
  }

  # If Z is present but BETA is not, calculate it based on sample size, allele frequency and Z score
  if('Z' %in% names(targ) & !('BETA' %in% names(targ))){
    if('FREQ' %in% names(targ)){
      frq_tmp<-targ$FREQ
      log_add(log_file = log_file, message = 'BETA and SE column inserted based on Z, FREQ, and N.')
    } else {
      frq_tmp<-targ$REF.FREQ
      log_add(log_file = log_file, message = 'BETA and SE column inserted based on Z, REF.FREQ, and N')
    }

    targ$SE <- 1/sqrt((2*frq_tmp)*(1-(frq_tmp))*(targ$N + (targ$Z^2)))
    targ$BETA <- targ$Z * targ$SE
  }
  return(targ)
}

#' Compute Standard Error from Beta and P-value
#'
#' This function calculates the standard error (SE) for each SNP's Beta coefficient, based on its P-value, if SE is not already present.
#'
#' @param targ A data frame containing SNP data with Beta coefficients and P-values.
#' @param log_file Optional; path to a log file where messages will be recorded.
#'
#' @return The input data frame with added `SE` column.
#' @export
#'
#' @import data.table
#' @import stats
#'
#' @examples
#' updated_data <- insert_se(targ = clean_sumstats_1)
insert_se <- function(targ, log_file = NULL){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  if(sum(names(targ) == 'SE') == 0){
    if(any(!(c('BETA','P') %in% names(targ)))){
      stop('BETA and P columns must be present to compute SE.\n')
    }

    z_tmp<-abs(qnorm(targ$P/2))
    targ$SE<-abs(targ$BETA/z_tmp)

    log_add(log_file = log_file, message = 'SE column inserted based on BETA and P.')
  }
  return(targ)
}

#' Adjust P-values to Avoid Genomic Control Artifacts
#'
#' This function checks for and adjusts P-values that may have been affected by genomic control. It recalculates P-values using available standard errors and Beta coefficients.
#'
#' @param targ A data frame containing SNP data with Beta coefficients, standard errors, and P-values.
#' @param log_file Optional; path to a log file where messages will be recorded.
#'
#' @return The input data frame with adjusted P-values.
#' @export
#'
#' @import data.table
#' @import stats
#'
#' @examples
#' # Read in example clean sumstats
#' adjusted_data <- avoid_gc(targ = clean_sumstats_1)
avoid_gc<-function(targ, log_file = NULL){
  # Convert data.frame to data.table if necessary
  if (is.data.frame(targ)) {
    targ <- data.table::as.data.table(targ)
  }

  if(sum(names(targ) == 'SE') == 1){
    targ$Z<-targ$BETA/targ$SE
    targ$P_check<-2*pnorm(-abs(targ$Z))
    targ$Z<-NULL

    if(abs(mean(targ$P[!is.na(targ$P_check)]) - mean(targ$P_check[!is.na(targ$P_check)])) > 0.01){
      targ$P<-targ$P_check
      targ$P_check<-NULL
      log_add(log_file = log_file, message = 'Genomic control detected. P-value recomputed using BETA and SE.')
    } else {
      log_add(log_file = log_file, message = 'Genomic control was not detected.')
      targ$P_check<-NULL
    }
  } else {
    log_add(log_file = log_file, message = 'SE column is not present, genomic control cannot be detected.')
  }
  return(targ)
}
