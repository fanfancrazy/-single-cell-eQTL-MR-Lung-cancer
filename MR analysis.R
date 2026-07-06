library(TwoSampleMR)
library(dplyr)
library(psych)
library(tools)

# ============================================================
# Helper functions for SNP-level R² calculation
# ============================================================

get_r_from_pn <- function(p, n) {
  z <- qnorm(p / 2, lower.tail = FALSE)
  r <- z / sqrt(z^2 + n - 2)
  abs(r)
}

get_r_from_lor <- function(lor, af, ncase, ncontrol) {
  n <- ncase + ncontrol
  abs(lor) * sqrt(af * (1 - af) / n)
}

# ============================================================
# Working directory and input files
# ============================================================

setwd("E:\\immene cell MR\\bulk eqtlgen")

outcome_file <- "GCST004744_buildGRCh37.tsv"
outcome_name <- gsub(".tsv", "", outcome_file)

# Alzheimer's disease GWAS sample size (Lambert et al., 2013)
N_CASE <- 11273
N_CONTROL <- 55483
N_OUTCOME <- N_CASE + N_CONTROL

exposure_files <- list.files(pattern = "\\.txt$", full.names = FALSE)

if (!dir.exists("output")) {
  dir.create("output")
}

# ============================================================
# Steiger filtering option
# TRUE  = strict SNP-level Steiger filtering
# FALSE = use harmonised SNPs only (default TwoSampleMR workflow)
# ============================================================

USE_STRICT_STEIGER <- FALSE

# ============================================================
# Main analysis
# ============================================================

for (exposure_file in exposure_files) {
  
  message("\n===== Processing: ", exposure_file, " =====")
  
  exposure_name <- file_path_sans_ext(exposure_file)
  
  # ----------------------------------------------------------
  # Read exposure data
  # ----------------------------------------------------------
  
  exposure_dat <- read_exposure_data(
    filename = exposure_file,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    pval_col = "pval.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    samplesize_col = "samplesize.exposure",
    chr_col = "chr",
    pos_col = "pos",
    phenotype_col = "gene",
    id_col = "gene.exposure",
    clump = FALSE
  )
  
  # ----------------------------------------------------------
  # Read outcome data
  # ----------------------------------------------------------
  
  outcome_dat <- read_outcome_data(
    filename = outcome_file,
    sep = "\t",
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
    eaf_col = "effect_allele_frequency"
  )
  
  outcome_dat$samplesize.outcome <- N_OUTCOME
  outcome_dat$outcome <- outcome_name
  
  # ----------------------------------------------------------
  # Harmonise exposure and outcome datasets
  # ----------------------------------------------------------
  
  dat <- harmonise_data(exposure_dat, outcome_dat)
  
  if (nrow(dat) == 0) {
    message("No overlapping SNPs. Skipping...")
    next
  }
  
  # ----------------------------------------------------------
  # Calculate SNP-level Steiger statistics
  # ----------------------------------------------------------
  
  dat$rsq.exposure <- get_r_from_pn(
    dat$pval.exposure,
    dat$samplesize.exposure
  )^2
  
  dat$rsq.outcome <- get_r_from_lor(
    dat$beta.outcome,
    dat$eaf.outcome,
    N_CASE,
    N_CONTROL
  )^2
  
  steiger_res <- psych::r.test(
    n = dat$samplesize.exposure,
    n2 = N_OUTCOME,
    r12 = sqrt(pmax(dat$rsq.exposure, 0)),
    r34 = sqrt(pmax(dat$rsq.outcome, 0))
  )
  
  dat$steiger_dir <- dat$rsq.exposure > dat$rsq.outcome
  dat$steiger_pval <- steiger_res$p
  
  # ----------------------------------------------------------
  # Select SNPs for MR analysis
  # ----------------------------------------------------------
  
  if (USE_STRICT_STEIGER) {
    
    mr_dat <- subset(
      dat,
      mr_keep & steiger_dir & steiger_pval < 0.05
    )
    
    message(
      "Using strict Steiger filtering: ",
      nrow(mr_dat),
      " SNPs retained."
    )
    
  } else {
    
    mr_dat <- subset(dat, mr_keep)
    
    message(
      "Using harmonised SNPs only: ",
      nrow(mr_dat),
      " SNPs retained."
    )
    
  }
  
  if (USE_STRICT_STEIGER && nrow(mr_dat) < 3) {
    message(
      "Warning: Fewer than 3 SNPs remained after strict Steiger filtering."
    )
  }
  
  # ----------------------------------------------------------
  # Export harmonised SNPs with Steiger statistics
  # ----------------------------------------------------------
  
  harmonised_all <- subset(dat, mr_keep) %>%
    select(
      SNP,
      exposure,
      outcome,
      chr.exposure,
      pos.exposure,
      beta.exposure,
      se.exposure,
      pval.exposure,
      samplesize.exposure,
      beta.outcome,
      se.outcome,
      pval.outcome,
      eaf.outcome,
      rsq.exposure,
      rsq.outcome,
      steiger_dir,
      steiger_pval
    )
  
  write.csv(
    harmonised_all,
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".harmonised_with_steiger.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # Mendelian randomization analysis
  # ----------------------------------------------------------
  
  mrResult <- mr(mr_dat)
  
  if (is.null(mrResult) || nrow(mrResult) == 0) {
    message("No valid MR results. Skipping remaining steps.")
    next
  }
  
  mrTab <- generate_odds_ratios(mrResult)
  
  write.csv(
    mrTab,
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".MR_results.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    directionality_test(mr_dat),
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".steiger_overall.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    mr_heterogeneity(mr_dat),
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".heterogeneity.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    mr_pleiotropy_test(mr_dat),
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".pleiotropy.csv"
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # Generate MR plots
  # ----------------------------------------------------------
  
  pdf(
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".scatter.pdf"
    ),
    7.5,
    7
  )
  
  print(mr_scatter_plot(mrResult, mr_dat))
  dev.off()
  
  sing <- mr_singlesnp(mr_dat)
  
  pdf(
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".forest.pdf"
    ),
    7,
    6
  )
  
  print(mr_forest_plot(sing))
  dev.off()
  
  pdf(
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".funnel.pdf"
    ),
    7,
    6.5
  )
  
  print(mr_funnel_plot(sing))
  dev.off()
  
  pdf(
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".leaveoneout.pdf"
    ),
    7,
    6
  )
  
  print(
    mr_leaveoneout_plot(
      mr_leaveoneout(mr_dat)
    )
  )
  
  dev.off()
  
  # ----------------------------------------------------------
  # FDR correction for the primary MR estimate
  # ----------------------------------------------------------
  
  main_res <- mrResult %>%
    group_by(exposure) %>%
    filter(
      method %in% c(
        "Inverse variance weighted",
        "Wald ratio"
      )
    ) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      pval_fdr = p.adjust(pval, method = "fdr")
    )
  
  write.csv(
    main_res,
    paste0(
      "output/",
      exposure_name,
      "_",
      outcome_name,
      ".MR_results_FDR.csv"
    ),
    row.names = FALSE
  )
  
  message(
    "Finished: ",
    exposure_name,
    " | SNPs used in MR: ",
    nrow(mr_dat),
    "\n"
  )
  
}

message("All analyses completed.")