# ============================================================
# Colocalization analysis using coloc.abf
# ============================================================

library(coloc)
library(data.table)
library(dplyr)
library(readxl)
library(TwoSampleMR)
library(tidyverse)

# ============================================================
# User-defined parameters
# ============================================================

WORK_DIR <- "D:/FJH Dynamic coloc"

INPUT_FILE <- file.path(WORK_DIR, "input.xlsx")
OUTCOME_FILE <- file.path(WORK_DIR, "GCST004748_buildGRCh38.tsv")
EQTL_DIR <- file.path(WORK_DIR, "Dynamic original data")
OUTPUT_DIR <- file.path(WORK_DIR, "results")

EQTL_SAMPLE_SIZE <- 119
GWAS_SAMPLE_SIZE <- 85716
CASE_PROPORTION <- 0.34
OUTCOME_NAME <- "Outcome"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ============================================================
# Read input file
# ============================================================

input <- read_excel(INPUT_FILE)

all_results <- list()

# ============================================================
# Perform colocalization analysis
# ============================================================

for (i in seq_len(nrow(input))) {
  
  message("Processing gene: ", input$exposure[i])
  
  gene_dir <- file.path(
    OUTPUT_DIR,
    input$exposure[i]
  )
  
  dir.create(
    gene_dir,
    showWarnings = FALSE
  )
  
  eqtl_file <- file.path(
    EQTL_DIR,
    paste0(
      input$type[i],
      "_500kb_processed.tsv.gz"
    )
  )
  
  gene_all <- fread(eqtl_file)
  
  gene <- gene_all %>%
    filter(
      exposure == input$exposure[i]
    )
  
  if (nrow(gene) == 0) {
    
    message(
      "No eQTL records found for ",
      input$exposure[i]
    )
    
    next
    
  }
  
  gene$samplesize.exposure <- EQTL_SAMPLE_SIZE
  
  gene <- format_data(
    as.data.frame(gene),
    type = "exposure",
    phenotype_col = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    eaf_col = "eaf.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure"
  )
  
  outcome_dat <- read_outcome_data(
    snps = gene$SNP,
    filename = OUTCOME_FILE,
    sep = "\t",
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
    eaf_col = "effect_allele_frequency"
  )
  
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    next
  }
  
  outcome_dat$samplesize.outcome <- GWAS_SAMPLE_SIZE
  outcome_dat$outcome <- OUTCOME_NAME
  
  dat <- harmonise_data(
    gene,
    outcome_dat,
    action = 1
  )
  
  dat <- dat[
    !duplicated(dat$SNP),
  ]
  
  result <- coloc.abf(
    
    dataset1 = list(
      snp = dat$SNP,
      pvalues = dat$pval.exposure,
      type = "quant",
      N = dat$samplesize.exposure,
      MAF = dat$eaf.exposure
    ),
    
    dataset2 = list(
      snp = dat$SNP,
      pvalues = dat$pval.outcome,
      type = "cc",
      s = CASE_PROPORTION,
      N = dat$samplesize.outcome,
      MAF = dat$eaf.outcome
    )
    
  )
  
  pp <- result$summary
  
  summary_result <- data.frame(
    
    type = input$type[i],
    exposure = input$exposure[i],
    nSNPs = pp["nsnps"],
    PP0 = pp["PP.H0.abf"],
    PP1 = pp["PP.H1.abf"],
    PP2 = pp["PP.H2.abf"],
    PP3 = pp["PP.H3.abf"],
    PP4 = pp["PP.H4.abf"]
    
  )
  
  write.table(
    summary_result,
    file = file.path(
      gene_dir,
      "coloc_summary.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  all_results[[length(all_results) + 1]] <- summary_result
  
}

# ============================================================
# Export summary of all genes
# ============================================================

final_results <- bind_rows(all_results)

write.table(
  
  final_results,
  
  file = file.path(
    OUTPUT_DIR,
    "all_coloc_results.txt"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)

message("Colocalization analysis completed.")