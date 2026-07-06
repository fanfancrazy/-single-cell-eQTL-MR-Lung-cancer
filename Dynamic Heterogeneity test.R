# ============================================================
# Time-point specific heterogeneity analysis across genes
# Outputs: I², Cochran's Q, Q p-value, and MR significance
# ============================================================

rm(list = ls())

library(dplyr)
library(metafor)

# ============================================================
# 1. Load data
# ============================================================

setwd("E:/immene cell MR/DSY Dynamic MR result/Time-specific analysis + Fisher test by outcome")

data <- read.csv(
  "Dynamic_with_timepoint.csv",
  stringsAsFactors = FALSE
)

# Standardize column order (optional but recommended)
data <- data %>%
  select(
    gene,
    exposure,
    beta,
    se,
    pval_fdr,
    outcome,
    gene_name,
    time_point,
    cell_type
  )

# ============================================================
# 2. Heterogeneity analysis per gene–outcome–cell type group
# ============================================================

results <- data %>%
  group_by(gene, gene_name, outcome, cell_type) %>%
  group_modify(~{
    
    dat <- .x
    n_timepoints <- nrow(dat)
    
    # --------------------------------------------------------
    # MR significance across time points
    # --------------------------------------------------------
    MR_significance <- ifelse(
      any(dat$pval_fdr < 0.05),
      "YES",
      "NO"
    )
    
    # --------------------------------------------------------
    # Heterogeneity analysis (fixed-effect model)
    # --------------------------------------------------------
    if (n_timepoints > 1) {
      
      res <- tryCatch(
        rma.uni(
          yi = dat$beta,
          sei = dat$se,
          method = "FE"
        ),
        error = function(e) NULL
      )
      
      if (!is.null(res)) {
        
        tibble(
          n_timepoints = n_timepoints,
          I2 = round(res$I2, 3),
          QE = round(res$QE, 3),
          QEp = signif(res$QEp, 3),
          MR_significance = MR_significance
        )
        
      } else {
        
        tibble(
          n_timepoints = n_timepoints,
          I2 = NA,
          QE = NA,
          QEp = NA,
          MR_significance = MR_significance
        )
        
      }
      
    } else {
      
      tibble(
        n_timepoints = n_timepoints,
        I2 = NA,
        QE = NA,
        QEp = NA,
        MR_significance = MR_significance
      )
      
    }
    
  }) %>%
  ungroup()

# ============================================================
# 3. Export results
# ============================================================

write.csv(
  results,
  "Dynamic_with_timepoint_heterogeneity_results.csv",
  row.names = FALSE
)

# Preview results
head(results)