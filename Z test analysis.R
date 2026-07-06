library(dplyr)
library(readxl)
library(writexl)

# ============================================================
# User-defined parameters
# ============================================================

WORK_DIR <- "E:/immene cell MR/DSY Dynamic MR result/Z-test"

INPUT_FILE <- "Dynamic_with_timepoint.xlsx"

FILTERED_FILE <- "filtered_results.xlsx"
REFERENCE_FILE <- "reference_timepoint.xlsx"
OUTPUT_FILE <- "Ztest_results.xlsx"

setwd(WORK_DIR)

# ============================================================
# Read input data
# ============================================================

df <- read_excel(INPUT_FILE)

df$time_point <- as.character(df$time_point)

activation_points <- c(
  "lowly active",
  "16h",
  "40h",
  "5d"
)

resting_point <- "0h"

# ============================================================
# Identify gene–outcome–cell type combinations for Z-test
# ============================================================

df_filtered <- df %>%
  group_by(gene, outcome, cell_type) %>%
  filter(
    any(
      time_point %in% activation_points &
        pval_fdr < 0.05,
      na.rm = TRUE
    ) &
      n() >= 2
  ) %>%
  ungroup()

write_xlsx(
  df_filtered,
  FILTERED_FILE
)

# ============================================================
# Select the reference time point
# (largest absolute MR effect)
# ============================================================

df_ref <- df_filtered %>%
  group_by(gene, outcome, cell_type) %>%
  mutate(abs_beta = abs(beta)) %>%
  slice(which.max(abs_beta)) %>%
  ungroup() %>%
  select(
    gene,
    outcome,
    cell_type,
    ref_time_point = time_point,
    beta_ref = beta,
    se_ref = se,
    pval_fdr_ref = pval_fdr
  )

write_xlsx(
  df_ref,
  REFERENCE_FILE
)

# ============================================================
# Perform Z-test
# Exclude the reference time point itself
# ============================================================

df_ztest <- df_filtered %>%
  left_join(
    df_ref,
    by = c(
      "gene",
      "outcome",
      "cell_type"
    )
  ) %>%
  filter(
    time_point != ref_time_point
  ) %>%
  mutate(
    Z = (beta - beta_ref) /
      sqrt(se^2 + se_ref^2),
    pval_Z = 2 * pnorm(-abs(Z))
  )

# ============================================================
# Export results
# ============================================================

write_xlsx(
  df_ztest,
  OUTPUT_FILE
)

message("Z-test analysis completed.")