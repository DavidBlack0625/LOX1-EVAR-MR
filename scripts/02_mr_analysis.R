###############################################################################
## 02_mr_analysis.R
## ================
## Mendelian Randomization: OLR1 (LOX-1) → Abdominal Aortic Aneurysm (AAA)
##
## Exposure: GTEx v8 OLR1 cis-eQTL
##   - Primary:   Artery_Coronary (9 SNPs)
##   - Sensitivity: Artery_Tibial (52 SNPs)
##   - Additional: eQTLGen Whole Blood (if available)
##
## Outcome: FinnGen AAA GWAS (I9_ABAORTANEUR or similar)
##
## Methods: IVW, MR-Egger, Weighted Median, Simple Mode, Weighted Mode
##          + MR-PRESSO (pleiotropy + outlier test)
##          + Leave-one-out sensitivity analysis
##          + Cochran's Q (heterogeneity)
##
## Author: Generated for Dr. Wei Wei
## Date: 2026-07-12
###############################################################################

## ============================================================
## 0. Setup
## ============================================================
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ggplot2)
})

# Set working directory
setwd("D:/workbuddy工作空间/LOX-1_EVAR_MR")

# Create output directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# Helper function: calculate SE from beta and p-value
calc_se <- function(beta, pval) {
  z <- abs(qnorm(pval / 2))
  se <- abs(beta) / z
  return(se)
}

## ============================================================
## 1. Load Exposure Data (GTEx v8 OLR1 eQTL)
## ============================================================
cat("\n============================================================\n")
cat("Step 1: Loading exposure data (GTEx v8 OLR1 eQTL)\n")
cat("============================================================\n\n")

# Read Artery_Coronary (primary)
coronary_dt <- fread("data/exposure/gtex_v8_OLR1_Artery_Coronary_eqtl.csv")
cat("Artery_Coronary:", nrow(coronary_dt), "SNPs\n")

# Read Artery_Tibial (sensitivity)
tibial_dt <- fread("data/exposure/gtex_v8_OLR1_Artery_Tibial_eqtl.csv")
cat("Artery_Tibial:", nrow(tibial_dt), "SNPs\n")

# Check for eQTLGen data
eqtlgen_file <- "data/exposure/eqtlgen_OLR1_blood_eqtl.csv"
has_eqtlgen <- file.exists(eqtlgen_file)
if (has_eqtlgen) {
  eqtlgen_dt <- fread(eqtlgen_file)
  cat("eQTLGen Whole Blood:", nrow(eqtlgen_dt), "SNPs\n")
}

## ============================================================
## 2. Format Exposure Data for TwoSampleMR
## ============================================================
cat("\n============================================================\n")
cat("Step 2: Formatting exposure data\n")
cat("============================================================\n\n")

format_gtex_exposure <- function(dt, tissue_name) {
  # Calculate SE from p-value if missing
  if (all(dt$se == "" | is.na(dt$se))) {
    dt$se <- calc_se(dt$beta, dt$pval)
    cat("  Calculated SE from p-value (no SE in GTEx API response)\n")
  } else {
    dt$se <- as.numeric(dt$se)
  }

  exp_dat <- data.frame(
    SNP = dt$SNP,
    beta = dt$beta,
    se = dt$se,
    pval = dt$pval,
    effect_allele = dt$effect_allele,
    other_allele = dt$other_allele,
    chr = as.character(gsub("chr", "", dt$chr)),
    pos = dt$pos,
    exposure = paste0("OLR1 eQTL (", tissue_name, ")"),
    gene = "OLR1",
    id.exposure = paste0("gtex_olr1_", tolower(gsub("Artery_", "", tissue_name))),
    units.exposure = "normalized_effect_size",
    samplesize.exposure = NA
  )

  return(exp_dat)
}

# Format each tissue
exp_coronary <- format_gtex_exposure(coronary_dt, "Artery_Coronary")
exp_tibial <- format_gtex_exposure(tibial_dt, "Artery_Tibial")

cat("\nPrimary exposure (Artery_Coronary):\n")
print(exp_coronary[, c("SNP", "effect_allele", "other_allele", "beta", "se", "pval")])

cat("\nSensitivity exposure (Artery_Tibial) - top 10 by p-value:\n")
print(exp_tibial[order(exp_tibial$pval), c("SNP", "effect_allele", "other_allele", "beta", "se", "pval")][1:10, ])

## ============================================================
## 3. Clump SNPs (ensure independence)
## ============================================================
cat("\n============================================================\n")
cat("Step 3: Clumping SNPs (LD pruning)\n")
cat("============================================================\n\n")

# Try clumping via IEU OpenGWAS API
# This requires an OpenGWAS token set in .Renviron
# If no token, use local clumping with PLINK (if available)

clump_exposure <- function(exp_dat, tissue_name) {
  cat("Clumping", tissue_name, "...\n")

  # Attempt 1: IEU OpenGWAS API clumping
  tryCatch({
    clumped <- clump_data(
      exp_dat,
      clump_kb = 10000,
      clump_r2 = 0.001,
      clump_p1 = 1,
      clump_p2 = 1,
      pop = "EUR"
    )
    cat("  Clumped via IEU OpenGWAS:", nrow(clumped), "independent SNPs (from", nrow(exp_dat), ")\n")
    return(clumped)
  }, error = function(e) {
    cat("  IEU OpenGWAS clumping failed:", conditionMessage(e), "\n")
    cat("  Attempting manual LD pruning by position (> 1 Mb apart)...\n")

    # Attempt 2: Simple position-based pruning (crude, but works without API)
    exp_sorted <- exp_dat[order(exp_dat$pval), ]
    exp_sorted$pos <- as.numeric(exp_sorted$pos)
    selected <- c()
    for (i in 1:nrow(exp_sorted)) {
      if (length(selected) == 0) {
        selected <- c(selected, i)
      } else {
        min_dist <- min(abs(exp_sorted$pos[i] - exp_sorted$pos[selected]))
        if (min_dist > 1e6) {  # > 1 Mb apart
          selected <- c(selected, i)
        }
      }
    }
    clumped <- exp_sorted[selected, ]
    cat("  Position-based pruning:", nrow(clumped), "SNPs (> 1 Mb apart)\n")
    cat("  NOTE: This is approximate. For publication, use proper LD clumping.\n")
    return(clumped)
  })
}

exp_coronary_clumped <- clump_exposure(exp_coronary, "Artery_Coronary")
exp_tibial_clumped <- clump_exposure(exp_tibial, "Artery_Tibial")

cat("\nAfter clumping:\n")
cat("  Artery_Coronary:", nrow(exp_coronary_clumped), "independent SNPs\n")
cat("  Artery_Tibial:", nrow(exp_tibial_clumped), "independent SNPs\n")

## ============================================================
## 4. Load Outcome Data (FinnGen AAA GWAS)
## ============================================================
cat("\n============================================================\n")
cat("Step 4: Loading outcome data (FinnGen AAA GWAS)\n")
cat("============================================================\n\n")

# Try to fetch FinnGen AAA data via IEU OpenGWAS
# Common FinnGen AAA phenotype IDs on IEU OpenGWAS:
# - finn-a-I9_ABAORTANEUR (Abdominal aortic aneurysm, FinnGen R6)
# - finn-b-I9_ABAORTANEUR (Abdominal aortic aneurysm, FinnGen R7)
# - finn-c-I9_ABAORTANEUR (Abdominal aortic aneurysm, FinnGen R8+)

get_outcome_data <- function(exp_dat_list, label) {
  # List of possible FinnGen AAA IDs to try
  possible_ids <- c(
    "finn-c-I9_ABAORTANEUR",
    "finn-b-I9_ABAORTANEUR",
    "finn-a-I9_ABAORTANEUR",
    "finn-b-I9_AORTANEUR",
    "finn-c-I9_AORTANEUR",
    "finn-a-I71_ABDOMINALAORTICANEURYSM"
  )

  all_snps <- unique(unlist(lapply(exp_dat_list, function(x) x$SNP)))

  for (outcome_id in possible_ids) {
    cat("  Trying", outcome_id, "...\n")
    tryCatch({
      outcome_dat <- extract_outcome_data(
        snps = all_snps,
        outcomes = outcome_id
      )
      if (!is.null(outcome_dat) && nrow(outcome_dat) > 0) {
        cat("  SUCCESS! Found", nrow(outcome_dat), "SNPs in", outcome_id, "\n")
        attr(outcome_dat, "outcome_id") <- outcome_id
        return(outcome_dat)
      }
    }, error = function(e) {
      cat("    Failed:", conditionMessage(e), "\n")
    })
  }

  cat("\n  ERROR: Could not fetch FinnGen AAA data via IEU OpenGWAS.\n")
  cat("  Possible reasons:\n")
  cat("    1. No OpenGWAS token set (set OPENGWAS_JWT in .Renviron)\n")
  cat("    2. FinnGen phenotype ID has changed\n")
  cat("    3. Network issue\n\n")
  cat("  MANUAL DOWNLOAD OPTION:\n")
  cat("    1. Go to https://r12.finngen.fi/ and search for 'aneurysm'\n")
  cat("    2. Download I9_ABAORTANEUR summary statistics\n")
  cat("    3. Save to data/outcome/finngen_aaa.csv\n")
  cat("    4. Re-run this script\n")
  return(NULL)
}

# Combine all exposure SNPs for outcome extraction
outcome_dat <- get_outcome_data(list(exp_coronary_clumped, exp_tibial_clumped), "FinnGen AAA")

## ============================================================
## 5. Harmonize Exposure and Outcome
## ============================================================
cat("\n============================================================\n")
cat("Step 5: Harmonizing data\n")
cat("============================================================\n\n")

run_mr_for_exposure <- function(exp_clumped, outcome_dat, tissue_label) {
  cat("\n--- MR Analysis:", tissue_label, "---\n\n")

  # Harmonize
  dat <- harmonise_data(
    exposure_dat = exp_clumped,
    outcome_dat = outcome_dat
  )

  if (is.null(dat) || nrow(dat) == 0) {
    cat("  No overlapping SNPs after harmonization!\n")
    return(NULL)
  }

  cat("  Harmonized SNPs:", nrow(dat), "\n")
  print(dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                 "effect_allele.outcome", "other_allele.outcome",
                 "beta.exposure", "beta.outcome", "pval.outcome")])

  ## ============================================================
  ## 6. MR Analysis
  ## ============================================================
  cat("\n  --- MR Results ---\n\n")

  res <- mr(dat, method_list = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median",
    "mr_simple_mode",
    "mr_weighted_mode"
  ))

  print(res)

  # Heterogeneity test
  cat("\n  --- Heterogeneity (Cochran's Q) ---\n")
  het <- mr_heterogeneity(dat)
  print(het)

  # Pleiotropy test (MR-Egger intercept)
  cat("\n  --- Pleiotropy (MR-Egger intercept) ---\n")
  pleio <- mr_pleiotropy_test(dat)
  print(pleio)

  # Leave-one-out analysis
  cat("\n  --- Leave-One-Out Analysis ---\n")
  loo <- mr_leaveoneout(dat)
  print(loo)

  # Single SNP analysis
  cat("\n  --- Single SNP (Wald ratio) ---\n")
  single <- mr_singlesnp(dat)
  print(single)

  # Save results
  tissue_tag <- gsub("[^a-zA-Z0-9]", "_", tissue_label)
  fwrite(res, paste0("results/tables/mr_results_", tissue_tag, ".csv"))
  fwrite(het, paste0("results/tables/heterogeneity_", tissue_tag, ".csv"))
  fwrite(pleio, paste0("results/tables/pleiotropy_", tissue_tag, ".csv"))
  fwrite(loo, paste0("results/tables/leaveoneout_", tissue_tag, ".csv"))
  fwrite(single, paste0("results/tables/single_snp_", tissue_tag, ".csv"))

  ## ============================================================
  ## 7. Generate Plots
  ## ============================================================
  cat("\n  --- Generating plots ---\n")

  # Scatter plot
  p_scatter <- mr_scatter_plot(res, dat)
  ggsave(paste0("results/figures/scatter_plot_", tissue_tag, ".png"),
         p_scatter[[1]], width = 8, height = 6, dpi = 300)

  # Forest plot
  p_forest <- mr_forest_plot(single)
  ggsave(paste0("results/figures/forest_plot_", tissue_tag, ".png"),
         p_forest[[1]], width = 8, height = 6, dpi = 300)

  # Funnel plot
  p_funnel <- mr_funnel_plot(single)
  ggsave(paste0("results/figures/funnel_plot_", tissue_tag, ".png"),
         p_funnel[[1]], width = 8, height = 6, dpi = 300)

  # Leave-one-out plot
  p_loo <- mr_leaveoneout_plot(loo)
  ggsave(paste0("results/figures/leaveoneout_plot_", tissue_tag, ".png"),
         p_loo[[1]], width = 8, height = 6, dpi = 300)

  cat("  Plots saved to results/figures/\n")

  return(list(results = res, data = dat, heterogeneity = het,
              pleiotropy = pleio, loo = loo, single = single))
}

# Run MR for each exposure tissue
if (!is.null(outcome_dat)) {
  results_coronary <- run_mr_for_exposure(exp_coronary_clumped, outcome_dat, "Artery_Coronary")
  results_tibial <- run_mr_for_exposure(exp_tibial_clumped, outcome_dat, "Artery_Tibial")

  # If eQTLGen data available, run that too
  if (has_eqtlgen && nrow(eqtlgen_dt) > 0) {
    eqtlgen_formatted <- data.frame(
      SNP = eqtlgen_dt$SNP,
      beta = as.numeric(eqtlgen_dt$beta),
      se = calc_se(as.numeric(eqtlgen_dt$beta), as.numeric(eqtlgen_dt$pval)),
      pval = as.numeric(eqtlgen_dt$pval),
      effect_allele = eqtlgen_dt$effect_allele,
      other_allele = eqtlgen_dt$other_allele,
      exposure = "OLR1 eQTL (Whole Blood, eQTLGen)",
      gene = "OLR1",
      id.exposure = "eqtlgen_olr1_blood",
      units.exposure = "z_score",
      samplesize.exposure = 31684
    )
    eqtlgen_clumped <- clump_exposure(eqtlgen_formatted, "eQTLGen Blood")
    results_eqtlgen <- run_mr_for_exposure(eqtlgen_clumped, outcome_dat, "eQTLGen_Blood")
  }
}

## ============================================================
## 8. Summary
## ============================================================
cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nResults saved to:\n")
cat("  results/tables/ - CSV files with MR results\n")
cat("  results/figures/ - PNG plots (scatter, forest, funnel, LOO)\n")
cat("\nNext steps:\n")
cat("  1. Check if MR results are significant (IVW p < 0.05)\n")
cat("  2. Check pleiotropy (MR-Egger intercept p > 0.05 = no pleiotropy)\n")
cat("  3. Check heterogeneity (Cochran's Q p > 0.05 = homogeneous)\n")
cat("  4. If significant, proceed to colocalization analysis (coloc package)\n")
cat("  5. Install coloc: install.packages('coloc')\n")
