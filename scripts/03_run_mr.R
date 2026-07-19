###############################################################################
## 03_run_mr.R  (v2 - fixed)
## Complete MR Analysis: OLR1 eQTL -> AAA (FinnGen R12)
###############################################################################

rm(list = ls())
BASE_DIR <- "D:/LOX1_EVAR_MR"

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

FIG_DIR <- file.path(BASE_DIR, "results/figures")
TBL_DIR <- file.path(BASE_DIR, "results/tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

calc_se <- function(beta, pval) {
  z <- abs(qnorm(pval / 2))
  abs(beta) / ifelse(z == 0, 1e-10, z)
}

## ============================================================
## 1. Load & Format Exposure Data
## ============================================================
cat("\n============================================================\n")
cat("Step 1: Loading and formatting exposure data\n")
cat("============================================================\n\n")

coronary <- fread(file.path(BASE_DIR, "data/exposure/gtex_v8_OLR1_Artery_Coronary_eqtl.csv"))
tibial   <- fread(file.path(BASE_DIR, "data/exposure/gtex_v8_OLR1_Artery_Tibial_eqtl.csv"))

cat("Artery_Coronary:", nrow(coronary), "SNPs | positions:", min(coronary$pos), "-", max(coronary$pos), "\n")
cat("Artery_Tibial:", nrow(tibial), "SNPs | positions:", min(tibial$pos), "-", max(tibial$pos), "\n")

# Use TwoSampleMR::format_data to ensure correct column names
format_exp <- function(dt, tissue) {
  dt$se <- calc_se(dt$beta, dt$pval)
  dt <- as.data.frame(dt)
  format_data(
    dt,
    type            = "exposure",
    snp_col         = "SNP",
    beta_col        = "beta",
    se_col          = "se",
    pval_col        = "pval",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    chr_col         = "chr",
    pos_col         = "pos"
  ) |> mutate(
    exposure = paste0("OLR1 eQTL (", tissue, ")"),
    gene = "OLR1",
    id.exposure = paste0("olr1_", tolower(gsub("Artery_", "", tissue))),
    units.exposure = "normalized_effect_size"
  )
}

exp_coronary <- format_exp(coronary, "Artery_Coronary")
exp_tibial   <- format_exp(tibial, "Artery_Tibial")

cat("\nArtery_Coronary formatted (", nrow(exp_coronary), "SNPs):\n")
print(exp_coronary[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                        "beta.exposure", "se.exposure", "pval.exposure")])

## ============================================================
## 2. Load & Format Outcome Data
## ============================================================
cat("\n============================================================\n")
cat("Step 2: Loading outcome data (FinnGen R12 AAA)\n")
cat("============================================================\n\n")

outcome_raw <- fread(file.path(BASE_DIR, "data/outcome/finngen_aaa_matched.csv"))
cat("Matched SNPs:", nrow(outcome_raw), "\n")

outcome_dat <- format_data(
  as.data.frame(outcome_raw),
  type             = "outcome",
  snp_col          = "matched_snp",
  beta_col         = "beta",
  se_col           = "sebeta",
  pval_col         = "pval",
  effect_allele_col = "alt",
  other_allele_col  = "ref",
  chr_col          = "#chrom",
  pos_col          = "pos"
) |> mutate(
  outcome = "Abdominal Aortic Aneurysm (FinnGen R12)",
  id.outcome = "finngen_r12_aaa",
  units.outcome = "log_odds_ratio",
  ncase.outcome = 4439,
  ncontrol.outcome = 463106,
  samplesize.outcome = 4439 + 463106
)

cat("Outcome formatted (", nrow(outcome_dat), "SNPs):\n")
print(head(outcome_dat[, c("SNP", "effect_allele.outcome", "other_allele.outcome",
                             "beta.outcome", "se.outcome", "pval.outcome")], 10))

## ============================================================
## 3. LD Pruning
## ============================================================
cat("\n============================================================\n")
cat("Step 3: LD pruning\n")
cat("============================================================\n\n")

# For cis-eQTLs, all SNPs are in a small region (~20-50 kb)
# Use r2 < 0.1 based on distance (< 10 kb = likely in LD)
# This is a crude proxy; proper clumping requires PLINK or API token

prune_local <- function(exp_dat, min_dist_kb = 10) {
  df <- exp_dat[order(exp_dat$pval.exposure), ]
  selected <- c()
  for (i in 1:nrow(df)) {
    if (length(selected) == 0) {
      selected <- c(selected, i)
    } else {
      min_dist <- min(abs(df$pos.exposure[i] - df$pos.exposure[selected]))
      if (min_dist >= min_dist_kb * 1000) selected <- c(selected, i)
    }
  }
  df[selected, ]
}

# Try API clumping first
try_api <- function(exp_dat) {
  tryCatch({
    clumped <- clump_data(exp_dat, clump_kb = 10000, clump_r2 = 0.001,
                          clump_p1 = 1, clump_p2 = 1, pop = "EUR")
    cat("  Clumped via API:", nrow(clumped), "SNPs\n")
    clumped
  }, error = function(e) {
    cat("  API clumping failed, using position-based (10 kb threshold)\n")
    NULL
  })
}

cat("Artery_Coronary:\n")
coronary_clumped <- try_api(exp_coronary)
if (is.null(coronary_clumped)) coronary_clumped <- prune_local(exp_coronary, 10)
cat("  ", nrow(coronary_clumped), "independent SNPs after pruning\n")
print(coronary_clumped[, c("SNP", "pos.exposure", "beta.exposure", "pval.exposure")])

cat("\nArtery_Tibial:\n")
tibial_clumped <- try_api(exp_tibial)
if (is.null(tibial_clumped)) tibial_clumped <- prune_local(exp_tibial, 10)
cat("  ", nrow(tibial_clumped), "independent SNPs after pruning\n")
print(tibial_clumped[, c("SNP", "pos.exposure", "beta.exposure", "pval.exposure")])

## ============================================================
## 4. Harmonize and Run MR
## ============================================================
cat("\n============================================================\n")
cat("Step 4: Harmonization and MR Analysis\n")
cat("============================================================\n")

run_mr <- function(exp_clumped, outcome_dat, tissue_label) {
  cat("\n========== MR:", tissue_label, "==========\n\n")

  dat <- harmonise_data(exp_clumped, outcome_dat)
  if (is.null(dat) || nrow(dat) == 0) {
    cat("  No overlapping SNPs!\n")
    return(NULL)
  }

  cat("Harmonized SNPs:", nrow(dat), "\n\n")
  print(dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                 "effect_allele.outcome", "other_allele.outcome",
                 "beta.exposure", "beta.outcome", "se.outcome", "pval.outcome")])

  # Determine which MR methods to use based on SNP count
  n_snps <- nrow(dat)
  if (n_snps == 1) {
    methods <- "mr_wald_ratio"
    cat("\n  Only 1 SNP available - using Wald ratio only\n")
  } else if (n_snps == 2) {
    methods <- c("mr_ivw", "mr_wald_ratio")
    cat("\n  2 SNPs available - using IVW and Wald ratio\n")
  } else {
    methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                 "mr_simple_mode", "mr_weighted_mode")
  }

  cat("\n--- MR Results ---\n")
  res <- mr(dat, method_list = methods)
  print(res)

  if (n_snps >= 3) {
    cat("\n--- Heterogeneity (Cochran's Q) ---\n")
    het <- mr_heterogeneity(dat)
    print(het)

    cat("\n--- Pleiotropy (MR-Egger intercept) ---\n")
    pleio <- mr_pleiotropy_test(dat)
    print(pleio)

    cat("\n--- Leave-One-Out ---\n")
    loo <- mr_leaveoneout(dat)
    print(loo)
  } else {
    het <- NULL; pleio <- NULL; loo <- NULL
  }

  cat("\n--- Single SNP ---\n")
  single <- mr_singlesnp(dat)
  print(single)

  # Save
  tag <- gsub("[^a-zA-Z0-9]", "_", tissue_label)
  fwrite(res, file.path(TBL_DIR, paste0("mr_results_", tag, ".csv")))
  fwrite(single, file.path(TBL_DIR, paste0("single_snp_", tag, ".csv")))
  if (!is.null(het)) fwrite(het, file.path(TBL_DIR, paste0("heterogeneity_", tag, ".csv")))
  if (!is.null(pleio)) fwrite(pleio, file.path(TBL_DIR, paste0("pleiotropy_", tag, ".csv")))
  if (!is.null(loo)) fwrite(loo, file.path(TBL_DIR, paste0("leaveoneout_", tag, ".csv")))

  # Plots
  cat("\n--- Generating plots ---\n")
  p1 <- mr_scatter_plot(res, dat)
  if (length(p1) > 0) ggsave(file.path(FIG_DIR, paste0("scatter_", tag, ".png")), p1[[1]], width = 8, height = 6, dpi = 300)
  p2 <- mr_forest_plot(single)
  if (length(p2) > 0) ggsave(file.path(FIG_DIR, paste0("forest_", tag, ".png")), p2[[1]], width = 8, height = 6, dpi = 300)
  p3 <- mr_funnel_plot(single)
  if (length(p3) > 0) ggsave(file.path(FIG_DIR, paste0("funnel_", tag, ".png")), p3[[1]], width = 8, height = 6, dpi = 300)
  if (!is.null(loo)) {
    p4 <- mr_leaveoneout_plot(loo)
    if (length(p4) > 0) ggsave(file.path(FIG_DIR, paste0("loo_", tag, ".png")), p4[[1]], width = 8, height = 6, dpi = 300)
  }
  cat("  Plots saved\n")

  return(list(res = res, dat = dat, het = het, pleio = pleio, loo = loo))
}

res_coronary <- run_mr(coronary_clumped, outcome_dat, "Artery_Coronary")
res_tibial   <- run_mr(tibial_clumped, outcome_dat, "Artery_Tibial")

## ============================================================
## 5. Summary
## ============================================================
cat("\n\n============================================================\n")
cat("              ANALYSIS COMPLETE - SUMMARY\n")
cat("============================================================\n\n")

cat("Exposure: OLR1 (LOX-1) cis-eQTL from GTEx v8\n")
cat("Outcome:  Abdominal Aortic Aneurysm (FinnGen R12, 4439 cases / 463106 controls)\n\n")

for (tissue in c("Artery_Coronary", "Artery_Tibial")) {
  tag <- gsub("[^a-zA-Z0-9]", "_", tissue)
  res_file <- file.path(TBL_DIR, paste0("mr_results_", tag, ".csv"))
  if (file.exists(res_file)) {
    res <- fread(res_file)
    cat(paste0("--- ", tissue, " ---\n"))
    for (i in 1:nrow(res)) {
      cat(sprintf("  %-25s  OR = %.4f  95%%CI = [%.4f, %.4f]  p = %.4e\n",
                  res$method[i],
                  exp(res$b[i]),
                  exp(res$b[i] - 1.96 * res$se[i]),
                  exp(res$b[i] + 1.96 * res$se[i]),
                  res$pval[i]))
    }
    cat("\n")
  }
}

cat("Tables:", TBL_DIR, "\n")
cat("Plots: ", FIG_DIR, "\n")
cat("\nNOTE: Position-based LD pruning (10 kb threshold) was used.\n")
cat("      For publication, use proper LD clumping (register at api.opengwas.io).\n")
