###############################################################################
## 04_mr_eqtlgen.R
## MR Analysis: OLR1 eQTL (eQTLGen, blood n=31684) -> AAA (FinnGen R12)
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

## ============================================================
## 1. Load eQTLGen Exposure Data
## ============================================================
cat("\n============================================================\n")
cat("MR: OLR1 eQTL (eQTLGen Blood) -> AAA (FinnGen R12)\n")
cat("============================================================\n\n")

eqtlgen <- fread(file.path(BASE_DIR, "data/exposure/eqtlgen_OLR1_blood_eqtl.csv"))
cat("eQTLGen OLR1 SNPs:", nrow(eqtlgen), "\n")
cat("Columns:", colnames(eqtlgen), "\n")

# eQTLGen provides Z-scores; use Z as beta, SE = 1
# (Z = beta/SE, so if we set SE=1, beta=Z preserves the ratio)
eqtlgen$beta <- as.numeric(eqtlgen$beta)  # This is actually Zscore
eqtlgen$pval <- as.numeric(eqtlgen$pval)

# Format for TwoSampleMR
exp_dat <- as.data.frame(eqtlgen) %>%
  mutate(
    se = 1,  # Z-score based, SE = 1
    chr = as.character(NA),
    pos = as.integer(NA)
  ) %>%
  format_data(
    .,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  ) %>%
  mutate(
    exposure = "OLR1 eQTL (Whole Blood, eQTLGen)",
    gene = "OLR1",
    id.exposure = "eqtlgen_olr1_blood",
    units.exposure = "z_score",
    samplesize.exposure = 31684
  )

cat("\neQTLGen exposure formatted:", nrow(exp_dat), "SNPs\n")
print(exp_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                   "beta.exposure", "se.exposure", "pval.exposure")])

## ============================================================
## 2. Load Outcome Data
## ============================================================
cat("\n--- Loading outcome data ---\n")
outcome_raw <- fread(file.path(BASE_DIR, "data/outcome/finngen_aaa_eqtlgen_matched.csv"))
cat("Matched SNPs:", nrow(outcome_raw), "\n")

outcome_dat <- as.data.frame(outcome_raw) %>%
  format_data(
    .,
    type = "outcome",
    snp_col = "matched_snp",
    beta_col = "beta",
    se_col = "sebeta",
    pval_col = "pval",
    effect_allele_col = "alt",
    other_allele_col = "ref"
  ) %>%
  mutate(
    outcome = "Abdominal Aortic Aneurysm (FinnGen R12)",
    id.outcome = "finngen_r12_aaa",
    units.outcome = "log_odds_ratio",
    ncase.outcome = 4439,
    ncontrol.outcome = 463106,
    samplesize.outcome = 467545
  )

## ============================================================
## 3. LD Pruning (position not available, use SNP-based)
## ============================================================
cat("\n--- LD Pruning ---\n")
# eQTLGen doesn't have positions, so try API clumping
exp_clumped <- tryCatch({
  clump_data(exp_dat, clump_kb = 10000, clump_r2 = 0.001,
             clump_p1 = 1, clump_p2 = 1, pop = "EUR")
}, error = function(e) {
  cat("  API clumping failed, using all SNPs (may have LD)\n")
  exp_dat
})
cat("  SNPs after clumping:", nrow(exp_clumped), "\n")

## ============================================================
## 4. Harmonize and Run MR
## ============================================================
cat("\n--- Harmonization ---\n")
dat <- harmonise_data(exp_clumped, outcome_dat)
cat("Harmonized SNPs:", nrow(dat), "\n\n")
print(dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
              "effect_allele.outcome", "other_allele.outcome",
              "beta.exposure", "beta.outcome", "se.outcome", "pval.outcome")])

cat("\n--- MR Results ---\n")
res <- mr(dat, method_list = c(
  "mr_ivw", "mr_egger_regression", "mr_weighted_median",
  "mr_simple_mode", "mr_weighted_mode"
))
print(res)

cat("\n--- Heterogeneity ---\n")
het <- mr_heterogeneity(dat)
print(het)

cat("\n--- Pleiotropy ---\n")
pleio <- mr_pleiotropy_test(dat)
print(pleio)

cat("\n--- Leave-One-Out ---\n")
loo <- mr_leaveoneout(dat)
print(loo)

cat("\n--- Single SNP ---\n")
single <- mr_singlesnp(dat)
print(single)

## ============================================================
## 5. Save Results and Plots
## ============================================================
fwrite(res,    file.path(TBL_DIR, "mr_results_eQTLGen_Blood.csv"))
fwrite(het,    file.path(TBL_DIR, "heterogeneity_eQTLGen_Blood.csv"))
fwrite(pleio,  file.path(TBL_DIR, "pleiotropy_eQTLGen_Blood.csv"))
fwrite(loo,    file.path(TBL_DIR, "leaveoneout_eQTLGen_Blood.csv"))
fwrite(single, file.path(TBL_DIR, "single_snp_eQTLGen_Blood.csv"))

p1 <- mr_scatter_plot(res, dat)
if (length(p1) > 0) ggsave(file.path(FIG_DIR, "scatter_eQTLGen_Blood.png"), p1[[1]], width = 8, height = 6, dpi = 300)
p2 <- mr_forest_plot(single)
if (length(p2) > 0) ggsave(file.path(FIG_DIR, "forest_eQTLGen_Blood.png"), p2[[1]], width = 8, height = 6, dpi = 300)
p3 <- mr_funnel_plot(single)
if (length(p3) > 0) ggsave(file.path(FIG_DIR, "funnel_eQTLGen_Blood.png"), p3[[1]], width = 8, height = 6, dpi = 300)
p4 <- mr_leaveoneout_plot(loo)
if (length(p4) > 0) ggsave(file.path(FIG_DIR, "loo_eQTLGen_Blood.png"), p4[[1]], width = 8, height = 6, dpi = 300)

## ============================================================
## 6. Summary
## ============================================================
cat("\n\n============================================================\n")
cat("eQTLGen Blood MR - SUMMARY\n")
cat("============================================================\n\n")

for (i in 1:nrow(res)) {
  cat(sprintf("  %-25s  OR = %.4f  95%%CI = [%.4f, %.4f]  p = %.4e\n",
              res$method[i],
              exp(res$b[i]),
              exp(res$b[i] - 1.96 * res$se[i]),
              exp(res$b[i] + 1.96 * res$se[i]),
              res$pval[i]))
}
