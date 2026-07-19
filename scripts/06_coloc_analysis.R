# ============================================================
# Colocalization Analysis: OLR1 eQTL vs AAA GWAS
# Tests whether OLR1 eQTL and AAA GWAS share the same causal variant
# ============================================================

rm(list = ls())
BASE_DIR <- "D:/LOX1_EVAR_MR"
COLOC_DIR <- file.path(BASE_DIR, "data/coloc")
RESULTS_DIR <- file.path(BASE_DIR, "results")
dir.create(file.path(RESULTS_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

cat("coloc version:", as.character(packageVersion("coloc")), "\n\n")

# ============================================================
# Function to run coloc.abf for one eQTL source
# ============================================================
run_coloc <- function(matched_file, tissue_name) {
  cat("=" , rep("=", 58), "\n", sep = "")
  cat(sprintf("Colocalization: OLR1 eQTL (%s) vs AAA GWAS\n", tissue_name))
  cat("=" , rep("=", 58), "\n", sep = "")

  # Read matched data
  dat <- fread(matched_file)
  cat(sprintf("Input: %d matched SNPs\n", nrow(dat)))

  # Filter: remove SNPs with invalid values
  dat <- dat[eqtl_varbeta > 0 & gwas_varbeta > 0]
  dat <- dat[eqtl_maf > 0 & eqtl_maf < 1]
  dat <- dat[gwas_maf > 0 & gwas_maf < 1]
  dat <- dat[!is.na(eqtl_beta) & !is.na(gwas_beta)]
  cat(sprintf("After filtering: %d SNPs\n", nrow(dat)))

  if (nrow(dat) < 3) {
    cat("  WARNING: Too few SNPs for meaningful coloc analysis. Skipping.\n\n")
    return(NULL)
  }

  # Prepare dataset1 (eQTL, quantitative trait)
  dataset1 <- list(
    snp      = dat$snp,
    beta     = dat$eqtl_beta,
    varbeta  = dat$eqtl_varbeta,
    MAF      = dat$eqtl_maf,
    N        = dat$eqtl_N[1],
    type     = "quant",
    sdY      = 1
  )

  # Prepare dataset2 (GWAS, case-control)
  dataset2 <- list(
    snp      = dat$snp,
    beta     = dat$gwas_beta,
    varbeta  = dat$gwas_varbeta,
    MAF      = dat$gwas_maf,
    N        = dat$gwas_N[1],
    type     = "cc",
    s        = dat$gwas_s[1]
  )

  # Run coloc.abf
  cat("Running coloc.abf...\n")
  result <- tryCatch({
    coloc.abf(dataset1 = dataset1, dataset2 = dataset2)
  }, error = function(e) {
    cat("  ERROR in coloc.abf:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(result)) return(NULL)

  # Extract posterior probabilities
  pp <- result$summary
  cat("\n--- Posterior Probabilities ---\n")
  cat(sprintf("  PP.H0 (no association):           %.4f\n", pp["PP.H0.abf"]))
  cat(sprintf("  PP.H1 (eQTL only):                %.4f\n", pp["PP.H1.abf"]))
  cat(sprintf("  PP.H2 (GWAS only):                %.4f\n", pp["PP.H2.abf"]))
  cat(sprintf("  PP.H3 (both, diff causal vars):   %.4f\n", pp["PP.H3.abf"]))
  cat(sprintf("  PP.H4 (colocalization):           %.4f\n", pp["PP.H4.abf"]))
  cat("\n")

  if (pp["PP.H4.abf"] > 0.7) {
    cat("  >>> STRONG evidence for colocalization (PP.H4 > 0.7) <<<\n")
  } else if (pp["PP.H4.abf"] > 0.5) {
    cat("  >> Moderate evidence for colocalization (0.5 < PP.H4 < 0.7) <<\n")
  } else {
    cat("  > No strong evidence for colocalization (PP.H4 < 0.5) <\n")
  }

  # Save summary
  summary_df <- data.frame(
    tissue = tissue_name,
    n_snps = nrow(dat),
    PP.H0 = as.numeric(pp["PP.H0.abf"]),
    PP.H1 = as.numeric(pp["PP.H1.abf"]),
    PP.H2 = as.numeric(pp["PP.H2.abf"]),
    PP.H3 = as.numeric(pp["PP.H3.abf"]),
    PP.H4 = as.numeric(pp["PP.H4.abf"])
  )
  fwrite(summary_df,
         file.path(RESULTS_DIR, "tables", paste0("coloc_summary_", tissue_name, ".csv")))

  # SNP-level results
  if (!is.null(result$results)) {
    snp_results <- as.data.frame(result$results)
    fwrite(snp_results,
           file.path(RESULTS_DIR, "tables", paste0("coloc_snp_", tissue_name, ".csv")))
    cat(sprintf("\n  SNP-level results saved (%d SNPs)\n", nrow(snp_results)))

    # --- Sensitivity plot ---
    png(file.path(RESULTS_DIR, "figures",
                  paste0("coloc_sensitivity_", tissue_name, ".png")),
        width = 2400, height = 1800, res = 300)
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

    # Merge position from original data
    sr <- merge(snp_results, as.data.frame(dat[, .(snp, pos)]), by = "snp", all.x = TRUE)
    sr <- sr[order(sr$pos), ]

    # Compute p-values from Z-scores (coloc v5 uses z.df1/z.df2)
    if ("z.df1" %in% names(sr)) sr$pval_eqtl <- 2 * pnorm(-abs(as.numeric(sr$z.df1)))
    if ("z.df2" %in% names(sr)) sr$pval_gwas <- 2 * pnorm(-abs(as.numeric(sr$z.df2)))

    # Plot 1: Regional association plot (-log10 p)
    ylim_max <- max(c(-log10(sr$pval_eqtl), -log10(sr$pval_gwas)), na.rm = TRUE) * 1.1
    if (ylim_max < 1) ylim_max <- 1

    plot(sr$pos, -log10(sr$pval_eqtl),
         pch = 20, col = "blue",
         xlab = "Position (chr12)", ylab = "-log10(p)",
         main = paste0("Regional Association (", tissue_name, ")"),
         cex = 0.8, ylim = c(0, ylim_max))
    points(sr$pos, -log10(sr$pval_gwas),
           pch = 18, col = "red", cex = 0.8)
    legend("topright", c("eQTL", "GWAS"), col = c("blue", "red"), pch = c(20, 18))

    # Plot 2: Per-SNP PP.H4
    if ("SNP.PP.H4" %in% names(sr)) {
      plot(sr$pos, as.numeric(sr$SNP.PP.H4),
           pch = 20, col = "darkgreen",
           xlab = "Position (chr12)", ylab = "SNP PP.H4",
           main = "Per-SNP Colocalization Probability",
           cex = 0.8)
      abline(h = 0.01, lty = 2, col = "gray")
    }

    # Plot 3: Z-scores comparison
    plot(as.numeric(sr$z.df1), as.numeric(sr$z.df2),
         pch = 20, col = "purple", cex = 0.8,
         xlab = "eQTL Z-score", ylab = "GWAS Z-score",
         main = "Z-score Comparison")
    abline(h = 0, lty = 2, col = "gray")
    abline(v = 0, lty = 2, col = "gray")

    # Plot 4: Posterior probability bar chart
    pp_vals <- c(pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"],
                 pp["PP.H3.abf"], pp["PP.H4.abf"])
    names(pp_vals) <- c("H0", "H1", "H2", "H3", "H4")
    bp <- barplot(pp_vals,
            col = c("gray80", "lightblue", "lightsalmon", "khaki", "forestgreen"),
            ylim = c(0, 1.15),
            ylab = "Posterior Probability",
            main = "Colocalization Posterior Probabilities")
    abline(h = 0.7, lty = 2, col = "red")
    text(bp[5], 1.08, sprintf("PP.H4=%.3f", pp["PP.H4.abf"]), cex = 0.7, font = 2)

    dev.off()
    cat("  Sensitivity plot saved\n")
  }

  cat("\n")
  return(summary_df)
}

# ============================================================
# Run coloc for each eQTL source
# ============================================================
results_list <- list()

# 1. GTEx Artery_Tibial (49 SNPs)
file1 <- file.path(COLOC_DIR, "coloc_input_Artery_Tibial.csv")
if (file.exists(file1)) {
  results_list[["Artery_Tibial"]] <- run_coloc(file1, "Artery_Tibial")
}

# 2. GTEx Artery_Coronary (9 SNPs)
file2 <- file.path(COLOC_DIR, "coloc_input_Artery_Coronary.csv")
if (file.exists(file2)) {
  results_list[["Artery_Coronary"]] <- run_coloc(file2, "Artery_Coronary")
}

# 3. eQTLGen Blood (31 SNPs)
file3 <- file.path(COLOC_DIR, "coloc_input_eQTLGen_Blood.csv")
if (file.exists(file3)) {
  results_list[["eQTLGen_Blood"]] <- run_coloc(file3, "eQTLGen_Blood")
}

# ============================================================
# Combined summary
# ============================================================
cat("\n")
cat("=" , rep("=", 58), "\n", sep = "")
cat("COMBINED COLOCALIZATION SUMMARY")
cat("\n")
cat("=" , rep("=", 58), "\n", sep = "")

combined <- do.call(rbind, results_list[!sapply(results_list, is.null)])
if (nrow(combined) > 0) {
  print(combined)
  fwrite(combined, file.path(RESULTS_DIR, "tables", "coloc_summary_all.csv"))
  cat("\nSaved to: results/tables/coloc_summary_all.csv\n")
}

cat("\nDone!\n")
