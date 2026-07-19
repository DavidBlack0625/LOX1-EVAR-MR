# Debug coloc results structure
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
})

dat <- fread("D:/LOX1_EVAR_MR/data/coloc/coloc_input_Artery_Tibial.csv")
dat <- dat[eqtl_varbeta > 0 & gwas_varbeta > 0 & eqtl_maf > 0 & eqtl_maf < 1 & gwas_maf > 0 & gwas_maf < 1]

dataset1 <- list(
  snp = dat$snp, beta = dat$eqtl_beta, varbeta = dat$eqtl_varbeta,
  MAF = dat$eqtl_maf, N = dat$eqtl_N[1], type = "quant", sdY = 1
)
dataset2 <- list(
  snp = dat$snp, beta = dat$gwas_beta, varbeta = dat$gwas_varbeta,
  MAF = dat$gwas_maf, N = dat$gwas_N[1], type = "cc", s = dat$gwas_s[1]
)

result <- coloc.abf(dataset1 = dataset1, dataset2 = dataset2)

cat("=== Result names ===\n")
cat(names(result), "\n")

cat("\n=== Summary ===\n")
print(result$summary)

if (!is.null(result$results)) {
  sr <- as.data.frame(result$results)
  cat("\n=== Results data frame ===\n")
  cat("Column names:", names(sr), "\n")
  cat("Column classes:", sapply(sr, class), "\n\n")
  cat("Head (first 3 rows):\n")
  print(head(sr, 3))
  cat("\nStr:\n")
  str(sr)
}
