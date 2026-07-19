# Check R packages
pkgs <- c('TwoSampleMR', 'ieugwasr', 'coloc', 'MRPRESSO', 'ggplot2', 'data.table')
for (p in pkgs) {
  cat(p, ':', ifelse(requireNamespace(p, quietly=TRUE), 'INSTALLED', 'NOT FOUND'), '\n')
}

# Try OpenGWAS connection
cat('\n--- Testing OpenGWAS API ---\n')
tryCatch({
  res <- ieugwasr::gwasinfo('ieu-b-30')
  if (nrow(res) > 0) {
    cat('OpenGWAS accessible! ieu-b-30:', res$trait, '\n')
    cat('n:', res$ncase + res$ncontrol, '\n')
  } else {
    cat('ieu-b-30 not found\n')
  }
}, error = function(e) {
  cat('OpenGWAS error:', conditionMessage(e), '\n')
  # Try with httr SSL config
  cat('\nTrying with SSL workaround...\n')
  tryCatch({
    options(httr_config = httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))
    httr::set_config(httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))
    res2 <- ieugwasr::gwasinfo('ieu-b-30')
    if (nrow(res2) > 0) {
      cat('OpenGWAS accessible with SSL workaround! ieu-b-30:', res2$trait, '\n')
    }
  }, error = function(e2) {
    cat('Still failed:', conditionMessage(e2), '\n')
  })
})
