#' Test case for the situation where the --window provided to `apex store` was too small to capture
#' all possible covariance values between variants. The dataset used here is purely synthetic; it has
#' two signals, one at each end of the "gene", and there is no covariance information available between them.
#' In this situation, we can either assume the missing covariance is 0, and carry on and fit the model, or
#' truly treat the covariance as missing and set the result to a sentinel value (p-value of 1.) If your --window
#' is quite large, it may be acceptable to assume long distance covariance is 0.
#'
#' This test case function tests the situation where covariance is assumed to be 0 if it is missing.
test_that("running with missing long range covariance set to 0", {
  # Temporarily required while Apex2R is not an actual R package
  withr::local_dir("../../")
  source("Apex2R.r")
  source("cond_p.R")

  # This is an unfortunate hack until the code can be refactored & cleaned up into a package
  # Ideally this would be done with a command line argument and passed all the way down rather than
  # relying on global variables like this
  withr::local_options(assumeMissingCovZero=TRUE)

  ssm = metaSumStats(
    "tests/testthat/fixtures/cov-window-too-small/cov_window_too_small_s1",
    "tests/testthat/fixtures/cov-window-too-small/cov_window_too_small_s2"
  )

  meta_stepwise = "tests/testthat/fixtures/cov-window-too-small/cov_window_too_small_meta.cis_meta.stepwise_het.tsv"

  res = data.frame(fread(meta_stepwise, stringsAsFactors=F))
  colnames(res)[1] = "gene"

  res_list = split(res$variant, res$gene)
  num_sig = lapply(1:length(res_list), function(x) {
    length(res_list[[x]])
  })
  res_list = res_list[num_sig > 1]

  # will result in error: x$.self$finalize(): attempt to apply non-function - Li advised to ignore this issue
  out = lapply(1:length(res_list), loop_gene_tryCatch, rlist=res_list, meta=ssm)

  qtl = do.call(rbind,out)
  qtl$variant = rownames(qtl)

  # Even though there is no covariance available between 1_2_A_G and 1_12_A_G, because we are assuming missing cov
  # is set to 0, the p-value can be calculated
  expect_equal(
    -log10(subset(qtl, (not_c == "1_12_A_G") & (variant == "1_12_A_G"))$pval),
    262.799,
    tolerance = 1e-3
  )
})

#' Similar to the test case above, except that missing covariance is set to NA, and tests with missing covariance
#' then are not performed and set to have a sentinel p-value of 1.
test_that("running with missing long range covariance set to NA", {
  withr::local_dir("../../")
  source("Apex2R.r")
  source("cond_p.R")

  withr::local_options(assumeMissingCovZero=FALSE)

  ssm = metaSumStats(
    "tests/testthat/fixtures/cov-window-too-small/cov_window_too_small_s1",
    "tests/testthat/fixtures/cov-window-too-small/cov_window_too_small_s2"
  )

  meta_stepwise = "tests/testthat/fixtures/cov-window-too-small/cov_window_too_small_meta.cis_meta.stepwise_het.tsv"

  res = data.frame(fread(meta_stepwise, stringsAsFactors=F))
  colnames(res)[1] = "gene"

  res_list = split(res$variant, res$gene)
  num_sig = lapply(1:length(res_list), function(x) {
    length(res_list[[x]])
  })
  res_list = res_list[num_sig > 1]

  # will result in error: x$.self$finalize(): attempt to apply non-function - Li advised to ignore this issue
  out = lapply(1:length(res_list), loop_gene_tryCatch, rlist=res_list, meta=ssm)

  qtl = do.call(rbind,out)
  qtl$variant = rownames(qtl)

  # In the test data, because the vcov was generated with a --window that is too small, there should be
  # no covariance available between 1_2_A_G and 1_12_A_G, so this p-value should be forced to 1
  expect_equal(
    subset(qtl, (not_c == "1_12_A_G") & (variant == "1_12_A_G"))$pval,
    1
  )
})