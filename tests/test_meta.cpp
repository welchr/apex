#include <gtest/gtest.h>
#include "setOptions.hpp"
#include "metaAnalysis.hpp"
#include <string>
#include <vector>
#include "meta_stepwise_reader.h"
#include "meta_svar_reader.h"
using namespace std;

/**
 * "Very small" meaning it requires long double precision or handling in log-scale.
 */
TEST(MetaTest, StepwiseConditionalWithVerySmallPvalue) {
  global_opts::reset();
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    "data/test_output.pval0",           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    0.05,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    false,        // use_ds (dosage)
    false,        // trim gene ids
    0.05,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,          // stepwise_marginal_thresh
    true           // write out log p-values
  );

  vector<string> meta_prefixes = {"data/pval0"};
  string region = "";
  cis_meta_data meta_dt(meta_prefixes, region);
  meta_dt.conditional_analysis();
  auto reader_truth = StepwiseReader("data/pval0.cis_meta.stepwise.tsv");
  auto reader_test = StepwiseReader("data/test_output.pval0.cis_meta.stepwise.tsv");
  ASSERT_TRUE(reader_truth == reader_test);
}

/**
 * Test the conditional_het() function with a very significant p-value (beyond double precision)
 */
TEST(MetaTest, StepwiseHetConditionalWithVerySmallPvalue) {
  global_opts::reset();
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    "data/test_output.pval0",           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    0.05,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    false,        // use_ds (dosage)
    false,        // trim gene ids
    0.05,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,          // stepwise_marginal_thresh
    true           // write out log p-values
  );

  vector<string> meta_prefixes = {"data/pval0"};
  string region = "";
  cis_meta_data meta_dt(meta_prefixes, region);
  meta_dt.conditional_analysis_het();
  auto reader_truth = StepwiseReader("data/pval0.cis_meta.stepwise_het.tsv");
  auto reader_test = StepwiseReader("data/test_output.pval0.cis_meta.stepwise_het.tsv");
  ASSERT_TRUE(reader_truth == reader_test);
}

/**
 * This test verifies that stepwise conditional analysis works correctly with simple input data (no wild p-values.)
 */
TEST(MetaTest, StepwiseConditionalRegressionTest) {
  global_opts::reset();
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    "data/test_output.simple_pvalue",           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    0.05,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    false,        // use_ds (dosage)
    false,        // trim gene ids
    0.05,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,         // stepwise_marginal_thresh
    false         // write out log p-values
  );

  vector<string> meta_prefixes = {"data/simple_pvalue"};
  string region = "";
  cis_meta_data meta_dt(meta_prefixes, region);
  meta_dt.conditional_analysis();

  // This is generated from latest stable version of APEX, commit 63b605e00766141693a69eee6d630ad8ec3b9de0
  auto reader_truth = StepwiseReader("data/simple_pvalue.cis_meta.stepwise.tsv");

  // This is the file generated by the test case
  auto reader_test = StepwiseReader("data/test_output.simple_pvalue.cis_meta.stepwise.tsv");
  ASSERT_TRUE(reader_truth == reader_test);
}

TEST(MetaTest, StepwiseConditionalHetRegressionTest) {
  global_opts::reset();
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    "data/test_output.simple_pvalue",           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    0.05,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    false,        // use_ds (dosage)
    false,        // trim gene ids
    0.05,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,         // stepwise_marginal_thresh
    false         // write out log p-values
  );

  vector<string> meta_prefixes = {"data/simple_pvalue"};
  string region = "";
  cis_meta_data meta_dt(meta_prefixes, region);
  meta_dt.conditional_analysis_het();

  // This is generated from latest stable version of APEX, commit 63b605e00766141693a69eee6d630ad8ec3b9de0
  auto reader_truth = StepwiseReader("data/simple_pvalue.cis_meta.stepwise_het.tsv");

  // This is the file generated by the test case
  auto reader_test = StepwiseReader("data/test_output.simple_pvalue.cis_meta.stepwise_het.tsv");
  ASSERT_TRUE(reader_truth == reader_test);
}

TEST(MetaTest, LogToString) {
  // Easy case
  double val1 = -16.8112428;
  string result1 = log_to_string(val1);
  ASSERT_EQ(result1, "5.000000e-8");

  // NaN
  double val2 = NAN;
  string result2 = log_to_string(val2);
  ASSERT_EQ(result2, "nan");

  // Infs
  double pos_inf = numeric_limits<double>::infinity();
  double neg_inf = -numeric_limits<double>::infinity();
  string result_posinf = log_to_string(pos_inf);
  string result_neginf = log_to_string(neg_inf);
  ASSERT_EQ(result_posinf, "inf");
  ASSERT_EQ(result_neginf, "-inf");
}

TEST(MetaTest, StringToLog) {
  // Easy case
  std::string test1 = "3.12e-3";
  double result1 = string_to_log(test1);
  ASSERT_NEAR(result1, log(3.12e-3), 1e-6);

  // Empty string
  std::string test2 = "";
  double result2 = string_to_log(test2);
  ASSERT_TRUE(std::isnan(result2));

  // Simple number
  std::string test3 = "0.02397462";
  double result3 = string_to_log(test3);
  ASSERT_NEAR(result3, log(0.02397462), 1e-6);

  // No decimal
  std::string test4 = "4e-10";
  double result4 = string_to_log(test4);
  ASSERT_NEAR(result4, log(4e-10), 1e-6);
}

/**
 * Regression test to check for nan pvalue caused by wrong covariance matrix used in stepwise procedure
 */
TEST(MetaTest, RegressionStepwiseCovFlip) {
  global_opts::reset();
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    "data/test_output.check_joint_cov_flips",           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    1e-6,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    true,        // use_ds (dosage)
    false,        // trim gene ids
    1e-6,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,         // stepwise_marginal_thresh
    false         // write out log p-values
  );

  vector<string> meta_prefixes = {"data/check_joint_cov_flips"};
  string region = "";
  cis_meta_data meta_dt(meta_prefixes, region);
  meta_dt.conditional_analysis();

  // Generated after fix in 909e885841b
  auto reader_truth = StepwiseReader("data/check_joint_cov_flips.cis_meta.stepwise.tsv");

  // This is the file generated by the test case
  auto reader_test = StepwiseReader("data/test_output.check_joint_cov_flips.cis_meta.stepwise.tsv");
  ASSERT_TRUE(reader_truth == reader_test);
}

/**
 * Test case that covers the backward selection code in meta --stepwise.
 */
TEST(MetaTest, StepwiseBackward) {
  global_opts::reset();
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    "data/test_output.bkward",           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    0.05,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    true,        // use_ds (dosage)
    false,        // trim gene ids
    0.05,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,         // stepwise_marginal_thresh
    false         // write out log p-values
  );

  vector<string> meta_prefixes = {"data/bkward"};
  string region = "";
  cis_meta_data meta_dt(meta_prefixes, region);
  meta_dt.conditional_analysis();

  auto reader_truth = StepwiseReader("data/bkward.cis_meta.stepwise.tsv");
  auto reader_test = StepwiseReader("data/test_output.bkward.cis_meta.stepwise.tsv");
  ASSERT_TRUE(reader_truth == reader_test);
}

TEST(MetaTest, SingleVarMetaVerySmallPvalue) {
  global_opts::reset();
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    "data/test_output.pval0",           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    0.05,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    false,        // use_ds (dosage)
    false,        // trim gene ids
    0.05,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,          // stepwise_marginal_thresh
    true           // write out log p-values
  );

  vector<string> meta_prefixes = {"data/pval0"};
  string region = "";
  cis_meta_data meta_dt(meta_prefixes, region);
  meta_dt.meta_analyze();
  auto reader_truth = SvarReader("data/pval0.cis_meta.svar.tsv");
  auto reader_test = SvarReader("data/test_output.pval0.cis_meta.svar.tsv");
  ASSERT_TRUE(reader_truth == reader_test);
}