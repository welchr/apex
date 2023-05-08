#!/bin/bash
set -euxo pipefail

APEX="../cmake-build-release/apex"

jupyter nbconvert --to notebook --execute --allow-errors cov_window_too_small.ipynb

${APEX} cis --vcf cov_window_too_small_s1.vcf.gz --field DS --bed cov_window_too_small_s1.bed.gz --out cov_window_too_small_s1 --long
${APEX} store --vcf cov_window_too_small_s1.vcf.gz --field DS --bed cov_window_too_small_s1.bed.gz --prefix cov_window_too_small_s1 --window 4
${APEX} meta --stepwise --het --sumstats cov_window_too_small_s1 --pvalue 0.05 --backward 0.05 --use-marginal-pval --prefix cov_window_too_small_s1

${APEX} cis --vcf cov_window_too_small_s2.vcf.gz --field DS --bed cov_window_too_small_s2.bed.gz --out cov_window_too_small_s2 --long
${APEX} store --vcf cov_window_too_small_s2.vcf.gz --field DS --bed cov_window_too_small_s2.bed.gz --prefix cov_window_too_small_s2 --window 4
${APEX} meta --stepwise --het --sumstats cov_window_too_small_s2 --pvalue 0.05 --backward 0.05 --use-marginal-pval --prefix cov_window_too_small_s2

${APEX} meta --sumstats cov_window_too_small_s1,cov_window_too_small_s2 --prefix cov_window_too_small_meta
${APEX} meta --stepwise --het --sumstats cov_window_too_small_s1,cov_window_too_small_s2 --pvalue 0.05 --backward 0.05 --use-marginal-pval --prefix cov_window_too_small_meta
