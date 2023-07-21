#!/usr/bin/env Rscript
source("Apex2R.r")
source("cond_p.R")
library(optparse)
library(data.table)

# Single study results:
ssm <- metaSumStats(
  "/path/to/studies/study1_prefix",
  "/path/to/studies/study2_prefix"
)

# Conditional results
f=paste0("/path/to/results/your_results.cis_meta.stepwise_het.tsv")

res = data.frame(fread(f,stringsAsFactors=F))
colnames(res)[1] = "gene"

res_list=split(res$variant,res$gene)
num_sig=lapply(1:length(res_list),function(x){length(res_list[[x]])})
res_list=res_list[num_sig>1]

# will result in error: x$.self$finalize(): attempt to apply non-function - Li advised to ignore this issue
ptm <- proc.time()
out=lapply(1:length(res_list),loop_gene_tryCatch)
proc.time() - ptm

save(out,file="cond_heter.RData")

qtl=do.call(rbind,out)
qtl$variant = rownames(qtl)
write.table(qtl,file="cond_heter.txt", quote = F, sep = "\t", row.names =F)
