logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

input <- snakemake@input
output <- snakemake@output

library(dplyr)
library(magrittr)

np_ids <- readRDS(input$np_ids)

cnames <- c("fid", "iid", "father", "mother", "sex", "pheno")
ids <- read.table(input$fam, sep = " ", col.names = cnames) %>%
    mutate(fid_iid = paste0(fid, "_", iid))
saveRDS(ids$iid, output$ids)

# HWE calculations should use a subset of monoethnic, unrelated controls,
# but I'm not sure if we have a way to define controls from the BPS variable.
# For now, use NHW samples without BPS measured.

hwe_ids <- filter(ids, !(fid_iid %in% np_ids))$iid

saveRDS(hwe_ids, output$hwe_ids)
cat("\n\nn samples: ")
cat(nrow(ids))
cat("\n\nn samples without BPS: ")
cat(length(hwe_ids))
cat("\n\nn samples with BPS: ")
cat(sum(np_ids %in% ids$fid_iid))
cat("\n\n")
