logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

input <- snakemake@input
output <- snakemake@output

library(dplyr)
library(magrittr)

np_ids <- readRDS(input$np_ids)

cnames <- c("fid", "iid", "father", "mother", "sex", "pheno")
fam <- read.table(input$fam, sep = " ", col.names = cnames) %>%
    mutate(fid_iid = paste0(fid, "_", iid))
saveRDS(fam$iid, output$ids)

# HWE calculations should use a subset of monoethnic, unrelated controls
# We'll use AD status to define cases and controls.

hwe_ids <- filter(fam, pheno == 1)$iid

saveRDS(hwe_ids, output$hwe_ids)
cat("\n\nn samples: ")
cat(nrow(fam))
cat("\n\nn samples with BPS: ")
cat(sum(np_ids %in% fam$fid_iid))
cat("\n\n")
