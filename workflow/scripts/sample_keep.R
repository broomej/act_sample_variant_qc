logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

library(dplyr)

gt <- lapply(input$gt, readRDS) %>% bind_rows() %>% filter(pass)
imputed <- lapply(input$imputed, readRDS) %>% bind_rows() %>% filter(pass)

keep_all <- intersect(gt$id, imputed$id)
keep_gwas <- filter(gt, np_measured, cohort %in% c("ACT1", "ACT3_NHW"),
                    id %in% keep_all)$id

cat("\nn samples: ")
unique(gt$id) %>% length()
cat("\n\nn NHW samples: ")
filter(gt, cohort %in% c("ACT1", "ACT3_NHW"))$id %>% unique() %>% length()
cat("\n\nn samples passing QC: ")
length(keep_all)
cat("\n\nn samples for GWAS (NP measured, pass QC, NHW): ")
length(keep_gwas)
cat("\n")

saveRDS(keep_all, output$keep_all)
saveRDS(keep_gwas, output$keep_gwas)
saveRDS(gt, output$sampmet_typed)
saveRDS(imputed, output$sampmet_imputed)
