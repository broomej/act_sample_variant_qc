logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

library(magrittr)
library(purrr)
library(dplyr)

sampmet <- map(input$sampmet, readRDS) %>%
    bind_rows() %>%
    filter(pass & np_measured)
if (any(is.na(c(sampmet$pass, sampmet$np_measured)))) {
    stop("all samples should have `pass` and `np_measured`")
}

varmet <- map(input$varmet, readRDS) %>% bind_rows()
select(varmet, hwe_rm, ms_rm, pass) %>% table(useNA = "ifany")

saveRDS(varmet, output$varmet)

# For multi-ancestry PCA, use variants that pass across all cohorts
fail_all <- filter(varmet, !pass)$id
pass_all <- filter(varmet, !(id %in% fail_all))$id %>% unique()

# For GWAS, keep variants that pass in both ACT1 and ACT3_NHW
fail_nhw <- filter(varmet, cohort %in% c("ACT1", "ACT3_NHW")) %>%
    filter(!pass) %$%
    id
pass_nhw <- filter(varmet, !(id %in% fail_nhw))$id %>% unique()

cat("\n\nvariants before filtering: ")
unique(varmet$id) %>% length() %>% cat()
cat("\n")

cat("\n\nn variants kept, all cohorts: ")
length(pass_all) %>% cat()
cat("\n")
saveRDS(pass_all, output$keep_all)

cat("\n\nn variants kept, NHW cohorts: ")
length(pass_nhw) %>% cat()
cat("\n")
saveRDS(pass_nhw, output$keep_nhw)
