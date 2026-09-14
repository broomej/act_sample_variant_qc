logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")
input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

library(dplyr)

var_qc <- readRDS(input$var_qc) %>%
    # exclude variants that are entirely missing from this cohort
    filter(missing_rate < 1L) %>%
    mutate(
        hwe_rm = p < params$hwe_p,
        ms_rm = missing_rate > params$var_ms,
        pass = !(hwe_rm | ms_rm)
    )
n_fail <- sum(!var_qc$pass, na.rm = TRUE)
pass <- filter(var_qc, pass | is.na(pass))
if (nrow(var_qc) - nrow(pass) != n_fail) {
    stop("Inconsistent number of failing variants. Check NA filtering")
}
saveRDS(pass$id, output$keep)
saveRDS(var_qc, output$varmet)
