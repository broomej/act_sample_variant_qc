logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

library(magrittr)
library(purrr)
library(dplyr)

sampmet <- readRDS(input$sampmet) %>% filter(pass & np_measured)
if (any(is.na(c(sampmet$pass, sampmet$np_measured)))) {
    stop("all samples should have `pass` and `np_measured`")
}
glimpse(sampmet)
cohorts_keep <- unique(sampmet$cohort)

cohort_labels <- params$cohort
names(cohort_labels) <- cohort_labels
varmet_files <- cohort_labels
for (cohort in cohort_labels) {
    idx <- grep(cohort, input$varmet)
    varmet_files[cohort] <- input$varmet[idx]
}
varmet <- map(varmet_files, readRDS)

if (length(cohorts_keep) == 1L) {
    var_summ <- varmet[[cohorts_keep]]
} else if (length(cohorts_keep) > 1L) {
    var_summ <- varmet[cohorts_keep] %>%
        map(select, id, pass, hwe_rm, ms_rm) %>%
        imap(~ rename_with(.x, .f = function(col) paste0(.y, "_", col), .cols = -id)) %>%
        reduce(inner_join, by = "id")
    pass <- select(var_summ, matches("pass")) %>%
        rowSums(na.rm = TRUE) %>%
        as.logical()
    var_summ %<>% select(-matches("pass"))
    var_summ$pass <- pass
    varmet$summary <- var_summ
} else {
    stop("Expecting one or more cohorts to get variant keep list from")
}

saveRDS(varmet, output$varmet)

cat("\n\nvariants before filtering: ")
nrow(var_summ) %>% cat()
cat("\n")

ids_keep <- filter(var_summ, pass | is.na(pass))$id
cat("\n\nn variants kept for analysis: ")
length(ids_keep) %>% cat()
cat("\n")
saveRDS(ids_keep, output$keep)
