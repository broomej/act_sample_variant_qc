logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

library(magrittr)
library(dplyr)
library(purrr)
library(ggplot2)

np_ids <- readRDS(input$np_ids) %>% gsub("^\\d+_", "", .)
samp_ms <- readRDS(input$samp_ms) %>% mutate(ms_rm = missing_rate > params$samp_ms)
samp_htz <- readRDS(input$samp_htz)

 htz_mean <- mean(samp_htz$heterozygosity, na.rm = TRUE)
 htz_sd <- sd(samp_htz$heterozygosity, na.rm = TRUE)
 samp_htz %<>%
     mutate(z_htz = (heterozygosity - htz_mean) / htz_sd) %>%
     mutate(htz_rm = abs(z_htz) > params$htz_sd) %>%
     select(!z_htz)


sexcheck <- readRDS(input$sexcheck) %>% mutate(sexcheck_rm = STATUS == "PROBLEM")
if (!all(sexcheck$STATUS %in% c("PROBLEM", "OK"))) {
    warning("Unexpected values in sexcheck STATUS column")
}
sexcheck %<>% select(id = IID, sexcheck_rm)

sampmet <- full_join(samp_htz, samp_ms, by = "id") %>%
    full_join(sexcheck, by = "id") %>%
    mutate(
        pass = !(ms_rm | htz_rm | sexcheck_rm),
        np_measured = id %in% np_ids,
        cohort = snakemake@wildcards$cohort
    ) %>%
    select(id, cohort, pass, np_measured, missing_rate, heterozygosity, ms_rm,
           htz_rm, sexcheck_rm, everything())
saveRDS(sampmet, output$sampmet)
select(sampmet, pass, np_measured) %>%
    ftable() %>%
    print()
samp_pass <- filter(sampmet, pass, np_measured)$id
saveRDS(samp_pass, output$keep)
