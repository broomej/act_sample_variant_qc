logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

input <- snakemake@input
output <- snakemake@output
library(magrittr)
library(dplyr)
library(SeqArray)
library(stringr)

cnames <- c("chr", "vid", "cmg", "bps", "a1", "a2")
bim <- read.table(input$bim,
    header = FALSE, stringsAsFactors = FALSE,
    sep = "\t", col.names = cnames
) %>%
    mutate(id = paste(chr, bps, sep = ":"))

gds <- seqOpen(input$gds)

variant <- data.frame(
    id = seqGetData(gds, "variant.id"),
    chr = seqGetData(gds, "chromosome"),
    pos = seqGetData(gds, "position"),
    annotation.id = seqGetData(gds, "annotation/id")
) %>%
    mutate(
        chr_pos = paste(chr, pos, sep = ":"),
        id_extracted = str_extract(annotation.id, "^\\d{1,2}:\\d+")
    ) %>%
    filter(!(is.na(id_extracted) & chr_pos == "0:0"))

filter(variant, chr_pos %in% bim$id)$id %>% saveRDS(output$var_ids)
