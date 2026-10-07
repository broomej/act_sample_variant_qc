logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

input <- snakemake@input
output <- snakemake@output
library(magrittr)
library(dplyr)
library(SeqArray)
library(stringr)

gds <- seqOpen(input$gds)
annot_ids <- readLines(input$annot_ids)

variant <- data.frame(
    id = seqGetData(gds, "variant.id"),
    annotation.id = seqGetData(gds, "annotation/id")
)

cat("\n")
cat("n variants in GDS: ")
cat(nrow(variant))
cat("\n\n")

variants_filtered <- filter(variant, annotation.id %in% annot_ids)
cat("n variants present in cohort: ")
cat(nrow(variants_filtered))
cat("\n\n")
variants_filtered$id %>% saveRDS(output[[1]])
