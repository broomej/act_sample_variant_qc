logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")
input <- snakemake@input
output <- snakemake@output

library(SeqArray)
var_include <- readRDS(input$var_include)
gds <- seqOpen(input$gds)
var_excl <- data.frame(
        id = seqGetData(gds, "variant.id"),
        annotation_id = seqGetData(gds, "annotation/id")
    ) |>
    dplyr::filter(!(id %in% var_include))
seqClose(gds)
writeLines(var_excl$annotation_id, output[[1]])
