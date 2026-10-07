logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

library(dplyr)
library(magrittr)

input <- snakemake@input
output <- snakemake@output

cnames <- c("chr", "vid", "cmg", "bps", "a1", "a2")
bim <- read.table(input$bim, header = FALSE, stringsAsFactors = FALSE,
                  sep = "\t", col.names = cnames)

cat("\nn variants in BIM: ")
cat(nrow(bim))
bim$vid %>% writeLines(output[[1]])
