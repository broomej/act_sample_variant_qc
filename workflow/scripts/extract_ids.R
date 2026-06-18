logcon <- file(snakemake@log[[1]], open = "wt")
sink(logcon)
sink(logcon, type = "message")

output <- snakemake@output
input <- snakemake@input
library(magrittr)
library(dplyr)
library(tidyr)
library(purrr)
bps_ids <- read.csv(input$bps, header = TRUE, stringsAsFactors = FALSE) %>%
    select(subject)
if (any(duplicated(bps_ids$subject))) {
    stop("Duplicate subject IDs found in the BPS file.")
}

# There are 5 duplicate rows in the xwalk file, remove them before proceding.
xwalk <- read.csv(input$xwalk, header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(across(ACT_ID, as.character))
dup_ids <- xwalk$ACT_ID[duplicated(xwalk$ACT_ID)]
dups <- map(xwalk, duplicated) %>% data.frame()
dup_row <- apply(dups, 1, any)
xwalk <- xwalk[!dup_row, ]

dup_check <- map(xwalk, duplicated) %>%
    unlist() %>%
    any()
if (dup_check) stop("There are still duplicates in the xwalk file")
if (!any(dup_ids %in% xwalk$ACT_ID)) {
    stop("Too many rows were removed from the xwalk when attempting to remove duplicates")
}

typed_ids <- readLines(input$typed_ids) %>%
    # For completion, this includes every ID an each VCF file. We'd expect the
    # same IDs in each VCF per cohort, but include all in case any aren't
    # in all files.
    unique() %>%
    data.frame(id = .) %>%
    separate_wider_delim(id, "_",
        too_many = "error", too_few = "error",
        cols_remove = FALSE, names = c("FID", "IID")
    ) %>%
    mutate(ID_type = case_when(
        grepl("^ACT", IID) ~ "ACT_ID",
        grepl("^A0|^B0|^C0|^D0", IID) ~ "IndNo",
        grepl("^\\d{8}$", IID) ~ "NumID",
        TRUE ~ "Unknown"
    ))

if (any(duplicated(typed_ids$id))) {
    stop("Duplicate IDs found in the genotype data.")
}

if (any(typed_ids$ID_type == "Unknown")) {
    print(sum(typed_ids$ID_type == "Unknown"))
    stop("Some IDs in the typed IDs file could not be classified:")
}

indno <- typed_ids %>%
    filter(ID_type == "IndNo") %>%
    left_join(xwalk, by = c("IID" = "IndNo")) %>%
    mutate(IndNo = IID)
if (any(is.na(indno$ACT_ID))) {
    warning("Some IndNo IDs could not be matched to ACT_IDs in the xwalk:")
    print(sum(is.na(indno$ACT_ID)))
}

numid <- typed_ids %>%
    filter(ID_type == "NumID") %>%
    left_join(xwalk, by = c("IID" = "IndNo")) %>%
    mutate(IndNo = IID)
if (any(is.na(numid$ACT_ID))) {
    warning("Some numeric IDs could not be matched to ACT_IDs in the xwalk.")
    print(sum(is.na(numid$ACT_ID)))
}
actid <- filter(typed_ids, ID_type == "ACT_ID") %>%
    mutate(ACT_ID = gsub("ACT", "", IID)) %>%
    left_join(xwalk, by = "ACT_ID")
if (any(is.na(actid$IndNo))) {
    warning("Some ACT_IDs could not be matched to ACT_IDs in the xwalk.")
    print(sum(is.na(actid$IndNo)))
}

ids_matched <- rbind(indno, numid, actid) %>%
    select(id, FID, IID, ID_type, ACT_ID, IndNo, DNANumber, everything())
if (any(duplicated(ids_matched$id))) {
    stop("Duplicate IDs found in the matched IDs data frame.")
}
ids_matched_na_rm <- filter(ids_matched, !is.na(ACT_ID))
ids_not_matched <- filter(ids_matched, is.na(ACT_ID))

np_typed <- filter(ids_matched_na_rm, ACT_ID %in% bps_ids$subject)$id

n_typed <- nrow(ids_matched)
n_np <- nrow(bps_ids)
n_np_typed <- length(np_typed)
cat("\nGenotyped samples:\n")
cat(n_typed)
cat("\nNP samples:\n")
cat(n_np)
cat("\nGenotyped NP samples:\n")
cat(n_np_typed)
cat("\nNP samples not genotyped:\n")
cat(n_np - n_np_typed)
cat("\n\n")

gsub("_", " ", np_typed) %>% writeLines(output$txt)
saveRDS(np_typed, output$rds)
write.table(ids_matched, output$ids_matched, row.names = FALSE, sep = ",", quote = FALSE)
