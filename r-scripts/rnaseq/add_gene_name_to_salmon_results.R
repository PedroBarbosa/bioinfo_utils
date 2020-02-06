library(tximport)
library(readr)
library(tidyverse)
mart_table <- read_tsv("/Users/pbarbosa/analysis/genome_utilities/hg38/mart_hg38.txt")
setwd("/Users/pbarbosa/analysis/marta_hcm/sandra_3_lines_rna_seq")

files <- list.files(path=".", pattern = "*sf")
names(files) <- c("DF6", "GIBCO", "TCLAB")

add_gene_names <- lapply(files, 
                         function(df) {
                           sample_name <- strsplit(df, "_")[[1]][1]
                           read_tsv(df) %>%
                             inner_join(mart_table, by = c("Name" = "Transcript stable ID version")) %>%
                             select(`Gene stable ID`, Name, `Gene name`, `Gene description`, Length, EffectiveLength, TPM, NumReads, ) %>%
                             arrange(`Gene stable ID`) %>%
                             rename( !!paste0(sample_name, "_TPM") := TPM, !!paste0(sample_name, "_NumReads") := NumReads, !!paste0(sample_name, "_EffectiveLength") := EffectiveLength) 
                         }
)

merged <- add_gene_names %>%
  reduce(left_join, by=c("Gene stable ID", "Name", "Gene name", "Gene description", "Length")) %>%
  write_excel_csv("transcript_expression.csv")
