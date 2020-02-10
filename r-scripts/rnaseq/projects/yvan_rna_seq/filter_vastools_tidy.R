library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(tidyverse)
library(writexl)
library("tidylog", warn.conflicts = FALSE)
args = commandArgs(trailingOnly=TRUE)

average_per_group=FALSE
if (length(args)<2) {
  stop("Input file and number of samples/groups should be specified.\n
       --3rd argument, if specified, refers to the delta PSI value to consider. Default:20
       --4th argument, if specified ,refers to whether PSI means per groups should be calculated. 
       It should be a file with the group name (1st col) and the column names (2nd col) referring to such group.
       One group per line ", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  tidy_file = args[1]
  n_samples = args[2]
  dPSI = 20
  outfile = "TIDY_INCLUSION_dPSI_20.txt"
} else if (length(args)>2) {
  tidy_file = args[1]
  n_samples = args[2]
  dPSI = as.numeric(args[3])
  outfile = paste0("TIDY_INCLUSION_dPSI_", dPSI, ".txt")
  if(file.exists(args[4])){
    average_per_group=TRUE
    groups = read_tsv(args[4], col_names = FALSE)
  }
}

#inclusion_table
#groups
#inclusion_table = read_tsv("../../cajal_bodies/vasttools/INCLUSION_LEVELS_TIDY_minFr0.5_noVLOW-pIR.tsv")
#groups = read_tsv("../../cajal_bodies/vasttools/groups_to_filter_vast.txt", col_names = c("group", "cols_to_avg")) 

inclusion_table = read_tsv(tidy_file)
inclusion_table <- inclusion_table %>%
  separate(col="EVENT", into=c("gene", "event"), sep = "=")

#no_na <- inclusion_table[rowSums(is.na(inclusion_table))==0,]
#with_na <- inclusion_table %>% 
#              filter(is.na(Ctrl) | is.na(ICM) | is.na(DCM))

avg_per_group <- function(group_df, df_inclusion) {
    print(group_df)
    group_name = as.name(group_df$group)
    colnames = unlist(strsplit(group_df$cols_to_avg, ","))
    df_inclusion %>%
      mutate(!! group_name := rowMeans(select(., !!colnames), na.rm = TRUE))
}

if (average_per_group){
  print("Groups were provided!")
  df_inclusion_new <- groups %>% 
    rowwise() %>%
    avg_per_group(., inclusion_table)
 # by(groups, 1:nrow(groups), avg_per_group(row),inclusion_table)
  
  #groups %>% 
   # t() %>% colwise() %>%
    #avg_per_group(inclusion_table)
}


# inclusion_table$t0 <- round(rowMeans(subset(inclusion_table, select = c(3,4,5)), na.rm = FALSE), digits = 2)
#  inclusion_table$t12 <- rowMeans(subset(inclusion_table, select = c(6,7,8)), na.rm = FALSE)
#  inclusion_table$t24 <- rowMeans(subset(inclusion_table, select = c(9,10,11)), na.rm = FALSE)
# inclusion_table$t48 <- rowMeans(subset(inclusion_table, select = c(12,13,14)), na.rm = FALSE)

filtered <- do.call(rbind,
                    lapply(combn(tail(names(inclusion_table), n_samples), 2, simplify = FALSE),
                           function(x) {
                             test <- abs(inclusion_table[, x[1]] - inclusion_table[, x[2]]) >= dPSI
                             inclusion_table[test,] })
)


unique_filtered <- unique(filtered[order(filtered$gene),])
unique_filtered %>% readr::write_delim(outfile, delim="\t", na="", col_names=T)
show(dim(unique_filtered))

#write_xlsx(as.data.frame(format(unique_filtered, digits=2)), path  = outfile, col_names=T, format_headers = T)
#write.table(format(unique_filtered, digits=2), sep="\t", 
#            row.names = F, quote = F, file = outfile)
#sessionInfo()
