library(dplyr)
library(tidyr)

inclusion_table = read.table("~/Desktop/INCLUSION_LEVELS_TIDY_minFr0.5_noVLOW-pIR.tsv", header=T)
inclusion_table <- inclusion_table %>%
  separate(col="EVENT", into=c("gene", "event"), sep = "=")


inclusion_table$t0 <- round(rowMeans(subset(inclusion_table, select = c(3,4,5)), na.rm = FALSE), digits = 2)
inclusion_table$t12 <- rowMeans(subset(inclusion_table, select = c(6,7,8)), na.rm = FALSE)
inclusion_table$t24 <- rowMeans(subset(inclusion_table, select = c(9,10,11)), na.rm = FALSE)
inclusion_table$t48 <- rowMeans(subset(inclusion_table, select = c(12,13,14)), na.rm = FALSE)

no_na <- inclusion_table[rowSums(is.na(inclusion_table))==0,]

filtered <- do.call(rbind,
        lapply(combn(tail(names(no_na), 4), 2, simplify = FALSE),
               function(x) {
                 test <- abs(no_na[, x[1]] - no_na[, x[2]]) >= 20
                 no_na[test,]
               })
)

unique_filtered <- unique(filtered[order(filtered$gene),])
dim(unique_filtered)
write.table(format(unique_filtered[,c("gene", "event", "t0", "t12", "t24", "t48")], digits=2), sep="\t", 
            row.names = F, quote = F, file = "~/Desktop/vasttools_deltaPSI0.2.tsv")
