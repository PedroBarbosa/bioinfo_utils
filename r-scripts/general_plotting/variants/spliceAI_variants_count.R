library(readxl)
library(dplyr)
spliceAI_all <- read_tsv("~/Desktop/2_all_spliceai_table.tsv")
spliceAI_all_maf <- read_tsv("~/Desktop/3_all_maf_spliceai_table.tsv")               

###FREQ HIST###
freq <- as.numeric(spliceAI_all$gnomADg_AF) 
freq[is.na(freq)]<- 0
ggplot() + aes(freq)+ geom_histogram()


###TYPE OF VARIANT####
counts <- spliceAI_all %>% summarize(donor_gain = sum(DS_DG > 0.2),
                           donor_lost = sum(DS_DL > 0.2),
                           acceptor_gain = sum(DS_AG > 0.2),
                           acceptor_lost = sum(DS_AL > 0.2))


spliceAI_all_maf %>%
  count(Gene_SpliceAI) %>%
  ggplot(aes(x=reorder(Gene_SpliceAI,n), y=n)) +
  geom_bar(stat="identity", width = 0.4, color = "black") +
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_light() +
  ylab("SpliceAI variants") +
  xlab("") +
  #theme(aspect.ratio = 1/1.5) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 22))
   

#%>% summarize_all(n_distinct)
barplot(spliceAI_all_maf$SYMBOL)
barplot(as.matrix(counts), cex.names = 2, cex.axis = 2 )
