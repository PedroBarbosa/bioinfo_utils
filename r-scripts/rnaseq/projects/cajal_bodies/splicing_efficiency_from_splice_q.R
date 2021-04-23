library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rstatix)
library(matrixStats)
library(tibble)
library(stringr)

setwd("~/Desktop/NFS_Carmo/mirek_CajalBodies/new_data_spike-ins/splicing/splice_q/")
IER_files <- list.files(".", pattern = "*IER*")
SPI_files <- list.files(".", pattern = "*SE_splice*")
SPI_level2_files <- list.files(".", pattern = "level2")

samples <- c("t0_rep1", "t0_rep2", "t0_rep3", "t12_rep1", "t12_rep2", 
             "t12_rep3", "t24_rep1", "t24_rep2", "t24_rep3", "t48_rep1", "t48_rep2", "t48_rep3")

do_pca <- function(df){
  rv <- rowVars(as.matrix(df))
  select <- order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
  pc <- prcomp(t(df[select,]))
  percentVar <- pc$sdev^2 / sum( pc$sdev^2 )
  
  pca_matrix <- as.data.frame(pc$x)
  pca_matrix <- pca_matrix %>% rownames_to_column(var="sample") %>% separate(sample, c("timepoints", "replicate"), sep = "_", remove = F) 
  
  ggplot(pca_matrix,aes(x=PC1, y=PC2, color=timepoints, shape=replicate)) +
    geom_point(size=5) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() +
    theme(legend.title=element_blank(), text = element_text(size = 20)) 
}

violin <- function(df){
  colors <- c("#999999", "#E69F00", "#56B4E9")
  ggplot(df, aes(x=timepoints, y=efficiency,fill=replicate)) + 
    geom_violin(width=1) + 
    geom_boxplot(width = 0.15, position = position_dodge(0.99), outlier.shape = NA) +
    facet_wrap(~timepoints, ncol = 4, scales = "free_x") + 
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0.8, 1)) + 
    theme(legend.title=element_blank(), text = element_text(size = 20))
}

###########################
###########################
##### IER processing ######
###########################
###########################
all_dfs <- lapply(IER_files, read_tsv)
for (i in seq_along(all_dfs)){
  colnames(all_dfs[[i]])[15]=samples[[i]]
  colnames_to_keep=c("#chr", "IStart", "IEnd", "strand", "intron_ID", "gene_ID", "transcript_ID", samples[[i]])
  all_dfs[[i]] = all_dfs[[i]] %>% dplyr::select(all_of(colnames_to_keep))
}

IER_df <- Reduce(function(x, y) full_join(x, y, by = c("#chr", "IStart", "IEnd", "strand", 
                                          "intron_ID", "gene_ID", "transcript_ID")), all_dfs)
IER_df <- drop_na(IER_df)

#################################
##### Intron-based analysis #####
#################################
# Many introns are redundant. To get an overall overview of splicing efficiencies, we will count the ocurrences of unique
# introns to see if their efficiencies values are altered.

# N = 317533
IER_df_subset <- dplyr::select(IER_df, -c(intron_ID, strand, transcript_ID))

# N = 119337
IER_distinct_efficiencies <- IER_df_subset %>% group_by(`#chr`, IStart, IEnd, gene_ID) %>% 
  distinct() %>% ungroup()

### PCA
to_pca <- IER_distinct_efficiencies %>% mutate_at("IStart", as.character) %>% mutate_at("IEnd", as.character) %>%
  unite("id", c("#chr", "IStart", "IEnd", "gene_ID"), sep="_")
to_pca$id <- seq(1, length(to_pca$id))
do_pca(to_pca %>% column_to_rownames("id"))


### Violin
long_IER <- IER_distinct_efficiencies %>% pivot_longer(!c(`#chr`, IStart, IEnd, gene_ID),
                                                      names_to = "samples", 
                                                      values_to = "efficiency") %>%
  separate(samples, c("timepoints", "replicate"), sep = "_", remove = F) %>%
  unite("intron", c("#chr", "IStart", "IEnd", "gene_ID"), sep="_") %>% 
  mutate(timepoints = str_replace(timepoints, "t","")) 

long_IER$timepoints <- as.numeric(long_IER$timepoints)

violin(long_IER)


### Wilcox
wilcox <- long_IER %>% group_by(replicate) %>%
  wilcox_test(efficiency ~ timepoints, 
              ref.group = "0", 
              exact=T, 
              p.adjust.method = "bonferroni", 
              paired=T,
              detailed = T,
              alternative = "two.sided")

## Linear model with mixed effects
long_IER <- mutate_at(long_IER, 'timepoints', as.integer)
lm_IER <- lmer(efficiency ~ timepoints + (1 | replicate) + (1 | intron),
             data = long_IER)
summary(lm_IER)
qqnorm(resid(lm_IER))
qqline(resid(lm_IER))

## GLM with gamma distribution (1-IER) - inneficiency
long_IER_inneficiency <- long_IER
long_IER_inneficiency$efficiency <- 1 - long_IER_inneficiency$efficiency
glm_IER <- glmer(efficiency ~ timepoints + (1 | replicate) + (1 | intron),
                 data = long_IER_inneficiency, family = Gamma)
summary(glm_IER)

#####################################
##### Transcript-based analysis #####
#####################################
IER_tx_df <- IER_df %>% group_by(transcript_ID) %>% summarise_at(vars(-`#chr`, -IStart, -IEnd, -strand, -intron_ID, -gene_ID), funs(mean(., na.rm=TRUE))) %>% 
  column_to_rownames("transcript_ID")

### PCA
do_pca(IER_tx_df)

### Violin
long_IER_tx <- IER_tx_df %>% rownames_to_column(var="transcript_ID") %>% 
  pivot_longer(!c(transcript_ID), names_to = "samples", values_to = "efficiency") %>%
  separate(samples, c("timepoints", "replicate"), sep = "_", remove = F) 
violin(long_IER_tx)

### Wilcox
wilcox <- long_IER_tx %>% group_by(replicate) %>%
  wilcox_test(efficiency ~ timepoints, 
              ref.group = "t0", 
              exact=T, 
              p.adjust.method = "bonferroni", 
              paired=T,
              detailed = T,
              alternative = "greater")

#####################################
######## Gene-based analysis ########
#####################################
tx_gene_map <- distinct(dplyr::select(IER_df, c(transcript_ID, gene_ID)))
IER_gene_df <- IER_tx_df %>% rownames_to_column("transcript_ID") %>% 
  left_join(tx_gene_map) %>% 
  group_by(gene_ID)%>% 
  summarise_at(vars(-transcript_ID), funs(mean(., na.rm=TRUE))) %>% 
  column_to_rownames("gene_ID")

## PCA
do_pca(IER_gene_df)

## Violin
long_IER_gene <- IER_gene_df %>% rownames_to_column(var="gene_ID") %>% 
  pivot_longer(!c(gene_ID), names_to = "samples", values_to = "efficiency") %>%
  separate(samples, c("timepoints", "replicate"), sep = "_", remove = F) 
violin(long_IER_gene)

## Wilcox
### Wilcox
wilcox <- long_IER_gene %>% group_by(replicate) %>%
  wilcox_test(efficiency ~ timepoints, 
              ref.group = "t0", 
              exact=T, 
              p.adjust.method = "bonferroni", 
              paired=T,
              detailed = T,
              alternative = "greater")

# Generalized linear mixed model
library(lmer4)
gm1 <- glmer(efficiency ~ timepoints + (1 | replicate) + (1 | gene_ID),
             data = long_IER_gene, family = binomial, nAGQ = 2)
summary(gm1)
library(car)
Anova(gm1, type=3)


##########################
##########################
#### SPI processing ######
##########################
##########################
all_dfs <- lapply(SPI_files, read_tsv)
for (i in seq_along(all_dfs)){
  colnames(all_dfs[[i]])[14]=samples[[i]]
  all_dfs[[i]] = all_dfs[[i]] %>% mutate("{samples[[i]]}_spliced" := sj5_cov_split+sj5_cov_split, 
                                         "{samples[[i]]}_unspliced" := sj5_cov_nonsplit+sj3_cov_nonsplit)
  
  colnames_to_remove=c("sj5_cov_split", "sj3_cov_split", "sj5_cov_nonsplit", "sj3_cov_nonsplit")
  all_dfs[[i]] = all_dfs[[i]] %>% dplyr::select(-all_of(colnames_to_remove))
}

SPI_df <- Reduce(function(x, y) full_join(x, y, by = c("#chr", "sj5start", "sj5end", "sj3start", "sj3end", "strand", 
                                                      "intron_ID", "gene_ID", "transcript_ID")), all_dfs) %>% drop_na()


#################################
##### Intron-based analysis #####
#################################
# N = 317685
SPI_df_subset <- dplyr::select(SPI_df, -c(intron_ID, strand, transcript_ID, gene_ID))

# N = 67410
SPI_distinct_efficiencies <- SPI_df_subset %>% group_by(`#chr`, sj5start, sj5end, sj3start, sj3end) %>% 
  distinct() %>% ungroup()

### PCA
to_pca <- SPI_distinct_efficiencies %>% mutate_at("sj5start", as.character) %>% mutate_at("sj5end", as.character) %>%
  mutate_at("sj3start", as.character) %>%  mutate_at("sj3start", as.character) %>%
  unite("id", c("#chr", "sj5start", "sj5end", "sj3start", "sj3end"), sep="_") %>% column_to_rownames("id")

do_pca(to_pca %>% dplyr::select(-contains("spliced"))) 

### Violin
long_SPI_efficiency <- SPI_distinct_efficiencies %>% 
  dplyr::select(-contains("spliced")) %>%
  pivot_longer(!c(`#chr`, sj5start, sj5end, sj3start, sj3end), names_to = "samples", values_to = "efficiency") %>%
  separate(samples, c("timepoints", "replicate"), sep = "_", remove = F) 

violin(long_SPI_efficiency)

# Generalized linear mixed model
to_remove <- c(as.character(c("#chr", "sj5start", "sj5end", "sj3start", "sj3end")), as.character(samples))
long_SPI_from_counts <- SPI_distinct_efficiencies %>%
  dplyr::select(!samples) %>%
  pivot_longer(!c(`#chr`, sj5start, sj5end, sj3start, sj3end), 
               names_to = c("samples", ".value"), 
               names_pattern = "(.+)_(.+$)") %>%
  separate(samples, c("timepoints", "replicate"), sep = "_", remove = F) %>%
  mutate(timepoints = str_replace(timepoints, "t","")) %>%
  mutate_at('timepoints', as.integer) %>% 
  unite("intron", c("#chr", "sj5start", "sj5end", "sj3start", "sj3end"), sep="_")


library(lme4)
glm_SPI <- glmer(cbind(spliced,unspliced) ~ timepoints + (1 | replicate) + (1 | intron),
             data = long_SPI_from_counts, family = binomial)

summary(glm_SPI)
# The obtained coefficient for time (-6.488e-03) is saying 
# that for every additional hour, the logit of p (the probability of getting a spliced read)
# decreases by 6.488e-03. 

plot(glm_SPI)
predict(glm_SPI, type="links", se.fit = T)

# Predict splicing efficiency with new timepoints
new_timepoints <- c(0:60)

# Plug each value into the regression formula using the intercept and coefficients
predict_logit <- fixef(glm_SPI)[1] + new_timepoints * fixef(glm_SPI)[2]

# Convert logits to p
inv_logit = function(input_logit) {1/(1 + exp(-input_logit))}
pred_p <- inv_logit(predict_logit)

plot(pred_p, ylim = c(0.90,1), cex.lab=1.5, cex.axis=1.5,xlab = "Time (hours)", ylab = "Predicted prob of splicing to have occurred")

#####################################
##### Transcript-based analysis #####
#####################################
SE_tx_df <- SE_df %>% group_by(transcript_ID) %>% summarise_at(vars(-`#chr`,  -sj5start, -sj5end, -sj3start, -sj3end, -strand, -intron_ID, -gene_ID), 
                                                                list(mean)) %>% 
  column_to_rownames("transcript_ID")

### PCA
do_pca(SE_tx_df)

### Violin
long_SE_tx <- SE_tx_df %>% rownames_to_column(var="transcript_ID") %>% 
  pivot_longer(!c(transcript_ID), names_to = "samples", values_to = "efficiency") %>%
  separate(samples, c("timepoints", "replicate"), sep = "_", remove = F) 
violin(long_SE_tx)



#####################################
######## Gene-based analysis ########
#####################################
tx_gene_map <- distinct(dplyr::select(SE_df, c(transcript_ID, gene_ID)))
SE_gene_df <- SE_tx_df %>% rownames_to_column("transcript_ID") %>% 
  left_join(tx_gene_map) %>% 
  group_by(gene_ID)%>% 
  summarise_at(vars(-transcript_ID), funs(mean(., na.rm=TRUE))) %>% 
  column_to_rownames("gene_ID")

## PCA
do_pca(SE_gene_df)

## Violin
long_SE_gene <- SE_gene_df %>% rownames_to_column(var="gene_ID") %>% 
  pivot_longer(!c(gene_ID), names_to = "samples", values_to = "efficiency") %>%
  separate(samples, c("timepoints", "replicate"), sep = "_", remove = F) 
violin(long_SE_gene)

## Wilcox
### Wilcox
wilcox <- long_IER_gene %>% group_by(replicate) %>%
  wilcox_test(efficiency ~ timepoints, 
              ref.group = "t0", 
              exact=T, 
              p.adjust.method = "bonferroni", 
              paired=T,
              detailed = T,
              alternative = "greater")

##################################
##################################
##### SE Level 2 processing ######
##################################
##################################
all_dfs <- lapply(SE_level2_files, read_tsv)
for (i in seq_along(all_dfs)){
  colnames(all_dfs[[i]])[14]=samples[[i]]
  colnames_to_keep=c("#chr", "sj5start", "sj5end", "sj3start", "sj3end", "strand", "intron_ID", "gene_ID", "transcript_ID", samples[[i]])
  all_dfs[[i]] = all_dfs[[i]] %>% dplyr::select(all_of(colnames_to_keep))
}

SE_level2_df <- Reduce(function(x, y) full_join(x, y, by = c("#chr", "sj5start", "sj5end", "sj3start", "sj3end", "strand", 
                                                             "intron_ID", "gene_ID", "transcript_ID")), all_dfs) %>% drop_na()



############################
############################
##### coSI processing ######
############################
############################
setwd("~/Desktop/NFS_Carmo/mirek_CajalBodies/new_data_spike-ins/splicing/ipsa/")
coSI <- read_tsv(file = "all.cosi.tsv") %>% column_to_rownames("id") %>% drop_na()

do_pca(coSI)

long_coSI <- coSI %>% rownames_to_column(var="id") %>% 
  pivot_longer(!c(id), names_to = "samples", values_to = "efficiency") %>%
  separate(samples, c("timepoints", "replicate"), sep = "_", remove = F) 
violin(long_coSI)
summary(coSI)

hist(coSI$t48_rep1)
wilcox <- long_coSI %>% group_by(replicate) %>%
  wilcox_test(efficiency ~ timepoints, 
              ref.group = "t0", 
              exact=T, 
              p.adjust.method = "bonferroni", 
              paired=T,
              detailed = T,
              alternative = "greater")
