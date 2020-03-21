setwd("/Users/pbarbosa/analysis/sandra_hcm_rnaseq")
library(tidyverse)

######TTN###########
ttn_data <- read_tsv(file = "TTN_transcript_analysis.tsv", locale = locale(decimal_mark = ",")) %>%
  group_by(Cell_line) %>%
  filter(row_number() <= 3L | stringr::str_detect(Info, 'N2BA_N1a|N2BA_N1b'))


ttn_data <- ttn_data %>% 
  mutate(Isoform_name_concat = paste0(Isoform_name, " (", Info, ")"))

ttn_data$Isoform_name_concat <- factor(ttn_data$Isoform_name_concat, 
                                       levels = c('TTN 212 (N2B; adult heart)','TTN 204 (Novex3; heart and skeletal muscle)', 'TTN 213 (Intron retained; no protein)',
                                                  'TTN 215 (N2BA_N1a; fetal heart)', 'TTN 216 (N2BA_N1b; fetal heart)'))
p<-ggplot(data=ttn_data, aes(x=Isoform_name_concat, y=TPM, fill=Cell_line)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Accent") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  xlab("Isoform name") + ylab("Transcripts per million (TPM)") 
p

#######TNNT2#######
tnnt2_data <- read_tsv(file = "TNNT2_transcript_analysis.tsv", locale = locale(decimal_mark = ",")) 
names(tnnt2_data) <- sub(" ", "", names(tnnt2_data))
tnnt2_data <- tnnt2_data %>%
  mutate(Isoform_name_concat = ifelse(is.na(Info), tnnt2_data$Isoform_name, paste0(Isoform_name, " (", Info, ")"))) %>%
  mutate(Isoform_name_concat = str_replace(Isoform_name_concat," (NA)", ""))

tnnt2_data$Isoform_name_concat <- factor(tnnt2_data$Isoform_name_concat, 
                                         levels = c('TNNT2-205 (hTNNT2 adult)','TNNT2-207 (hTNNT2 adult)', 'TNNT2-208 (hTNNT2 adult)',
                                                    'TNNT2-227 (hTNNT2 fetal)', 'TNNT2-204', 'TNNT2-221','TNNT2-226'))
p<-ggplot(data=tnnt2_data, aes(x=Isoform_name_concat, y=TPM, fill=Cell_line)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Accent") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  xlab("Isoform name") + ylab("Transcripts per million (TPM)") 
p

view(test)

tnnt2_data$Isoform_name_concat

