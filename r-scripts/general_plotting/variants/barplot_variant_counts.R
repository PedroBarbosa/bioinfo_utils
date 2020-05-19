library(ggplot2)
library(RColorBrewer)

deep_intronic=132
in_gnmoad=20
in_clinvar=59
in_clinvar_and_patho_or_likely=50


v<-c(deep_intronic,in_gnmoad,in_clinvar,in_clinvar_and_patho_or_likely)

cols <- brewer.pal(n = 4, name = "PuBuGn")
names(v)<-c("deep intronic","in gnomAD",
            "in Clinvar","in Clinvar & Pathogenic/Likely_pathogenic")
          
barplot(v, ylim = c(0, 140),col = cols, xaxt='n')
legend(locator(1), legend = names(v),fill = cols)

