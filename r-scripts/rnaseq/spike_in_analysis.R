spike_in_table=read.table("spike_in_concentrations.txt", sep="\t",header = TRUE)
mix1=subset(spike_in_table,select=c("ERCC.ID","concentration.in.Mix.1..attomoles.ul."))
mix1[3] = log(mix1[2],2)
colnames(mix1)=c("ERCC","Mix1_concentration","Log2_concentration")



raw_all_ERCC_counts
raw_all_counts
fraction_ERCC=raw_all_ERCC_counts/raw_all_counts
