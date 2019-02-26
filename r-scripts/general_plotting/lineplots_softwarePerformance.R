library(ggplot2)
library(reshape2)
library(RColorBrewer)
timedata=read.table(file = "Desktop/time.txt",
              header=TRUE,check.names = FALSE, dec = ",")
df=melt(timedata,id.vars = "Samples")
p<-ggplot(df, aes(x=variable,y=value, group=Samples, color=factor(Samples))) +
      geom_line()+
      geom_point()+
      ylab(label="Running time (h)")+
      xlab(label="CPUs")+    
      ylim(0,80)+
      labs(color= "Number of samples")+
      theme(legend.justification = c(1, 1), legend.position = c(1, 1))
p


##################
###Speedup plot###
##################
speedupdata=read.table(file = "Desktop/speedup.txt",
                       header=TRUE,check.names = FALSE, dec = ",")
df=melt(speedupdata,id.vars = "Samples")
cols  <- brewer.pal(n=4, "Set1")
p<-ggplot(df, aes(x=variable,y=value, group=Samples, color=factor(Samples))) +
  geom_line()+
  geom_point()+
  scale_color_manual(values = cols) +
  ylab(label="Speedup (times)")+
  xlab(label="CPUs")+    
  ylim(0,7)+
  
  labs(color= "Number of samples")+
  theme(legend.justification = c(1, 0), legend.position = c(1, 0))
p

