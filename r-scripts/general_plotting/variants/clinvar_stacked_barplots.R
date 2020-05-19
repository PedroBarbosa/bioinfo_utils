library(ggplot2)


data <- expand.grid('location' = c("exonic","splicesite","deepintronic","unknown"),
                    'class'= c("Benign","Pathogenic"),
                    'df' = c("1s","2s","3s"))
#sure
data$counts <-  c(18571,5914,1148,3110,
                  13979,2918,79,15267, 
                  6454,1909,869,644,
                  4100,571,26,5014,
                  460,110,803,209,
                  1413,119,3,3489)


#with likely
data$counts <- c(86185,22710,1905,12656,
                 13979,2918,79,15267,
                 16895,3635,911,1031,
                 4100,571,26,5014,
                 1831,239,803,238,
                 1413,119,3,3489)

ggplot() +
  geom_bar(data=data, aes(y = counts, x = class, fill = location), stat="identity",
           position='stack') +
  theme_bw(base_size = 22) + 
  facet_grid( ~ df) +
  scale_fill_manual(values = c("darkgrey","burlywood4","darkslateblue","brown4"))



                
