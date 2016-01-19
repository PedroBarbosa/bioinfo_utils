library(sleuth)
help(package = 'sleuth')
base_dir <- "C:/Users/User/Desktop/sleuth/stem"
sample_id <- dir(file.path(base_dir),"C")
sample_id
kal_dirs <- sapply(sample_id, function(sample_id) file.path(base_dir,sample_id))
kal_dirs                  
s2c <- read.table(file.path(base_dir,"info.txt"),header=TRUE,stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path=kal_dirs)
s2c


#Now, we can create the sleuth object


so <- sleuth_prep(s2c, ~ condition)

#Input our own design matrix
experiment_design <- design_matrix(so)
experiment_design[, 2] <- as.numeric(!experiment_design[, 2])
so <- sleuth_prep(s2c, experiment_design)
so <- sleuth_fit(so)


models(so)
so <- sleuth_wt(so, 'conditionc2')
sleuth_live(so)

results_table <- sleuth_results(so, 'conditionc2')
results_table

                  