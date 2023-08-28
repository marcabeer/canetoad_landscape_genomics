setwd("~/analyses/pop_struct/construct")


library(conStruct)
library(geodist)
library(parallel)
library(foreach)
library(doParallel)

###
#formatting data
##do it by sampling locality, as this reduces the matrix to 59 rows (localities) instead of 932 individuals

#read in locality allele frequencies
##summarize individual genetic data into their corresponding sampling localities
##this reduces the genetic data matrix to 59 rows (localities) instead of 932 individuals
allele_freq <- as.matrix(read.table("be_freq", sep="\t", header=TRUE))

#read in locality data
locality_data <- as.matrix(read.table("env_filtered", header=TRUE, sep="\t"))
latlong <- locality_data[,c(2,3)]

#calculate geodesic distances between localities
##divide by 1000 to get values in km
geodist <- geodist::geodist(x=latlong, measure="geodesic")/1000

###
#run conStruct with cross-validation
##want 2000 iters, 10 replicates, for K=1:20, with training.prop=0.75
##kamiak storfer partition has 28 CPUs per node
##parallelization happens across replicates, so it is not clear whether having n.nodes > n.reps is beneficial

cl <- makeCluster(20, type="FORK")
registerDoParallel(cl)

my.xvals <- x.validation(train.prop=0.90,
                         n.reps=20,
                         n.iter=10000,
                         K=1:10,
                         freqs=allele_freq[,-1],
                         geoDist=geodist,
                         coords=latlong,
                         save.files=TRUE,
                         make.figs=FALSE,
                         prefix="cv",
                         parallel=TRUE,
                         n.nodes=20)

stopCluster(cl)

saveRDS(my.xvals, "my.xvals")
