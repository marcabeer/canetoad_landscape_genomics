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
#run conStruct final models for K=3:5

final_k3 <- conStruct::conStruct(spatial=TRUE,
                     K=3,
                     freqs=allele_freq[,-1],
                     geoDist=geodist,
                     coords=latlong,
                     prefix="final_k3",
                     n.chains=3,
                     n.iter = 10000,
                     save.files=TRUE,
                     make.figs=TRUE)

final_k4 <- conStruct::conStruct(spatial=TRUE,
                                 K=4,
                                 freqs=allele_freq[,-1],
                                 geoDist=geodist,
                                 coords=latlong,
                                 prefix="final_k4",
                                 n.chains=3,
                                 n.iter = 10000,
                                 save.files=TRUE,
                                 make.figs=TRUE)

final_k5 <- conStruct::conStruct(spatial=TRUE,
                                 K=5,
                                 freqs=allele_freq[,-1],
                                 geoDist=geodist,
                                 coords=latlong,
                                 prefix="final_k5",
                                 n.chains=3,
                                 n.iter = 10000,
                                 save.files=TRUE,
                                 make.figs=TRUE)


saveRDS(final_k3, "final_k3")
saveRDS(final_k4, "final_k4")
saveRDS(final_k5, "final_k5")
