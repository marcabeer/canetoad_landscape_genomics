setwd("~/analyses/pop_struct/construct")

library(conStruct)

##############################
#Vignettes to explain plotting
##############################
#https://github.com/gbradburd/conStruct

# formatting data
vignette(topic="format-data",package="conStruct")

# how to run a conStruct analysis
vignette(topic="run-conStruct",package="conStruct")

# how to visualize the output of a conStruct model
vignette(topic="visualize-results",package="conStruct")

# how to compare and select between different conStruct models
vignette(topic="model-comparison",package="conStruct")


##############################
#Read in outputs of running cross-validation
##############################

#path to crossvalidation outputs
cv_subdir <- "construct_crossvalidation_best_runs"

#read in cross-validation results
my.xvals <- readRDS(paste(getwd(), cv_subdir, "my.xvals", sep="/"))

sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)

# first, get the 95% confidence intervals for the spatial and nonspatial
#   models over values of K (mean +/- 1.96 the standard error)

sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:10 with 20 replicates
par(mfrow=c(2,1))
plot(rowMeans(sp.results),
     pch=19,col="blue", cex=1.5,
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.results,nsp.results),
     main="cross-validation results")
points(rowMeans(nsp.results),col="green",pch=19, cex=1.5,)

# finally, visualize results for the spatial model
#   separately with its confidence interval bars
#
# note that you could do the same with the spatial model, 
#   but the confidence intervals don't really show up 
#   because the differences between predictive accuracies
#   across values of K are so large.

plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.CIs),
     main="spatial cross-validation results", cex=1.5,)
segments(x0 = 1:nrow(sp.results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp.results),
         y1 = sp.CIs[2,],
         col = "blue",lwd=2)


#find best spatial model replicate for each value of K to compare layer contributions
best_rep <- data.frame(K=1:10, rep=NA)
for (i in 1:nrow(sp.results)){
  best_rep[i,2] <- which.max(sp.results[i,])
}


##############################
#Compare layer contributions
##############################
layer.contributions <- matrix(NA,nrow=10,ncol=10)

rep_path <- paste(getwd(), cv_subdir, sep="/")

load(paste(getwd(), cv_subdir,"cv_sp_rep1K1_conStruct.results.Robj", sep="/"))
load(paste(getwd(), cv_subdir,"cv_sp_rep1K1_data.block.Robj", sep="/"))

tmp <- conStruct.results[[1]]$MAP$admix.proportions
layer.contributions[,1] <- c(calculate.layer.contribution(conStruct.results[[1]],data.block),rep(0,9))


for(i in 2:10){
  # load the conStruct.results.Robj and data.block.Robj
  #   files saved at the end of a conStruct run
  load(paste(rep_path, "/cv_sp_rep", best_rep[i,2], "K", best_rep[i,1], "_conStruct.results.Robj",  sep=""))
  load(paste(rep_path, "/cv_sp_rep", best_rep[i,2], "K", best_rep[i,1], "_data.block.Robj",  sep=""))
  
  
  # match layers up across runs to keep plotting colors consistent
  #   for the same layers in different runs
  tmp.order <- match.layers.x.runs(tmp,conStruct.results[[1]]$MAP$admix.proportions)  
  
  # calculate layer contributions
  layer.contributions[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                                            data.block=data.block,
                                                            layer.order=tmp.order),
                               rep(0,10-i))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order]
}


par(mfrow=c(1,1))
barplot(layer.contributions,
        xlab="",
        ylab="layer contributions",
        names.arg=paste0("K=",1:10))



##############################
#Make STRUCTURE-esque ancestry plots
##############################

#read in locality metadata
md <- read.table(paste(getwd(),"env_filtered", sep="/"), header=TRUE, sep="\t")


###
#load in best K=3 crossvalidation run
load(paste(rep_path, "/cv_sp_rep15K3_conStruct.results.Robj", sep="/"))
load(paste(rep_path, "cv_sp_rep15K3_data.block.Robj", sep="/"))

# assign model results to new variable names
k3sp_cr <- conStruct.results
k3sp_db <- data.block

#get admixture/ancestry proportions
admix.props_k3 <- k3sp_cr$chain_1$MAP$admix.proportions


###
#load in best K=4 crossvalidation run
load(paste(rep_path, "cv_sp_rep7K4_conStruct.results.Robj", sep="/"))
load(paste(rep_path, "cv_sp_rep7K4_data.block.Robj", sep="/"))

# assign model results to new variable names
k4sp_cr <- conStruct.results
k4sp_db <- data.block

#get admixture/ancestry proportions
admix.props_k4 <- k4sp_cr$chain_1$MAP$admix.proportions


###
#load in best K=5 crossvalidation run
load(paste(rep_path, "cv_sp_rep19K5_conStruct.results.Robj", sep="/"))
load(paste(rep_path, "cv_sp_rep19K5_data.block.Robj", sep="/"))

# assign model results to new variable names
k5sp_cr <- conStruct.results
k5sp_db <- data.block

admix.props_k5 <- k5sp_cr$chain_1$MAP$admix.proportions


###
#visual inspection
par(mfrow=c(3,2))
make.structure.plot(admix.proportions = admix.props_k3, sample.order=order(md$long))
make.admix.pie.plot(admix.proportions = admix.props_k3, coords = data.block$coords, radii=2)

make.structure.plot(admix.proportions = admix.props_k4, sample.order=order(md$long))
make.admix.pie.plot(admix.proportions = admix.props_k4, coords = data.block$coords, radii=2)

make.structure.plot(admix.proportions = admix.props_k5, sample.order=order(md$long))
make.admix.pie.plot(admix.proportions = admix.props_k5, coords = data.block$coords, radii=2)

