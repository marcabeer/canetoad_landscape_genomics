#read in Bayenv2 GEA results
##five replicate runs
gea_r1<-read.table("gea_r1_cat", sep="\t", header=FALSE)[,-c(1,14)]
gea_r2<-read.table("gea_r2_cat", sep="\t", header=FALSE)[,-c(1,14)]
gea_r3<-read.table("gea_r3_cat", sep="\t", header=FALSE)[,-c(1,14)]
gea_r4<-read.table("gea_r4_cat", sep="\t", header=FALSE)[,-c(1,14)]
gea_r5<-read.table("gea_r5_cat", sep="\t", header=FALSE)[,-c(1,14)]

#calculate median values for each SNP/test statistic combination
gea_medians<-as.data.frame(matrix(data=NA, nrow=nrow(gea_r1), ncol=ncol(gea_r1)))
for (i in 1:nrow(gea_r1)){
  for (j in 1:ncol(gea_r1)){
    gea_medians[i,j]<-median(c(gea_r1[i,j],gea_r2[i,j],gea_r3[i,j],gea_r4[i,j], gea_r5[i,j]))
  }
}
colnames(gea_medians)<-c(paste("Bio2", c("_bf", "_rho", "_pearson"), sep=""), paste("Bio4", c("_bf", "_rho", "_pearson"), sep=""), paste("Bio12", c("_bf", "_rho", "_pearson"), sep=""), paste("Elev", c("_bf", "_rho", "_pearson"), sep="")) 

#remove Pearson's r statistic (more prone to outliers than Spearman's rho)
gea_medians_nopearson<-gea_medians[,-c(3,6,9,12)]

#for each environmental factor
##identify SNPs with >=10 Bayes Factor (bf) and >=95% quantile absolute Spearman's rho

be_gea_sig_kf<-list()

index<-c(1,3,5,7)
for (i in 1:4){
  top_bf<-which(gea_medians_nopearson[,index[i]]>=10)
  
  rho_thresh<-quantile(abs(gea_medians_nopearson[,index[i]+1]), 0.95)
  top_rho<-which(abs(gea_medians_nopearson[,index[i]+1])>=rho_thresh)
  
  be_gea_sig_kf[[i]]<-intersect(top_bf, top_rho)
  
}
names(be_gea_sig_kf)<-c("Bio2", "Bio4", "Bio12", "Elev")

#save output
saveRDS(be_gea_sig_kf, "be_gea_sig_kf_r95")
