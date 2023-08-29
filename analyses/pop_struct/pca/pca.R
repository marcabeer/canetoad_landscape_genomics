

library(LEA)
#provides read.matrix() function
library(tseries)
library(adegenet)
library(vcfR)

#provide plotting functions for PCA results
library(factoextra)
library(ggplot2)
library(ggnewscale)
library(cowplot)


##############################
#Load cane toad metadata
##############################
md <- read.table("metadata_individuals", sep="\t", header=TRUE)
latlong <- subset(ct_metadata, select=c(long,lat))

#load VCF
##need to gunzip vcf first
lfmm_genotypes<-vcf2lfmm("ct.vcf")


##############################
#Impute missing genotypes
##############################

#estimate admixture coefficients assuming different values of K supported by TESS3
k_vec<-c(seq(from=1, to=40, by=1))
lfmm_snmf<-snmf(input.file=lfmm_genotypes, K=k_vec, project="new", entropy=TRUE, repetitions=3, CPU=5)

#identify run with lowest cross-entropy (best genotype prediction accuracy)
k_min_cross<-as.data.frame(matrix(NA, nrow=length(k_vec), ncol=4))
colnames(k_min_cross)<-c("k", "min_cross", "med_cross", "max_cross")

for (i in 1:length(k_vec)){
  k_min<-min(cross.entropy(lfmm_snmf, K = k_vec[i]))
  k_med<-median(cross.entropy(lfmm_snmf, K = k_vec[i]))
  k_max<-max(cross.entropy(lfmm_snmf, K = k_vec[i]))
  k_min_cross$k[i]<-k_vec[i]
  k_min_cross$min_cross[i]<-k_min
  k_min_cross$med_cross[i]<-k_med
  k_min_cross$max_cross[i]<-k_max
}

#output crossvalidation results
#write.table(k_min_cross, "k_cross_entropy", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#if not re-running preceding code, read in cross entropy results
k_min_cross <- read.table("k_cross_entropy", header=TRUE, sep="\t")
plot(x=k_min_cross$k, y=k_min_cross$med_cross, col="black")

#check delta cross entropies to determine when change in cross entropy decreases
k_vec_offset <- k_vec[-1]
k_delta_ce <- as.data.frame(matrix(NA, nrow=length(k_vec_offset), ncol=2))
colnames(k_delta_ce) <- c("k", "delta_ce")

for (i in 1:length(k_vec_offset)){
  k_delta_ce$k[i] <- k_vec_offset[i]
  k_delta_ce$delta_ce[i] <- k_min_cross$med_cross[k_vec_offset[i]] - k_min_cross$med_cross[k_vec_offset[i]-1]
}

#K ~ 12 is where cross entropy asymptotes
plot(x=k_delta_ce$k, y=k_delta_ce$delta_ce)

#isolate best replicate of K = 12
best<-which.min(cross.entropy(lfmm_snmf, K=12))

#impute missing genotypes
lfmm_imputed<-impute(lfmm_snmf, input.file=lfmm_genotypes, method="mode", K=12, run=best)

#load in imputed genotype matrix
lfmm_imputed_genotype_matrix<-read.matrix("ct.lfmm_imputed.lfmm", header=FALSE, sep=" ")

#save copy of imputed genotype matrix for manipulation
imputed_geno<-lfmm_imputed_genotype_matrix


##############################
#Run PCA using imputed genotypes
##############################

#imputed_geno[imputed_geno==0]<-20
#imputed_geno[imputed_geno==2]<-60
#imputed_geno[imputed_geno==20]<-2
#imputed_geno[imputed_geno==60]<-0

#create genotype matrix in proper format for genind
geno_matrix <- data.frame(matrix(ncol=dim(imputed_geno)[2]*2, nrow=dim(imputed_geno)[1]))
loci <- seq(1, dim(geno_matrix)[2], by=2)
loci2 <- seq(2, dim(geno_matrix)[2], by=2)
geno_matrix [,loci] <- imputed_geno
geno_matrix[,loci2] <- abs(geno_matrix[,loci]-2)

#convert to integer
geno_matrix_unlist <- as.integer(unlist(geno_matrix))
geno_matrix_int <- data.frame(matrix(ncol=ncol(geno_matrix), nrow=nrow(geno_matrix), data=geno_matrix_unlist))

#create genind object out of original VCF
ct_vcf <- read.vcfR("ct.vcf")
ct_genind <- vcfR2genind(ct_vcf)

#insert imputed genotype matrix into copy of genind
ct_genind_imp <- ct_genind
ct_genind_imp$tab <- as.matrix(geno_matrix_int)

#run PCA
ct_int_pca<-dudi.pca(ct_genind_imp, scannf=FALSE, nf=10)
#saveRDS(ct_int_pca, "ct_pca_Rfile")

#if not re-running preceding code, read in PCA results
ct_int_pca <- readRDS("ct_pca_Rfile")

#make screeplot (visualize % explained variance for each PC)
##PC1 and PC2 explain quite a bit of variance; PC1 most important for explaining pop structure
pvariance <- data.frame(PC=seq(from=1, to=length(ct_int_pca$eig), by=1), pvar=ct_int_pca$eig/sum(ct_int_pca$eig)*100)
scree <- ggplot()+
  geom_line(data=pvariance[1:20,], aes(x=PC, y=pvar), size=1, colour="gray50")+
  geom_point(data=pvariance[1:20,], aes(x=PC, y=pvar), size=2, colour="gray50")+
  scale_x_discrete(limits=c(seq(from=1, to=20, by=1)))+
  xlab("PC")+
  ylab("Percent variance explained")+
  theme_bw()
scree


##############################
#Plot PCA
##############################

#create PCA df
ct_pca_df <- data.frame(samp_id=md$sample, pop=md$pop, region_name=md$region_name, long=md$long, lat=md$lat, pc1=ct_int_pca$li[,1], pc2=ct_int_pca$li[,2], pc3=ct_int_pca$li[,3])
#write.table(pc_df, "ct_pca_df", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#plot parameters
plot_marg <- 15
pal <- c("tan1", "turquoise3", "slateblue3")

p_pca<-ggplot()+
  geom_point(data=ct_pca_df, aes(x=pc1, y=pc2, colour=as.factor(region_name)), size=2.5, alpha=0.5)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  scale_colour_manual(values=pal)+
  scale_x_continuous(limits=c(-50, 50), breaks=seq(from=-50, to=50, by=10), expand=c(0,0))+
  scale_y_continuous(limits=c(-45, 45), breaks=seq(from=-50, to=50, by=10), expand=c(0,0))+
  labs(x="PC1 (7.18%)", y="PC2 (2.53%)", colour="Region")+
  theme_bw()+
  theme(legend.position=c(0.15,0.15))+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
p_pca
ggsave("pca_plot.svg", width=8, height=8)

