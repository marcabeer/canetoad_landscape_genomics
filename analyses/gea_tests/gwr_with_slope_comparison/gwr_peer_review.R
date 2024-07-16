
#PCA packages
library(vegan)

#GWR-related packages
library(GWmodel)
library(sp)

#parallel computing
library(parallel)

#data manipulation
library(dplyr)
library(plotrix)

#plotting packages
library(ggplot2)
library(scales)
library(patchwork)
library(ggpubr)

#remotes::install_version("dplyr", version="1.0.10")



#geographic data plotting
library("rnaturalearth")
library("rnaturalearthhires")
library("rnaturalearthdata")

#resources for GWR
#https://rstudio-pubs-static.s3.amazonaws.com/44975_0342ec49f925426fa16ebcdc28210118.html
#https://statswithr.github.io/book/bayesian-model-selection.html


##############################
#Load cane toad metadata
##############################
md <- read.table("metadata_individuals", sep="\t", header=TRUE)
md_pop <- md[!duplicated(md$pop),][,-1]


##############################
#read in environmental data
##############################
env <- read.table("env_filtered", header=TRUE, sep="\t")


##############################
#read in locality allele frequency data
##############################
freq <- read.table("allele_freq", header=TRUE, sep="\t")


##############################
#PCA using localitiy allele frequencies
##############################

#save copy of allele frequency matrix
freq_pca <- freq

#there are 20 instances where a population is missing allele frequency data at a SNP
##replace with median allele frequency for that SNP
sum(is.na(freq_pca))
for (i in 1:ncol(freq_pca)){
  freq_pca[which(is.na(freq_pca[,i])),i] <- median(freq_pca[-which(is.na(freq_pca[,i])),i], na.rm=TRUE)
}

#run PCA
pca <- rda(freq_pca[,-1], scale=T) # PCA in vegan uses the rda() call without any predictors

#make screeplot
screeplot(pca, type = "barplot", npcs=15, main="PCA Eigenvalues")

#get % variance explained by first two eigenvalues
pca$CA$eig[1]/sum(pca$CA$eig)
pca$CA$eig[2]/sum(pca$CA$eig)

#get locality PC scores for first two axes
##represent neutral population structure for GWR
PCs <- scores(pca, choices=c(1:2), display="sites", scaling=0)
PopStruct <- data.frame(pop=freq[,1], PCs)
colnames(PopStruct) <- c("pop", "PC1", "PC2")


##############################
#Geographically weighted regression (GWR)
##############################

#format data to input into GWR
env_df <- data.frame(env[,-c(1:3)], PopStruct[,-1], freq_pca[,-1])
xy <- env[,2:3]
spdf <- SpatialPointsDataFrame(coords = xy, data = env_df,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

#calculate distance matrix for GWR
gw.dist <- gw.dist(dp.locat=as.matrix(env[,2:3]), longlat=TRUE)

#establish possible models
model_vars_univ <- list(
  c("PC1", "PC2", "year"),
  c("PC1", "PC2", "year", "bio2"),
  c("PC1", "PC2", "year", "bio4"),
  c("PC1", "PC2", "year", "bio12"),
  c("PC1", "PC2", "year", "elev")
)
paste(colnames(spdf@data)[8], " ~ ", paste0(model_vars_univ[[1]], collapse=" + "), sep="")

#establish function that runs GWR and stores models
evaluate_gwr <- function(locus_index, model_vars, spatialdf, geodist){
  covariates <- model_vars
  bw_formula <- paste(colnames(spatialdf@data)[locus_index], " ~ ", paste0(model_vars, collapse=" + "), sep="")
  bw_calc <- bw.gwr(formula=bw_formula, data=spatialdf, dMat=geodist, adaptive=FALSE)
  bw_mod <- gwr.basic(formula=bw_formula, data=spatialdf, dMat=geodist, bw=bw_calc)
  
  out <- c(formula=bw_formula, ncovariates=length(model_vars), unlist(bw_mod$GW.diagnostic))
  return(out)
}

#establish function that compares GWR models for a given SNP using an approximate Bayes Factor
choose_top_gwr <- function(x){
  modcomp <- as.data.frame(t(sapply(X=model_vars_univ, evaluate_gwr, locus_index=(x+7), spatialdf=spdf, geodist=gw.dist, simplify="array")))
  
  bf_exp <- (as.numeric(modcomp$BIC[which(modcomp$ncovariates==3)]) - as.numeric(modcomp$BIC))/2
  bf_exp_inv <- (as.numeric(modcomp$BIC) - as.numeric(modcomp$BIC[which(modcomp$ncovariates==3)]))/2
  modcomp$bf <- exp(bf_exp)
  modcomp$bf_inv <- exp(bf_exp_inv)
  
  modcomp_ordered <- modcomp[order(modcomp$bf, decreasing=TRUE),]
  modcomp_ordered$snp <- x
  
  return(modcomp_ordered[1,])
}

###
#run GWR
x <- 1:ncol(freq[,-1])

#set up cluster
cl <- makeCluster(8)
clusterExport(cl, c("evaluate_gwr", "model_vars_univ", "spdf", "gw.dist"))
clusterEvalQ(cl, library("GWmodel"))

#run GWR
gwr_run <- parLapply(X=x, fun=choose_top_gwr, cl)
stopCluster(cl)

#save raw results
#saveRDS(gwr_run, "gwr_output_03032023")
gwr_result <- readRDS("gwr_output_03032023")


###
#reformat GWR output
gwr_result_df <- as.data.frame(do.call(rbind, gwr_result))
#write.table(gwr_result_df, "GWR_results_03032023_df", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

#if not re-running GWR, load in our data
setwd("C:/Users/icepi/Documents/GitHub/canetoad_landscape_genomics/analyses/gea_tests/gwr_with_slope_comparison")
gwr_result_df <- read.table("GWR_results_03032023_df", header=TRUE, sep="\t")

#identify SNPs with Bayes Factors in the top 0.05 quantile
gwr_result_df_sig <- gwr_result_df[which(gwr_result_df$bf>quantile(gwr_result_df$bf, 0.95)),]
min(gwr_result_df_sig$bf)


##############################
#overlap significant SNPs identified by GWR with those identified by Bayenv2
##############################
#load in significant BayEnv GEAs
be_sig <- sort(unique(unlist(readRDS("be_gea_sig_kf_r95"))))

#keep only the SNPs in the top 5% GWR results that ALSO overlap with BayEnv detections
sig_intersect <- gwr_result_df[intersect(gwr_result_df_sig$snp, be_sig),]


##############################
#save SNPs showing signatures of selection in order to evaluate genetic diversity
##############################

#load locus information
##uses locus missingness output from vcftools
lmiss <- read.table("ct_lmiss.lmiss", sep="\t", header=TRUE)

#intersect GWR and Bayenv2 significant SNPs
be_gwr_snp_intersect <- lmiss[sig_intersect$snp,1:2]

#obtain the total set of significant SNPs identified by EITHER GWR or Bayenv2
be_gwr_snp_all <- lmiss[sort(unique(c(gwr_result_df_sig$snp, be_sig))), 1:2]

#save SNPs without significant GEAs to a dataset of putatively neutral loci
snp_neutral <- lmiss[-sort(unique(c(gwr_result_df_sig$snp, be_sig))), 1:2]

#write outputs
#write.table(be_gwr_snp_intersect, "be_gwr_snp_intersect", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
#write.table(be_gwr_snp_all, "be_gwr_snp_all", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
#write.table(snp_neutral, "snp_neutral", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


##############################
#run permutation test to evaluate significance of observed overlap between GWR and BayEnv2
##############################

#represent SNPs based on their indices
snp_indices <- 1:5723

#establish permutation test function
gea_perm <- function(x){
  gwr_sim <- sample(x=x, size=nrow(gwr_result_df_sig), replace=FALSE)
  be_sim <- sample(x=x, size=length(be_sig), replace=FALSE)
  overlap <- length(intersect(gwr_sim, be_sim))
  return(overlap)
}

#run 500,000 permutation replicates
gea_perm_res <- replicate(n=500000, expr=gea_perm(x=snp_indices))

#calculate p-value of observed overlap
sum(gea_perm_res>=nrow(sig_intersect))/length(gea_perm_res)

#plot permutation replicates and observed overlap
perm_test_overlap <- data.frame(overlap=gea_perm_res)

#set plot margins
plot_marg <- 20

p_permtest<-ggplot(perm_test_overlap, mapping=aes(x=overlap))+
  geom_histogram(color="black", fill="gray70", binwidth=1)+
  geom_vline(xintercept=nrow(sig_intersect), colour="red", size=1)+
  scale_x_continuous(limits=c(-5,30), breaks=seq(from=-5, to=30, by=5))+
  theme_bw()+
  labs(x="Overlap between GWR and BayEnv2\nsignificant SNPs", y="Frequency")+
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18))+
  annotate(x=nrow(sig_intersect)+3,y=+Inf,label="Observed\noverlap", vjust=2, geom="label", size=6, color="red", label.size=NA)+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
p_permtest

ggsave(filename="permtest_overlap.svg", width=12, height=8)


##############################
#for overlapping SNPs, refit GWR models with center-scaled environmental data and extract beta coefficients
##############################

#format data for GWR
env_df_refit <- data.frame(scale(env[,-c(1:3)], center=TRUE, scale=TRUE), PopStruct[,-1], freq_pca[,-1])
xy <- env[,2:3]
spdf_refit <- SpatialPointsDataFrame(coords = xy, data = env_df_refit,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

###
#refit models
refit_models <- list()
for (i in 1:nrow(sig_intersect)){
  refit_gwr <- bw.gwr(formula=sig_intersect$formula[i], data=spdf_refit, dMat=gw.dist, adaptive=FALSE)
  refit_mod <- gwr.basic(formula=sig_intersect$formula[i], data=spdf_refit, dMat=gw.dist, bw=refit_gwr)
  
  refit_out <- list(gwr=refit_gwr, mod=refit_mod)
  refit_models[[i]] <- refit_out
  
}

#extract beta coefficients
refit_betas <- as.data.frame(matrix(data=NA, nrow=1, ncol=2))
colnames(refit_betas) <- c("env", "year")

for (i in 1:length(refit_models)){
  betas_temp <- data.frame(refit_models[[i]]$mod[["SDF"]]@data[,5], refit_models[[i]]$mod[["SDF"]]@data[,4])
  colnames(betas_temp) <- c("env", "year")
  
  refit_betas <- rbind(refit_betas, betas_temp)
}
refit_betas <- refit_betas[-1,]

#add snp ids
refit_betas$snp <- rep(sig_intersect$snp, each=59)

#add regions
refit_betas$region <- rep(md_pop$region_name, nrow(sig_intersect))

#add latlongs
refit_betas$long <- rep(md_pop$long, nrow(sig_intersect))
refit_betas$lat <- rep(md_pop$lat, nrow(sig_intersect))


##############################
#determine directions of GWR slopes for each SNP in each region
##############################
library(gghalves)

#get average slopes per SNP per region
refit_betas_summ_regions <- refit_betas %>%
  group_by(region, snp) %>%
  summarise(env_beta_mean=mean(env), year_beta_mean=mean(year), env_beta_se=std.error(abs(env)), year_beta_se=std.error(abs(year)))

#for each SNP, check whether the sign of NW region is the same as sign of S region
snp_list = unique(refit_betas_summ_regions$snp)
snp_same_sign = rep(NA, length(snp_list))

for (i in 1:length(snp_list)){
  data_temp = refit_betas_summ_regions[which(refit_betas_summ_regions$snp == snp_list[i]),]
  
  snp_same_sign[i] = sign(data_temp$env_beta_mean[which(data_temp$region=="NW")]) == sign(data_temp$env_beta_mean[which(data_temp$region=="S")])
  
}

#plot locality env slopes with region mean slopes as vertical lines
pal <- c("tan1", "turquoise3", "slateblue3")

ggplot(data=refit_betas, mapping=aes(x=env, y=region, fill=region, colour=region))+
  geom_point()+
  geom_vline(xintercept=0)+
  geom_vline(data=refit_betas_summ_regions, mapping=aes(xintercept=env_beta_mean, colour=as.factor(region)))+
  scale_x_continuous(limits=c(-0.5, 0.3), breaks=seq(from=-0.5, to=0.3, by=0.1))+
  facet_wrap(~snp)+
  scale_colour_manual(values=pal)+
  scale_colour_manual(values=pal)+
  theme_bw()

##############################
#plot maps of locality-specific GWR slopes for the overlapping SNPs
##############################

#get map data
world <- ne_countries(scale = "medium", returnclass = "sf")

#get ids of overlapping SNPs
overlapping_snp_ids <- unique(refit_betas$snp)

#make one map per SNP
overlap_maps <- list()
for (i in 1:length(overlapping_snp_ids)){
  overlap_maps[[i]] <- ggplot(data = world) +
    geom_sf(fill="gray80", colour="gray70") +
    coord_sf(xlim = c(127, 155), ylim = c(-30, -10))+
    labs(x="Longitude", y="Latitude", fill="Slope")+
    geom_point(data = refit_betas[which(refit_betas$snp==overlapping_snp_ids[i]),], aes(x = long, y = lat, fill=env), shape=21, colour="black", size=2.5)+
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, labels=scientific)+
    theme_bw()+
    theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
    theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
    theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
    theme(legend.key.height = unit(0.4, "cm"))+
    theme(legend.position=c(0.2, 0.25))+
    theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
  
  names(overlap_maps)[i] <- paste(paste("snp_", overlapping_snp_ids[i], sep=""))
}

#combine all maps into one plot
overlap_maps_plot <- (overlap_maps[[1]] | overlap_maps[[2]] | overlap_maps[[3]]) / (overlap_maps[[4]] | overlap_maps[[5]] | overlap_maps[[6]]) / 
  (overlap_maps[[7]] | overlap_maps[[8]] | overlap_maps[[9]]) / (overlap_maps[[10]] | overlap_maps[[11]] | overlap_maps[[12]]) / 
  (overlap_maps[[13]] | overlap_maps[[14]] | overlap_maps[[15]]) / (overlap_maps[[16]] | overlap_maps[[17]] | overlap_maps[[18]]) +
  plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 25))
overlap_maps_plot


#output plot parameters
sample_map_width <- 18
scale_factor <- 20/28
ggsave(filename="gwr_overlap_maps_midres.svg", width=sample_map_width, height=(scale_factor * 6)*(sample_map_width/3))



##############################
#Regional comparison of locality environmental beta coefficients
##############################

#summarize beta coefficients by region
refit_betas_summ <- refit_betas %>%
  group_by(region) %>%
  summarise(env_beta_mean=mean(abs(env)), year_beta_mean=mean(abs(year)), env_beta_se=std.error(abs(env)), year_beta_se=std.error(abs(year)))

#instead of summarizing everything together, first summarize to get mean beta for each SNP within each region
refit_betas_summ_snp <- refit_betas %>%
  group_by(region, snp) %>%
  summarise(env_beta_mean=mean(abs(env)), year_beta_mean=mean(abs(year)), env_beta_se=std.error(abs(env)), year_beta_se=std.error(abs(year)))

#now summarize across SNPs within regions
refit_betas_summ_snpregion <- refit_betas_summ_snp %>%
  group_by(region) %>%
  summarise(env_beta_mean_final=mean(env_beta_mean), year_beta_mean_final=mean(year_beta_mean), env_beta_se=std.error(env_beta_mean), year_beta_se=std.error(year_beta_mean))


#make plots of regional beta coefficients paired by SNP
pal <- c("tan1", "turquoise3", "slateblue3")

my_comparisons <- list( c("NE", "NW"), c("S", "NW"), c("NE", "S") )

sig_symbols <- list(cutpoints = c(0, 0.00001, 0.0001, 0.001, 0.01, 1), symbols = c("****", "***", "**", "*", "ns"))

overlap_pairplot <- ggscatter(refit_betas_summ_snp, x = "region", y = "env_beta_mean",
                            color = "white", alpha=0)+ 
                    stat_compare_means(method="wilcox.test", paired=TRUE, comparisons = my_comparisons, label="p.format", symnum.args=sig_symbols, size=7)+
                    geom_line(refit_betas_summ_snp, mapping=aes(x=region, y=env_beta_mean, group = snp), colour="gray50", size=0.67)+
                    geom_point(refit_betas_summ_snp, mapping=aes(x = region, y = env_beta_mean, fill=region), shape=21, colour="black", size=3) +
                                      scale_fill_manual(values=pal)+
                    #geom_line(refit_betas_summ_snp, mapping=aes(x=region, y=env_beta_mean, group = snp), colour="gray50", size=0.67)+
                    labs(fill="Region", colour="Region", x="Region", y="Absolute GWR\nenvironmental beta coefficient")+
                    theme_bw()+
                    theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
                    theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),
                          axis.title.x=element_text(size=18), axis.title.y=element_text(size=18),
                          legend.title=element_text(size=18), legend.text=element_text(size=14))+
                    theme(legend.position = c(0.9, 0.75))
overlap_pairplot
ggsave(filename="gwr_overlap_pairplot.svg", width=12, height=8)


#paired statistical tests; summarize by both SNP and Region
refit_betas_summ_snp <- refit_betas %>%
  group_by(region, snp) %>%
  summarise(env_beta_mean=mean(abs(env)), year_beta_mean=mean(abs(year)), env_beta_se=std.error(abs(env)), year_beta_se=std.error(abs(year)))


##########################################
##relevant df is refit_betas_summ_snp

#Keep only unlinked SNPs and re-evaluate significance


#define some functions
#merge intrachrom and interchrom dfs
merge_ld_dfs <- function(intrachrom, interchrom){
  
  #reformat intrachrom df to look like interchrom df
  intrachrom_v2 <- data.frame(CHR1=intrachrom$CHR, POS1=intrachrom$POS1, CHR2=intrachrom$CHR, POS2=intrachrom$POS2, N_INDV=intrachrom$N_INDV, R.2=intrachrom$R.2)
  
  #merge
  merged_ld <- rbind(intrachrom_v2, interchrom)
  
  #rename loci using a combined chrom-pos id
  merged_ld_chrompos <- data.frame(LOC1=paste(merged_ld$CHR1, merged_ld$POS1, sep=";"), LOC2=paste(merged_ld$CHR2, merged_ld$POS2, sep=";"), R.2=merged_ld$R.2)
  
  return(merged_ld_chrompos)
}

check_linkage <- function(ld_data, LD_threshold=0.1){
  
  #save unique locus names
  loci <- unique(c(ld_data$LOC1, ld_data$LOC2))
  loci_hiLD <- data.frame(locus=loci, count=rep(0, length(loci)))
  
  #count times a locus appears in a high-LD pair
  for (i in 1:nrow(loci_hiLD)){
    focal_locus <- ld_data[which(ld_data$LOC1==loci_hiLD$locus[i] | ld_data$LOC2==loci_hiLD$locus[i]),]
    focal_locus <- focal_locus[which(focal_locus$R.2 > LD_threshold),]
    loci_hiLD$count[i] <- nrow(focal_locus)
  }
  
  return(loci_hiLD)
  
}

###
#read in LD data produced by VCFtools
nw_ld1 <- read.table("ct_nw_gea18snps.geno.ld", sep="\t", header=TRUE)
nw_ld2 <- read.table("ct_nw_gea18snps.interchrom.geno.ld", sep="\t", header=TRUE)

ne_ld1 <- read.table("ct_ne_gea18snps.geno.ld", sep="\t", header=TRUE)
ne_ld2 <- read.table("ct_ne_gea18snps.interchrom.geno.ld", sep="\t", header=TRUE)

s_ld1 <- read.table("ct_s_gea18snps.geno.ld", sep="\t", header=TRUE)
s_ld2 <- read.table("ct_s_gea18snps.interchrom.geno.ld", sep="\t", header=TRUE)


#merged dfs
nw_ld <- merge_ld_dfs(intrachrom=nw_ld1, interchrom=nw_ld2)
ne_ld <- merge_ld_dfs(intrachrom=ne_ld1, interchrom=ne_ld2)
s_ld <- merge_ld_dfs(intrachrom=s_ld1, interchrom=s_ld2)


#successively remove loci with highest number of high-LD pairs
LD_threshold <- 0.1
nw_df <- nw_ld
ne_df <- ne_ld
s_df <- s_ld
bad_locus <- NA


for(i in 1:100){
  
  print(paste(i))
  
  nw_linkage <- check_linkage(ld_data=nw_df, LD_threshold=LD_threshold)
  ne_linkage <- check_linkage(ld_data=ne_df, LD_threshold=LD_threshold)
  s_linkage <- check_linkage(ld_data=s_df, LD_threshold=LD_threshold)
  all_linkage <- data.frame(locus=nw_linkage$locus, count_sum=(nw_linkage$count+ne_linkage$count+s_linkage$count))
  
  if(max(all_linkage$count_sum)[1] > 0){
    
    bad_locus <- all_linkage$locus[which(all_linkage$count_sum == max(all_linkage$count_sum))[1]]
    
    nw_df <- nw_df[-which(nw_df$LOC1==bad_locus| nw_df$LOC2==bad_locus),]
    ne_df <- ne_df[-which(ne_df$LOC1==bad_locus| ne_df$LOC2==bad_locus),]
    s_df <- s_df[-which(s_df$LOC1==bad_locus| s_df$LOC2==bad_locus),]
    
    bad_locus <- NA
    
  }else{
    break
  }
  
}

#snps to keep for paired test
library(stringr)

unlinked_loci_id <- str_split_fixed(nw_linkage$locus, fixed(";"), n=2)
colnames(unlinked_loci_id) <- c("CHROM", "POS")


#match unlinked SNP info to the significant SNPs
ct_loci <- read.table("ct_lmiss.lmiss", header=TRUE, sep="\t")

locus_unlinked_index <- rep(NA, nrow(unlinked_loci_id))
for(i in 1:nrow(unlinked_loci_id)){
  locus_unlinked_index[i] <- which(paste(unlinked_loci_id[i,1], unlinked_loci_id[i,2], sep=";") == 
                                     paste(ct_loci[,1], ct_loci[,2], sep=";"))
}

refit_betas_summ_snp_unlinked <- refit_betas_summ_snp
refit_betas_summ_snp_unlinked <- refit_betas_summ_snp_unlinked[which(refit_betas_summ_snp_unlinked$snp%in%locus_unlinked_index),]

#test for normality of differences
#NW - NE comparison just barely fails shapiro-wilk test, so must use nonparametric test
shapiro.test(refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="NW")] - refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="NE")])
shapiro.test(refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="NW")] - refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="S")])
shapiro.test(refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="NE")] - refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="S")])


#Wilcoxon signed-rank test
wilcox_p_adj <- 0.05/3 #Bonferroni correction

wilcox_nw_ne <- wilcox.test(refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="NW")],
                            refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="NE")],
                            paired = TRUE, alternative = "two.sided")

wilcox_nw_s <- wilcox.test(refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="NW")],
                           refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="S")],
                           paired = TRUE, alternative = "two.sided")

wilcox_ne_s <- wilcox.test(refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="NE")],
                           refit_betas_summ_snp_unlinked$env_beta_mean[which(refit_betas_summ_snp_unlinked$region=="S")],
                           paired = TRUE, alternative = "two.sided")

paste("NW:NE", "signed-rank sum T = ", wilcox_nw_ne$statistic, ";", "p =", wilcox_nw_ne$p.value, sep=" ")
paste("NW:S", "signed-rank sum T = ", wilcox_nw_s$statistic, ";", "p =", wilcox_nw_s$p.value, sep=" ")
paste("NE:S", "signed-rank sum T = ", wilcox_ne_s$statistic, ";", "p =", wilcox_ne_s$p.value, sep=" ")

#make plots of regional beta coefficients paired by SNP
pal <- c("tan1", "turquoise3", "slateblue3")

my_comparisons <- list( c("NE", "NW"), c("S", "NW"), c("NE", "S") )

sig_symbols <- list(cutpoints = c(0, 0.00001, 0.0001, 0.001, 0.01, 0.016667), symbols = c("****", "***", "**", "*", "ns"))

overlap_pairplot <- ggscatter(refit_betas_summ_snp_unlinked, x = "region", y = "env_beta_mean",
                              color = "white", alpha=0)+ 
  stat_compare_means(method="wilcox.test", paired=TRUE, comparisons = my_comparisons, label="p.format", symnum.args=sig_symbols, size=7)+
  geom_line(refit_betas_summ_snp_unlinked, mapping=aes(x=region, y=env_beta_mean, group = snp), colour="gray50", size=0.67)+
  geom_point(refit_betas_summ_snp_unlinked, mapping=aes(x = region, y = env_beta_mean, fill=region), shape=21, colour="black", size=3) +
  scale_y_continuous(limits=c(0,0.4))+
  scale_fill_manual(values=pal)+
  #geom_line(refit_betas_summ_snp_unlinked, mapping=aes(x=region, y=env_beta_mean, group = snp), colour="gray50", size=0.67)+
  labs(fill="Region", colour="Region", x="Region", y="Absolute GWR\nenvironmental beta coefficient")+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=14))+
  theme(legend.position = c(0.9, 0.75))
overlap_pairplot
ggsave(filename="gwr_overlap_unlinked_pairplot.svg", width=8, height=8)


