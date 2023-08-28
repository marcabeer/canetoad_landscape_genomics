
library(diveRsity)
library(geodist)

library(dplyr)
library(SciViews)
library(ggplot2)
library(ggpubr)
library(patchwork)

############################################################
#FST AND ISOLATION BY DISTANCE (IBD) BASED ON LOCALITIES (TESS3-BASED ANALYSIS IN SUBSEQUENT SECTION)
############################################################

##############################
#Load cane toad metadata
##############################
md <- read.table("metadata_individuals", header=TRUE, sep="\t")
md_pop <- md[!duplicated(md$pop),][,-1]

#each pair of localities, determine which regions the two localities are in
##e.g. for a pair of localities where both localities are in the NW region, the value will be "NW_NW"
reg_pairs<-t(combn(md_pop$region_name, 2))
reg_pairs<-t(combn(md_pop$region_name, 2, simplify=FALSE))
reg_pairs_comb<-paste(md_pop[,1], reg_pairs[,2], sep="_")

reg_pairs<-expand.grid(md_pop$region_name, md_pop$region_name)
reg_pairs_comb<-paste(reg_pairs[,1], reg_pairs[,2], sep="_")


##############################
#Calculate pairwise geodesic distances between localities
##############################

#division by 1000 converts to km
geodist<-geodist(x=data.frame(long=md_pop$long, lat=md_pop$lat), measure="geodesic")/1000
geodist[upper.tri(geodist, diag=TRUE)]<-NA
geodist_long<-c(unlist(geodist))


##############################
#Estimate pairwise FST between localities
##############################

#code to re-estimate FST using the R package diveRsity
#fst<-diveRsity::diffCalc(infile="localities_genepop.txt", fst=TRUE, pairwise=TRUE)
#saveRDS(fst, "fst")

#load in FST estimates used in publication
##need to unzip fst.7z first
fst<-readRDS("fst")

#isolate pairwise FST matrix
fst_pair<-fst$pairwise$Fst

#convert FST matrix to vector
fst_long<-unlist(c(fst_pair))

#calculate linearized FST for correlation with geographic distance
fst_long_lin<-fst_long/(1-fst_long)

#combine variables into df
fst_df<-data.frame(reg_pair=reg_pairs_comb, dist=geodist_long, dist_ln=ln(geodist_long), fst=fst_long, fst_lin=fst_long_lin)

#convert negative FST to 0
fst_df$fst[which(fst_df$fst<0)]<-0
fst_df$fst_lin[which(fst_df$fst_lin<0)]<-0

#remove rows with NA values
fst_df_trim<-fst_df[-which(is.na(fst_df$dist)==1),]

##############################
#Characterize isolation by distance (IBD)
##############################

#first isolate FST values corresponding to pairs of localities within regions
##will be used for linear models
within_reg_index<-c(which(fst_df_trim$reg_pair=="NW_NW"), which(fst_df_trim$reg_pair=="NE_NE"), which(fst_df_trim$reg_pair=="S_S"))

#make plot
pal2 <- c("NE_NE"="tan1", "NW_NW"="turquoise3", "S_S"="slateblue3", "NE_NW"="firebrick2", "S_NE"="dodgerblue", "S_NW"="darkolivegreen3")
pal3 <- c("NE_NE"="tan4", "NW_NW"="turquoise4", "S_S"="slateblue4")
pal_shape <- c("NE_NE"=19, "NW_NW"=19, "S_S"=19, "NE_NW"=17, "S_NE"=15, "S_NW"=18)

p_ibd <- ggplot()+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_point(data=fst_df_trim, aes(x=dist_ln, y=fst_lin, colour=reg_pair, shape=reg_pair), size=2.5, alpha=0.5)+
  scale_colour_manual(values=pal2)+
  scale_shape_manual(values=pal_shape)+
  labs(x= "Ln (geographic distance)", y="Linearized FST (FST / 1-FST)", colour="Region pair", fill="Region pair", shape="Region pair")+
  new_scale_colour()+
  new_scale_fill()+
  geom_smooth(data=fst_df_trim[within_reg_index,], aes(x=dist_ln, y=fst_lin, colour=reg_pair, fill=reg_pair), method='lm')+
  scale_colour_manual(values=pal3[1:3])+
  scale_fill_manual(values=pal2[1:3])+
  scale_y_continuous(limits=c(-0.025, 0.3), breaks=seq(from=-0.05, to=0.3, by=0.05))+
  scale_x_continuous(limits=c(3,8))+
  labs(x= "Geographic distance (Ln[km])", y="Linearized FST (FST / 1-FST)", colour="Region pair", fill="Region pair", shape="Region pair")+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
p_ibd

ggsave("ibd_localities.svg", width=12, height=8)


############################################################
#TESS3-BASED ANALYSIS OF FST AND IBD
############################################################

##############################
#load in K-specific popmaps (pm) based on TESS3 outputs
##############################

pm_k4 <- read.table("tess3_popmap_k4.txt", sep="\t", header=FALSE)
pm_k5 <- read.table("tess3_popmap_k5.txt", sep="\t", header=FALSE)
pm_k6 <- read.table("tess3_popmap_k6.txt", sep="\t", header=FALSE)
pm_k7 <- read.table("tess3_popmap_k7.txt", sep="\t", header=FALSE)
pm_k8 <- read.table("tess3_popmap_k8.txt", sep="\t", header=FALSE)
pm_k9 <- read.table("tess3_popmap_k9.txt", sep="\t", header=FALSE)

#rename columns
##painfully hard-coded -  I know, I know
pm_colnames <- c("sample", "tess3_cluster")
colnames(pm_k4) <- pm_colnames
colnames(pm_k5) <- pm_colnames
colnames(pm_k6) <- pm_colnames
colnames(pm_k7) <- pm_colnames
colnames(pm_k8) <- pm_colnames
colnames(pm_k9) <- pm_colnames

###
#establish plot margins for later code
plot_marg <- 20


##############################
#Characterize IBD for K = 4
##############################

#create metadata
pm_k4_md <- merge.data.frame(md, pm_k4, by="sample", sort=FALSE)

#get mean coordinates of each of K genetic clusters
pm_k4_popsumm <- pm_k4_md %>%
  group_by(pop) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=round(mean(tess3_cluster), 0)) %>%
  group_by(tess3_cluster) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=mean(tess3_cluster))

#calculate geodesic distance between each of K genetic clusters
k4_geodist<-geodist(x=data.frame(long=pm_k4_popsumm$long, lat=pm_k4_popsumm$lat), measure="geodesic")/1000
#k4_geodist[upper.tri(k4_geodist, diag=TRUE)]<-NA
k4_geodist_long<-c(unlist(k4_geodist))

k4_dist_pops <- pm_k4_popsumm$tess3_cluster
k4_dist_popcomp <- matrix(data=NA, nrow=length(k4_dist_pops), ncol=length(k4_dist_pops))
for (i in 1:length(k4_dist_pops)){
  for (j in 1:length(k4_dist_pops)){
    k4_dist_popcomp[i,j] <- paste(k4_dist_pops[i], k4_dist_pops[j], sep="_")
  }
}

#retrieve pairwise FST values
#k4_fst <- diveRsity::diffCalc(infile="tess3_k4_genepop.txt", fst=TRUE, pairwise=TRUE)
k4_fst<-readRDS("fst_k4")

#isolate pairwise FST matrix
k4_fst_pair<-k4_fst$pairwise$Fst

k4_fst_repsamp <- substring(rownames(k4_fst_pair),1, nchar(rownames(k4_fst_pair))-1)
k4_fst_pops <- c()
for (i in 1:length(k4_fst_repsamp)){
  k4_fst_pops[i] <- pm_k4_md$tess3_cluster[which(pm_k4_md$sample==k4_fst_repsamp[i])]
}

k4_fst_popcomp <- matrix(data=NA, nrow=length(k4_fst_pops), ncol=length(k4_fst_pops))
for (i in 1:length(k4_fst_pops)){
  for (j in 1:length(k4_fst_pops)){
    k4_fst_popcomp[i,j] <- paste(k4_fst_pops[i], k4_fst_pops[j], sep="_")
  }
}

#convert FST matrix to vector
k4_fst_long<-data.frame(comp=unlist(c(k4_fst_popcomp)), fst=unlist(c(k4_fst_pair)))

#match FST values to distance values
k4_dist_long<-data.frame(comp=unlist(c(k4_dist_popcomp)), dist=k4_geodist_long)
k4_fst_dist <- merge.data.frame(k4_fst_long, k4_dist_long, by="comp", sort=FALSE)

#calculate linearized FST for correlation with geographic distance
k4_fst_dist$fst_lin<-k4_fst_dist$fst/(1-k4_fst_dist$fst)

#calculate ln geographic distance
k4_fst_dist$dist_ln<-ln(k4_fst_dist$dist)

#convert negative FST to 0
k4_fst_dist$fst[which(k4_fst_dist$fst<0)]<-0
k4_fst_dist$fst_lin[which(k4_fst_dist$fst_lin<0)]<-0

#remove rows with NA values
k4_fst_df_trim<-k4_fst_dist[-which(is.na(k4_fst_dist$fst)==1),]

#plot linearized FST vs geodesic distance for K = 4
ibd_k4 <- ggscatter(data=k4_fst_df_trim, x="dist_ln", y="fst_lin", add="reg.line", conf.int=TRUE, ellipse.alpha=0.25)+
  stat_cor(method="pearson", label.x=4.1, label.y=0.18, show.legend=FALSE, size=6)+
  labs(x="Ln(km)", y="FST / (1 - FST)")+
  xlim(4,8)+
  ylim(0,0.2)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#view plot
ibd_k4

##############################
#Characterize IBD for K = 5
##############################

#create metadata
pm_k5_md <- merge.data.frame(md, pm_k5, by="sample", sort=FALSE)

#get mean coordinates of each of K genetic clusters
pm_k5_popsumm <- pm_k5_md %>%
  group_by(pop) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=round(mean(tess3_cluster), 0)) %>%
  group_by(tess3_cluster) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=mean(tess3_cluster))

#calculate geodesic distance between each of K genetic clusters
k5_geodist<-geodist(x=data.frame(long=pm_k5_popsumm$long, lat=pm_k5_popsumm$lat), measure="geodesic")/1000
#k5_geodist[upper.tri(k5_geodist, diag=TRUE)]<-NA
k5_geodist_long<-c(unlist(k5_geodist))

k5_dist_pops <- pm_k5_popsumm$tess3_cluster
k5_dist_popcomp <- matrix(data=NA, nrow=length(k5_dist_pops), ncol=length(k5_dist_pops))
for (i in 1:length(k5_dist_pops)){
  for (j in 1:length(k5_dist_pops)){
    k5_dist_popcomp[i,j] <- paste(k5_dist_pops[i], k5_dist_pops[j], sep="_")
  }
}

#retrieve pairwise FST values
#k5_fst <- diveRsity::diffCalc(infile="tess3_k5_genepop.txt", fst=TRUE, pairwise=TRUE)
k5_fst<-readRDS("fst_k5")

#isolate pairwise FST matrix
k5_fst_pair<-k5_fst$pairwise$Fst

k5_fst_repsamp <- substring(rownames(k5_fst_pair),1, nchar(rownames(k5_fst_pair))-1)
k5_fst_pops <- c()
for (i in 1:length(k5_fst_repsamp)){
  k5_fst_pops[i] <- pm_k5_md$tess3_cluster[which(pm_k5_md$sample==k5_fst_repsamp[i])]
}

k5_fst_popcomp <- matrix(data=NA, nrow=length(k5_fst_pops), ncol=length(k5_fst_pops))
for (i in 1:length(k5_fst_pops)){
  for (j in 1:length(k5_fst_pops)){
    k5_fst_popcomp[i,j] <- paste(k5_fst_pops[i], k5_fst_pops[j], sep="_")
  }
}

#convert FST matrix to vector
k5_fst_long<-data.frame(comp=unlist(c(k5_fst_popcomp)), fst=unlist(c(k5_fst_pair)))

#match FST values to distance values
k5_dist_long<-data.frame(comp=unlist(c(k5_dist_popcomp)), dist=k5_geodist_long)
k5_fst_dist <- merge.data.frame(k5_fst_long, k5_dist_long, by="comp", sort=FALSE)

#calculate linearized FST for correlation with geographic distance
k5_fst_dist$fst_lin<-k5_fst_dist$fst/(1-k5_fst_dist$fst)

#calculate ln geographic distance
k5_fst_dist$dist_ln<-ln(k5_fst_dist$dist)

#convert negative FST to 0
k5_fst_dist$fst[which(k5_fst_dist$fst<0)]<-0
k5_fst_dist$fst_lin[which(k5_fst_dist$fst_lin<0)]<-0

#remove rows with NA values
k5_fst_df_trim<-k5_fst_dist[-which(is.na(k5_fst_dist$fst)==1),]

#plot linearized FST vs geodesic distance for K = 4
ibd_k5 <- ggscatter(data=k5_fst_df_trim, x="dist_ln", y="fst_lin", add="reg.line", conf.int=TRUE, ellipse.alpha=0.25)+
  stat_cor(method="pearson", label.x=4.1, label.y=0.18, show.legend=FALSE, size=6)+
  labs(x="Ln(km)", y="FST / (1 - FST)")+
  xlim(4,8)+
  ylim(0,0.2)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#view plot
ibd_k5

##############################
#Characterize IBD for K = 6
##############################

#create metadata
pm_k6_md <- merge.data.frame(md, pm_k6, by="sample", sort=FALSE)

#get mean coordinates of each of K genetic clusters
pm_k6_popsumm <- pm_k6_md %>%
  group_by(pop) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=round(mean(tess3_cluster), 0)) %>%
  group_by(tess3_cluster) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=mean(tess3_cluster))

#calculate geodesic distance between each of K genetic clusters
k6_geodist<-geodist(x=data.frame(long=pm_k6_popsumm$long, lat=pm_k6_popsumm$lat), measure="geodesic")/1000
#k6_geodist[upper.tri(k6_geodist, diag=TRUE)]<-NA
k6_geodist_long<-c(unlist(k6_geodist))

k6_dist_pops <- pm_k6_popsumm$tess3_cluster
k6_dist_popcomp <- matrix(data=NA, nrow=length(k6_dist_pops), ncol=length(k6_dist_pops))
for (i in 1:length(k6_dist_pops)){
  for (j in 1:length(k6_dist_pops)){
    k6_dist_popcomp[i,j] <- paste(k6_dist_pops[i], k6_dist_pops[j], sep="_")
  }
}

#retrieve pairwise FST values
#k6_fst <- diveRsity::diffCalc(infile="tess3_k6_genepop.txt", fst=TRUE, pairwise=TRUE)
k6_fst<-readRDS("fst_k6")

#isolate pairwise FST matrix
k6_fst_pair<-k6_fst$pairwise$Fst

k6_fst_repsamp <- substring(rownames(k6_fst_pair),1, nchar(rownames(k6_fst_pair))-1)
k6_fst_pops <- c()
for (i in 1:length(k6_fst_repsamp)){
  k6_fst_pops[i] <- pm_k6_md$tess3_cluster[which(pm_k6_md$sample==k6_fst_repsamp[i])]
}

k6_fst_popcomp <- matrix(data=NA, nrow=length(k6_fst_pops), ncol=length(k6_fst_pops))
for (i in 1:length(k6_fst_pops)){
  for (j in 1:length(k6_fst_pops)){
    k6_fst_popcomp[i,j] <- paste(k6_fst_pops[i], k6_fst_pops[j], sep="_")
  }
}

#convert FST matrix to vector
k6_fst_long<-data.frame(comp=unlist(c(k6_fst_popcomp)), fst=unlist(c(k6_fst_pair)))

#match FST values to distance values
k6_dist_long<-data.frame(comp=unlist(c(k6_dist_popcomp)), dist=k6_geodist_long)
k6_fst_dist <- merge.data.frame(k6_fst_long, k6_dist_long, by="comp", sort=FALSE)

#calculate linearized FST for correlation with geographic distance
k6_fst_dist$fst_lin<-k6_fst_dist$fst/(1-k6_fst_dist$fst)

#calculate ln geographic distance
k6_fst_dist$dist_ln<-ln(k6_fst_dist$dist)

#convert negative FST to 0
k6_fst_dist$fst[which(k6_fst_dist$fst<0)]<-0
k6_fst_dist$fst_lin[which(k6_fst_dist$fst_lin<0)]<-0

#remove rows with NA values
k6_fst_df_trim<-k6_fst_dist[-which(is.na(k6_fst_dist$fst)==1),]

#plot linearized FST vs geodesic distance for K = 4
ibd_k6 <- ggscatter(data=k6_fst_df_trim, x="dist_ln", y="fst_lin", add="reg.line", conf.int=TRUE, ellipse.alpha=0.25)+
  stat_cor(method="pearson", label.x=4.1, label.y=0.18, show.legend=FALSE, size=6)+
  labs(x="Ln(km)", y="FST / (1 - FST)")+
  xlim(4,8)+
  ylim(0,0.2)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#view plot
ibd_k6

##############################
#Characterize IBD for K = 7
##############################

#create metadata
pm_k7_md <- merge.data.frame(md, pm_k7, by="sample", sort=FALSE)

#get mean coordinates of each of K genetic clusters
pm_k7_popsumm <- pm_k7_md %>%
  group_by(pop) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=round(mean(tess3_cluster), 0)) %>%
  group_by(tess3_cluster) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=mean(tess3_cluster))

#calculate geodesic distance between each of K genetic clusters
k7_geodist<-geodist(x=data.frame(long=pm_k7_popsumm$long, lat=pm_k7_popsumm$lat), measure="geodesic")/1000
#k7_geodist[upper.tri(k7_geodist, diag=TRUE)]<-NA
k7_geodist_long<-c(unlist(k7_geodist))

k7_dist_pops <- pm_k7_popsumm$tess3_cluster
k7_dist_popcomp <- matrix(data=NA, nrow=length(k7_dist_pops), ncol=length(k7_dist_pops))
for (i in 1:length(k7_dist_pops)){
  for (j in 1:length(k7_dist_pops)){
    k7_dist_popcomp[i,j] <- paste(k7_dist_pops[i], k7_dist_pops[j], sep="_")
  }
}

#retrieve pairwise FST values
#k7_fst <- diveRsity::diffCalc(infile="tess3_k7_genepop.txt", fst=TRUE, pairwise=TRUE)
k7_fst<-readRDS("fst_k7")

#isolate pairwise FST matrix
k7_fst_pair<-k7_fst$pairwise$Fst

k7_fst_repsamp <- substring(rownames(k7_fst_pair),1, nchar(rownames(k7_fst_pair))-1)
k7_fst_pops <- c()
for (i in 1:length(k7_fst_repsamp)){
  k7_fst_pops[i] <- pm_k7_md$tess3_cluster[which(pm_k7_md$sample==k7_fst_repsamp[i])]
}

k7_fst_popcomp <- matrix(data=NA, nrow=length(k7_fst_pops), ncol=length(k7_fst_pops))
for (i in 1:length(k7_fst_pops)){
  for (j in 1:length(k7_fst_pops)){
    k7_fst_popcomp[i,j] <- paste(k7_fst_pops[i], k7_fst_pops[j], sep="_")
  }
}

#convert FST matrix to vector
k7_fst_long<-data.frame(comp=unlist(c(k7_fst_popcomp)), fst=unlist(c(k7_fst_pair)))

#match FST values to distance values
k7_dist_long<-data.frame(comp=unlist(c(k7_dist_popcomp)), dist=k7_geodist_long)
k7_fst_dist <- merge.data.frame(k7_fst_long, k7_dist_long, by="comp", sort=FALSE)

#calculate linearized FST for correlation with geographic distance
k7_fst_dist$fst_lin<-k7_fst_dist$fst/(1-k7_fst_dist$fst)

#calculate ln geographic distance
k7_fst_dist$dist_ln<-ln(k7_fst_dist$dist)

#convert negative FST to 0
k7_fst_dist$fst[which(k7_fst_dist$fst<0)]<-0
k7_fst_dist$fst_lin[which(k7_fst_dist$fst_lin<0)]<-0

#remove rows with NA values
k7_fst_df_trim<-k7_fst_dist[-which(is.na(k7_fst_dist$fst)==1),]

#plot linearized FST vs geodesic distance for K = 4
ibd_k7 <- ggscatter(data=k7_fst_df_trim, x="dist_ln", y="fst_lin", add="reg.line", conf.int=TRUE, ellipse.alpha=0.25)+
  stat_cor(method="pearson", label.x=4.1, label.y=0.18, show.legend=FALSE, size=6)+
  labs(x="Ln(km)", y="FST / (1 - FST)")+
  xlim(4,8)+
  ylim(0,0.2)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#view plot
ibd_k7

##############################
#Characterize IBD for K = 8
##############################

#create metadata
pm_k8_md <- merge.data.frame(md, pm_k8, by="sample", sort=FALSE)

#get mean coordinates of each of K genetic clusters
pm_k8_popsumm <- pm_k8_md %>%
  group_by(pop) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=round(mean(tess3_cluster), 0)) %>%
  group_by(tess3_cluster) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=mean(tess3_cluster))

#calculate geodesic distance between each of K genetic clusters
k8_geodist<-geodist(x=data.frame(long=pm_k8_popsumm$long, lat=pm_k8_popsumm$lat), measure="geodesic")/1000
#k8_geodist[upper.tri(k8_geodist, diag=TRUE)]<-NA
k8_geodist_long<-c(unlist(k8_geodist))

k8_dist_pops <- pm_k8_popsumm$tess3_cluster
k8_dist_popcomp <- matrix(data=NA, nrow=length(k8_dist_pops), ncol=length(k8_dist_pops))
for (i in 1:length(k8_dist_pops)){
  for (j in 1:length(k8_dist_pops)){
    k8_dist_popcomp[i,j] <- paste(k8_dist_pops[i], k8_dist_pops[j], sep="_")
  }
}

#retrieve pairwise FST values
#k8_fst <- diveRsity::diffCalc(infile="tess3_k8_genepop.txt", fst=TRUE, pairwise=TRUE)
k8_fst<-readRDS("fst_k8")

#isolate pairwise FST matrix
k8_fst_pair<-k8_fst$pairwise$Fst

k8_fst_repsamp <- substring(rownames(k8_fst_pair),1, nchar(rownames(k8_fst_pair))-1)
k8_fst_pops <- c()
for (i in 1:length(k8_fst_repsamp)){
  k8_fst_pops[i] <- pm_k8_md$tess3_cluster[which(pm_k8_md$sample==k8_fst_repsamp[i])]
}

k8_fst_popcomp <- matrix(data=NA, nrow=length(k8_fst_pops), ncol=length(k8_fst_pops))
for (i in 1:length(k8_fst_pops)){
  for (j in 1:length(k8_fst_pops)){
    k8_fst_popcomp[i,j] <- paste(k8_fst_pops[i], k8_fst_pops[j], sep="_")
  }
}

#convert FST matrix to vector
k8_fst_long<-data.frame(comp=unlist(c(k8_fst_popcomp)), fst=unlist(c(k8_fst_pair)))

#match FST values to distance values
k8_dist_long<-data.frame(comp=unlist(c(k8_dist_popcomp)), dist=k8_geodist_long)
k8_fst_dist <- merge.data.frame(k8_fst_long, k8_dist_long, by="comp", sort=FALSE)

#calculate linearized FST for correlation with geographic distance
k8_fst_dist$fst_lin<-k8_fst_dist$fst/(1-k8_fst_dist$fst)

#calculate ln geographic distance
k8_fst_dist$dist_ln<-ln(k8_fst_dist$dist)

#convert negative FST to 0
k8_fst_dist$fst[which(k8_fst_dist$fst<0)]<-0
k8_fst_dist$fst_lin[which(k8_fst_dist$fst_lin<0)]<-0

#remove rows with NA values
k8_fst_df_trim<-k8_fst_dist[-which(is.na(k8_fst_dist$fst)==1),]

#plot linearized FST vs geodesic distance for K = 4
ibd_k8 <- ggscatter(data=k8_fst_df_trim, x="dist_ln", y="fst_lin", add="reg.line", conf.int=TRUE, ellipse.alpha=0.25)+
  stat_cor(method="pearson", label.x=4.1, label.y=0.18, show.legend=FALSE, size=6)+
  labs(x="Ln(km)", y="FST / (1 - FST)")+
  xlim(4,8)+
  ylim(0,0.2)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#view plot
ibd_k8

##############################
#Characterize IBD for K = 9
##############################

#create metadata
pm_k9_md <- merge.data.frame(md, pm_k9, by="sample", sort=FALSE)

#get mean coordinates of each of K genetic clusters
pm_k9_popsumm <- pm_k9_md %>%
  group_by(pop) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=round(mean(tess3_cluster), 0)) %>%
  group_by(tess3_cluster) %>%
  summarise(long=mean(long), lat=mean(lat), tess3_cluster=mean(tess3_cluster))

#calculate geodesic distance between each of K genetic clusters
k9_geodist<-geodist(x=data.frame(long=pm_k9_popsumm$long, lat=pm_k9_popsumm$lat), measure="geodesic")/1000
#k9_geodist[upper.tri(k9_geodist, diag=TRUE)]<-NA
k9_geodist_long<-c(unlist(k9_geodist))

k9_dist_pops <- pm_k9_popsumm$tess3_cluster
k9_dist_popcomp <- matrix(data=NA, nrow=length(k9_dist_pops), ncol=length(k9_dist_pops))
for (i in 1:length(k9_dist_pops)){
  for (j in 1:length(k9_dist_pops)){
    k9_dist_popcomp[i,j] <- paste(k9_dist_pops[i], k9_dist_pops[j], sep="_")
  }
}

#retrieve pairwise FST values
#k9_fst <- diveRsity::diffCalc(infile="tess3_k9_genepop.txt", fst=TRUE, pairwise=TRUE)
k9_fst<-readRDS("fst_k9")

#isolate pairwise FST matrix
k9_fst_pair<-k9_fst$pairwise$Fst

k9_fst_repsamp <- substring(rownames(k9_fst_pair),1, nchar(rownames(k9_fst_pair))-1)
k9_fst_pops <- c()
for (i in 1:length(k9_fst_repsamp)){
  k9_fst_pops[i] <- pm_k9_md$tess3_cluster[which(pm_k9_md$sample==k9_fst_repsamp[i])]
}

k9_fst_popcomp <- matrix(data=NA, nrow=length(k9_fst_pops), ncol=length(k9_fst_pops))
for (i in 1:length(k9_fst_pops)){
  for (j in 1:length(k9_fst_pops)){
    k9_fst_popcomp[i,j] <- paste(k9_fst_pops[i], k9_fst_pops[j], sep="_")
  }
}

#convert FST matrix to vector
k9_fst_long<-data.frame(comp=unlist(c(k9_fst_popcomp)), fst=unlist(c(k9_fst_pair)))

#match FST values to distance values
k9_dist_long<-data.frame(comp=unlist(c(k9_dist_popcomp)), dist=k9_geodist_long)
k9_fst_dist <- merge.data.frame(k9_fst_long, k9_dist_long, by="comp", sort=FALSE)

#calculate linearized FST for correlation with geographic distance
k9_fst_dist$fst_lin<-k9_fst_dist$fst/(1-k9_fst_dist$fst)

#calculate ln geographic distance
k9_fst_dist$dist_ln<-ln(k9_fst_dist$dist)

#convert negative FST to 0
k9_fst_dist$fst[which(k9_fst_dist$fst<0)]<-0
k9_fst_dist$fst_lin[which(k9_fst_dist$fst_lin<0)]<-0

#remove rows with NA values
k9_fst_df_trim<-k9_fst_dist[-which(is.na(k9_fst_dist$fst)==1),]

#plot linearized FST vs geodesic distance for K = 4
ibd_k9 <- ggscatter(data=k9_fst_df_trim, x="dist_ln", y="fst_lin", add="reg.line", conf.int=TRUE, ellipse.alpha=0.25)+
  stat_cor(method="pearson", label.x=4.1, label.y=0.18, show.legend=FALSE, size=6)+
  labs(x="Ln(km)", y="FST / (1 - FST)")+
  xlim(4,8)+
  ylim(0,0.2)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#view plot
ibd_k9

##############################
#Create multi-panel figure of IBD for TESS3 K = 4 - 9
##############################

ibd_plots <- (ibd_k4 | ibd_k5) / (ibd_k6 | ibd_k7) / (ibd_k8 | ibd_k9) + plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 25))
ibd_plots

plot_width <- 12
ggsave(filename="ibd_tess3_plots.svg", width=plot_width, height=plot_width*1)
