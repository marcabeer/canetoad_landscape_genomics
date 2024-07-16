
library(diveRsity)
library(geodist)

library(Matrix)
library(corMLPE)
library(MuMIn)

library(tibble)
library(tidyr)
library(dplyr)
library(SciViews)
library(ggplot2)
library(ggpubr)
library(patchwork)

setwd("C:/Users/icepi/Documents/Github/canetoad_landscape_genomics/analyses/pop_struct/fst")

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
fst<-readRDS("./fst/fst")


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
#Plot isolation by distance (IBD)
## NOTE this is a first pass at visualizing the data, not a proper analysis (see below)
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
#Prepare geographic data
############################################################

###
#calculate geodesic distances
nw_latlong = data.frame(long=md_pop$long[which(md_pop$region_name=="NW")], lat=md_pop$lat[which(md_pop$region_name=="NW")])
ne_latlong = data.frame(long=md_pop$long[which(md_pop$region_name=="NE")], lat=md_pop$lat[which(md_pop$region_name=="NE")])
s_latlong = data.frame(long=md_pop$long[which(md_pop$region_name=="S")], lat=md_pop$lat[which(md_pop$region_name=="S")])

#natural log for ibd analysis
nw_geodist = ln(geodist(nw_latlong, measure="geodesic")/1000)
ne_geodist = ln(geodist(ne_latlong, measure="geodesic")/1000)
s_geodist = ln(geodist(s_latlong, measure="geodesic")/1000)

#replace -Inf diagonal with 0s
diag(nw_geodist) = 0
diag(ne_geodist) = 0
diag(s_geodist) = 0

###
#isolate pairwise FST matrix
fst_pair = fst$pairwise$Fst

nw_fst = fst_pair[which(md_pop$region_name=="NW"),which(md_pop$region_name=="NW")]
ne_fst = fst_pair[which(md_pop$region_name=="NE"),which(md_pop$region_name=="NE")]
s_fst = fst_pair[which(md_pop$region_name=="S"),which(md_pop$region_name=="S")]

diag(nw_fst) <- 0
diag(ne_fst) <- 0
diag(s_fst) <- 0


#make fst matrices symmetric

nw_fst_symm <- Matrix::forceSymmetric(nw_fst, uplo="L")
nw_fst_symm_mat <- matrix(nw_fst_symm, nrow=nw_fst_symm@Dim[1], ncol=nw_fst_symm@Dim[1])
nw_fst_symm_mat[which(nw_fst_symm_mat < 0)] = 0

ne_fst_symm <- Matrix::forceSymmetric(ne_fst, uplo="L")
ne_fst_symm_mat <- matrix(ne_fst_symm, nrow=ne_fst_symm@Dim[1], ncol=ne_fst_symm@Dim[1])
ne_fst_symm_mat[which(ne_fst_symm_mat < 0)] = 0

s_fst_symm <- Matrix::forceSymmetric(s_fst, uplo="L")
s_fst_symm_mat <- matrix(s_fst_symm, nrow=s_fst_symm@Dim[1], ncol=s_fst_symm@Dim[1])
s_fst_symm_mat[which(s_fst_symm_mat < 0)] = 0


#convert fst to linearized form (fst/[1 - fst)])
nw_fst_lin <- nw_fst_symm_mat/(1-nw_fst_symm_mat)
ne_fst_lin <- ne_fst_symm_mat/(1-ne_fst_symm_mat)
s_fst_lin <- s_fst_symm_mat/(1-s_fst_symm_mat)


############################################################
#MLPE GLS test of IBD among localities
############################################################


###
#NW region data prep
nw_fst_lin_df <- data.frame(gendist=as.vector(nw_fst_lin), geodist=as.vector(nw_geodist), pop1=rep(1:nrow(nw_fst_lin), each=nrow(nw_fst_lin)), pop2=rep(1:nrow(nw_fst_lin), nrow(nw_fst_lin)))
nw_fst_lin_df <- nw_fst_lin_df[-which(nw_fst_lin_df$pop1==nw_fst_lin_df$pop2),]

#remove duplicate entries
nw_fst_lin_df$pop1pop2 <- NA

#creates a two-population index with the smaller value first
##Both entries in a square symmetric matrix will have the same two-population index, allowing filtering out duplicates 
for(i in 1:nrow(nw_fst_lin_df)){
  nw_fst_lin_df$pop1pop2[i] <- paste(min(c(nw_fst_lin_df$pop1[i], nw_fst_lin_df$pop2[i])),
                                     max(c(nw_fst_lin_df$pop1[i], nw_fst_lin_df$pop2[i])),
                                     sep="_")
}

#remove rows with duplicated two-population indices
nw_fst_lin_df_nodupl <- nw_fst_lin_df[!duplicated(nw_fst_lin_df$pop1pop2),]
nw_fst_lin_df_nodupl$species <- 1


###
#NE region data prep
ne_fst_lin_df <- data.frame(gendist=as.vector(ne_fst_lin), geodist=as.vector(ne_geodist), pop1=rep(1:nrow(ne_fst_lin), each=nrow(ne_fst_lin)), pop2=rep(1:nrow(ne_fst_lin), nrow(ne_fst_lin)))
ne_fst_lin_df <- ne_fst_lin_df[-which(ne_fst_lin_df$pop1==ne_fst_lin_df$pop2),]

#remove duplicate entries
ne_fst_lin_df$pop1pop2 <- NA

#creates a two-population index with the smaller value first
##Both entries in a square symmetric matrix will have the same two-population index, allowing filtering out duplicates 
for(i in 1:nrow(ne_fst_lin_df)){
  ne_fst_lin_df$pop1pop2[i] <- paste(min(c(ne_fst_lin_df$pop1[i], ne_fst_lin_df$pop2[i])),
                                     max(c(ne_fst_lin_df$pop1[i], ne_fst_lin_df$pop2[i])),
                                     sep="_")
}

#remove rows with duplicated two-population indices
ne_fst_lin_df_nodupl <- ne_fst_lin_df[!duplicated(ne_fst_lin_df$pop1pop2),]
ne_fst_lin_df_nodupl$species <- 1


###
#S region data prep
s_fst_lin_df <- data.frame(gendist=as.vector(s_fst_lin), geodist=as.vector(s_geodist), pop1=rep(1:nrow(s_fst_lin), each=nrow(s_fst_lin)), pop2=rep(1:nrow(s_fst_lin), nrow(s_fst_lin)))
s_fst_lin_df <- s_fst_lin_df[-which(s_fst_lin_df$pop1==s_fst_lin_df$pop2),]

#remove duplicate entries
s_fst_lin_df$pop1pop2 <- NA

#creates a two-population index with the smaller value first
##Both entries in a square symmetric matrix will have the same two-population index, allowing filtering out duplicates 
for(i in 1:nrow(s_fst_lin_df)){
  s_fst_lin_df$pop1pop2[i] <- paste(min(c(s_fst_lin_df$pop1[i], s_fst_lin_df$pop2[i])),
                                     max(c(s_fst_lin_df$pop1[i], s_fst_lin_df$pop2[i])),
                                     sep="_")
}

#remove rows with duplicated two-population indices
s_fst_lin_df_nodupl <- s_fst_lin_df[!duplicated(s_fst_lin_df$pop1pop2),]
s_fst_lin_df_nodupl$species <- 1


###
#run MLPE GLS models (with geographic distance and intercept-only null model for comparison)

mlpe_mod_stats <- data.frame(region=c("NW", "NE", "S"), AICc_ibd=NA, AICc_null=NA, dAICc_vs_null=NA, effect_geodist_ibd=NA, effect_intercept_ibd=NA, effect_intercept_null=NA)

#loop through each region (pretty useless loop, I know, but helps me organize the script)
for(i in 1:nrow(mlpe_mod_stats)){
  
  if(mlpe_mod_stats$region[i]=="NW"){
    
    #run geographic distance model
    nw_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=nw_fst_lin_df_nodupl, method="ML")
    
    #run null model
    nw_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=nw_fst_lin_df_nodupl, method="ML")

    #calculate AICc
    mlpe_mod_stats$AICc_ibd[i] = MuMIn::AICc(nw_gls_mlpe_dist)
    mlpe_mod_stats$AICc_null[i] = MuMIn::AICc(nw_gls_mlpe_null)
    
    #compare AICc of geographic distance to null model; negative means geographic distance model is better
    mlpe_mod_stats$dAICc_vs_null[i] = mlpe_mod_stats$AICc_ibd[i] - mlpe_mod_stats$AICc_null[i]
    
    #grab model coefficients
    mlpe_mod_stats$effect_geodist_ibd[i] = nw_gls_mlpe_dist$coefficients[2]
    mlpe_mod_stats$effect_intercept_ibd[i] = nw_gls_mlpe_dist$coefficients[1]
    mlpe_mod_stats$effect_intercept_null[i] = nw_gls_mlpe_null$coefficients[1]
    
  }else{
    if(mlpe_mod_stats$region[i]=="NE"){
      
      #run geographic distance model
      ne_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=ne_fst_lin_df_nodupl, method="ML")
      
      #run null model
      ne_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=ne_fst_lin_df_nodupl, method="ML")
      
      #calculate AICc
      mlpe_mod_stats$AICc_ibd[i] = MuMIn::AICc(ne_gls_mlpe_dist)
      mlpe_mod_stats$AICc_null[i] = MuMIn::AICc(ne_gls_mlpe_null)
      
      #compare AICc of geographic distance to null model; negative means geographic distance model is better
      mlpe_mod_stats$dAICc_vs_null[i] = mlpe_mod_stats$AICc_ibd[i] - mlpe_mod_stats$AICc_null[i]
      
      #grab model coefficients
      mlpe_mod_stats$effect_geodist_ibd[i] = ne_gls_mlpe_dist$coefficients[2]
      mlpe_mod_stats$effect_intercept_ibd[i] = ne_gls_mlpe_dist$coefficients[1]
      mlpe_mod_stats$effect_intercept_null[i] = ne_gls_mlpe_null$coefficients[1]
      
      
    }else{
      if(mlpe_mod_stats$region[i]=="S"){
        
        #run geographic distance model
        s_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=s_fst_lin_df_nodupl, method="ML")
        
        #run null model
        s_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=s_fst_lin_df_nodupl, method="ML")
        
        #calculate AICc
        mlpe_mod_stats$AICc_ibd[i] = MuMIn::AICc(s_gls_mlpe_dist)
        mlpe_mod_stats$AICc_null[i] = MuMIn::AICc(s_gls_mlpe_null)
        
        #compare AICc of geographic distance to null model; negative means geographic distance model is better
        mlpe_mod_stats$dAICc_vs_null[i] = mlpe_mod_stats$AICc_ibd[i] - mlpe_mod_stats$AICc_null[i]
        
        #grab model coefficients
        mlpe_mod_stats$effect_geodist_ibd[i] = s_gls_mlpe_dist$coefficients[2]
        mlpe_mod_stats$effect_intercept_ibd[i] = s_gls_mlpe_dist$coefficients[1]
        mlpe_mod_stats$effect_intercept_null[i] = s_gls_mlpe_null$coefficients[1]
        
        
      }
    }
    
  }
  
}

#view and save output
print(tibble::as_tibble(mlpe_mod_stats))
#write.csv(mlpe_mod_stats, "mlpe_mod_stats.csv", row.names=FALSE)


###
#graph models

#predict data to plot lines
nw_newdata = data.frame(geodist=seq(from=min(nw_fst_lin_df_nodupl$geodist), to=max(nw_fst_lin_df_nodupl$geodist), by=0.01))
ne_newdata = data.frame(geodist=seq(from=min(ne_fst_lin_df_nodupl$geodist), to=max(ne_fst_lin_df_nodupl$geodist), by=0.01))
s_newdata = data.frame(geodist=seq(from=min(s_fst_lin_df_nodupl$geodist), to=max(s_fst_lin_df_nodupl$geodist), by=0.01))

nw_trend = data.frame(geodist=nw_newdata, gendist=predict(object=nw_gls_mlpe_dist, newdata=nw_newdata, se=TRUE), region="NW")
ne_trend = data.frame(geodist=ne_newdata, gendist=predict(object=ne_gls_mlpe_dist, newdata=ne_newdata, se=TRUE), region="NE")
s_trend = data.frame(geodist=s_newdata, gendist=predict(object=s_gls_mlpe_dist, newdata=s_newdata, se=TRUE), region="S")

intraregion_trends = rbind(nw_trend, ne_trend, s_trend)

#make big df of all intra-region inter-locality FST values
nw_fst_lin_df_nodupl$region = "NW"
ne_fst_lin_df_nodupl$region = "NE"
s_fst_lin_df_nodupl$region = "S"

intraregion_fst_lin_df_nodupl = rbind(nw_fst_lin_df_nodupl, ne_fst_lin_df_nodupl, s_fst_lin_df_nodupl)

pal <- c("NE"="tan1", "NW"="turquoise3", "S"="slateblue3")
plot_marg = 4

ggplot()+
  geom_line(data=intraregion_trends, mapping=aes(x=geodist, y=gendist.fit, colour=region), size=1.5)+
  scale_colour_manual(values=pal)+
  geom_ribbon(data = intraregion_trends, aes(x=geodist, ymin = gendist.fit-1.96*gendist.se.fit, ymax = gendist.fit+1.96*gendist.se.fit, fill=region), alpha=0.25)+
  geom_point(data=intraregion_fst_lin_df_nodupl, aes(x=geodist, y=gendist, fill=region), shape=21, size=3, alpha=0.75)+
  scale_fill_manual(values=pal)+
  theme_bw()+
  scale_y_continuous(limits=c(-0.025, 0.3), breaks=seq(from=-0.05, to=0.25, by=0.05), expand=c(0,0))+
  scale_x_continuous(limits=c(3.25,6.75), breaks=seq(from=3, to=6.75, by=0.5), expand=c(0,0))+
  labs(x= "Ln (geographic distance)", y="FST / 1-FST", fill="Region", colour="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))


ggsave("mlpe_ibd_localities.svg", width=12, height=8)

#plot for seminar
ggplot()+
  geom_line(data=intraregion_trends, mapping=aes(x=geodist, y=gendist.fit, colour=region), size=1.5)+
  scale_colour_manual(values=pal)+
  geom_ribbon(data = intraregion_trends, aes(x=geodist, ymin = gendist.fit-1.96*gendist.se.fit, ymax = gendist.fit+1.96*gendist.se.fit, fill=region), alpha=0.25)+
  geom_point(data=intraregion_fst_lin_df_nodupl, aes(x=geodist, y=gendist, fill=region), shape=21, size=2.5, alpha=0.667)+
  scale_fill_manual(values=pal)+
  theme_bw()+
  scale_y_continuous(limits=c(-0.025, 0.27), breaks=seq(from=-0.05, to=0.25, by=0.05), expand=c(0,0))+
  scale_x_continuous(limits=c(3.25,6.75), breaks=seq(from=3, to=6.75, by=0.5), expand=c(0,0))+
  labs(x= "Ln (geographic distance)", y="FST / 1-FST\n", fill="Region", colour="Region")+
  theme(axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))


ggsave("mlpe_ibd_localities_seminar.svg", width=10, height=4)


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

##############################
#K = 4
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
##natural log for IBD analysis
k4_geodist = ln(geodist(x=data.frame(long=pm_k4_popsumm$long, lat=pm_k4_popsumm$lat), measure="geodesic")/1000)
#replace -Inf diagonal with 0s
diag(k4_geodist) = 0

#reorganize geodesic distance matrix with correct population order
k4_dist_pops <- pm_k4_popsumm$tess3_cluster
k4_dist_popcomp <- matrix(data=NA, nrow=length(k4_dist_pops), ncol=length(k4_dist_pops))
for (i in 1:length(k4_dist_pops)){
  for (j in 1:length(k4_dist_pops)){
    k4_dist_popcomp[i,j] <- paste(min(c(k4_dist_pops[i], k4_dist_pops[j])), max(c(k4_dist_pops[i], k4_dist_pops[j])), sep="_")
  }
}

k4_geodist_df = data.frame(geodist=as.vector(k4_geodist), popcomp=as.vector(k4_dist_popcomp))

#retrieve pairwise FST values
#k4_fst <- diveRsity::diffCalc(infile="tess3_k4_genepop.txt", fst=TRUE, pairwise=TRUE)
k4_fst<-readRDS("fst_k4")

#isolate pairwise FST matrix
k4_fst_pair<-k4_fst$pairwise$Fst
diag(k4_fst_pair) <- 0

#make fst matrices symmetric
k4_fst_symm <- Matrix::forceSymmetric(k4_fst_pair, uplo="L")
k4_fst_symm_mat <- matrix(k4_fst_symm, nrow=k4_fst_symm@Dim[1], ncol=k4_fst_symm@Dim[1])
k4_fst_symm_mat[which(k4_fst_symm_mat < 0)] = 0

#convert fst to linearized form (fst/[1 - fst)])
k4_fst_lin <- k4_fst_symm_mat/(1-k4_fst_symm_mat)

#assign genetic cluster based on sample name used to represent the cluster in the genepop
k4_fst_repsamp <- substring(rownames(k4_fst_pair),1, nchar(rownames(k4_fst_pair))-1)
k4_fst_pops <- c()
for (i in 1:length(k4_fst_repsamp)){
  k4_fst_pops[i] <- pm_k4_md$tess3_cluster[which(pm_k4_md$sample==k4_fst_repsamp[i])]
}

#make pairwise pop IDs
k4_fst_popcomp <- matrix(data=NA, nrow=length(k4_fst_pops), ncol=length(k4_fst_pops))
for (i in 1:length(k4_fst_pops)){
  for (j in 1:length(k4_fst_pops)){
    k4_fst_popcomp[i,j] <- paste(min(c(k4_fst_pops[i], k4_fst_pops[j])), max(c(k4_fst_pops[i], k4_fst_pops[j])), sep="_")
  }
}

#k=4 region data prep
k4_fst_lin_df = data.frame(gendist=as.vector(k4_fst_lin), pop1pop2=as.vector(k4_fst_popcomp))

#attach geodist data
k4_fst_lin_df = merge.data.frame(x=k4_fst_lin_df, y=k4_geodist_df, by.x="pop1pop2", by.y="popcomp")

#split composite pop ID
k4_fst_lin_df = tidyr::separate(k4_fst_lin_df, col="pop1pop2", into=c("pop1", "pop2"), remove=FALSE)

#remove self comparisons
k4_fst_lin_df <- k4_fst_lin_df[-which(k4_fst_lin_df$pop1==k4_fst_lin_df$pop2),]

#remove rows with duplicated two-population indices
k4_fst_lin_df_nodupl <- k4_fst_lin_df[!duplicated(k4_fst_lin_df$pop1pop2),]


###
#run MLPE models for K=4

#create df of model stats
k4_mlpe_mod_stats <- data.frame(k_value=4, AICc_ibd=NA, AICc_null=NA, dAICc_vs_null=NA, effect_geodist_ibd=NA, effect_intercept_ibd=NA, effect_intercept_null=NA)

#run geographic distance model
k4_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=k4_fst_lin_df_nodupl, method="ML")

#run null model
k4_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=k4_fst_lin_df_nodupl, method="ML")

#calculate AICc
k4_mlpe_mod_stats$AICc_ibd = MuMIn::AICc(k4_gls_mlpe_dist)
k4_mlpe_mod_stats$AICc_null = MuMIn::AICc(k4_gls_mlpe_null)

#compare AICc of geographic distance to null model; negative means geographic distance model is better
k4_mlpe_mod_stats$dAICc_vs_null = k4_mlpe_mod_stats$AICc_ibd - k4_mlpe_mod_stats$AICc_null

#grab model coefficients
k4_mlpe_mod_stats$effect_geodist_ibd = k4_gls_mlpe_dist$coefficients[2]
k4_mlpe_mod_stats$effect_intercept_ibd = k4_gls_mlpe_dist$coefficients[1]
k4_mlpe_mod_stats$effect_intercept_null = k4_gls_mlpe_null$coefficients[1]

#predict data to plot lines
k4_newdata = data.frame(geodist=seq(from=min(k4_fst_lin_df_nodupl$geodist), to=max(k4_fst_lin_df_nodupl$geodist), by=0.01))
k4_trend = data.frame(k_value=4, geodist=k4_newdata, gendist=predict(object=k4_gls_mlpe_dist, newdata=k4_newdata, se=TRUE))

ibd_k4 = ggplot()+
  geom_line(data=k4_trend, mapping=aes(x=geodist, y=gendist.fit), size=1.5)+
  geom_ribbon(data = k4_trend, aes(x=geodist, ymin = gendist.fit-1.96*gendist.se.fit, ymax = gendist.fit+1.96*gendist.se.fit), alpha=0.25)+
  geom_point(data=k4_fst_lin_df_nodupl, aes(x=geodist, y=gendist), size=3, alpha=1)+
  theme_bw()+
  scale_x_continuous(limits=c(4,8), breaks=seq(from=4, to=8, by=0.5))+
  scale_y_continuous(limits=c(0,0.2), breaks=seq(from=0, to=0.2, by=0.05))+
  labs(x= "Ln (geographic distance)", y="FST / (1 - FST)", fill="Region", colour="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))



##############################
#K = 5
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
##natural log for IBD analysis
k5_geodist = ln(geodist(x=data.frame(long=pm_k5_popsumm$long, lat=pm_k5_popsumm$lat), measure="geodesic")/1000)
#replace -Inf diagonal with 0s
diag(k5_geodist) = 0

#reorganize geodesic distance matrix with correct population order
k5_dist_pops <- pm_k5_popsumm$tess3_cluster
k5_dist_popcomp <- matrix(data=NA, nrow=length(k5_dist_pops), ncol=length(k5_dist_pops))
for (i in 1:length(k5_dist_pops)){
  for (j in 1:length(k5_dist_pops)){
    k5_dist_popcomp[i,j] <- paste(min(c(k5_dist_pops[i], k5_dist_pops[j])), max(c(k5_dist_pops[i], k5_dist_pops[j])), sep="_")
  }
}

k5_geodist_df = data.frame(geodist=as.vector(k5_geodist), popcomp=as.vector(k5_dist_popcomp))

#retrieve pairwise FST values
#k5_fst <- diveRsity::diffCalc(infile="tess3_k5_genepop.txt", fst=TRUE, pairwise=TRUE)
k5_fst<-readRDS("fst_k5")

#isolate pairwise FST matrix
k5_fst_pair<-k5_fst$pairwise$Fst
diag(k5_fst_pair) <- 0

#make fst matrices symmetric
k5_fst_symm <- Matrix::forceSymmetric(k5_fst_pair, uplo="L")
k5_fst_symm_mat <- matrix(k5_fst_symm, nrow=k5_fst_symm@Dim[1], ncol=k5_fst_symm@Dim[1])
k5_fst_symm_mat[which(k5_fst_symm_mat < 0)] = 0

#convert fst to linearized form (fst/[1 - fst)])
k5_fst_lin <- k5_fst_symm_mat/(1-k5_fst_symm_mat)

#assign genetic cluster based on sample name used to represent the cluster in the genepop
k5_fst_repsamp <- substring(rownames(k5_fst_pair),1, nchar(rownames(k5_fst_pair))-1)
k5_fst_pops <- c()
for (i in 1:length(k5_fst_repsamp)){
  k5_fst_pops[i] <- pm_k5_md$tess3_cluster[which(pm_k5_md$sample==k5_fst_repsamp[i])]
}

#make pairwise pop IDs
k5_fst_popcomp <- matrix(data=NA, nrow=length(k5_fst_pops), ncol=length(k5_fst_pops))
for (i in 1:length(k5_fst_pops)){
  for (j in 1:length(k5_fst_pops)){
    k5_fst_popcomp[i,j] <- paste(min(c(k5_fst_pops[i], k5_fst_pops[j])), max(c(k5_fst_pops[i], k5_fst_pops[j])), sep="_")
  }
}

#k=5 data prep
k5_fst_lin_df = data.frame(gendist=as.vector(k5_fst_lin), pop1pop2=as.vector(k5_fst_popcomp))

#attach geodist data
k5_fst_lin_df = merge.data.frame(x=k5_fst_lin_df, y=k5_geodist_df, by.x="pop1pop2", by.y="popcomp")

#split composite pop ID
k5_fst_lin_df = tidyr::separate(k5_fst_lin_df, col="pop1pop2", into=c("pop1", "pop2"), remove=FALSE)

#remove self comparisons
k5_fst_lin_df <- k5_fst_lin_df[-which(k5_fst_lin_df$pop1==k5_fst_lin_df$pop2),]

#remove rows with duplicated two-population indices
k5_fst_lin_df_nodupl <- k5_fst_lin_df[!duplicated(k5_fst_lin_df$pop1pop2),]


###
#run MLPE models for K=5

#create df of model stats
k5_mlpe_mod_stats <- data.frame(k_value=5, AICc_ibd=NA, AICc_null=NA, dAICc_vs_null=NA, effect_geodist_ibd=NA, effect_intercept_ibd=NA, effect_intercept_null=NA)

#run geographic distance model
k5_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=k5_fst_lin_df_nodupl, method="ML")

#run null model
k5_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=k5_fst_lin_df_nodupl, method="ML")

#calculate AICc
k5_mlpe_mod_stats$AICc_ibd = MuMIn::AICc(k5_gls_mlpe_dist)
k5_mlpe_mod_stats$AICc_null = MuMIn::AICc(k5_gls_mlpe_null)

#compare AICc of geographic distance to null model; negative means geographic distance model is better
k5_mlpe_mod_stats$dAICc_vs_null = k5_mlpe_mod_stats$AICc_ibd - k5_mlpe_mod_stats$AICc_null

#grab model coefficients
k5_mlpe_mod_stats$effect_geodist_ibd = k5_gls_mlpe_dist$coefficients[2]
k5_mlpe_mod_stats$effect_intercept_ibd = k5_gls_mlpe_dist$coefficients[1]
k5_mlpe_mod_stats$effect_intercept_null = k5_gls_mlpe_null$coefficients[1]

#predict data to plot lines
k5_newdata = data.frame(geodist=seq(from=min(k5_fst_lin_df_nodupl$geodist), to=max(k5_fst_lin_df_nodupl$geodist), by=0.01))
k5_trend = data.frame(k_value=5, geodist=k5_newdata, gendist=predict(object=k5_gls_mlpe_dist, newdata=k5_newdata, se=TRUE))

ibd_k5 = ggplot()+
  geom_line(data=k5_trend, mapping=aes(x=geodist, y=gendist.fit), size=1.5)+
  geom_ribbon(data = k5_trend, aes(x=geodist, ymin = gendist.fit-1.96*gendist.se.fit, ymax = gendist.fit+1.96*gendist.se.fit), alpha=0.25)+
  geom_point(data=k5_fst_lin_df_nodupl, aes(x=geodist, y=gendist), size=3, alpha=1)+
  theme_bw()+
  scale_x_continuous(limits=c(4,8), breaks=seq(from=4, to=8, by=0.5))+
  scale_y_continuous(limits=c(0,0.2), breaks=seq(from=0, to=0.2, by=0.05))+
  labs(x= "Ln (geographic distance)", y="FST / (1 - FST)", fill="Region", colour="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))



##############################
#K = 6
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
##natural log for IBD analysis
k6_geodist = ln(geodist(x=data.frame(long=pm_k6_popsumm$long, lat=pm_k6_popsumm$lat), measure="geodesic")/1000)
#replace -Inf diagonal with 0s
diag(k6_geodist) = 0

#reorganize geodesic distance matrix with correct population order
k6_dist_pops <- pm_k6_popsumm$tess3_cluster
k6_dist_popcomp <- matrix(data=NA, nrow=length(k6_dist_pops), ncol=length(k6_dist_pops))
for (i in 1:length(k6_dist_pops)){
  for (j in 1:length(k6_dist_pops)){
    k6_dist_popcomp[i,j] <- paste(min(c(k6_dist_pops[i], k6_dist_pops[j])), max(c(k6_dist_pops[i], k6_dist_pops[j])), sep="_")
  }
}

k6_geodist_df = data.frame(geodist=as.vector(k6_geodist), popcomp=as.vector(k6_dist_popcomp))

#retrieve pairwise FST values
#k6_fst <- diveRsity::diffCalc(infile="tess3_k6_genepop.txt", fst=TRUE, pairwise=TRUE)
k6_fst<-readRDS("fst_k6")

#isolate pairwise FST matrix
k6_fst_pair<-k6_fst$pairwise$Fst
diag(k6_fst_pair) <- 0

#make fst matrices symmetric
k6_fst_symm <- Matrix::forceSymmetric(k6_fst_pair, uplo="L")
k6_fst_symm_mat <- matrix(k6_fst_symm, nrow=k6_fst_symm@Dim[1], ncol=k6_fst_symm@Dim[1])
k6_fst_symm_mat[which(k6_fst_symm_mat < 0)] = 0

#convert fst to linearized form (fst/[1 - fst)])
k6_fst_lin <- k6_fst_symm_mat/(1-k6_fst_symm_mat)

#assign genetic cluster based on sample name used to represent the cluster in the genepop
k6_fst_repsamp <- substring(rownames(k6_fst_pair),1, nchar(rownames(k6_fst_pair))-1)
k6_fst_pops <- c()
for (i in 1:length(k6_fst_repsamp)){
  k6_fst_pops[i] <- pm_k6_md$tess3_cluster[which(pm_k6_md$sample==k6_fst_repsamp[i])]
}

#make pairwise pop IDs
k6_fst_popcomp <- matrix(data=NA, nrow=length(k6_fst_pops), ncol=length(k6_fst_pops))
for (i in 1:length(k6_fst_pops)){
  for (j in 1:length(k6_fst_pops)){
    k6_fst_popcomp[i,j] <- paste(min(c(k6_fst_pops[i], k6_fst_pops[j])), max(c(k6_fst_pops[i], k6_fst_pops[j])), sep="_")
  }
}

#k=6 data prep
k6_fst_lin_df = data.frame(gendist=as.vector(k6_fst_lin), pop1pop2=as.vector(k6_fst_popcomp))

#attach geodist data
k6_fst_lin_df = merge.data.frame(x=k6_fst_lin_df, y=k6_geodist_df, by.x="pop1pop2", by.y="popcomp")

#split composite pop ID
k6_fst_lin_df = tidyr::separate(k6_fst_lin_df, col="pop1pop2", into=c("pop1", "pop2"), remove=FALSE)

#remove self comparisons
k6_fst_lin_df <- k6_fst_lin_df[-which(k6_fst_lin_df$pop1==k6_fst_lin_df$pop2),]

#remove rows with duplicated two-population indices
k6_fst_lin_df_nodupl <- k6_fst_lin_df[!duplicated(k6_fst_lin_df$pop1pop2),]


###
#run MLPE models for K=6

#create df of model stats
k6_mlpe_mod_stats <- data.frame(k_value=6, AICc_ibd=NA, AICc_null=NA, dAICc_vs_null=NA, effect_geodist_ibd=NA, effect_intercept_ibd=NA, effect_intercept_null=NA)

#run geographic distance model
k6_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=k6_fst_lin_df_nodupl, method="ML")

#run null model
k6_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=k6_fst_lin_df_nodupl, method="ML")

#calculate AICc
k6_mlpe_mod_stats$AICc_ibd = MuMIn::AICc(k6_gls_mlpe_dist)
k6_mlpe_mod_stats$AICc_null = MuMIn::AICc(k6_gls_mlpe_null)

#compare AICc of geographic distance to null model; negative means geographic distance model is better
k6_mlpe_mod_stats$dAICc_vs_null = k6_mlpe_mod_stats$AICc_ibd - k6_mlpe_mod_stats$AICc_null

#grab model coefficients
k6_mlpe_mod_stats$effect_geodist_ibd = k6_gls_mlpe_dist$coefficients[2]
k6_mlpe_mod_stats$effect_intercept_ibd = k6_gls_mlpe_dist$coefficients[1]
k6_mlpe_mod_stats$effect_intercept_null = k6_gls_mlpe_null$coefficients[1]

#predict data to plot lines
k6_newdata = data.frame(geodist=seq(from=min(k6_fst_lin_df_nodupl$geodist), to=max(k6_fst_lin_df_nodupl$geodist), by=0.01))
k6_trend = data.frame(k_value=6, geodist=k6_newdata, gendist=predict(object=k6_gls_mlpe_dist, newdata=k6_newdata, se=TRUE))

ibd_k6 = ggplot()+
  geom_line(data=k6_trend, mapping=aes(x=geodist, y=gendist.fit), size=1.5)+
  geom_ribbon(data = k6_trend, aes(x=geodist, ymin = gendist.fit-1.96*gendist.se.fit, ymax = gendist.fit+1.96*gendist.se.fit), alpha=0.25)+
  geom_point(data=k6_fst_lin_df_nodupl, aes(x=geodist, y=gendist), size=3, alpha=1)+
  theme_bw()+
  scale_x_continuous(limits=c(4,8), breaks=seq(from=4, to=8, by=0.5))+
  scale_y_continuous(limits=c(0,0.2), breaks=seq(from=0, to=0.2, by=0.05))+
  labs(x= "Ln (geographic distance)", y="FST / (1 - FST)", fill="Region", colour="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))



##############################
#K = 7
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
##natural log for IBD analysis
k7_geodist = ln(geodist(x=data.frame(long=pm_k7_popsumm$long, lat=pm_k7_popsumm$lat), measure="geodesic")/1000)
#replace -Inf diagonal with 0s
diag(k7_geodist) = 0

#reorganize geodesic distance matrix with correct population order
k7_dist_pops <- pm_k7_popsumm$tess3_cluster
k7_dist_popcomp <- matrix(data=NA, nrow=length(k7_dist_pops), ncol=length(k7_dist_pops))
for (i in 1:length(k7_dist_pops)){
  for (j in 1:length(k7_dist_pops)){
    k7_dist_popcomp[i,j] <- paste(min(c(k7_dist_pops[i], k7_dist_pops[j])), max(c(k7_dist_pops[i], k7_dist_pops[j])), sep="_")
  }
}

k7_geodist_df = data.frame(geodist=as.vector(k7_geodist), popcomp=as.vector(k7_dist_popcomp))

#retrieve pairwise FST values
#k7_fst <- diveRsity::diffCalc(infile="tess3_k7_genepop.txt", fst=TRUE, pairwise=TRUE)
k7_fst<-readRDS("fst_k7")

#isolate pairwise FST matrix
k7_fst_pair<-k7_fst$pairwise$Fst
diag(k7_fst_pair) <- 0

#make fst matrices symmetric
k7_fst_symm <- Matrix::forceSymmetric(k7_fst_pair, uplo="L")
k7_fst_symm_mat <- matrix(k7_fst_symm, nrow=k7_fst_symm@Dim[1], ncol=k7_fst_symm@Dim[1])
k7_fst_symm_mat[which(k7_fst_symm_mat < 0)] = 0

#convert fst to linearized form (fst/[1 - fst)])
k7_fst_lin <- k7_fst_symm_mat/(1-k7_fst_symm_mat)

#assign genetic cluster based on sample name used to represent the cluster in the genepop
k7_fst_repsamp <- substring(rownames(k7_fst_pair),1, nchar(rownames(k7_fst_pair))-1)
k7_fst_pops <- c()
for (i in 1:length(k7_fst_repsamp)){
  k7_fst_pops[i] <- pm_k7_md$tess3_cluster[which(pm_k7_md$sample==k7_fst_repsamp[i])]
}

#make pairwise pop IDs
k7_fst_popcomp <- matrix(data=NA, nrow=length(k7_fst_pops), ncol=length(k7_fst_pops))
for (i in 1:length(k7_fst_pops)){
  for (j in 1:length(k7_fst_pops)){
    k7_fst_popcomp[i,j] <- paste(min(c(k7_fst_pops[i], k7_fst_pops[j])), max(c(k7_fst_pops[i], k7_fst_pops[j])), sep="_")
  }
}

#k=7 data prep
k7_fst_lin_df = data.frame(gendist=as.vector(k7_fst_lin), pop1pop2=as.vector(k7_fst_popcomp))

#attach geodist data
k7_fst_lin_df = merge.data.frame(x=k7_fst_lin_df, y=k7_geodist_df, by.x="pop1pop2", by.y="popcomp")

#split composite pop ID
k7_fst_lin_df = tidyr::separate(k7_fst_lin_df, col="pop1pop2", into=c("pop1", "pop2"), remove=FALSE)

#remove self comparisons
k7_fst_lin_df <- k7_fst_lin_df[-which(k7_fst_lin_df$pop1==k7_fst_lin_df$pop2),]

#remove rows with duplicated two-population indices
k7_fst_lin_df_nodupl <- k7_fst_lin_df[!duplicated(k7_fst_lin_df$pop1pop2),]


###
#run MLPE models for K=7

#create df of model stats
k7_mlpe_mod_stats <- data.frame(k_value=7, AICc_ibd=NA, AICc_null=NA, dAICc_vs_null=NA, effect_geodist_ibd=NA, effect_intercept_ibd=NA, effect_intercept_null=NA)

#run geographic distance model
k7_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=k7_fst_lin_df_nodupl, method="ML")

#run null model
k7_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=k7_fst_lin_df_nodupl, method="ML")

#calculate AICc
k7_mlpe_mod_stats$AICc_ibd = MuMIn::AICc(k7_gls_mlpe_dist)
k7_mlpe_mod_stats$AICc_null = MuMIn::AICc(k7_gls_mlpe_null)

#compare AICc of geographic distance to null model; negative means geographic distance model is better
k7_mlpe_mod_stats$dAICc_vs_null = k7_mlpe_mod_stats$AICc_ibd - k7_mlpe_mod_stats$AICc_null

#grab model coefficients
k7_mlpe_mod_stats$effect_geodist_ibd = k7_gls_mlpe_dist$coefficients[2]
k7_mlpe_mod_stats$effect_intercept_ibd = k7_gls_mlpe_dist$coefficients[1]
k7_mlpe_mod_stats$effect_intercept_null = k7_gls_mlpe_null$coefficients[1]

#predict data to plot lines
k7_newdata = data.frame(geodist=seq(from=min(k7_fst_lin_df_nodupl$geodist), to=max(k7_fst_lin_df_nodupl$geodist), by=0.01))
k7_trend = data.frame(k_value=7, geodist=k7_newdata, gendist=predict(object=k7_gls_mlpe_dist, newdata=k7_newdata, se=TRUE))

ibd_k7 = ggplot()+
  geom_line(data=k7_trend, mapping=aes(x=geodist, y=gendist.fit), size=1.5)+
  geom_ribbon(data = k7_trend, aes(x=geodist, ymin = gendist.fit-1.96*gendist.se.fit, ymax = gendist.fit+1.96*gendist.se.fit), alpha=0.25)+
  geom_point(data=k7_fst_lin_df_nodupl, aes(x=geodist, y=gendist), size=3, alpha=1)+
  theme_bw()+
  scale_x_continuous(limits=c(4,8), breaks=seq(from=4, to=8, by=0.5))+
  scale_y_continuous(limits=c(0,0.2), breaks=seq(from=0, to=0.2, by=0.05))+
  labs(x= "Ln (geographic distance)", y="FST / (1 - FST)", fill="Region", colour="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))



##############################
#K = 8
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
##natural log for IBD analysis
k8_geodist = ln(geodist(x=data.frame(long=pm_k8_popsumm$long, lat=pm_k8_popsumm$lat), measure="geodesic")/1000)
#replace -Inf diagonal with 0s
diag(k8_geodist) = 0

#reorganize geodesic distance matrix with correct population order
k8_dist_pops <- pm_k8_popsumm$tess3_cluster
k8_dist_popcomp <- matrix(data=NA, nrow=length(k8_dist_pops), ncol=length(k8_dist_pops))
for (i in 1:length(k8_dist_pops)){
  for (j in 1:length(k8_dist_pops)){
    k8_dist_popcomp[i,j] <- paste(min(c(k8_dist_pops[i], k8_dist_pops[j])), max(c(k8_dist_pops[i], k8_dist_pops[j])), sep="_")
  }
}

k8_geodist_df = data.frame(geodist=as.vector(k8_geodist), popcomp=as.vector(k8_dist_popcomp))

#retrieve pairwise FST values
#k8_fst <- diveRsity::diffCalc(infile="tess3_k8_genepop.txt", fst=TRUE, pairwise=TRUE)
k8_fst<-readRDS("fst_k8")

#isolate pairwise FST matrix
k8_fst_pair<-k8_fst$pairwise$Fst
diag(k8_fst_pair) <- 0

#make fst matrices symmetric
k8_fst_symm <- Matrix::forceSymmetric(k8_fst_pair, uplo="L")
k8_fst_symm_mat <- matrix(k8_fst_symm, nrow=k8_fst_symm@Dim[1], ncol=k8_fst_symm@Dim[1])
k8_fst_symm_mat[which(k8_fst_symm_mat < 0)] = 0

#convert fst to linearized form (fst/[1 - fst)])
k8_fst_lin <- k8_fst_symm_mat/(1-k8_fst_symm_mat)

#assign genetic cluster based on sample name used to represent the cluster in the genepop
k8_fst_repsamp <- substring(rownames(k8_fst_pair),1, nchar(rownames(k8_fst_pair))-1)
k8_fst_pops <- c()
for (i in 1:length(k8_fst_repsamp)){
  k8_fst_pops[i] <- pm_k8_md$tess3_cluster[which(pm_k8_md$sample==k8_fst_repsamp[i])]
}

#make pairwise pop IDs
k8_fst_popcomp <- matrix(data=NA, nrow=length(k8_fst_pops), ncol=length(k8_fst_pops))
for (i in 1:length(k8_fst_pops)){
  for (j in 1:length(k8_fst_pops)){
    k8_fst_popcomp[i,j] <- paste(min(c(k8_fst_pops[i], k8_fst_pops[j])), max(c(k8_fst_pops[i], k8_fst_pops[j])), sep="_")
  }
}

#k=8 data prep
k8_fst_lin_df = data.frame(gendist=as.vector(k8_fst_lin), pop1pop2=as.vector(k8_fst_popcomp))

#attach geodist data
k8_fst_lin_df = merge.data.frame(x=k8_fst_lin_df, y=k8_geodist_df, by.x="pop1pop2", by.y="popcomp")

#split composite pop ID
k8_fst_lin_df = tidyr::separate(k8_fst_lin_df, col="pop1pop2", into=c("pop1", "pop2"), remove=FALSE)

#remove self comparisons
k8_fst_lin_df <- k8_fst_lin_df[-which(k8_fst_lin_df$pop1==k8_fst_lin_df$pop2),]

#remove rows with duplicated two-population indices
k8_fst_lin_df_nodupl <- k8_fst_lin_df[!duplicated(k8_fst_lin_df$pop1pop2),]


###
#run MLPE models for K=8

#create df of model stats
k8_mlpe_mod_stats <- data.frame(k_value=8, AICc_ibd=NA, AICc_null=NA, dAICc_vs_null=NA, effect_geodist_ibd=NA, effect_intercept_ibd=NA, effect_intercept_null=NA)

#run geographic distance model
k8_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=k8_fst_lin_df_nodupl, method="ML")

#run null model
k8_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=k8_fst_lin_df_nodupl, method="ML")

#calculate AICc
k8_mlpe_mod_stats$AICc_ibd = MuMIn::AICc(k8_gls_mlpe_dist)
k8_mlpe_mod_stats$AICc_null = MuMIn::AICc(k8_gls_mlpe_null)

#compare AICc of geographic distance to null model; negative means geographic distance model is better
k8_mlpe_mod_stats$dAICc_vs_null = k8_mlpe_mod_stats$AICc_ibd - k8_mlpe_mod_stats$AICc_null

#grab model coefficients
k8_mlpe_mod_stats$effect_geodist_ibd = k8_gls_mlpe_dist$coefficients[2]
k8_mlpe_mod_stats$effect_intercept_ibd = k8_gls_mlpe_dist$coefficients[1]
k8_mlpe_mod_stats$effect_intercept_null = k8_gls_mlpe_null$coefficients[1]

#predict data to plot lines
k8_newdata = data.frame(geodist=seq(from=min(k8_fst_lin_df_nodupl$geodist), to=max(k8_fst_lin_df_nodupl$geodist), by=0.01))
k8_trend = data.frame(k_value=8, geodist=k8_newdata, gendist=predict(object=k8_gls_mlpe_dist, newdata=k8_newdata, se=TRUE))

ibd_k8 = ggplot()+
  geom_line(data=k8_trend, mapping=aes(x=geodist, y=gendist.fit), size=1.5)+
  geom_ribbon(data = k8_trend, aes(x=geodist, ymin = gendist.fit-1.96*gendist.se.fit, ymax = gendist.fit+1.96*gendist.se.fit), alpha=0.25)+
  geom_point(data=k8_fst_lin_df_nodupl, aes(x=geodist, y=gendist), size=3, alpha=1)+
  theme_bw()+
  scale_x_continuous(limits=c(4,8), breaks=seq(from=4, to=8, by=0.5))+
  scale_y_continuous(limits=c(0,0.2), breaks=seq(from=0, to=0.2, by=0.05))+
  labs(x= "Ln (geographic distance)", y="FST / (1 - FST)", fill="Region", colour="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))



##############################
#K = 9
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
##natural log for IBD analysis
k9_geodist = ln(geodist(x=data.frame(long=pm_k9_popsumm$long, lat=pm_k9_popsumm$lat), measure="geodesic")/1000)
#replace -Inf diagonal with 0s
diag(k9_geodist) = 0

#reorganize geodesic distance matrix with correct population order
k9_dist_pops <- pm_k9_popsumm$tess3_cluster
k9_dist_popcomp <- matrix(data=NA, nrow=length(k9_dist_pops), ncol=length(k9_dist_pops))
for (i in 1:length(k9_dist_pops)){
  for (j in 1:length(k9_dist_pops)){
    k9_dist_popcomp[i,j] <- paste(min(c(k9_dist_pops[i], k9_dist_pops[j])), max(c(k9_dist_pops[i], k9_dist_pops[j])), sep="_")
  }
}

k9_geodist_df = data.frame(geodist=as.vector(k9_geodist), popcomp=as.vector(k9_dist_popcomp))

#retrieve pairwise FST values
#k9_fst <- diveRsity::diffCalc(infile="tess3_k9_genepop.txt", fst=TRUE, pairwise=TRUE)
k9_fst<-readRDS("fst_k9")

#isolate pairwise FST matrix
k9_fst_pair<-k9_fst$pairwise$Fst
diag(k9_fst_pair) <- 0

#make fst matrices symmetric
k9_fst_symm <- Matrix::forceSymmetric(k9_fst_pair, uplo="L")
k9_fst_symm_mat <- matrix(k9_fst_symm, nrow=k9_fst_symm@Dim[1], ncol=k9_fst_symm@Dim[1])
k9_fst_symm_mat[which(k9_fst_symm_mat < 0)] = 0

#convert fst to linearized form (fst/[1 - fst)])
k9_fst_lin <- k9_fst_symm_mat/(1-k9_fst_symm_mat)

#assign genetic cluster based on sample name used to represent the cluster in the genepop
k9_fst_repsamp <- substring(rownames(k9_fst_pair),1, nchar(rownames(k9_fst_pair))-1)
k9_fst_pops <- c()
for (i in 1:length(k9_fst_repsamp)){
  k9_fst_pops[i] <- pm_k9_md$tess3_cluster[which(pm_k9_md$sample==k9_fst_repsamp[i])]
}

#make pairwise pop IDs
k9_fst_popcomp <- matrix(data=NA, nrow=length(k9_fst_pops), ncol=length(k9_fst_pops))
for (i in 1:length(k9_fst_pops)){
  for (j in 1:length(k9_fst_pops)){
    k9_fst_popcomp[i,j] <- paste(min(c(k9_fst_pops[i], k9_fst_pops[j])), max(c(k9_fst_pops[i], k9_fst_pops[j])), sep="_")
  }
}

#k=9 data prep
k9_fst_lin_df = data.frame(gendist=as.vector(k9_fst_lin), pop1pop2=as.vector(k9_fst_popcomp))

#attach geodist data
k9_fst_lin_df = merge.data.frame(x=k9_fst_lin_df, y=k9_geodist_df, by.x="pop1pop2", by.y="popcomp")

#split composite pop ID
k9_fst_lin_df = tidyr::separate(k9_fst_lin_df, col="pop1pop2", into=c("pop1", "pop2"), remove=FALSE)

#remove self comparisons
k9_fst_lin_df <- k9_fst_lin_df[-which(k9_fst_lin_df$pop1==k9_fst_lin_df$pop2),]

#remove rows with duplicated two-population indices
k9_fst_lin_df_nodupl <- k9_fst_lin_df[!duplicated(k9_fst_lin_df$pop1pop2),]


###
#run MLPE models for K=9

#create df of model stats
k9_mlpe_mod_stats <- data.frame(k_value=9, AICc_ibd=NA, AICc_null=NA, dAICc_vs_null=NA, effect_geodist_ibd=NA, effect_intercept_ibd=NA, effect_intercept_null=NA)

#run geographic distance model
k9_gls_mlpe_dist = gls(gendist ~ geodist, correlation=corMLPE(form=~pop1+pop2), data=k9_fst_lin_df_nodupl, method="ML")

#run null model
k9_gls_mlpe_null = gls(gendist ~ 1, correlation=corMLPE(form=~pop1+pop2), data=k9_fst_lin_df_nodupl, method="ML")

#calculate AICc
k9_mlpe_mod_stats$AICc_ibd = MuMIn::AICc(k9_gls_mlpe_dist)
k9_mlpe_mod_stats$AICc_null = MuMIn::AICc(k9_gls_mlpe_null)

#compare AICc of geographic distance to null model; negative means geographic distance model is better
k9_mlpe_mod_stats$dAICc_vs_null = k9_mlpe_mod_stats$AICc_ibd - k9_mlpe_mod_stats$AICc_null

#grab model coefficients
k9_mlpe_mod_stats$effect_geodist_ibd = k9_gls_mlpe_dist$coefficients[2]
k9_mlpe_mod_stats$effect_intercept_ibd = k9_gls_mlpe_dist$coefficients[1]
k9_mlpe_mod_stats$effect_intercept_null = k9_gls_mlpe_null$coefficients[1]

#predict data to plot lines
k9_newdata = data.frame(geodist=seq(from=min(k9_fst_lin_df_nodupl$geodist), to=max(k9_fst_lin_df_nodupl$geodist), by=0.01))
k9_trend = data.frame(k_value=9, geodist=k9_newdata, gendist=predict(object=k9_gls_mlpe_dist, newdata=k9_newdata, se=TRUE))

ibd_k9 = ggplot()+
  geom_line(data=k9_trend, mapping=aes(x=geodist, y=gendist.fit), size=1.5)+
  geom_ribbon(data = k9_trend, aes(x=geodist, ymin = gendist.fit-1.96*gendist.se.fit, ymax = gendist.fit+1.96*gendist.se.fit), alpha=0.25)+
  geom_point(data=k9_fst_lin_df_nodupl, aes(x=geodist, y=gendist), size=3, alpha=1)+
  theme_bw()+
  scale_x_continuous(limits=c(4,8), breaks=seq(from=4, to=8, by=0.5))+
  scale_y_continuous(limits=c(0,0.2), breaks=seq(from=0, to=0.2, by=0.05))+
  labs(x= "Ln (geographic distance)", y="FST / (1 - FST)", fill="Region", colour="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))


##############################
#Create multi-panel figure of IBD for TESS3 K = 4 - 9
##############################

ibd_plots <- (ibd_k4 | ibd_k5) / (ibd_k6 | ibd_k7) / (ibd_k8 | ibd_k9) + plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 25))
ibd_plots

plot_width <- 12
ggsave(filename="mlpe_ibd_tess3_plots.svg", width=plot_width, height=plot_width*1)


##############################
#Create output table of IBD MLPE models for TESS3 K = 4 - 9
##############################

all_k_mlpe_mod_stats = rbind(k4_mlpe_mod_stats, k5_mlpe_mod_stats, k6_mlpe_mod_stats, k7_mlpe_mod_stats, k8_mlpe_mod_stats, k9_mlpe_mod_stats)
print(all_k_mlpe_mod_stats)
#write.csv(all_k_mlpe_mod_stats, "all_k_mlpe_mod_stats.csv", row.names=FALSE)
