library(tess3r)
library(LEA)
library(maps)

library(dplyr)
library(plotrix)
library(ggplot2)
library(patchwork)
library(rworldmap)
library(RColorBrewer)

##############################
#Load data
##############################

#read in metadata
md <- read.table("metadata_individuals", sep="\t", header=TRUE)
latlong <- as.matrix(data.frame(long=md$long, lat=md$lat))

#read in SNP data and convert formats
##need to gunzip vcf first
ct_lfmm <- vcf2lfmm("ct.vcf")
ct_geno <- lfmm2geno(ct_lfmm)

#load in genotype matrix
##NA values stored as value 9; recode back to NA
ct_geno <- read.geno("ct.geno")
ct_geno[which(ct_geno==9)] <- NA


##############################
#Run TESS3
##############################

#run TESS3 for K = 1:40
##need to run this through multiple times for replication (I did 10; see below)
tess3.obj_v3 <- tess3(X=ct_geno, coord=latlong, K=1:40, method="projected.ls", ploidy=2, openMP.core.num = 8)

#save TESS3 run
saveRDS(tess3.obj_v3, file="tess3_to40_v3") #can be loaded with readRDS()

#load in TESS3 replicate runs
##these files are prohibitively large for GitHub, so skip to section below
##to read in cross-validation results based on the replicates
tess3.obj <- readRDS("tess3_to40.gz")
tess3.obj_v2 <- readRDS("tess3_to40_v2.gz")
tess3.obj_v3 <- readRDS("tess3_to40_v3.gz")
tess3.obj_v4 <- readRDS("tess3_to40_v4")
tess3.obj_v5 <- readRDS("tess3_to40_v5")
tess3.obj_v6 <- readRDS("tess3_to40_v6")
tess3.obj_v7 <- readRDS("tess3_to40_v7")
tess3.obj_v8 <- readRDS("tess3_to40_v8")
tess3.obj_v9 <- readRDS("tess3_to40_v9")
tess3.obj_v10 <- readRDS("tess3_to40_v10")


##############################
#Evaluate cross-validation results
##############################

#list each replicate run
tess3_run_lists <- list(tess3.obj, tess3.obj_v2, tess3.obj_v3, tess3.obj_v4, tess3.obj_v5, tess3.obj_v6, tess3.obj_v7, tess3.obj_v8, tess3.obj_v9, tess3.obj_v10)
tess3_run_df <- as.data.frame(matrix(data=NA, nrow=length(tess3_run_lists)*length(tess3_run_lists[[1]]), ncol=4))
tess3_run_df <- as.data.frame(matrix(data=NA, nrow=1, ncol=4))
colnames(tess3_run_df) <- c("run", "k", "rmse", "ce")

#for each run and value of K, obtain relevant performance statistics
for (i in 1:length(tess3_run_lists)){
  for (j in 1:length(tess3_run_lists[[i]])){
    run <- i
    k <- j
    rmse <- tess3_run_lists[[i]][[j]][[3]]
    ce <- tess3_run_lists[[i]][[j]][[4]]
    runk_stats <- c(run, k, rmse, ce)
    tess3_run_df <- rbind(tess3_run_df, runk_stats)
  }
}
tess3_run_df <- tess3_run_df[-1,]

#get summary statistics (based on replicate runs) for each value of K
tess3_run_df_summ <- tess3_run_df %>%
  group_by(k) %>%
  summarise(rmse_mean = mean(rmse), rmse_se = std.error(rmse), rmse_min = min(rmse), ce_mean = mean(ce), ce_se = std.error(ce), ce_min = min(ce))

#if not re-running previous code, read in per-replicate and summarized cross-validation results
tess3_run_df <- read.csv("tess3_run_df.csv")
tess3_run_df_summ <- read.csv("tess3_run_df_summary.csv")

##############################
#Plot summarized cross-validation results
##############################

#plot parameters
plot_marg <- 20

#RMSE
tess3_rmse <- ggplot()+
  geom_line(tess3_run_df_summ, mapping=aes(x=k, y=rmse_min, colour="Minimum"), size=0.75)+
  geom_point(tess3_run_df_summ, mapping=aes(x=k, y=rmse_min, colour="Minimum"), size=2.5, shape=15)+
  geom_line(tess3_run_df_summ, mapping=aes(x=k, y=rmse_mean, colour="Mean"), size=0.75)+
  geom_point(tess3_run_df_summ, mapping=aes(x=k, y=rmse_mean, colour="Mean"))+
  geom_errorbar(tess3_run_df_summ, mapping=aes(x=k, ymin=rmse_mean-rmse_se, ymax=rmse_mean+rmse_se), width=.5, position=position_dodge(.9), inherit.aes=FALSE)+
  scale_colour_manual(limits=c("Mean", "Minimum"), values=c(Minimum="red", Mean="black"))+
  scale_x_continuous(breaks=seq(from=0, to=40, by=4))+
  scale_y_continuous(breaks=seq(from=0.35, to=0.40, by=0.005))+
  labs(x="K (No. genetic clusters)", y="Root mean square error", colour="Value")+
  theme_bw()+
  theme(panel.grid.major = element_line(color="gray80", size=0.5, linetype=1),
        panel.grid.minor = element_line(color="gray95", size=0.5, linetype=1),
        legend.position=c(0.9, 0.85))+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=14))+
  theme(legend.position = c(0.75, 0.75))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#Cross entropy
tess3_ce <- ggplot()+
  geom_line(tess3_run_df_summ, mapping=aes(x=k, y=ce_min, colour="Minimum"), size=0.75)+
  geom_point(tess3_run_df_summ, mapping=aes(x=k, y=ce_min, colour="Minimum"), size=2.5, shape=15)+
  geom_line(tess3_run_df_summ, mapping=aes(x=k, y=ce_mean, colour="Mean"), size=0.75)+
  geom_point(tess3_run_df_summ, mapping=aes(x=k, y=ce_mean, colour="Mean"))+
  geom_errorbar(tess3_run_df_summ, mapping=aes(x=k, ymin=ce_mean-ce_se, ymax=ce_mean+ce_se), width=.5, position=position_dodge(.9), inherit.aes=FALSE)+
  scale_colour_manual(limits=c("Mean", "Minimum"), values=c(Minimum="red", Mean="black"))+
  scale_x_continuous(breaks=seq(from=0, to=40, by=4))+
  scale_y_continuous(breaks=seq(from=0.47, to=0.54, by=0.01))+
  labs(x="K (No. genetic clusters)", y="Cross entropy", colour="Value")+
  theme_bw()+
  theme(panel.grid.major = element_line(color="gray80", size=0.5, linetype=1),
        panel.grid.minor = element_line(color="gray95", size=0.5, linetype=1),
        legend.position=c(0.9, 0.85))+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=14))+
  theme(legend.position = c(0.75, 0.75))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#combined plot
tess3_rmse + tess3_ce + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 25))
ggsave("tess3_rmse_ce.svg", width=12, height=5)


##############################
#Plot TESS3 genetic cluster surfaces
##############################

#obtain Q matrix for K clusters
##taking the replicate with the minimum value of RMSE
##K = 4:9 had similar cross validation support, so we'll plot all six values
k9_q <- qmatrix(tess3_run_lists[[which.min(tess3_run_df[which(tess3_run_df$k==9),c("run", "rmse")][,2])]]
                , K = 9)
k8_q <- qmatrix(tess3_run_lists[[which.min(tess3_run_df[which(tess3_run_df$k==8),c("run", "rmse")][,2])]]
                , K = 8)
k7_q <- qmatrix(tess3_run_lists[[which.min(tess3_run_df[which(tess3_run_df$k==7),c("run", "rmse")][,2])]]
                , K = 7)
k6_q <- qmatrix(tess3_run_lists[[which.min(tess3_run_df[which(tess3_run_df$k==6),c("run", "rmse")][,2])]]
                , K = 6)
k5_q <- qmatrix(tess3_run_lists[[which.min(tess3_run_df[which(tess3_run_df$k==5),c("run", "rmse")][,2])]]
                , K = 5)
k4_q <- qmatrix(tess3_run_lists[[which.min(tess3_run_df[which(tess3_run_df$k==4),c("run", "rmse")][,2])]]
                , K = 4)


#plot parameters
plot_marg <- 10
map.polygon <- getMap(resolution="low")

#create color palette
##difficult because geographic location of a given genetic cluster varies across runs
##and I wanted consistent colors across values of K

#red = #FB8072
#pink = #FCCDE5
#orange = #FDB462
#yellow = #FFFFB3
#lime = #B3DE69
#green = #8DD3C7
#blue = #80B1D3
#purple = #BEBADA
#gray = #D9D9D9

k4_colors <- c("#FDB462", "#BEBADA", "#8DD3C7", "#80B1D3")
k4_palette <- CreatePalette(k4_colors, palette.length = 9)

k5_colors <- c("#FB8072", "#8DD3C7", "#BEBADA", "#FDB462", "#80B1D3")
k5_palette <- CreatePalette(k5_colors, palette.length = 9)

k6_colors <- c("#80B1D3", "#FFFFB3", "#FB8072", "#FDB462", "#BEBADA", "#8DD3C7")
k6_palette <- CreatePalette(k6_colors, palette.length = 9)

k7_colors <- c("#FFFFB3", "#FB8072", "#8DD3C7", "#B3DE69", "#BEBADA", "#80B1D3", "#FDB462")
k7_palette <- CreatePalette(k7_colors, palette.length = 9)

k8_colors <- c("#BEBADA", "#80B1D3", "#8DD3C7", "#FFFFB3", "#B3DE69", "#FDB462", "#FCCDE5", "#FB8072")
k8_palette <- CreatePalette(k8_colors, palette.length = 9)

k9_colors <- c("#FB8072", "#FDB462", "#B3DE69", "#BEBADA", "#80B1D3", "#FCCDE5", "#FFFFB3", "#8DD3C7", "#D9D9D9")
k9_palette <- CreatePalette(k9_colors, palette.length = 9)


#maps
pl_k4 <- ggtess3Q(k4_q, latlong, map.polygon=map.polygon, col.palette=k4_palette, interpol=FieldsKrigModel(10))
t3_k4 <- pl_k4 +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), colour="gray50") +
  scale_x_continuous(limits=c(127, 155), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-30, -10), expand=c(0,0)) + 
  coord_equal() + 
  geom_point(data = as.data.frame(latlong), aes(x = long, y = lat), size = 0.5, shape=19, fill="white", colour="gray30") + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
t3_k4


pl_k5 <- ggtess3Q(k5_q, latlong, map.polygon=map.polygon, col.palette=k5_palette, interpol=FieldsKrigModel(10))
t3_k5 <- pl_k5 +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), colour="gray50") +
  scale_x_continuous(limits=c(127, 155), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-30, -10), expand=c(0,0)) + 
  coord_equal() + 
  geom_point(data = as.data.frame(latlong), aes(x = long, y = lat), size = 0.5, shape=19, fill="white", colour="gray30") + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
t3_k5

pl_k6 <- ggtess3Q(k6_q, latlong, map.polygon=map.polygon, col.palette=k6_palette, interpol=FieldsKrigModel(10))
t3_k6 <- pl_k6 +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), colour="gray50") +
  scale_x_continuous(limits=c(127, 155), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-30, -10), expand=c(0,0)) + 
  coord_equal() + 
  geom_point(data = as.data.frame(latlong), aes(x = long, y = lat), size = 0.5, shape=19, fill="white", colour="gray30") + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
t3_k6

pl_k7 <- ggtess3Q(k7_q, latlong, map.polygon=map.polygon, col.palette=k7_palette, interpol=FieldsKrigModel(10))
t3_k7 <- pl_k7 +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), colour="gray50") +
  scale_x_continuous(limits=c(127, 155), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-30, -10), expand=c(0,0)) + 
  coord_equal() + 
  geom_point(data = as.data.frame(latlong), aes(x = long, y = lat), size = 0.5, shape=19, fill="white", colour="gray30") + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
t3_k7

pl_k8 <- ggtess3Q(k8_q, latlong, map.polygon=map.polygon, col.palette=k8_palette, interpol=FieldsKrigModel(10))
t3_k8 <- pl_k8 +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), colour="gray50") +
  scale_x_continuous(limits=c(127, 155), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-30, -10), expand=c(0,0)) +  
  coord_equal() + 
  geom_point(data = as.data.frame(latlong), aes(x = long, y = lat), size = 0.5, shape=19, fill="white", colour="gray30") + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
t3_k8

pl_k9 <- ggtess3Q(k9_q, latlong, map.polygon=map.polygon, col.palette=k9_palette, interpol=FieldsKrigModel(10))
t3_k9 <- pl_k9 +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), colour="gray50") +
  scale_x_continuous(limits=c(127, 155), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-30, -10), expand=c(0,0)) +  
  coord_equal() + 
  geom_point(data = as.data.frame(latlong), aes(x = long, y = lat), size = 0.5, shape=19, fill="white", colour="gray30") + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
t3_k9

(t3_k4 | t3_k5) / (t3_k6 | t3_k7) / (t3_k8 | t3_k9) + plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 25))

map_width <- 12/2
scale_factor <- 20/28
ggsave(filename="tess3_maps.svg", width=map_width*2, height=(scale_factor * 3)*(map_width))

