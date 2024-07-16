setwd("~/github_prep/analyses/env_data")

library(raster)
library(geobuffer)

library(cowplot)
library(corrplot)
library(ggplot2)
library(ggcorrplot)
library(patchwork)

##############################
#Download environmental data
##############################

#environmental data were downloaded from WorldClim database v2 at 1km resolution
##pre-processing took place to match extents
##of colonization year raster and bioclimatic variables
##code below documents extraction of environmental conditions
##at sampling localities


##############################
#Load environmental data (if not re-processing)
##############################

#data found in rasterstack
env_all<-stack("env_all.tif")
names(env_all)<-c("bio1", "bio2", "bio3", "bio4", "bio7", "bio12", "bio15", "elev", "year")


##############################
#Load cane toad metadata
##############################
md <- read.table("metadata_individual", header=TRUE, sep="\t")

#keep one entry per sampling site ("pop")
md_sites <- md[!duplicated(md$pop),]

#isolate latlongs
md_latlong <- subset(md_sites, select=c(long, lat))


##############################
#Extract environmental data at localities
##############################

#buffer latlongs
##5km radius
buf_rad<-5
latlong_buffer<-geobuffer::geobuffer_pts(xy=md_latlong, dist_m=(buf_rad*1000), step_dg=10)

#extract environmental data within buffers
env_mean<-raster::extract(x=env_all, y=latlong_buffer, fun=mean, na.rm=TRUE, df=TRUE, layer=1, nl=9, sp=FALSE)

#turns out one set of cells in arrival year raster are missing values
na_index<-which(is.na(env_mean[,10]))
na_latlongs<-ct_latlong[na_index,]

#convert raster to df, omitting NA values
##yields all cells with non-NA values
r_df<-raster::as.data.frame(x=ct_year, xy=TRUE, na.rm=TRUE)

#calculate distance between sampling site and non-NA cells
na_dist<-geodist::geodist(x=na_latlongs, y=r_df[,1:2], measure="geodesic")

#identify minimum distance between sampling site and non-NA cells
##minimum distance is 1.45km, which is negligible
na_dist_min_index<-c()
for (i in 1:nrow(na_latlongs)){
  na_dist_vec<-as.vector(na_dist[i,])
  na_dist_min_index[i]<-which(na_dist_vec==min(na_dist_vec))
}

#obtain latlongs of nearest non-NA points, buffer them, and extract data
na_new_pts<-r_df[na_dist_min_index,1:2]
na_new_buffer<-geobuffer::geobuffer_pts(xy=na_new_pts, dist_m=(buf_rad*1000), step_dg=10)
na_env_mean<-raster::extract(x=env_all, y=na_new_buffer, fun=mean, na.rm=TRUE, df=TRUE, layer=1, nl=9, sp=FALSE)

#input new non-na values into original env_mean df
for (i in 1:length(na_index)){
  env_mean[na_index[i], 2:ncol(na_env_mean)]<-na_env_mean[i,-1]
}

#replace colnames with simpler ones
colnames(env_mean)<-c("pop", "bio1", "bio2", "bio3", "bio4", "bio7", "bio12", "bio15", "elev", "year")

#bio3 and bio4 are multiplied by 100 in the WorldClim database
##largely matters for reporting values (e.g., in table)
env_mean$bio3<-env_mean$bio3/100
env_mean$bio4<-env_mean$bio4/100

#pop ids were replaced by latlong indices (row numbers)
## put them back
env_mean$pop<-ct_sites$pop

#export environmental data for all variables
#write.table(env_mean, "env_mean", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


##############################
#Filter environmental data to remove excessively collinear variables
##############################

#read in locality environmental data from previous section of this script
env_mean<-read.table("env_mean", header=TRUE, sep="\t")

#assess correlations between environmental factors and filter those that are highly multicollinear
#attach region IDs to environmental data
env_mean$region<-md_sites$region_name

#assess overall correlations (all sampling localities together)

#subset data into different regions
env_nw<-env_mean[which(env_mean$region == "NW"),]
env_ne<-env_mean[which(env_mean$region == "NE"),]
env_s<-env_mean[which(env_mean$region == "S"),]

#make plots
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

#make separate correlation plots for:
##1) overall correlations across all localities together
##2) correlations for only NW localities
##3) correlations for only NE localities
##4) correlations for only S localities
par(mfrow=c(2,2))
corrplot(cor(env_mean[,-c(1,11)]), method="color", col=col(200), type="upper", addgrid.col="gray70", addCoef.col="gray10", tl.cex=1, tl.col="black", number.cex=0.9, diag=TRUE)
corrplot(cor(env_nw[,-c(1,11)]), method="color", col=col(200), type="upper", addgrid.col="gray70", addCoef.col="gray10", tl.cex=1, tl.col="black", number.cex=0.9, diag=TRUE)
corrplot(cor(env_ne[,-c(1,11)]), method="color", col=col(200), type="upper", addgrid.col="gray70", addCoef.col="gray10", tl.cex=1, tl.col="black", number.cex=0.9, diag=TRUE)
corrplot(cor(env_s[,-c(1,11)]), method="color", col=col(200), type="upper", addgrid.col="gray70", addCoef.col="gray10", tl.cex=1, tl.col="black", number.cex=0.9, diag=TRUE)
par(mfrow=c(1,1))

#remove environmental factors that correlate with colonization year at |r| > 0.70 across all localities
##retains only bio2, bio4, bio7, bio12, and elev
env_v2<-subset(env_mean, select=c(pop, region, bio2, bio4, bio7, bio12, elev, year))
corrplot(cor(env_v2[,-c(1,2)]), method="color", col=col(200), type="upper", addgrid.col="gray70", addCoef.col="gray10", tl.cex=1, tl.col="black", number.cex=0.9, diag=TRUE)

#now remove environmental factors that correlate with one another at |r| > 0.70 across all localities
#bio7 is strongly correlated with many other variables, so remove it
##end up with three variables: bio2, bio4, bio12, and elev (plus colonization year)
env_v3<-subset(env_v2, select=-c(bio7))
corrplot(cor(env_v3[,-c(1,2)]), method="color", col=col(200), type="upper", addgrid.col="gray70", addCoef.col="gray10", tl.cex=1, tl.col="black", number.cex=0.9, diag=TRUE)

#export filtered set of environmental variables
env_filtered<-data.frame(pop=env_v3$pop, long=md_sites$long, lat=md_sites$lat, env_v3[,-c(1,2)])
#write.table(env_filtered, "env_filtered", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


##############################
#Produce correlation plots for supplementary material
##############################

#calculate correlations
global_corr_prefilter = round(cor(env_mean[,c("bio2", "bio4", "bio12", "elev", "bio1", "bio3", "bio7", "bio15", "year")]),2)
global_corr <- round(cor(env_filtered[,c("bio2", "bio4", "bio12", "elev", "year")]),2)
nw_corr <- round(cor(env_filtered[which(md_sites$region_name=="NW"),c("bio2", "bio4", "bio12", "elev", "year")]),2)
ne_corr <- round(cor(env_filtered[which(md_sites$region_name=="NE"),c("bio2", "bio4", "bio12", "elev", "year")]), 2)
s_corr <- round(cor(env_filtered[which(md_sites$region_name=="S"),c("bio2", "bio4", "bio12", "elev", "year")]), 2)

#rename variables with full names
colnames(global_corr_prefilter) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation",
                                     "Mean annual\ntemperature", "Isothermality", "Temperature\nannual range", "Precipitation\nseasonality", "Colonization\nyear")

rownames(global_corr_prefilter) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation",
                                     "Mean annual\ntemperature", "Isothermality", "Temperature\nannual range", "Precipitation\nseasonality", "Colonization\nyear")

colnames(global_corr) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation", "Colonization\nyear")
rownames(global_corr) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation", "Colonization\nyear")

colnames(nw_corr) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation", "Colonization\nyear")
rownames(nw_corr) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation", "Colonization\nyear")

colnames(ne_corr) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation", "Colonization\nyear")
rownames(ne_corr) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation", "Colonization\nyear")

colnames(s_corr) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation", "Colonization\nyear")
rownames(s_corr) <- c("Mean temperature\ndiurnal range", "Temperature\nseasonality", "Annual\nprecipitation", "Elevation", "Colonization\nyear")


#make correlation plots
plot_marg <- 5
plot_marg_right = 10

corrplot_global_prefilter <- ggcorrplot(global_corr_prefilter, method="square", lab=TRUE, lab_size=5)+
  labs(x="", y="")+
  theme_minimal()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14, angle=90, vjust=0.5), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg_right, b = plot_marg, l = plot_marg))+
  scale_fill_gradient2(low = "royalblue2", mid="white", high = "firebrick3", breaks=c(-1, 0, 1), limit=c(-1, 1))+
  labs(x="", y="", fill="Pearson's\nr")+
  theme(legend.position=c(-0.2, -0.25), legend.key.height=unit(0.33, "cm"))
corrplot_global_prefilter

corrplot_global <- ggcorrplot(global_corr, method="square", lab=TRUE, lab_size=5)+
  labs(x="", y="")+
  theme_minimal()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14, angle=90, vjust=0.5), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg_right, b = plot_marg, l = plot_marg))+
  scale_fill_gradient2(low = "royalblue2", mid="white", high = "firebrick3", breaks=c(-1, 0, 1), limit=c(-1, 1))+
  labs(x="", y="", fill="Pearson's\nr")+
  theme(legend.position=c(-0.2, -0.25), legend.key.height=unit(0.33, "cm"))

corrplot_nw <- ggcorrplot(nw_corr, method="square", lab=TRUE, lab_size=5)+
  labs(x="", y="")+
  theme_minimal()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14, angle=90, vjust=0.5), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg_right, b = plot_marg, l = plot_marg))+
  scale_fill_gradient2(low = "royalblue2", mid="white", high = "firebrick3", breaks=c(-1, 0, 1), limit=c(-1, 1))+
  labs(x="", y="", fill="Pearson's\nr")+
  theme(legend.position=c(-0.2, -0.25), legend.key.height=unit(0.33, "cm"))

corrplot_ne <- ggcorrplot(ne_corr, method="square", lab=TRUE, lab_size=5)+
  labs(x="", y="")+
  theme_minimal()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14, angle=90, vjust=0.5), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg_right, b = plot_marg, l = plot_marg))+
  scale_fill_gradient2(low = "royalblue2", mid="white", high = "firebrick3", breaks=c(-1, 0, 1), limit=c(-1, 1))+
  labs(x="", y="", fill="Pearson's\nr")+
  theme(legend.position=c(-0.2, -0.25), legend.key.height=unit(0.33, "cm"))

corrplot_s <- ggcorrplot(s_corr, method="square", lab=TRUE, lab_size=5)+
  theme_minimal()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14, angle=90, vjust=0.5), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg_right, b = plot_marg, l = plot_marg))+
  scale_fill_gradient2(low = "royalblue2", mid="white", high = "firebrick3", breaks=c(-1, 0, 1), limit=c(-1, 1))+
  labs(x="", y="", fill="Pearson's\nr")+
  theme(legend.position=c(-0.2, -0.25), legend.key.height=unit(0.33, "cm"))

(corrplot_global | corrplot_nw) / (corrplot_ne | corrplot_s) + plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 25))

ggsave("corrplot_filtered.svg", width=16, height=12)


#plot including all starting environmental factors
corrplot_global_prefilter / (corrplot_global + corrplot_nw) / (corrplot_ne + corrplot_s) + plot_annotation(tag_levels="A") + plot_layout(heights=c(2,1,1), guides="collect") & theme(plot.tag = element_text(size = 25), legend.position="right")
ggsave("corrplot_filtered.png", width=15, height=20)
ggsave("corrplot_filtered.svg", width=15, height=20)


16*1.25
12*1.25
