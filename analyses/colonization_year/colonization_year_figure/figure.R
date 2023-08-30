#basic data manipulation and plotting
library(plyr)
library(ggplot2)
library(gghalves)
library(patchwork)
library(ggnewscale)
library(viridis)

#geographic data plotting
library("rnaturalearth")
library("rnaturalearthhires")
library("rnaturalearthdata")

library(stringr)
library(raster)
library(rangeBuilder)


##############################
#Load filtered set of occurrence records
##used in creation of colonization year raster surface
##############################

#load in filtered occurrence records
occurrence <- read.csv("occurrence_filtered.csv")

#load sampling localities
localities <- read.csv("ct_sites.csv")


##############################
#Load colonization year raster
##############################

#load colonization year raster
colyear <- raster("saga_gauss_20cell_v7.tif")

#aggregate to lower resolution for plotting
colyear <- aggregate(colyear, fact=3)

#store as df
colyear_df <- as.data.frame(colyear, xy=TRUE)
colnames(colyear_df) <- c("x", "y", "year")
colyear_df_nona <- colyear_df[-which(is.na(colyear_df$year)),]


##############################
#Create plot/map
##############################

#specify color palette
colorpal <- viridis::mako(n=20, begin=0.3, end=0.9)

#specify color breaks
colorbreaks <- seq(from=1930, to=2010, by=10)

#load world geographic data
world <- ne_countries(scale = "large", returnclass = "sf")

#make map
ggplot(data = world) +
  geom_sf(fill="gray75", colour="gray70") +
  coord_sf(xlim = c(127, 155), ylim = c(-30, -10))+
  ggnewscale::new_scale_fill()+
  geom_tile(data=colyear_df_nona, mapping=aes(x=x, y=y, fill=year))+
  scale_fill_gradientn(colours=colorpal, breaks=colorbreaks)+
  labs(fill="Colonization\nyear")+
  ggnewscale::new_scale_fill()+
  geom_point(data = occurrence[-c(which(occurrence$year > 2011)),], aes(x = long, y = lat), shape=19, colour="black", size=1, alpha=0.3)+
  geom_point(data = localities, aes(x = long, y = lat), shape=21, colour="black", fill="white", size=3)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=12))+
  labs(x="Longitude", y="Latitude")+
  theme(legend.position=c(0.90, 0.65), legend.key.height=unit(1.25, "cm"))


sample_map_width <- 12
scale_factor <- 20/28
ggsave(filename="invasion_map.svg", width=sample_map_width, height=scale_factor*sample_map_width)



##############################
#Create plot/map showing decadal progression of cane toads
##############################

#generate isolines
##painfully hard-coded; I know, I know
isoline_df <- colyear_df_nona
isoline_df$year[which(isoline_df$year<1940)] <- 0
isoline_df$year[which(isoline_df$year>=1940 & isoline_df$year<1950)] <- 1
isoline_df$year[which(isoline_df$year>=1950 & isoline_df$year<1960)] <- 2
isoline_df$year[which(isoline_df$year>=1960 & isoline_df$year<1970)] <- 3
isoline_df$year[which(isoline_df$year>=1970 & isoline_df$year<1980)] <- 4
isoline_df$year[which(isoline_df$year>=1980 & isoline_df$year<1990)] <- 5
isoline_df$year[which(isoline_df$year>=1990 & isoline_df$year<2000)] <- 6
isoline_df$year[which(isoline_df$year>=2000)] <- 7

isoline_df[isoline_df==0] <- 1940
isoline_df[isoline_df==1] <- 1950
isoline_df[isoline_df==2] <- 1960
isoline_df[isoline_df==3] <- 1970
isoline_df[isoline_df==4] <- 1980
isoline_df[isoline_df==5] <- 1990
isoline_df[isoline_df==6] <- 2000
isoline_df[isoline_df==7] <- 2010

#make map
colorpal2 <- viridis::mako(n=8, begin=0.3, end=0.9) 

ggplot(data = world) +
  geom_sf(fill="gray80", colour="gray70") +
  coord_sf(xlim = c(127, 155), ylim = c(-30, -10))+
  ggnewscale::new_scale_fill()+
  geom_tile(data=isoline_df, mapping=aes(x=x, y=y, fill=as.factor(year)))+
  scale_fill_manual(values=colorpal2, breaks=seq(from=1940, to=2010, by=10), labels=seq(from=1940, to=2010, by=10), guide=guide_legend(reverse=TRUE))+
  labs(fill="Colonization\nyear")+
  ggnewscale::new_scale_fill()+
  geom_point(data = localities, aes(x = long, y = lat), shape=21, colour="black", fill="white", size=3)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=12))+
  labs(x="Longitude", y="Latitude")+
  theme(legend.position=c(0.90, 0.65), legend.key.height=unit(1.25, "cm"))

ggsave(filename="invasion_map_isolines.svg", width=sample_map_width, height=scale_factor*sample_map_width)

