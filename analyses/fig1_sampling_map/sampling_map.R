
#basic data manipulation and plotting
library(plyr)
library(ggplot2)
library(gghalves)
library(patchwork)

#geographic data plotting
library(raster)
library("rnaturalearth")
library("rnaturalearthhires")
library("rnaturalearthdata")
library(ggspatial)
#basic data manipulation and plotting
library(plyr)
library(ggplot2)
library(gghalves)
library(patchwork)
library(ggnewscale)


##############################
#Load metadata
##############################

#read metadata
md <- read.table("metadata_individuals", sep="\t", header=TRUE)
md_pops <- md[!duplicated(md$pop),][,-1]

#calculate sample size per locality
sample_size <- c()
for (i in 1:nrow(md_pops)){
  sample_size[i] <- length(which(md$pop==md_pops$pop[i]))
}
md_pops$sample_size <- sample_size

#read in environmental data
env <- read.table("env_filtered", header=TRUE, sep="\t")

#merged md_pops env data
df <- merge.data.frame(md_pops, env, by.x=c("pop", "long", "lat"), by.y=c("pop", "long", "lat"))

##############################
#read in cane toad introduction locations
##############################

#load introduction sites
intro_sites <- read.csv("introduction_sites.csv")
colnames(intro_sites)[1] <- "site"

#combine Giru and Ayr locations into one since they are too close to map separately
intro_sites <- as.data.frame(rbind(intro_sites, c("Giru & Ayr", median(intro_sites$latitude[which(intro_sites$site == "Giru" | intro_sites$site == "Ayr")]),
                                                  median(intro_sites$longitude[which(intro_sites$site == "Giru" | intro_sites$site == "Ayr")]), NA)))
intro_sites$longitude <- as.numeric(intro_sites$longitude)
intro_sites$latitude <- as.numeric(intro_sites$latitude)

#combine Bundaberg and Isis locations into one since they are too close to map separately
intro_sites <- as.data.frame(rbind(intro_sites, c("Bundaberg\n& Isis", median(intro_sites$latitude[which(intro_sites$site == "Bundaberg" | intro_sites$site == "Isis")]),
                                                  median(intro_sites$longitude[which(intro_sites$site == "Bundaberg" | intro_sites$site == "Isis")]), NA)))
intro_sites$longitude <- as.numeric(intro_sites$longitude)
intro_sites$latitude <- as.numeric(intro_sites$latitude)

intro_sites <- intro_sites[-c(which(intro_sites$site == "Giru" | intro_sites$site == "Ayr")),]
intro_sites <- intro_sites[-c(which(intro_sites$site == "Bundaberg" | intro_sites$site == "Isis")),]

##############################
#Load colonization year raster and habitat suitability raster
##colonization year raster comes from output of ~/analyses/colonization_year
##############################

#load colonization year raster
colyear <- raster("colonization_year.tif")

#load habitat suitability raster
##need to gunzip suitability_median.asc first
enm <- raster("suitability_median.asc")
projection(enm) <- "+proj=longlat +datum=WGS84 +no_defs"

#resample habitat suitability raster
##so that it stacks with colonization year raster
enm_resamp <- resample(x=enm, y=colyear, method="bilinear")

#stack the two rasters
colyear_enm <- stack(colyear, enm_resamp)
plot(colyear_enm)

#aggregate to lower resolution for plotting
colyear_enm_aggfact3 <- aggregate(colyear_enm, fact=3)

#convert to df
colyear_enm_df <- as.data.frame(colyear_enm_aggfact3, xy=TRUE)
colnames(colyear_enm_df) <- c("x", "y", "year", "suitability")

#remove NA values
colyear_enm_df_nona <- colyear_enm_df[-which(is.na(colyear_enm_df$year)),]
colyear_enm_df_nona$suitability[which(is.na(colyear_enm_df_nona$suitability))] <- 0

#remove values where median habitat suitability is < 0.20
##avoids misrepresenting spatial extend of cane toad
colyear_enm_df_nona_suit025 <- colyear_enm_df_nona[-which(colyear_enm_df_nona$suitability < 0.2),]

#plotting invasion/colonization years as-is ultimately makes the figure less clear
##visualize broad decadal progression of cane toad instead
isoline_df <- colyear_enm_df_nona_suit025
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


##############################
#create sampling map figure
##############################

#define color palette for sampling localities
pal <- c("tan1", "turquoise3", "slateblue3")

#specify locations of region labels (i.e., NW, NE, and S)
region_labels <- data.frame(long=c(134.5, 144, 149.5), lat=c(-13.5, -18, -24), region=c("NW", "NE", "S"))

#read in global geographic data
world <- ne_countries(scale = "large", returnclass = "sf")

#define color palette for colonization year surface
colorpal <- viridis::mako(n=8, begin=0.3, end=0.9) 
color_ramp <- colorRampPalette(colors=colorpal)(8)

#create preliminary sampling map
sample_map <- ggplot(data = world) +
  geom_sf(fill="gray75", colour="gray70") +
  coord_sf(xlim = c(127, 155), ylim = c(-30, -10))+
  ggnewscale::new_scale_fill()+
  geom_tile(data=isoline_df, mapping=aes(x=x, y=y, fill=as.factor(year)))+
  scale_fill_manual(values=colorpal, breaks=seq(from=1940, to=2010, by=10), labels=seq(from=1940, to=2010, by=10), guide=guide_legend(reverse=TRUE))+
  labs(fill="Colonization\nyear")+
  ggnewscale::new_scale_fill()+
  geom_star(data = df, aes(x = long, y = lat, size=sample_size, fill=region_name, colour=region_name), starshape="hexagon")+
  scale_colour_manual(values=c("black", "black", "white"))+
  scale_size_continuous(range = c(2.75, 4.25))+
  scale_fill_manual(values=pal)+
  geom_point(data=intro_sites[which(is.na(intro_sites$year)),], mapping=aes(x=longitude, y=latitude), shape=24, colour="black", fill="white", size=3)+
  geom_point(data=intro_sites[which(!is.na(intro_sites$year)),], mapping=aes(x=longitude, y=latitude), shape=24, colour="black", fill="firebrick1", size=3.5)+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=12))+
  labs(x="Longitude", y="Latitude", size="Sample size", fill="Region")+
  geom_text(region_labels, mapping=aes(x=long, y=lat, label=region), size=10, colour="white")+##
  geom_text(intro_sites, mapping=aes(x=longitude+1.5, y=latitude+c(0,0,-0.25,0,0,0,0), label=site), size=5)+
  theme(legend.position=c(0.90, 0.65))

sample_map

#add scale bar and north arrow
sample_map2 <- sample_map +
  ggspatial::annotation_north_arrow(location = "br", which_north = "true", style=north_arrow_fancy_orienteering())+
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white"),
    text_cex=1
  )
sample_map2

#add inset showing Australia as a whole
rectangle <- data.frame(x1=127, x2=155, y1=-30, y2=-10)
rectangle$x2 - rectangle$x1
rectangle$y2 - rectangle$y1
scale_factor <- 20/28

inset_map <- ggplot(data = world) +
  geom_sf(fill="gray90", colour="gray70") +
  coord_sf(xlim = c(115, 155), ylim = c(-38, -10))+
  labs(x="Longitude", y="Latitude")+
  geom_point(data = df, aes(x = long, y = lat), colour="black", size=1)+
  scale_fill_manual(values=pal)+
  scale_size_continuous(range = c(3, 6))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())+
  theme(plot.margin = margin(0,0,0,0))+
  theme(axis.ticks.length = unit(0, "mm"))+
  scale_x_continuous(limits=c(115,155), breaks=seq(from=115, to=155, by=5))+
  geom_rect(rectangle, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), colour="black", alpha=0, size=1, linetype=2)
inset_map

#combine sampling map with inset map
sample_map_final <- sample_map2 + inset_element(p=inset_map, left=0.0001, bottom=0.015, right = 0.5, top=0.5)
sample_map_final

#output
sample_map_width <- 12
#ggsave(filename="sample_map.svg", width=sample_map_width, height=scale_factor*sample_map_width)

sample_map_final_labelled <- sample_map_final + plot_annotation(tag_levels = list(c("A", ""))) & theme(plot.tag = element_text(size = 25))
sample_map_final_labelled

#save output
##note that this figure was processed further in graphic design software to move around labels
ggsave(filename="sample_map.svg", width=sample_map_width, height=scale_factor*sample_map_width)


##############################
#add additional panels showing environmental variation at sampling localities
##############################

#plotting parameters
plot_marg = 20

#temperature diurnal range
p_bio2 <- ggplot()+
  geom_half_violin(df, mapping=aes(y=bio2, x=region_name, fill=region_name, colour=region_name), alpha=0.25, size=1, side="r")+
  geom_half_boxplot(df, mapping=aes(y=bio2, x=region_name, fill=region_name), colour="black", width=0.2, size=0.5, side="r")+
  geom_half_point(df, mapping=aes(y=bio2, x=region_name, fill=region_name), colour="black", shape=21, size=3, side="l")+
  scale_fill_manual(values=pal)+
  scale_colour_manual(values=pal)+
  scale_y_continuous(limits = c(round_any(min(df$bio2), 1, f=floor), round_any(max(df$bio2), 1, f=ceiling)),
                     breaks=seq(from = round_any(min(df$bio2), 1, f=floor), to=round_any(max(df$bio2), 1, f=ceiling), by=1)
  )+
  theme_bw()+
  labs(y="Mean diurnal temperature range (degC)", x="Region", colour="Region", fill="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  coord_flip()

#temperature seasonality
p_bio4 <- ggplot()+
  geom_half_violin(df, mapping=aes(y=bio4, x=region_name, fill=region_name, colour=region_name), alpha=0.25, size=1, side="r")+
  geom_half_boxplot(df, mapping=aes(y=bio4, x=region_name, fill=region_name), colour="black", width=0.2, size=0.5, side="r")+
  geom_half_point(df, mapping=aes(y=bio4, x=region_name, fill=region_name), colour="black", shape=21, size=3, side="l")+
  scale_fill_manual(values=pal)+
  scale_colour_manual(values=pal)+
  scale_y_continuous(limits = c(round_any(min(df$bio4), 1, f=floor), round_any(max(df$bio4), 1, f=ceiling)),
                     breaks=seq(from = round_any(min(df$bio4), 1, f=floor), to=round_any(max(df$bio4), 1, f=ceiling), by=1)
  )+
  theme_bw()+
  labs(y="Temperature seasonality (degC)", x="Region", colour="Region", fill="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  coord_flip()

#precipitation
p_bio12 <- ggplot()+
  geom_half_violin(df, mapping=aes(y=bio12, x=region_name, fill=region_name, colour=region_name), alpha=0.25, size=1, side="r")+
  geom_half_boxplot(df, mapping=aes(y=bio12, x=region_name, fill=region_name), colour="black", width=0.2, size=0.5, side="r")+
  geom_half_point(df, mapping=aes(y=bio12, x=region_name, fill=region_name), colour="black", shape=21, size=3, side="l")+
  scale_fill_manual(values=pal)+
  scale_colour_manual(values=pal)+
  scale_y_continuous(limits = c(round_any(min(df$bio12), 100, f=floor), round_any(max(df$bio12), 10, f=ceiling)),
                     breaks=seq(from = round_any(min(df$bio12), 100, f=floor), to=round_any(max(df$bio12), 10, f=ceiling), by=200)
  )+
  theme_bw()+
  labs(y="Annual precipitation (mm)", x="Region", colour="Region", fill="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  coord_flip()

#elevation
p_elev <- ggplot()+
  geom_half_violin(df, mapping=aes(y=elev, x=region_name, fill=region_name, colour=region_name), alpha=0.25, size=1, side="r")+
  geom_half_boxplot(df, mapping=aes(y=elev, x=region_name, fill=region_name), colour="black", width=0.2, size=0.5, side="r")+
  geom_half_point(df, mapping=aes(y=elev, x=region_name, fill=region_name), colour="black", shape=21, size=3, side="l")+
  scale_fill_manual(values=pal)+
  scale_colour_manual(values=pal)+
  scale_y_continuous(limits = c(round_any(min(df$elev), 10, f=floor), round_any(max(df$elev), 10, f=ceiling)),
                     breaks=seq(from = round_any(min(df$elev), 10, f=floor), to=round_any(max(df$elev), 10, f=ceiling), by=100)
  )+
  theme_bw()+
  labs(y="Elevation (m)", x="Region", colour="Region", fill="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.position = "none")+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  coord_flip()

#combined plot

layout <- "
AAAAAA
AAAAAA
AAAAAA
AAAAAA
AAAAAA
AAAAAA
BBBCCC
BBBCCC
DDDEEE
DDDEEE
"

sample_env_var <- sample_map_final + p_bio2 + p_bio4 + p_bio12 + p_elev + plot_layout(design=layout) +  plot_annotation(tag_levels = list(c("A", "", "B", "C", "D", "E"))) & theme(plot.tag = element_text(size = 30))
sample_env_var

#save output
##note that this figure was processed further in graphic design software to move around labels
ggsave(filename="sample_map_env_var.svg", width=sample_map_width, height=2*scale_factor*sample_map_width)

