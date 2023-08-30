
library(ggplot2)
library(ggpubr)
library(gghalves)


##############################
#Load cane toad metadata
##############################

#read metadata
md <- read.table("metadata_individuals", sep="\t", header=TRUE)
md_pop <- md[!duplicated(md$pop),][,-1]

#read environmental data (includes colonization year)
year <- read.table("env_filtered", sep="\t", header=TRUE)$year

##############################
#Read NeEstimator LD-based Ne outputs
##############################

#read NeEstimator output
ne_maf001 <- read.table("neutral_genepop_LD_Ne.txt", sep="\t", header=FALSE)
colnames(ne_maf001) <- c("pop", "sampsize", "weighted_hmean", "n_ind_alleles", "r2", "exp_r2", "ne", "loci_par", "hici_par", "loci_jack", "hici_jack", "eff_df")

#attach metadata to Ne estimate
pop_ne <- data.frame(pop=md_pop$pop, ne_maf001=ne_maf001$ne)

#create df for plotting
ne_df <- data.frame(pop=md_pop$pop, region_name=md_pop$region_name, ne=ne_maf001$ne, year=year)
ne_df <- ne_df[-which(ne_df$ne<0),]

##############################
#Estimate Ne - colonization year correlations
##############################

#NW region correlation
cor.test(x=ne_df$year[which(ne_df$region_name=="NW")], y=ne_df$ne[which(ne_df$region_name=="NW")])

#S region correlation
cor.test(x=ne_df$year[which(ne_df$region_name=="S")], y=ne_df$ne[which(ne_df$region_name=="S")])

#NE region correlation
cor.test(x=ne_df$year[which(ne_df$region_name=="NE")], y=ne_df$ne[which(ne_df$region_name=="NE")])

#check whether northeast region correlation is still significant after removing outlier
ne_ne <- ne_df[which(ne_df$region_name=="NE"),]
ne_ne_no_outlier <- ne_ne[which(ne_ne$ne<1000),]
cor.test(x=ne_ne_no_outlier$year, y=ne_ne_no_outlier$ne)



##############################
#Create plots
##############################

#plotting parameters
pal <- c("tan1", "turquoise3", "slateblue3")
plot_marg <- 10

#create boxplot
ne_neutral_box<-ggplot()+
  geom_half_boxplot(data=ne_df, aes(x=region_name, y=ne, fill=region_name), outlier.shape=4, outlier.size=2, width=0.5)+
  scale_fill_manual(values=pal)+
  geom_half_point(data=ne_df, aes(x=region_name, y=ne, shape=region_name, color=region_name), alpha=0.5, size=2.5, position="dodge2", transformation=position_jitter(width=0.075, height=0))+
  scale_color_manual(values=pal)+
  scale_y_continuous(limits=c(-20,1200), breaks=seq(from=0, to=1200, by=100), expand=c(0,0))+
  xlab("Region")+
  ylab("Effective population\nsize (Ne)")+
  theme_bw()+
  labs(colour="Region", fill="Region", shape="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  theme(panel.grid.minor=element_blank())
ne_neutral_box

#create scatterplot
ne_neutral_scatter <- ggscatter(data=ne_df, x='year', y='ne',
                                shape="region_name", color="region_name", size=2.5,
                                palette=pal, add="reg.line", conf.int=TRUE)+
  stat_cor(method="pearson", aes(color = region_name), label.x=c(1935, 1960, 1985), label.y=rep(x=1150, 3))+
  scale_y_continuous(limits=c(-20,1200), breaks=seq(from=0, to=1200, by=100), expand=c(0,0), oob = scales::oob_squish)+
  scale_x_continuous(limits=c(1930,2015), breaks=seq(from=1930, to=2015, by=10), expand=c(0,0))+
  xlab("Colonization year")+
  ylab("Effective population\nsize (Ne)")+
  theme_bw()+
  labs(colour="Region", fill="Region", shape="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.key.height = unit(1, "cm"))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  theme(panel.grid.minor=element_blank())
ne_neutral_scatter

