library(diveRsity)

library(ggplot2)
library(ggpubr)
library(gghalves)
library(patchwork)

##############################
#Load cane toad metadata
##############################

#read metadata
md <- read.table("metadata_individuals", sep="\t", header=TRUE)
md_pop <- md[!duplicated(md$pop),][,-1]

#read environmental data (includes colonization year)
env <- read.table("env_filtered", sep="\t", header=TRUE)


##############################
#Run diveRsity to obtain He, Ho, FIS, and rarefied Allelic richness (Ar)
##############################

#specify input genepop
input <- "~/localities_genepop.txt"

#specify number of bootstraps for FIS estimation
fis_bootstraps <- 5000

#run diveRsity
#fis <- diveRsity::basicStats(infile=input, fis_ci=TRUE, fis_boots=fis_bootstraps, rarefaction=TRUE)
#saveRDS(object=fis, file="fis_localities")

#if not re-running above code, read in diveRsity output here:
fis <- readRDS("fis_localities")

#isolate overall (i.e., multilocus) estimates of diversity metrics
div_pars <- as.data.frame(matrix(data=NA, ncol=nrow(fis$main_tab[[1]]), nrow=length(fis$main_tab)))
colnames(div_pars) <- rownames(fis$main_tab[[1]])

for (i in 1:length(fis$main_tab)){
  div_pars[i,] <- fis$main_tab[[i]][,5724]
}
div_pops <- data.frame(md_pop, env[,-1], div_pars)


##############################
#Make plots
##############################

#plotting parameters
pal <- c("NE"="tan1", "NW"="turquoise3", "S"="slateblue3")
plot_marg <- 5

###
#disable scientific notation
options(scipen=999)


#diversity versus colonization year
fis_scatter<-ggscatter(data=div_pops, x='year', y='fis',
                       shape="region_name", color="region_name", size=2.5,
                       palette=pal, add="reg.line", conf.int=TRUE)+
  stat_cor(method="pearson", aes(color = region_name), label.x=c(1930, 1955, 1980), label.y=rep(x=0.1, 3), p.accuracy=0.01, r.accuracy=0.1)+
  scale_y_continuous(limits=c(-0.30, 0.1), breaks=round(seq(from=-0.3, to=0.1, by=0.1), 1))+
  scale_x_continuous(limits=c(1930,2015), breaks=seq(from=1930, to=2015, by=10))+
  xlab("Colonization year")+
  ylab("FIS")+
  theme_bw()+
  labs(colour="Region", fill="Region", shape="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
fis_scatter

ho_scatter<-ggscatter(data=div_pops, x='year', y='obs_het',
                      shape="region_name", color="region_name", size=2.5,
                      palette=pal, add="reg.line", conf.int=TRUE)+
  stat_cor(method="pearson", aes(color = region_name), label.x=c(1930, 1955, 1980), label.y=rep(x=0.4, 3), p.accuracy=0.01, r.accuracy=0.1)+
  scale_y_continuous(limits=c(0.25, 0.40), breaks=seq(from=0.25, to=0.40, by=0.05))+
  scale_x_continuous(limits=c(1930,2015), breaks=seq(from=1930, to=2015, by=10))+
  xlab("Colonization year")+
  ylab("Observed SNP\nheterozygosity")+
  theme_bw()+
  labs(colour="Region", fill="Region", shape="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
ho_scatter

ar_scatter<-ggscatter(data=div_pops, x='year', y='ar',
                      shape="region_name", color="region_name", size=2.5,
                      palette=pal, add="reg.line", conf.int=TRUE)+
  stat_cor(method="pearson", aes(color = region_name), label.x=c(1930, 1955, 1980), label.y=rep(x=1.8, 3), p.accuracy=0.01, r.accuracy=0.1)+
  scale_y_continuous(limits=c(1.6, 1.8), breaks=seq(from=1.6, to=1.8, by=0.05))+
  scale_x_continuous(limits=c(1930,2015), breaks=seq(from=1930, to=2015, by=10))+
  xlab("Colonization year")+
  ylab("Allelic richness")+
  theme_bw()+
  labs(colour="Region", fill="Region", shape="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
ar_scatter


#diversity boxplots
fis_box<-ggplot()+
  geom_half_boxplot(data=div_pops, aes(x=region_name, y=fis, fill=region_name), outlier.shape=4, outlier.size=2, width=0.5)+
  scale_fill_manual(values=pal)+
  geom_half_point(data=div_pops, aes(x=region_name, y=fis, shape=region_name, color=region_name), alpha=0.5, size=2.5, position="dodge2", transformation=position_jitter(width=0.075, height=0))+
  scale_color_manual(values=pal)+
  scale_y_continuous(limits=c(-0.3, 0.1), breaks=round(seq(from=-0.3, to=0.1, by=0.1), digits=2))+
  xlab("Region")+
  ylab("FIS")+
  theme_bw()+
  labs(colour="Region", fill="Region", shape="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
fis_box

ho_box<-ggplot()+
  geom_half_boxplot(data=div_pops, aes(x=region_name, y=obs_het, fill=region_name), outlier.shape=4, outlier.size=2, width=0.5)+
  scale_fill_manual(values=pal)+
  geom_half_point(data=div_pops, aes(x=region_name, y=obs_het, shape=region_name, color=region_name), alpha=0.5, size=2.5, position="dodge2", transformation=position_jitter(width=0.075, height=0))+
  scale_color_manual(values=pal)+
  scale_y_continuous(limits=c(0.25, 0.40), breaks=seq(from=0.25, to=0.40, by=0.05))+
  xlab("Region")+
  ylab("Observed SNP\nheterozygosity")+
  theme_bw()+
  labs(colour="Region", fill="Region", shape="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
ho_box

ar_box<-ggplot()+
  geom_half_boxplot(data=div_pops, aes(x=region_name, y=ar, fill=region_name), outlier.shape=4, outlier.size=2, width=0.5)+
  scale_fill_manual(values=pal)+
  geom_half_point(data=div_pops, aes(x=region_name, y=ar, shape=region_name, color=region_name), alpha=0.5, size=2.5, position="dodge2", transformation=position_jitter(width=0.075, height=0))+
  scale_color_manual(values=pal)+
  scale_y_continuous(limits=c(1.6, 1.8), breaks=seq(from=1.6, to=1.8, by=0.05))+
  xlab("Region")+
  ylab("Allelic richness")+
  theme_bw()+
  labs(colour="Region", fill="Region", shape="Region")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.key.height = unit(0.4, "cm"))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
ar_box

