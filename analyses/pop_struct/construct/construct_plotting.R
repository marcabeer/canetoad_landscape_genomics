wd_base <- getwd()

#########
#load packages

library(conStruct)

library(dplyr)
library(tidyr)
library(plotrix)

library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)


##########
###
#read in metadata
md <- read.table(paste(getwd(), "/env_filtered", sep=""), sep="\t", header=TRUE)

##########
###
#read in cross-validation results

setwd(paste(getwd(), "/construct_crossvalidation_best_runs", sep=""))

#read in cross-validation R object
my.xvals <- readRDS("my.xvals")

#isolate results for spatial (sp) and non-spatial (nsp) models
sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)


#find best spatial model replicate for each value of K to compare layer contributions
best_rep <- data.frame(K=1:10, rep=NA)
for (i in 1:nrow(sp.results)){
  best_rep[i,2] <- which.max(sp.results[i,])
}


###
#compare layer contributions
layer.contributions <- matrix(NA,nrow=10,ncol=10)


load("cv_sp_rep1K1_conStruct.results.Robj")
load("cv_sp_rep1K1_data.block.Robj")

tmp <- conStruct.results[[1]]$MAP$admix.proportions
layer.contributions[,1] <- c(calculate.layer.contribution(conStruct.results[[1]],data.block),rep(0,9))


for(i in 2:10){
  # load the conStruct.results.Robj and data.block.Robj
  #   files saved at the end of a conStruct run
  load(paste("cv_sp_rep", best_rep[i,2], "K", best_rep[i,1], "_conStruct.results.Robj",  sep=""))
  load(paste("cv_sp_rep", best_rep[i,2], "K", best_rep[i,1], "_data.block.Robj",  sep=""))
  
  
  # match layers up across runs to keep plotting colors consistent
  #   for the same layers in different runs
  tmp.order <- match.layers.x.runs(tmp,conStruct.results[[1]]$MAP$admix.proportions)  
  
  # calculate layer contributions
  layer.contributions[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                                            data.block=data.block,
                                                            layer.order=tmp.order),
                               rep(0,10-i))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order]
}


##############################
##########
###make cross-validation plots

###
#set some plotting parameters
plot_marg <- 20


###
#cross-validation model accuracy plots
cv_accuracy_long <- data.frame(k=rep(rep(1:10, 20),2), model=rep(c("spatial", "non-spatial"), each=20*10), accuracy=c(c(sp.results), c(nsp.results)))

cv_accuracy_summ <- cv_accuracy_long %>%
                    group_by(k, model) %>%
                    summarise(acc_mean = mean(accuracy), acc_95ci = 1.96*std.error(accuracy))

clrblind_pal <- c("#000000", "#D81B60")

#predictive accuracy for all values of K
p_acc <- ggplot()+
  geom_errorbar(data=cv_accuracy_summ, mapping=aes(ymin=acc_mean - acc_95ci, ymax=acc_mean + acc_95ci, x=as.factor(k), width=0.2))+
  geom_point(data=cv_accuracy_summ, mapping=aes(y=acc_mean, x=as.factor(k), fill=model), shape=21, colour="black", size=3)+
  scale_colour_manual(values=clrblind_pal)+
  scale_fill_manual(values=clrblind_pal)+
  theme_bw()+
  labs(x="K (No. geogenetic layers)", y="Predictive accuracy", fill="Model", colour="Model")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.position = c(0.8,0.2))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#predictive accuracy for K >= 2
p_acc_zoom <- ggplot()+
  geom_errorbar(data=cv_accuracy_summ[-which(cv_accuracy_summ$k==1),], mapping=aes(ymin=acc_mean - acc_95ci, ymax=acc_mean + acc_95ci, x=as.factor(k), width=0.2))+
  geom_point(data=cv_accuracy_summ[-which(cv_accuracy_summ$k==1),], mapping=aes(y=acc_mean, x=as.factor(k), fill=model), shape=21, colour="black", size=3.5)+
  scale_colour_manual(values=clrblind_pal)+
  scale_fill_manual(values=clrblind_pal)+
  theme_bw()+
  labs(x="K (No. geogenetic layers)", y="Predictive accuracy", fill="Model", colour="Model")+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.position = c(0.8,0.2))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))

#combine plots
(p_acc | p_acc_zoom) + plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 25))


##########
###
#cross-validation layer contribution plot
pal = c("deepskyblue4", "goldenrod1", "orangered", "darkviolet", "antiquewhite", "gray60", "chocolate1", "cadetblue1", "firebrick4", "deeppink")

lyrs <- data.frame(k=rep(1:10,10), layer=rep(1:10, each=10), contribution=c(t(layer.contributions)))

#stackplot
p_lyr_contribs <- ggplot(lyrs, aes(fill=as.factor(layer), y=contribution, x=k)) + 
  geom_bar(position="stack", stat="identity", width=1, colour="black", linewidth=0.75)+
  scale_fill_manual(values=pal)+
  labs(x = "K (No. geogenetic layers)", y="Layer contribution", fill="Layer")+
  scale_y_continuous(limits=c(0,1), breaks=seq(from=0, to=1, by=0.2), expand=c(0,0))+
  scale_x_discrete(limits=factor(1:10), expand=c(0,0))+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))
p_lyr_contribs


#minimum non-zero layer contribution
lyrs_min <- data.frame(k=1:10, min_contrib=NA)
for (i in 1:10){
  tmp_contrib <- lyrs$contribution[which(lyrs$k==i)]
  tmp_contrib[tmp_contrib==0] <- NA
  lyrs_min$min_contrib[i] <- min(tmp_contrib, na.rm=TRUE)
}

p_lyrs_min <- ggplot()+
  geom_hline(yintercept=0, size=1, colour="gray70")+
  geom_line(lyrs_min, mapping=aes(x=k, y=min_contrib), colour="black", linewidth=1)+
  geom_point(lyrs_min, mapping=aes(x=k, y=min_contrib), colour="black", size=3)+
  geom_hline(yintercept=0.05, linetype=2, size=1, colour="red")+ ###line indicates 5% layer contribution
  scale_y_continuous(limits=c(0,1), breaks=seq(from=0, to=1, by=0.2))+
  scale_x_continuous(limits=c(1,10), breaks=seq(from=1, to=10, by=1))+
  labs(x="K (No. geogenetic layers)", y="Minimum layer contribution")+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))


#combine plots
cv_layout <- "
AAABBB
AAABBB
CCCDDD
CCCDDD
"

p_acc + p_acc_zoom + p_lyr_contribs + p_lyrs_min + plot_layout(design=cv_layout) + plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 25))

ggsave(paste(wd_base, "/plots/crossval_plots.svg", sep=""), width=12, height=10)
  



##############################
##########
###
#locality ancestry proportions from final models
setwd(paste(wd_base, "/construct_final_models", sep=""))

###
#load final models
load("final_k3_conStruct.results.Robj")
load("final_k3_data.block.Robj")

###
# assign model results to new variable names
k3sp_cr <- conStruct.results
k3sp_db <- data.block

# load output files from a run with 
#   the spatial model and K=3
load("final_k4_conStruct.results.Robj")
load("final_k4_data.block.Robj")

# assign to new variable names
k4sp_cr <- conStruct.results
k4sp_db <- data.block

#extract locality admixture proportions 
admix.props_k3 <- k3sp_cr$chain_3$MAP$admix.proportions
admix.props_k4 <- k4sp_cr$chain_3$MAP$admix.proportions

#rearrange K=3 layers to match those of K=2
admix.props_k4 <- data.frame(admix.props_k4[,c(1,4,2,3)])

#add column names
colnames(admix.props_k3) <- c("Layer_1", "Layer_2", "Layer_3")
colnames(admix.props_k4) <- c("Layer_1", "Layer_2", "Layer_3", "Layer_4")

#attach metadata to locality admixture proportions
k3_admix_md <- cbind(admix.props_k3, md)
k4_admix_md <- cbind(admix.props_k4, md)

#reorder localities by longitude
k3_admix_ordered <- arrange(k3_admix_md, long)
k4_admix_ordered <- arrange(k4_admix_md, long)

#add numbers for ordering localities in plot
k3_admix_ordered$number <- seq(1:nrow(k3_admix_ordered))
k4_admix_ordered$number <- seq(1:nrow(k4_admix_ordered))

###
#save locality ancestry proportions in tables
#write.table(k3_admix_ordered[,c("SRR", "Longitude", "Latitude", "Layer_1", "Layer_2")], paste(wd_base, "/tables/k3_anc_prop", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#write.table(k4_admix_ordered[,c("SRR", "Longitude", "Latitude", "Layer_1", "Layer_2", "Layer_3")], paste(wd_base, "/tables/k4_anc_prop", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


###
#convert to long format for ggplot
#K = 3
keycol <- "cs_layer"
valuecol <- "anc_prop"
gathercols <- c("Layer_1", "Layer_2", "Layer_3")
k3_long <- gather_(k3_admix_ordered, keycol, valuecol, gathercols)

#k = 4
keycol <- "cs_layer"
valuecol <- "anc_prop"
gathercols <- c("Layer_1", "Layer_2", "Layer_3", "Layer_4")
k4_long <- gather_(k4_admix_ordered, keycol, valuecol, gathercols)

#new palettes
pal_k3 <- c("turquoise3", "tan1", "slateblue3")
pal_k4 <- c("turquoise3", "tan1", "slateblue3", "orangered")

###
#make stackplots / ancestry proportion plots
p_k3_stackplot <- ggplot(k3_long, aes(fill=as.factor(cs_layer), y=anc_prop, x=number)) + 
  geom_bar(position="stack", stat="identity", width=1, colour="black", size=0.33)+
  scale_fill_manual(values=pal_k3)+
  labs(x = "Locality\n(ordered by longitude)", y="Ancestry proportion", fill="Geogenetic\nlayer")+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))+
  scale_x_discrete(limits=factor(1:max(k3_long$number)), expand=c(0,0))+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.position = "none")
p_k3_stackplot

p_k4_stackplot <- ggplot(k4_long, aes(fill=as.factor(cs_layer), y=anc_prop, x=number)) + 
  geom_bar(position="stack", stat="identity", width=1, colour="black", size=0.33)+
  scale_fill_manual(values=pal_k4)+
  labs(x = "Locality\n(ordered by longitude)", y="Ancestry proportion", fill="Geogenetic\nlayer")+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))+
  scale_x_discrete(limits=factor(1:max(k3_long$number)), expand=c(0,0))+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.position = "none")
p_k4_stackplot

###
#make maps showing locality ancestry proportions
#maps
world <- ne_countries(scale = "medium", returnclass = "sf")

###
#create locality pies for map

plot_marg <- 10

#k=2
k3_pie.list <- k3_long %>% 
  tidyr::nest(cs_layer, anc_prop) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data,
                               function(d) ggplotGrob(ggplot(d, 
                                                             aes(x = 1, y = anc_prop, fill = cs_layer)) +
                                                        geom_col(color = "black",
                                                                 show.legend = FALSE) +
                                                        coord_polar(theta = "y") +
                                                        scale_fill_manual(values=pal_k3)+
                                                        theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(radius = 0.25) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = long - radius, xmax = long + radius,
                                          ymin = lat - radius, ymax = lat + radius)))
#k=3
k4_pie.list <- k4_long %>% 
  tidyr::nest(cs_layer, anc_prop) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data,
                               function(d) ggplotGrob(ggplot(d, 
                                                             aes(x = 1, y = anc_prop, fill = cs_layer)) +
                                                        geom_col(color = "black",
                                                                 show.legend = FALSE) +
                                                        coord_polar(theta = "y") +
                                                        scale_fill_manual(values=pal_k4)+
                                                        theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(radius = 0.25) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = long - radius, xmax = long + radius,
                                          ymin = lat - radius, ymax = lat + radius)))




#The dimensions of this figure have been carefully tailored
##plotting maps with other plots is challenging
##because the dimensions of maps are constrained to accurately reflect latitude and longitude

p_k3_map <- ggplot(data = world) +
  geom_sf(fill="gray80", colour="gray70") +
  coord_sf(xlim = c(min(md$long), max(md$long)+0.5), ylim = c(min(md$lat)-1, max(md$lat)))+
  labs(x="Longitude", y="Latitude", fill="Slope")+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(legend.position=c(0.2, 0.25))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  k3_pie.list$subgrob
  
p_k4_map <- ggplot(data = world) +
  geom_sf(fill="gray80", colour="gray70") +
  coord_sf(xlim = c(min(md$long), max(md$long)+0.5), ylim = c(min(md$lat)-1, max(md$lat)))+
  labs(x="Longitude", y="Latitude", fill="Slope")+
  theme_bw()+
  theme(axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = plot_marg, r = plot_marg, b = plot_marg, l = plot_marg))+
  k4_pie.list$subgrob+
  theme(legend.position=c(0.5,0.5))

layout <- "
AAAAABBBBB
AAAAABBBBB
AAAAABBBBB
AAAAABBBBB
CCCCCDDDDD
CCCCCDDDDD
"
#combine final model plots
p_k3_map + p_k4_map + p_k3_stackplot + p_k4_stackplot +
  plot_layout(design=layout) +
  plot_annotation(tag_levels=list(c("A", "B", "", ""))) & theme(plot.tag = element_text(size = 25)) 


p_k3_map_final <- p_k3_map + inset_element(p=p_k3_stackplot, left=0.015, bottom=0.015, right = 0.5, top=0.5)
p_k4_map_final <- p_k4_map + inset_element(p=p_k4_stackplot, left=0.015, bottom=0.015, right = 0.5, top=0.5)

p_k3_map_final / p_k4_map_final +
  plot_annotation(tag_levels=list(c("A", "", "B", ""))) & theme(plot.tag = element_text(size = 25)) 

ggsave(paste(wd_base, "/conStruct_final_models.svg", sep=""), width=10, height=14)


