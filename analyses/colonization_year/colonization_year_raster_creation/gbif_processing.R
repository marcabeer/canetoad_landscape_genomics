library(stringr)
library(raster)
library(lubridate)
library(rangeBuilder)


##############################
#Process cane toad occurrence records from GBIF
##############################

#load and manipulate occurrence records
gbif_raw <- readLines("ct_gbif.csv")

#test<-str_replace_all(gbif_raw[[2]], pattern="\t\t", replacement="\tNA\t")

#replace successive \t with \tNA\t in order to insert missing NAs
##then split strings at \t
gbif_fillna <- list()
gbif_fillna2 <- list()
for (i in 1:length(gbif_raw)){
  gbif_fillna[[i]] <- str_replace_all(gbif_raw[[i]], pattern="\t\t", replacement="\tNA\t")
  gbif_fillna2[[i]] <- strsplit(gbif_fillna[[i]], split="\t")
}

#keep only first 40 fields (some observations are missing later fields)
gbif_fillna_trunc <- list()
for (i in 1:length(gbif_fillna)){
  gbif_fillna_trunc[[i]] <- gbif_fillna2[[i]][[1]][c(1:40)]
  
}

#turn list into df
gbif <- t(as.data.frame(gbif_fillna_trunc))
colnames(gbif) <- c(gbif[1,])
rownames(gbif) <- 1:nrow(gbif)
gbif <- gbif[-1,]

gbif<-data.frame(gbif)

#turn occurrence observation dates in to decimal format
gbif_ymd <- ymd(paste(gbif$year, gbif$month, gbif$day, sep="-"))
gbif_yeardec <- decimal_date(gbif_ymd)

#isolate year and coordinates
gbif_all <- data.frame(long=as.numeric(gbif$decimalLongitude), lat=as.numeric(gbif$decimalLatitude), year=gbif_yeardec)
na_year_index <- which(is.na(gbif_all$year))
gbif_all$year[na_year_index] <- gbif$year[na_year_index]

#load in sampling metadata for genetic dataset
ct_sites<-read.table("ct_sites.csv", sep=",", header=TRUE)
ct_sites<-data.frame(long=ct_sites$long, lat=ct_sites$lat, year=2010)

#combine gbif observations and those from the genetic dataset
gbif_all <- data.frame(rbind(gbif_all, ct_sites))
gbif_all <- gbif_all[which(gbif_all$lat > -34.5),]
gbif_all <- gbif_all[which(gbif_all$long > 125),]

#remove unlikely occurrence locations (e.g., central Australia in the desert)
gbif_badlat<-gbif_all$lat < -22
gbif_badlong<-gbif_all$long < 140
gbif_badlat_badlong<-gbif_badlat+gbif_badlong
gbif_all<-gbif_all[-which(gbif_badlat_badlong == 2),]

#save filtered set of occurrence records
#write.csv(gbif_all, "occurrence_filtered.csv", row.names=FALSE)


##############################
#Generate colonization year surface using the rangeBuilder package
##############################

###
#load in raster of Australia (uses part of the elevation raster from the WorldClim database)
##this will be the reference raster that is overwritten by colonization year data
elev_aus <- raster("elev_aus.tif")
spplot(elev_aus)

#replace non-NA values with 1
##wipes original elevation data, creating a blank raster of the Australian continent's terrestrial landmass
r_aus <- elev_aus
r_aus[!is.na(r_aus)] <- 1
spplot(r_aus)

#plot to check geographic distribution of occurrence records
plot(r_aus)
points(gbif_all[which(gbif_all$year<=2011),1:2], pch=16, cex=0.75, col="darkred")

###
#generate alpha shape of occurence records for each year
##repeat process multiple times with different starting alphas
##will ultimately average all of the runs together

#start at 1936 because it is the first year with >3 unique coordinates (necessary for creating alpha hull)
##stop at 2011 because it is the latest year of sampling of genetic data
years <- seq(from=1936, to=2011, by=0.5)

#set up parameters for loop
alpha_hulls <- list()
alpha_vec <- seq(from=2, to=4, by=1)
base_reps <- list()

#use buffer of 10km in the generation of alpha shapes because it represents the 90% quantile of coordinate uncertainty in meters
p2011_coord_uncertainty<-as.numeric(gbif$coordinateUncertaintyInMeters[which(gbif_all$year<=2011)])
quantile(p2011_coord_uncertainty, 0.90, na.rm=TRUE)

#loop through starting alpha values and years
for (k in 1:length(alpha_vec)){
  
  alpha_hulls<-list()
  
  for (i in 1:length(years)){
    
    #keep only occurrence records through a given year
    gbif_subset <- gbif_all[which(gbif_all$year<=years[i]),]
    
    #remove any duplicate records within a given year
    gbif_subset <- gbif_subset[!duplicated(gbif_subset[,c(1:2)]),]
    
    #generate alpha shape [hull] of occurrence records in a given year
    ahull <- getDynamicAlphaHull(gbif_subset[,1:2], fraction = 0.975, partCount = 5, buff = 10000, 
                               initialAlpha = alpha_vec[k], coordHeaders = c('long', 'lat'), 
                               clipToCoast = 'terrestrial', proj = "+proj=longlat +datum=WGS84", 
                               alphaIncrement = 0.5, verbose = FALSE)
    
    #rasterize the alpha shape
    ahull_rast <- rasterize(ahull[[1]], r_aus, mask=TRUE)
    ahull_rast[!is.na(ahull_rast)] <- years[i]
    alpha_hulls[[i]] <- ahull_rast
  }
  names(alpha_hulls)<-years
  
  base<-alpha_hulls$`2011`
  
  index<-rev(seq(from=1, to=length(years), by=1))
  
  for (i in 1:length(index)){
    base[!is.na(alpha_hulls[[index[i]]])]<-years[index[i]]
  }
  
  base_reps[[k]]<-base
  print(paste("completed", k, "replicates of", length(alpha_vec), sep=" "))
  
}

#plot the results
base_stack<-stack( base_reps[[1]], base_reps[[2]], base_reps[[3]])
plot(base_stack)
spplot(base_stack)

#average all surfaces together
base_mean <- mean(base_stack, na.rm=TRUE)

#save results
writeRaster(base_reps[[1]], "ct_a2_v7", format="GTiff")
writeRaster(base_reps[[2]], "ct_a3_v7", format="GTiff")
writeRaster(base_reps[[3]], "ct_a4_v7", format="GTiff")
writeRaster(base_mean, "ct_mean_v7", format="GTiff")

#post-process using SAGA
##used Gaussian blurring at 20-cell distance to smooth the raster


