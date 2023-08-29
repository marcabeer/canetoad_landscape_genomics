#read in raw environmental data
env <- read.table("env_filtered", header=TRUE, sep="\t")

#remove population label, latlongs, and year
env_nopop <- env[,-c(1,2,3,8)]

#standardize environmental data
env_scaled<-scale(env_nopop, center=TRUE, scale=TRUE)

#transpose to meet matrix orientation needed by bayenv
env_scaled_t<-t(env_scaled)
write.table(env_scaled_t, "environfile.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
