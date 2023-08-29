
library(gsubfn)
library(stringr)

##############################
#create BED files for SNPs with significant GEAs
##############################

#read in all SNPs to make reference gene list (for Gene Ontolog enrichment)
ref_raw <- read.table("ct_lmiss.lmiss", sep="\t", header=TRUE)[,1:2]
ref_bed <- data.frame(chrom=ref_raw[,1], chromStart=ref_raw[,2]-1, chromEnd=ref_raw[,2])
#write.table(ref_bed, "ref_bedfile", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

###
#read in all SNPs with significant GEAs identified by either Bayenv2 or GWR
be_gwr_snp_all <- read.table("be_gwr_snp_all")
colnames(be_gwr_snp_all) <- c("contig", "pos")

#calculate 0-index start position for bedtools
be_gwr_snp_all_start <- be_gwr_snp_all$pos - 1

be_gwr_all_bed <- data.frame(chrom=be_gwr_snp_all$contig, chromStart=be_gwr_snp_all_start, chromEnd=be_gwr_snp_all$pos)
#write.table(be_gwr_all_bed, "be_gwr_all_bedfile", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


##############################
#Read in bedtools closest -d output
##############################

#gene_annotation.gff3 (found in this subdirectory) is one of the inputs for bedtools

#read in output based on gene.gff3
genes_raw <- read.table("bedtools_closest_genes", header=FALSE, sep="\t")

#reorder by contig name
genes_contigs <- matrix(unlist(stringr::str_split(genes_raw$V1, "_")), ncol=4, byrow=TRUE)
genes_contigs[,1] <- sub('...', '', genes_contigs[,1])
genes_ordered <- genes_raw[order(as.numeric(genes_contigs[,1])),]
genes_ordered[genes_ordered == "."] <- NA
genes_ordered[genes_ordered == "-1"] <- NA

#now extract just the ID of possible ortholog
gene_ids <- c()
for (i in 1:nrow(genes_ordered)){
  #get matching string
  stringmatch <- strapplyc(genes_ordered$V7[i], "Similar to (.*;)", simplify = TRUE)[[1]]
  
  #if no matching string, output NA
  if (length(stringmatch)==0){
    gene_ids[i] <- NA
  }else{
    gene_ids[i] <- strsplit(stringmatch, ";")[[1]][1]
  }
}

#organize SNP and gene genomic information
genes_ordered <- cbind(genes_ordered[,-8], gene_ids, genes_ordered[,8])
colnames(genes_ordered) <- c("snp_contig", "snp_start", "snp_end", "feature_contig", "feature_start", "feature_end", "feature_info", "ortholog_id", "distance")

###
#add environmental associations to the SNPs/candidate genes

#read in full locus list
##uses vcftools locus missingness output
ct_loci <- read.table("ct_lmiss.lmiss", header=TRUE)[,1:2]
ct_loci <- paste(ct_loci[,1], ct_loci[,2], sep="_")

#read in per-environmental factor significant SNPs for Bayenv2
##these are SNP indices that lack locus information
be_sig <- readRDS("be_gea_sig_kf_r95")
be_sig_bio2 <- be_sig$Bio2
be_sig_bio4 <- be_sig$Bio4
be_sig_bio12 <- be_sig$Bio12
be_sig_elev <- be_sig$Elev

#isolate locus names for Bayenv2 significant SNPs
be_sig_bio2_loci <- ct_loci[be_sig_bio2]
be_sig_bio4_loci <- ct_loci[be_sig_bio4]
be_sig_bio12_loci <- ct_loci[be_sig_bio12]
be_sig_elev_loci <- ct_loci[be_sig_elev]

#read in per-environmental factor significant SNPs for GWR
##these are SNP indices that lack locus information
gwr_sig <- readRDS("gwr_output_03032023")
gwr_par_df <- as.data.frame(do.call(rbind, gwr_sig))
gwr_par_df_sig <- gwr_par_df[which(gwr_par_df$bf>quantile(gwr_par_df$bf, 0.95)),]
gwr_sig_env <- as.data.frame(matrix(unlist(str_split(gwr_par_df_sig$formula, " ")), ncol=9, byrow=TRUE))
gwr_sig_env$snp <- gwr_par_df_sig$snp
gwr_sig_env <- gwr_sig_env[,c(10,9)]

gwr_sig_bio2 <- gwr_sig_env$snp[which(gwr_sig_env$V9=="bio2")]
gwr_sig_bio4 <- gwr_sig_env$snp[which(gwr_sig_env$V9=="bio4")]
gwr_sig_bio12 <- gwr_sig_env$snp[which(gwr_sig_env$V9=="bio12")]
gwr_sig_elev <- gwr_sig_env$snp[which(gwr_sig_env$V9=="elev")]

#isolate locus names for GWR significant SNPs
gwr_sig_bio2_loci <- ct_loci[gwr_sig_bio2]
gwr_sig_bio4_loci <- ct_loci[gwr_sig_bio4]
gwr_sig_bio12_loci <- ct_loci[gwr_sig_bio12]
gwr_sig_elev_loci <- ct_loci[gwr_sig_elev]

###
#construct final table
genes_ordered$gwr <- NA
genes_ordered$bayenv <- NA
genes_ordered$gwr_bio2 <- NA
genes_ordered$bayenv_bio2 <- NA
genes_ordered$gwr_bio4 <- NA
genes_ordered$bayenv_bio4 <- NA
genes_ordered$gwr_bio12 <- NA
genes_ordered$bayenv_bio12 <- NA
genes_ordered$gwr_elev <- NA
genes_ordered$bayenv_elev <- NA

gene_loci <- paste(genes_ordered[,1], genes_ordered[,3], sep="_")

for (i in 1:nrow(genes_ordered)){
  #GWR detection
  genes_ordered$gwr_bio2[i] <- !is.na(match(gene_loci[i], gwr_sig_bio2_loci))
  genes_ordered$gwr_bio4[i] <- !is.na(match(gene_loci[i], gwr_sig_bio4_loci))
  genes_ordered$gwr_bio12[i] <- !is.na(match(gene_loci[i], gwr_sig_bio12_loci))
  genes_ordered$gwr_elev[i] <- !is.na(match(gene_loci[i], gwr_sig_elev_loci))
  
  #Bayenv2 detection
  genes_ordered$bayenv_bio2[i] <- !is.na(match(gene_loci[i], be_sig_bio2_loci))
  genes_ordered$bayenv_bio4[i] <- !is.na(match(gene_loci[i], be_sig_bio4_loci))
  genes_ordered$bayenv_bio12[i] <- !is.na(match(gene_loci[i], be_sig_bio12_loci))
  genes_ordered$bayenv_elev[i] <- !is.na(match(gene_loci[i], be_sig_elev_loci))
  
  #detected at all by GWR or Bayenv2
  gwr_any <- sum(c(genes_ordered$gwr_bio2[i], genes_ordered$gwr_bio4[i], genes_ordered$gwr_bio12[i], genes_ordered$gwr_elev[i]), na.rm=TRUE)
  bayenv_any <- sum(c(genes_ordered$bayenv_bio2[i], genes_ordered$bayenv_bio4[i], genes_ordered$bayenv_bio12[i], genes_ordered$bayenv_elev[i]), na.rm=TRUE)
  
  if (gwr_any > 0){
    genes_ordered$gwr[i] <- TRUE
  }else{
    genes_ordered$gwr[i] <- FALSE
  }
  
  if (bayenv_any > 0){
    genes_ordered$bayenv[i] <- TRUE
  }else{
    genes_ordered$bayenv[i] <- FALSE
  }
  
}

#remove FALSE values to make table easier to read
genes_ordered[genes_ordered==FALSE] <- NA

#isolate candidate genes for only SNPs detected by both Bayenv and GWR
genes_ordered_intersect <- genes_ordered[which(genes_ordered$gwr==TRUE & genes_ordered$bayenv==TRUE),]

#output CSV files
##will be post-processed to remove NA values for environmental conditions (makes it easier to read)
##final versions of these tables are in supp_table_s3.xlsx
write.csv(genes_ordered, "candidate_genes_all.csv", row.names=FALSE)
write.csv(genes_ordered_intersect, "candidate_genes_intersect.csv", row.names=FALSE)
