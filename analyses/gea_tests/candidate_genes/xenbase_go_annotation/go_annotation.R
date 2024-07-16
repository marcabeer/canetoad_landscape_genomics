
library(openxlsx)

#read in Xenbase gene annotations
xenb_genes = read.table("C:/Users/marc_/Documents/wsu/research/cane_toad/peer_review_05-2024/go_annotation/GeneGoTerms_chd.txt", header=TRUE, sep="\t")

#capitalize xenbase gene symbols
xenb_genes$gene_symbol_caps = toupper(xenb_genes$gene_symbol)

#read in candidate SNPs/genes
cand_sheet1 = openxlsx::read.xlsx("table_s4-s5.xlsx", sheet=1)
cand_18 = openxlsx::read.xlsx("table_s4-s5.xlsx", sheet=2)
cand_all = openxlsx::read.xlsx("table_s4-s5.xlsx", sheet=3)

#attach GO terms to candidate SNPs/genes
#cand_18_go_merge = merge.data.frame(x=cand_18, y=xenb_genes, by.x="gene", by.y="gene_symbol_caps", all.x=TRUE, sort=FALSE)
#cand_all_go_merge = merge.data.frame(x=cand_all, y=xenb_genes, by.x="gene", by.y="gene_symbol_caps", all.x=TRUE)

#18 overlapping candidate genes
cand_18_go = data.frame(Xenbase_page=rep(NA, nrow(cand_18)), Xenbase_gene_id=rep(NA, nrow(cand_18)), GO_ids=rep(NA, nrow(cand_18)))

for (i in 1:nrow(cand_18)){
  xenb_gene_info = xenb_genes[which(xenb_genes$gene_symbol_caps == cand_18$gene[i]), c("Xenbase", "gene_ID", "GO_Ids")]
  if (nrow(xenb_gene_info)==0){
    cand_18_go[i,] = cand_18_go[i,]
  }else{
    cand_18_go[i,] = xenb_gene_info
  }
}

#all candidate genes
cand_all_go = data.frame(Xenbase_page=rep(NA, nrow(cand_all)), Xenbase_gene_id=rep(NA, nrow(cand_all)), GO_ids=rep(NA, nrow(cand_all)))

for (i in 1:nrow(cand_all)){
  xenb_gene_info = xenb_genes[which(xenb_genes$gene_symbol_caps == cand_all$gene[i]), c("Xenbase", "gene_ID", "GO_Ids")]
  if (nrow(xenb_gene_info)==0){
    cand_all_go[i,] = cand_all_go[i,]
  }else{
    cand_all_go[i,] = xenb_gene_info
  }
}

#replace NAs with character NAs
cand_18_go[is.na(cand_18_go)] = "NA"
cand_all_go[is.na(cand_all_go)] = "NA"

#combine go annotations with candidate gene info
cand_18_go_df = cbind(cand_18[,1:3], cand_18_go, cand_18[,4:ncol(cand_18)])
cand_all_go_df = cbind(cand_all[,1:3], cand_all_go, cand_all[,4:ncol(cand_all)])

#write output
write.xlsx(x=cand_18_go_df, file="table_s4_wGO.xlsx", append=TRUE, overwrite=FALSE)
write.xlsx(x=cand_all_go_df, file="table_s5_wGO.xlsx", append=TRUE, overwrite=FALSE)

