library(purrr)

mat1<-read.table("draw_matrix01", header=FALSE, sep="\t")[,-60]
mat2<-read.table("draw_matrix02", header=FALSE, sep="\t")[,-60]
mat3<-read.table("draw_matrix03", header=FALSE, sep="\t")[,-60]
mat4<-read.table("draw_matrix04", header=FALSE, sep="\t")[,-60]
mat5<-read.table("draw_matrix05", header=FALSE, sep="\t")[,-60]

mat_list = list(mat1, mat2, mat3, mat4, mat5)
mat_avg = reduce(mat_list, `+`) / length(mat_list)

write.table(mat_avg, "cov_mat_avg", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


