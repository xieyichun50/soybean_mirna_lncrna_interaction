####Load package
library(edgeR)

####Read in TMM expression matrix
rawcount<-read.delim("gene_count_matrix.csv", header = TRUE)
row.names(rawcount)<-rawcount$gene_id
TMM<-rawcount[,2:ncol(rawcount)]

data=as.matrix(TMM) 

####Read in sample group
trait<-read.delim("run/Gsoja/samples_n_reads_decribed.txt", header = FALSE)
group <- factor(trait[,2])
y<-DGEList(counts=data,group=group)
y <- calcNormFactors(y, method = "TMM")
TMM.norm<-cpm(y)

write.table(TMM.norm, file = "gene_count_TMMmatrix.txt", sep = "\t")

TMM.norm.log2<-log2(TMM.norm+1)
write.table(TMM.norm.log2, file = "gene_count_log2TMMmatrix.txt", sep = "\t")

##mean
TMM.norm<-read.delim("run/Gsoja/gene_count_matrix.TMM.xls", header = TRUE)
trait<-read.delim("run/Gsoja/samples_n_reads_decribed.txt", header = FALSE)
treatment<-as.data.frame(unique(trait[,1]))
names(treatment)[1]<-"treatment"
TMM.mean<-matrix(NA, nrow = nrow(TMM.norm), ncol = nrow(treatment),
                 dimnames = list(TMM.norm$X,treatment$treatment))
TMM.mean<-as.data.frame(TMM.mean)
for (i in 1:nrow(treatment)) {
  TMM.mean[i] <-(TMM.norm[(i*3)-1]+TMM.norm[i*3]+TMM.norm[(i*3)+1])/3
  
}

write.csv(TMM.mean, file = "Gsoja_mRNA_filter_gene_count_matrix.csv", quote = FALSE)
