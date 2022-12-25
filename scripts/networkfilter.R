setwd("D:/beetles/soybean_miRNA_lncRNA_interaction")
setwd("WGCNA_Gsoja/")
genemodule<-read.delim("turquoise.unsigned.edges.txt", header = TRUE)
summary(genemodule)
genemodule.sub<-subset(genemodule, weight >= 0.635)
a<-as.data.frame(unique(genemodule.sub$fromNode))
names(a)[1]="gene"
b<-as.data.frame(unique(genemodule.sub$toNode))
names(b)[1]="gene"
ab<-rbind(a,b)
nrow(unique(ab))
write.table(genemodule.sub,
            file = "turquoise.unsigned.edges.filter635.txt", 
            sep = '\t', row.names = FALSE,
            quote = FALSE)
