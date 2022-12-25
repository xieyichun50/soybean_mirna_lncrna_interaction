##WGCNA
library(WGCNA)
library(psych)
library(tidyr)
library(dplyr)

####read in expression
exp.mRNA<-read.csv(file = "run/Gsoja/Gsoja_mRNA_filter_gene_count_matrix.csv", 
                   header = TRUE)
DEG<-read.table("run/Gsoja/DEGlist.Gsoja.txt", header = TRUE)
names(exp.mRNA)[1]="gene_id"
exp.mRNA<-exp.mRNA[which(exp.mRNA$gene_id %in% DEG$gene_id),]
row.names(exp.mRNA)=exp.mRNA$gene_id
exp.mRNA<-exp.mRNA[,2:ncol(exp.mRNA)]
exp.mRNA$max = apply(exp.mRNA, 1, max)
exp.mRNA$mean = apply(exp.mRNA, 1, mean)
exp.mRNA<-exp.mRNA[which(exp.mRNA$mean > 10 & exp.mRNA$max >= 50), 1:(ncol(exp.mRNA)-2)]
exp.mRNA<-log2(exp.mRNA+1)
exp.mRNA<-t(exp.mRNA)

setwd("WGCNA_Gsoja/")

####read in trait
trait<-read.delim(file = "trait.txt", header = TRUE)
rownames(trait) = trait$ID
trait <- subset(trait, select = -c(ID))

WGCNA.data <- exp.mRNA

allowWGCNAThreads(4)
save.image("WGCNA.RData")

###########################################################
####Unsigned
###########################################################
####pickSoftThreshold
#similarity.unsigned<-abs(bicor(WGCNA.data))
powers.unsigned = c(c(1:10), seq(from = 1, to=30, by=1))
sft.unsigned = pickSoftThreshold(WGCNA.data, powerVector = powers.unsigned, networkType = "unsigned" , verbose = 5)
png(file = "SoftThreshold.unsigned.png", width = 6, height = 6*0.75, units = "in", res = 300)
SoftThreshold.unsigned<-plot(sft.unsigned$fitIndices[,1], 
                             -sign(sft.unsigned$fitIndices[,3])*sft.unsigned$fitIndices[,2], 
                             xlab='Soft Threshold (power)', 
                             ylab='Scale Free Topology Model Fit, unsigned R^2', 
                             type = 'n', main = paste('Scale independence'))+
  text(sft.unsigned$fitIndices[,1],
       -sign(sft.unsigned$fitIndices[,3])*sft.unsigned$fitIndices[,2],
       labels=powers.unsigned,cex=1,col='red')+
  abline(h=0.80,col='red')
dev.off()

####pick the first power that reach 0.90
#for (i in 1:nrow(sft.unsigned$fitIndices)){
#  if (-sign(sft.unsigned$fitIndices[i,3])*sft.unsigned$fitIndices[i,2] < 0.9){
#    beta.unsigned=sft.unsigned$fitIndices[i,1]+1
#  } else {beta.unsigned = beta.unsigned}
#}
beta.unsigned=18
cat(paste("beta.unsigned=",beta.unsigned, sep = ""))  

save.image("WGCNA.RData")

####Mean Connectivity
png(file = "meanconnectivity.unsigned.png", width = 6, height = 6*0.75, units = "in", res = 300)
meanconnectivity.unsigned<-plot(sft.unsigned$fitIndices[,1], sft.unsigned$fitIndices[,5], 
                                xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
                                type = 'n', main = paste("Mean connectivity, unsigned"))+
  text(sft.unsigned$fitIndices[,1], sft.unsigned$fitIndices[,5], 
       labels=powers.unsigned, cex=1, col="red")
dev.off()
save.image("WGCNA.RData")

####network construction
net.unsigned=blockwiseModules(WGCNA.data, power = beta.unsigned,
                              maxBlockSize = 20000, 
                              minModuleSize = 30,
                              mergeCutHeight = 0.20,
                              corType = "pearson", 
                              TOMType = "unsigned", saveTOMs = TRUE, verbose = 5)

table(net.unsigned$colors)
write.table(table(net.unsigned$colors), 
            file = "net.unsigned.colors.txt", 
            sep = '\t', row.names = FALSE,
            quote = FALSE)
moduleLabels.unsigned = net.unsigned$colors
moduleColors.unsigned = labels2colors(moduleLabels.unsigned)
png(file = "dendrograms.unsigned.png", width = 6, height = 6, units = "in", res = 300)
dendrograms.unsigned<-plotDendroAndColors(net.unsigned$dendrograms[[1]], 
                                          moduleColors.unsigned[net.unsigned$blockGenes[[1]]], 
                                          "Module colors", dendroLabels = FALSE, hang = 0.03, 
                                          addGuide = TRUE, guideHang = 0.05)
dev.off()
save.image("WGCNA.RData")

####correlation among modules
MEs.unsigned = net.unsigned$MEs
MEs_col.unsigned = MEs.unsigned
MEs_col.unsigned = orderMEs(MEs_col.unsigned)
png(file = "Eigengene.unsigned.png", width = 6, height = 6*2, units = "in", res = 300)
Eigengene.unsigned<-plotEigengeneNetworks(MEs_col.unsigned, "Eigengene adjacency heatmap", 
                                          marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), 
                                          plotDendrograms = T, xLabelsAngle = 90)
dev.off()
save.image("WGCNA.RData")

####Trait
moduleTraitCor.unsigned = cor(MEs_col.unsigned, trait, use = "p")
moduleTraitPvalue.unsigned = corPvalueStudent(moduleTraitCor.unsigned, nrow(MEs_col.unsigned))
textMatrix = paste(signif(moduleTraitCor.unsigned, 2), "\n(", signif(moduleTraitPvalue.unsigned, 1), ")", sep = ""); 

png(file = "Module-trait.relationships.png", width = 8, height = 6, units = "in", res = 300)
ModuleTrait.relationheatmap<-
  labeledHeatmap(Matrix = moduleTraitCor.unsigned,
                 xLabels = names(trait),
                 xSymbols = names(trait),
                 yLabels = names(MEs_col.unsigned),
                 ySymbols = names(MEs_col.unsigned),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = TRUE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
dev.off()

save.image("WGCNA.RData")

####TOMplot
load("WGCNA.RData")
load(net.unsigned$TOMFiles[1], verbose=T)
TOM.unsigned <- as.matrix(TOM)
#dissTOM.unsigned = 1-TOM.unsigned
#plotTOM.unsigned = dissTOM.unsigned^beta.unsigned

#png(file = "TOMplot.unsigned.large.png", width = 20, height = 20, units = "in", res = 300)
#TOMplot.unsigned<-TOMplot(plotTOM.unsigned, net.unsigned$dendrograms, 
#                          moduleColors.unsigned,  main = "Network heatmap plot, all genes")
#dev.off()

#geneTree = hclust(as.dist(dissTOM.unsigned), method = "average")
#dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM.unsigned,
#                             deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
#dynamicColors <- labels2colors(dynamicMods)
#table(dynamicColors)
#plot_sim <-dissTOM.unsigned
#diag(plot_sim) <- NA

####Export to cytoscape
probes = colnames(WGCNA.data)
cyt = exportNetworkToCytoscape(TOM.unsigned, edgeFile = paste("all", ".unsigned", ".edges.txt", sep=""), 
                               nodeFile = paste("all", ".unsigned", ".nodes.txt", sep=""), 
                               weighted = TRUE, threshold = 0.5, nodeNames = probes, 
                               nodeAttr = moduleColors.unsigned)
##by color
modulename<-as.data.frame(names(MEs.unsigned))
modulename$`names(MEs.unsigned)`<-gsub("ME", "", modulename$`names(MEs.unsigned)`)
dimnames(TOM.unsigned)<-list(probes, probes)
for (i in 1:nrow(modulename)) {
  modules= modulename[i,1]
  inmodules=is.finite(match(moduleLabels.unsigned, modules))
  modProbes=probes[inmodules]
  modTOM=TOM.unsigned[modProbes,modProbes]
  cyt = exportNetworkToCytoscape(modTOM, edgeFile = paste(modules, ".unsigned", ".edges.txt", sep=""), 
                                 nodeFile = paste(modules, ".unsigned", ".nodes.txt", sep=""), 
                                 weighted = TRUE, threshold = 0.2, nodeNames = modProbes, 
                                 nodeAttr = moduleLabels.unsigned[inmodules])
}

gene.name.fromNode <- gene.name
names(gene.name.fromNode)[1] = "fromNode"
names(gene.name.fromNode)[2] = "fromGene"
gene.name.toNode <- gene.name
names(gene.name.toNode)[1] = "toNode"
names(gene.name.toNode)[2] = "toGene"
edge<-merge(edge, gene.name.fromNode, by = c("fromNode"), all.x = TRUE)
edge<-merge(edge, gene.name.toNode, by = c("toNode"), all.x = TRUE)
edge<-subset(edge, select = c("fromGene", "toGene", "weight"))
write.table(edge, file = "black.unsigned.edges.named.txt")


##Annotation
enrich.input<-cbind(probes, moduleLabels.unsigned)
enrich.input<-as.data.frame(enrich.input)
write.table(enrich.input, file = "gene.color.txt", 
            sep = '\t', row.names = FALSE,
            quote = FALSE)
