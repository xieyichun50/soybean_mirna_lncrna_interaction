##WGCNA
library(WGCNA)
library(psych)
library(tidyr)
library(dplyr)

####read in expression
exp.mRNA<-read.csv(file = "run/Gmax_mRNA_filter_gene_count_matrix.csv", 
                   header = TRUE)
DEG<-read.table("run/DEGlist.Gmax.txt", header = TRUE)
exp.mRNA<-exp.mRNA[which(exp.mRNA$gene_id %in% DEG$gene_id),]
row.names(exp.mRNA)=exp.mRNA$gene_id
exp.mRNA<-exp.mRNA[,2:ncol(exp.mRNA)]
exp.mRNA$max = apply(exp.mRNA, 1, max)
exp.mRNA$mean = apply(exp.mRNA, 1, mean)
exp.mRNA<-exp.mRNA[which(exp.mRNA$mean > 10 & exp.mRNA$max >= 50), 1:(ncol(exp.mRNA)-2)]
exp.mRNA<-t(exp.mRNA)

setwd("WGCNA_Gmax/")

####read in trait
trait<-read.delim(file = "trait.txt", header = TRUE)
rownames(trait) = trait$ID
trait <- subset(trait, select = -c(ID))

WGCNA.data <- exp.mRNA

allowWGCNAThreads(4)
save.image("WGCNA.RData")

###########################################################
####Signed
###########################################################
####pickSoftThreshold
#similarity.signed<-abs(bicor(WGCNA.data)+1)/2
powers.signed = c(c(1:10), seq(from = 1, to=30, by=1))
sft.signed = pickSoftThreshold(WGCNA.data, powerVector = powers.signed, networkType = "signed" , verbose = 5)
png(file = "SoftThreshold.signed.png", width = 6, height = 6*0.75, units = "in", res = 300)
SoftThreshold.signed<-plot(sft.signed$fitIndices[,1], 
                             -sign(sft.signed$fitIndices[,3])*sft.signed$fitIndices[,2], 
                             xlab='Soft Threshold (power)', 
                             ylab='Scale Free Topology Model Fit, signed R^2', 
                             type = 'n', main = paste('Scale independence'))+
  text(sft.signed$fitIndices[,1],
       -sign(sft.signed$fitIndices[,3])*sft.signed$fitIndices[,2],
       labels=powers.signed,cex=1,col='red')+
  abline(h=0.90,col='red')
dev.off()

####pick the first power that reach 0.90
for (i in 1:nrow(sft.signed$fitIndices)){
  if (-sign(sft.signed$fitIndices[i,3])*sft.signed$fitIndices[i,2] < 0.9){
    beta.signed=sft.signed$fitIndices[i,1]+1
  } else {beta.signed = beta.signed}
}

cat(paste("beta.signed=",beta.signed, sep = ""))  

save.image("WGCNA.RData")

####Mean Connectivity
png(file = "meanconnectivity.signed.png", width = 6, height = 6*0.75, units = "in", res = 300)
meanconnectivity.signed<-plot(sft.signed$fitIndices[,1], sft.signed$fitIndices[,5], 
                                xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
                                type = 'n', main = paste("Mean connectivity, signed"))+
  text(sft.signed$fitIndices[,1], sft.signed$fitIndices[,5], labels=powers.signed, cex=1, col="red")
dev.off()
save.image("WGCNA.RData")

####network construction
net.unsigned=blockwiseModules(WGCNA.data, power = beta.unsigned, maxBlockSize = 20000, corType = "bicor", 
                              TOMType = "unsigned", saveTOMs = TRUE, verbose = 5)

table(net.unsigned$colors)
moduleLabels.unsigned = net.unsigned$colors
moduleColors.unsigned = labels2colors(moduleLabels.unsigned)
png(file = "dendrograms.unsigned.png", width = 6, height = 6, units = "in", res = 300)
dendrograms.unsigned<-plotDendroAndColors(net.unsigned$dendrograms[[1]], 
                                          moduleColors.unsigned[net.unsigned$blockGenes[[1]]], 
                                          "Module colors", dendroLabels = FALSE, hang = 0.03, 
                                          addGuide = TRUE, guideHang = 0.05)
dev.off()
save.image("~/huilab/drosophila/expression/mRNA_miRNA_filter/WGCNA.RData")

######merge small difference group
net.unsigned.merge=blockwiseModules(WGCNA.data, power = beta.unsigned, maxBlockSize = 20000, corType = "bicor", 
                                    TOMType = "unsigned", saveTOMs = TRUE, mergeCutHeight = 0.2, verbose = 3)

table(net.unsigned.merge$colors)
moduleLabels.unsigned.merge = net.unsigned.merge$colors
moduleColors.unsigned.merge = labels2colors(moduleLabels.unsigned.merge)
png(file = "dendrograms.unsigned.merge.png", width = 6, height = 6, units = "in", res = 300)
dendrograms.unsigned.merge<-plotDendroAndColors(net.unsigned.merge$dendrograms[[1]], 
                                                moduleColors.unsigned.merge[net.unsigned.merge$blockGenes[[1]]], 
                                                "Module colors", dendroLabels = FALSE, hang = 0.03, 
                                                addGuide = TRUE, guideHang = 0.05)
dev.off()
save.image("~/huilab/drosophila/expression/mRNA_miRNA_filter/WGCNA.RData")

####TOMplot
TOM.unsigned = TOMsimilarityFromExpr(WGCNA.data, power=beta.unsigned, corType="bicor", networkType="unsigned")
TOM.unsigned <- as.matrix(TOM.unsigned)
dissTOM.unsigned = 1-TOM.unsigned
plotTOM.unsigned = dissTOM.unsigned^beta.unsigned


png(file = "TOMplot.unsigned.large.png", width = 20, height = 20, units = "in", res = 300)
TOMplot.unsigned<-TOMplot(plotTOM.unsigned, net.unsigned$dendrograms, 
                          moduleColors.unsigned,  main = "Network heatmap plot, all genes")
dev.off()
save.image("~/huilab/drosophila/expression/mRNA_miRNA_filter/WGCNA.RData")

######Select part of the genes
minModuleSize <- 30  #模块基因数目
geneTree = hclust(as.dist(dissTOM.unsigned), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM.unsigned,
                             deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
plot_sim <-dissTOM.unsigned
diag(plot_sim) <- NA
#png(file = "TOMplot.unsigned.large.png", width = 20, height = 20, units = "in", res = 300)
#TOMplot.unsigned.selected<-TOMplot(plot_sim, WGCNA.data, dynamicColors, 
#                                   main = 'Network heatmap plot, selected genes')
#dev.off()

save.image("~/huilab/drosophila/expression/mRNA_miRNA_filter/WGCNA.RData")
####correlation among modules
MEs.unsigned = net.unsigned$MEs
MEs_col.unsigned = MEs.unsigned
MEs_col.unsigned = orderMEs(MEs_col.unsigned)
png(file = "Eigengene.unsigned.png", width = 6, height = 6*2, units = "in", res = 300)
Eigengene.unsigned<-plotEigengeneNetworks(MEs_col.unsigned, "Eigengene adjacency heatmap", 
                                          marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), 
                                          plotDendrograms = T, xLabelsAngle = 90)
dev.off()
save.image("~/huilab/drosophila/expression/mRNA_miRNA_filter/WGCNA.RData")

####correlation among modules.merge
MEs.unsigned.merge = net.unsigned.merge$MEs
MEs_col.unsigned.merge = MEs.unsigned.merge
MEs_col.unsigned.merge = orderMEs(MEs_col.unsigned.merge)

png(file = "Eigengene.unsigned.merge.png", width = 6, height = 6*2, units = "in", res = 300)
Eigengene.unsigned.merge<-plotEigengeneNetworks(MEs_col.unsigned.merge, "Eigengene adjacency heatmap", 
                                                marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), 
                                                plotDendrograms = T, xLabelsAngle = 90)
dev.off()
save.image("~/huilab/drosophila/expression/mRNA_miRNA_filter/WGCNA.RData")

####Trait

moduleTraitCor.unsigned.merge = cor(MEs_col.unsigned.merge, trait, use = "p")
moduleTraitPvalue.unsigned.merge = corPvalueStudent(moduleTraitCor.unsigned.merge, nrow(MEs_col.unsigned.merge))
textMatrix = paste(signif(moduleTraitCor.unsigned.merge, 2), "\n(", signif(moduleTraitPvalue.unsigned.merge, 1), ")", sep = ""); 

png(file = "Module-trait.relationships.png", width = 8, height = 6, units = "in", res = 300)
ModuleTrait.relationheatmap<-
  labeledHeatmap(Matrix = moduleTraitCor.unsigned.merge,
                 xLabels = names(trait),
                 xSymbols = names(trait),
                 yLabels = names(MEs_col.unsigned.merge),
                 ySymbols = names(MEs_col.unsigned.merge),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = TRUE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
dev.off()

save.image("~/huilab/drosophila/expression/mRNA_miRNA_filter/WGCNA.RData")

####Export to cytoscape

probes = colnames(WGCNA.data)
modulename<-as.data.frame(names(MEs.unsigned.merge))
modulename$`names(MEs.unsigned.merge)`<-gsub("ME", "", modulename$`names(MEs.unsigned.merge)`)
for (i in 1:nrow(modulename)) {
  modules= modulename[i,1]
  inmodules=is.finite(match(moduleColors.unsigned.merge, modules))
  modProbes=probes[inmodules]
  modTOM=TOM.unsigned[inmodules,inmodules]
  dimnames(modTOM) <- list(modProbes,modProbes)
  cyt = exportNetworkToCytoscape(modTOM, edgeFile = paste(modules, ".unsigned", ".edges.txt", sep=""), 
                                 nodeFile = paste(modules, ".unsigned", ".nodes.txt", sep=""), 
                                 weighted = TRUE, threshold = 0.3, nodeNames = modProbes, 
                                 nodeAttr = moduleColors.unsigned[inmodules])
  
}


modules= "blue"
inmodules=is.finite(match(moduleColors.unsigned.merge, modules))
modProbes=probes[inmodules]
modTOM=TOM.unsigned[inmodules,inmodules]
dimnames(modTOM) <- list(modProbes,modProbes)
cyt = exportNetworkToCytoscape(modTOM, edgeFile = paste(modules, ".unsigned", ".edges.txt", sep=""), 
                               nodeFile = paste(modules, ".unsigned", ".nodes.txt", sep=""), 
                               weighted = TRUE, threshold = 0.3, nodeNames = modProbes, 
                               nodeAttr = moduleColors.unsigned[inmodules])

edge<-read.delim(file = "~/huilab/drosophila/expression/black.unsigned.edges.txt")
gene.name.fromNode <- gene.name
names(gene.name.fromNode)[1] = "fromNode"
names(gene.name.fromNode)[2] = "fromGene"
gene.name.toNode <- gene.name
names(gene.name.toNode)[1] = "toNode"
names(gene.name.toNode)[2] = "toGene"
edge<-merge(edge, gene.name.fromNode, by = c("fromNode"), all.x = TRUE)
edge<-merge(edge, gene.name.toNode, by = c("toNode"), all.x = TRUE)
edge<-subset(edge, select = c("fromGene", "toGene", "weight"))
write.table(edge, file = "~/huilab/drosophila/expression/black.unsigned.edges.named.txt")