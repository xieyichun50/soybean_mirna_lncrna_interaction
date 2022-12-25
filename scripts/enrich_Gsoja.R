library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)

#Whole genome Annotation file
eggnog<-read.delim("eggnog-online/Gsoja.eggnog.emapper.annotations", 
                   header = TRUE, skip = 4)
#ID<-read.delim("eggnog/Gene_Protein_ID", header = FALSE)
#names(ID)[1]="X.query"
#names(ID)[2]="Genes"
#eggnog<-merge(eggnog, ID, by = "X.query", all.x = TRUE)
eggnog$X.query<-gsub("Glysoja.","",eggnog$X.query)
eggnog$X.query<-gsub(".\\d$", "",eggnog$X.query)
names(eggnog)[1]="Genes"
##KEGG from eggnog
{
  #read in kegg2name
  kegg2name <- read.delim("D:/3enrichment/kegg2name.txt", 
                          sep = "\t", colClasses = "character")
  
  #eggnog<-separate(eggnog, V1, c("Genes", "Transcript"), sep = "-", remove = FALSE)
  pathways<-eggnog[,c("Genes","KEGG_Pathway")]
  pathways<-separate(pathways, KEGG_Pathway, c("ko","map"), sep = ",map", remove = FALSE)
  pathways<-pathways[,c(1,3)]
  names(pathways)[1]="Genes"
  names(pathways)[2]="KEGG"
  pathways<-subset(pathways, KEGG != "-" & KEGG != "" & is.na(KEGG)==FALSE,
                   select = c("Genes", "KEGG"))
  pathways<-unique(pathways)
  #Format GenesKEGGpair
  GenesKEGGpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
  GenesKEGGpair.1v1<-as.data.frame(GenesKEGGpair.1v1)
  names(GenesKEGGpair.1v1)[1]="Genes"
  names(GenesKEGGpair.1v1)[2]="KEGG"
  
  for (i in 1:nrow(pathways)) {
    subtable<-pathways[i,]
    rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$KEGG[1], ',')[[1]]))
    pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
    pairtable<-as.data.frame(pairtable)
    pairtable$Genes<-rownames(pairtable)
    rownames(pairtable)<-1:nrow(pairtable)
    pairtable<-as.data.frame(pairtable)
    pairtable.new<-pairtable %>% gather(KEGG, pair, c(1:ncol(pairtable)-1))
    pairtable.new<-pairtable.new[,c(1:2)]
    GenesKEGGpair.1v1<-rbind(GenesKEGGpair.1v1, pairtable.new)
  }
  GenesKEGGpair.1v1<-subset(GenesKEGGpair.1v1, 
                            is.na(GenesKEGGpair.1v1$Genes)==FALSE, 
                            select = c("KEGG", "Genes"))
  GenesKEGGpair.1v1<-unique(GenesKEGGpair.1v1)
  
  write.table(GenesKEGGpair.1v1, 
              file = "eggnog-online/Gsoja.KEGG.1v1.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  rm(pairtable, pairtable.new, rcnames, subtable, pathways)
}
kegg2name <- read.delim("D:/3enrichment/kegg2name.txt", 
                        sep = "\t", colClasses = "character")
GenesKEGGpair.1v1<-read.delim(file = "eggnog-online/Gsoja.KEGG.1v1.txt",
                              header = TRUE)
##KOG from eggnog
{
  #read in kog2name
  kog2name<-read.delim("D:/3enrichment/kog2name.txt", 
                       sep = "\t", colClasses = "character")
  
  pathways<-eggnog[,c("Genes","COG_category")]
  names(pathways)[1]="Genes"
  names(pathways)[2]="KOG"
  pathways<-subset(pathways, KOG != "" & KOG != "-" & is.na(KOG)==FALSE)
  pathways<-unique(pathways)
  
  #Format GenesKOGpair
  GenesKOGpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
  GenesKOGpair.1v1<-as.data.frame(GenesKOGpair.1v1)
  names(GenesKOGpair.1v1)[1]="Genes"
  names(GenesKOGpair.1v1)[2]="KOG"
  
  for (i in 1:nrow(pathways)) {
    subtable<-pathways[i,]
    rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$KOG[1], '')[[1]]))
    pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
    pairtable<-as.data.frame(pairtable)
    pairtable$Genes<-rownames(pairtable)
    rownames(pairtable)<-1:nrow(pairtable)
    pairtable<-as.data.frame(pairtable)
    pairtable.new<-pairtable %>% gather(KOG, pair, c(1:ncol(pairtable)-1))
    pairtable.new<-pairtable.new[,c(1:2)]
    GenesKOGpair.1v1<-rbind(GenesKOGpair.1v1, pairtable.new)
  }
  GenesKOGpair.1v1<-subset(GenesKOGpair.1v1, 
                           is.na(GenesKOGpair.1v1$Genes)==FALSE, 
                           select = c("KOG", "Genes"))
  GenesKOGpair.1v1<-unique(GenesKOGpair.1v1)
  
  write.table(GenesKOGpair.1v1, 
              file = "eggnog-online/Gsoja.KOG.1v1.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  rm(pairtable, pairtable.new, rcnames, subtable, pathways)
  
}
kog2name<-read.delim("D:/3enrichment/kog2name.txt", 
                     sep = "\t", colClasses = "character")
GenesKOGpair.1v1<-read.delim(file = "eggnog-online/Gsoja.KOG.1v1.txt",
                             header = TRUE)
##GO from eggnog
{
  #read in go2name
  go2name<-read.delim("D:/3enrichment/go2name.txt", 
                      sep = "\t", colClasses = "character",
                      header = FALSE)
  names(go2name)[1]="goClass"
  names(go2name)[2]="goName"
  names(go2name)[3]="ONTOLOGY"
  
  pathways<-eggnog[,c("Genes", "GOs")]
  names(pathways)[1]="Genes"
  names(pathways)[2]="GO"
  pathways<-subset(pathways, GO != "" & GO != "-" & is.na(GO)==FALSE)
  pathways<-unique(pathways)
  
  #Format GenesGOGpair
  GenesGOpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
  GenesGOpair.1v1<-as.data.frame(GenesGOpair.1v1)
  names(GenesGOpair.1v1)[1]="Genes"
  names(GenesGOpair.1v1)[2]="GO"
  
  for (i in 1:nrow(pathways)) {
    subtable<-pathways[i,]
    rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$GO[1], ',')[[1]]))
    pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
    pairtable<-as.data.frame(pairtable)
    pairtable$Genes<-rownames(pairtable)
    rownames(pairtable)<-1:nrow(pairtable)
    pairtable<-as.data.frame(pairtable)
    pairtable.new<-pairtable %>% gather(GO, pair, c(1:ncol(pairtable)-1))
    pairtable.new<-pairtable.new[,c(1:2)]
    GenesGOpair.1v1<-rbind(GenesGOpair.1v1, pairtable.new)
  }
  GenesGOpair.1v1<-subset(GenesGOpair.1v1, 
                          is.na(GenesGOpair.1v1$Genes)==FALSE, 
                          select = c("GO", "Genes"))
  GenesGOpair.1v1<-unique(GenesGOpair.1v1)
  
  write.table(GenesGOpair.1v1, 
              file = "eggnog-online/Gsoja.GO.1v1.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  rm(pairtable, pairtable.new, rcnames, subtable, pathways)
}
go2name<-read.delim("D:/3enrichment/go2name.txt", 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
GenesGOpair.1v1<-read.delim(file = "eggnog-online/Gsoja.GO.1v1.txt",
                            header = TRUE)

##Output of gene color from WGCNA.R
{
setwd("WGCNA_Gsoja/")
DEG<-read.delim("gene.color.txt", 
                header = TRUE)
names(DEG)[1]="Genes"
names(DEG)[2]="Groups"
DEG$Genes<-gsub("Glysoja.", "", DEG$Genes)
}
##Output of target prediction
{
  target3p<-read.delim("armswitch_target/Gsoja_3p.txt", 
                       skip = 1, header = TRUE)
  target3p<-separate(target3p, "Target_Acc.", c("UTR", "mRNA"),
                     sep = "\\.", remove = FALSE, convert = TRUE)
  target3p$gene<-target3p$mRNA
  target3p$Groups<-paste(target3p$miRNA_Acc., "3'-UTR", sep = " to ")
  target3p<-target3p[,c("gene","miRNA_Acc.", "Groups")]
  target3p<-unique(target3p)
  
  target5p<-read.delim("armswitch_target/Gsoja_5p.txt", 
                       skip = 1, header = TRUE)
  target5p<-separate(target5p, "Target_Acc.", c("UTR", "mRNA"),
                     sep = "\\.", remove = FALSE, convert = TRUE)
  target5p$gene<-target5p$mRNA
  target5p$Groups<-paste(target5p$miRNA_Acc., "5'-UTR", sep = " to ")
  target5p<-target5p[,c("gene","miRNA_Acc.", "Groups")]
  target5p<-unique(target5p)
  
  targetcds<-read.delim("armswitch_target/Gsoja_cds.txt", 
                        skip = 1, header = TRUE)
  targetcds<-separate(targetcds, "Target_Acc.", c("UTR", "mRNA"),
                     sep = "\\.", remove = FALSE, convert = TRUE)
  targetcds$gene<-targetcds$mRNA
  targetcds$Groups<-paste(targetcds$miRNA_Acc., "cds", sep = " to ")
  targetcds<-targetcds[,c("gene","miRNA_Acc.", "Groups")]
  targetcds<-unique(targetcds)
  
  DEG<-rbind(target3p, target5p, targetcds)
  DEG1<-DEG
  DEG1$Groups<-paste(DEG1$miRNA_Acc., "gene", sep = " to ")
  DEG1<-unique(DEG1)
  DEG<-rbind(DEG,DEG1)
  names(DEG)[1]<-"Genes"
  
  setwd("armswitch_target/Gsoja/")
}

addline_format <- function(x,...){ 
  gsub('_','\n',x) 
} 

setwd("significantmodule/")
modulelist<-c("blue", "turquoise", "grey")
DEG<-DEG[which(DEG$Groups %in% modulelist),]
####KOG enrich
{
  KOG.all.1<-compareCluster(Genes ~ Groups, 
                            data = DEG, 
                            fun = 'enricher',
                            TERM2GENE = GenesKOGpair.1v1,
                            TERM2NAME = kog2name,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 1,
                            minGSSize = 1,
                            maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KOG.all.1)
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[3]<-"kogClass"
  kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)
  plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  plotdata$Description<-paste(plotdata$Description,"(", plotdata$BGnumerator, ")")
  
  write.table(plotdata, file = "KOGenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  a<-ggplot(plotdata, aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "KOG Enrichment", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
    #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("KOGp1.bybackground.tiff", width = 10, height = 8, units = "in", dpi = 300)
  ggsave("KOGp1.bybackground.png", width = 10, height = 8, units = "in", dpi = 300)
  
  a<-ggplot(plotdata, aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "KOG Enrichment", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
    #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("KOGp1.bygroup.tiff", width = 10, height = 8, units = "in", dpi = 300)
  ggsave("KOGp1.bygroup.png", width = 10, height = 8, units = "in", dpi = 300)
}

####KEGG enrich
{
  KEGG.all.1<-compareCluster(Genes ~ Groups, 
                             data = DEG, 
                             fun = 'enricher',
                             TERM2GENE = GenesKEGGpair.1v1,
                             TERM2NAME = kegg2name,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             minGSSize = 1,
                             maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KEGG.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  plotdata$Description<-paste(plotdata$Description,"(", plotdata$BGnumerator, ")")
  
  write.table(plotdata, file = "KEGG.enrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  a<-ggplot(subset(plotdata, ratio1 >= 0.05 & p.adjust < 0.2), 
            aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "KEGG Enrichment", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
    #                 labels = addline_format(c("Knot\n(117)", "Oidia\n(19)\n \nup", "Scl\n(286)", "Knot\n(151)", "Oidia\n(211)\n \ndown", "Scl\n(113)")))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("KEGGp1.bybackground.r05p20.tiff", width = 10, height = 20, units = "in", dpi = 300)
  ggsave("KEGGp1.bybackground.r05p20.png", width = 10, height = 20, units = "in", dpi = 300)
  
  a<-ggplot(subset(plotdata, ratio2 >= 0.05 & p.adjust < 1),
            aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "KEGG Enrichment", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
    #                 labels = addline_format(c("Knot\n(117)", "Oidia\n(19)\n \nup", "Scl\n(286)", "Knot\n(151)", "Oidia\n(211)\n \ndown", "Scl\n(113)")))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("KEGGp1.bygroup.r05.tiff", width = 10, height = 10, units = "in", dpi = 300)
  ggsave("KEGGp1.bygroup.r05.png", width = 10, height = 10, units = "in", dpi = 300)
}

####GO enrich
{
  GO.all.1<-compareCluster(Genes ~ Groups, 
                           data = DEG, 
                           fun = 'enricher',
                           TERM2GENE = GenesGOpair.1v1,
                           TERM2NAME = go2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 2000000)
  
  #Plot
  plotin<-as.data.frame(GO.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[3]<-"goClass"
  go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)
  
  plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  plotdata$Description<-paste(plotdata$Description,"(", plotdata$BGnumerator, ")")
  
  write.table(plotdata, file = "GOenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  ##all GO
  {
    a<-ggplot(subset(plotdata, ratio1 >=0.2 & p.adjust < 0.05), 
              aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.05))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GOp05.bybackground.r20p05.tiff", width = 15, height = 60, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("GOp05.bybackground.r20p05.png", width = 15, height = 60, units = "in", dpi = 300, limitsize = FALSE)
    
    a<-ggplot(subset(plotdata, ratio2 >=0.1 & p.adjust < 0.2), 
              aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GO.bygroup.r10p20.tiff", width = 15, height = 60, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("GO.bygroup.r10p20.png", width = 15, height = 60, units = "in", dpi = 300, limitsize = FALSE)
  }
  
  ##BP
  {
    a<-ggplot(subset(plotdata, ratio1 >=0.25 & p.adjust < 0.05 & ONTOLOGY == "Biological Process" ), 
              aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment: Biological Process", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.05))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GO.BP.bybackground.r25p05.tiff", width = 15, height = 20, units = "in", dpi = 300)
    ggsave("GO.BP.bybackground.r25p05.png", width = 15, height = 20, units = "in", dpi = 300)
    
    a<-ggplot(subset(plotdata, ratio2 >= 0.1 & p.adjust < 0.20 & ONTOLOGY == "Biological Process" ), 
              aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment: Biological Process", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.20))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GO.BP.bygroup.r10p20.tiff", width = 10, height = 8, units = "in", dpi = 300)
    ggsave("GO.BP.bygroup.r10p20.png", width = 10, height = 8, units = "in", dpi = 300)
  }
  
  ##MF
  {
    a<-ggplot(subset(plotdata, ratio1 >=0.05 & p.adjust < 0.05 & ONTOLOGY == "Molecular Function" ), 
              aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment: Molecular Function", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.05))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GO.MF.bybackground.r05p05.tiff", width = 10, height = 20, units = "in", dpi = 300)
    ggsave("GO.MF.bybackground.r05p05.png", width = 10, height = 20, units = "in", dpi = 300)
    
    a<-ggplot(subset(plotdata, ratio2 >= 0.03 & p.adjust < 0.20 & ONTOLOGY == "Molecular Function" ), 
              aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment: Molecular Function", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.20))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GO.MF.bygroup.r03p20.tiff", width = 12, height = 8, units = "in", dpi = 300)
    ggsave("GO.MF.bygroup.r03p20.png", width = 12, height = 8, units = "in", dpi = 300)
  }
  
  ##CC
  {
    a<-ggplot(subset(plotdata, ratio1 >=0.05 & p.adjust < 0.20 & ONTOLOGY == "Cellular Component" ), 
              aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment: Cellular Component", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.05))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GO.CC.bybackground.r0p20.tiff", width = 10, height = 20, units = "in", dpi = 300)
    ggsave("GO.CC.bybackground.r0p20.png", width = 10, height = 20, units = "in", dpi = 300)
    
    a<-ggplot(subset(plotdata, ratio2 >= 0.03 & p.adjust < 0.20 & ONTOLOGY == "Cellular Component" ), 
              aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment: Cellular Component", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.20))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GO.CC.bygroup.r03p20.tiff", width = 10, height = 15, units = "in", dpi = 300)
    ggsave("GO.CC.bygroup.r03p20.png", width = 10, height = 15, units = "in", dpi = 300)
  }
}
