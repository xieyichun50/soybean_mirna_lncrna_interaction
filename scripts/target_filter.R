library(tidyr)
library(stringr)

##mrna and mirna
cor.mrna2mirna<-read.table("run/Gsoja/selected_miRNA/PCC_Gsoja_mRNA_miRNA_corr.all.txt", header = TRUE)
#cor.mrna2mirna$gene<-gsub(".Wm82.a2.v1", "", cor.mrna2mirna$r.signif.mRNA)
cor.mrna2mirna$gene<-gsub("Glysoja.", "", cor.mrna2mirna$r.signif.mRNA)
cor.mrna2mirna$miRNA<-cor.mrna2mirna$r.signif.miRNA

#tar.mrna3p2mirna<-read.delim("PCC/psRNATarget_Gsoja_mRNA_miRNA_3p.txt", header = TRUE)
tar.mrna3p2mirna<-read.delim("PCC/psRNATarget_Gsoja_mRNA_miRNA_3p.txt", header = TRUE,skip = 1)
tar.mrna3p2mirna<-separate(tar.mrna3p2mirna, "Target_Acc.", c("UTR", "mRNA"),
                           sep = "\\.", remove = FALSE, convert = TRUE)
                           #sep = ";Parent=", remove = FALSE, convert = TRUE)
tar.mrna3p2mirna$gene<-gsub(".\\d.Wm82.a2.v1", "", tar.mrna3p2mirna$mRNA)
tar.mrna3p2mirna$miRNA<-tar.mrna3p2mirna$miRNA_Acc.
#tar.mrna3p2mirna<-tar.mrna3p2mirna[,-c(2)]
#tar.mrna3p2mirna<-unique(tar.mrna3p2mirna)

cor.tar.mrna3p2mirna<-merge(cor.mrna2mirna, tar.mrna3p2mirna,
                            by = c("gene", "miRNA"), all.x = TRUE)
cor.tar.mrna3p2mirna<-subset(cor.tar.mrna3p2mirna, is.na(Target_Acc.) ==FALSE)

#tar.mrna5p2mirna<-read.csv("PCC/psRNATarget_Gsoja_mRNA_miRNA_5p.csv", header = TRUE)
tar.mrna5p2mirna<-read.delim("PCC/psRNATarget_Gsoja_mRNA_miRNA_5p.txt", header = TRUE,skip = 1)
tar.mrna5p2mirna<-separate(tar.mrna5p2mirna, "Target_Acc.", c("UTR", "mRNA"),
                           sep = "\\.", remove = FALSE, convert = TRUE)
                           #sep = ";Parent=", remove = FALSE, convert = TRUE)
tar.mrna5p2mirna$gene<-gsub(".\\d.Wm82.a2.v1", "", tar.mrna5p2mirna$mRNA)
tar.mrna5p2mirna$miRNA<-tar.mrna5p2mirna$miRNA_Acc.

cor.tar.mrna5p2mirna<-merge(cor.mrna2mirna, tar.mrna5p2mirna,
                            by = c("gene", "miRNA"), all.x = TRUE)
cor.tar.mrna5p2mirna<-subset(cor.tar.mrna5p2mirna, is.na(Target_Acc.) ==FALSE)

cor.tar.mrna2mirna<-rbind(cor.tar.mrna3p2mirna, cor.tar.mrna5p2mirna)

tar.mrnacds2mirna<-read.delim("PCC/psRNATarget_Gsoja_mRNA_miRNA_cds.txt", header = TRUE,skip = 1)
tar.mrnacds2mirna<-separate(tar.mrnacds2mirna, "Target_Acc.", c("UTR", "mRNA"),
                           sep = "\\.", remove = FALSE, convert = TRUE)
                          #sep = ";Parent=", remove = FALSE, convert = TRUE)
tar.mrnacds2mirna$gene<-gsub(".\\d.Wm82.a2.v1", "", tar.mrnacds2mirna$mRNA)
tar.mrnacds2mirna$miRNA<-tar.mrnacds2mirna$miRNA_Acc.
tar.mrnacds2mirna$Multiplicity<-NA
tar.mrnacds2mirna$alignment<-NA
cor.tar.mrnacds2mirna<-merge(cor.mrna2mirna, tar.mrnacds2mirna,
                            by = c("gene", "miRNA"), all.x = TRUE)
cor.tar.mrnacds2mirna<-subset(cor.tar.mrnacds2mirna, is.na(Target_Acc.) ==FALSE)
cor.tar.mrnacds2mirna$UPE.<-cor.tar.mrnacds2mirna$UPE
cor.tar.mrnacds2mirna<-subset(cor.tar.mrnacds2mirna, 
       select = c("gene","miRNA","r.signif.mRNA", "r.signif.miRNA",
                  "r.signif", "p.signif", "p.adj.signif",
                  "miRNA_Acc.","Target_Acc.","UTR","mRNA", "Expectation", "UPE.", "miRNA_start",
                  "miRNA_end", "Target_start", "Target_end", "miRNA_aligned_fragment",
                  "alignment", "Target_aligned_fragment",
                  "Inhibition", "Target_Desc.", "Multiplicity"))
cor.tar.mrna2mirna<-rbind(cor.tar.mrna2mirna, cor.tar.mrnacds2mirna)
setwd("run/Gsoja/selected_miRNA/")
write.table(cor.tar.mrna2mirna, 
            #file = "run/Gsoja/Gsoja_cor.tar.mrna2mirna.txt",
            file = "Gsoja_cor.tar.mrna2mirna.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)

##enrich analysis

kegg2ont <- read.delim("D:/3enrichment/kegglevel.AC.txt", 
                       sep = "\t", colClasses = "character")

kegg2name <- read.delim("D:/3enrichment/kegg2name.txt", 
                        sep = "\t", colClasses = "character")

kog2name<-read.delim("D:/3enrichment/kog2name.txt", 
                     sep = "\t", colClasses = "character")

go2name<-read.delim("D:/3enrichment/go2name.txt", 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"


#read in GO levels
GOl3<-read.delim("D:/3enrichment/GOlevel3.txt", 
                 sep = "\t", colClasses = "character",
                 header = TRUE)
GOl4<-read.delim("D:/3enrichment/GOlevel4.txt", 
                 sep = "\t", colClasses = "character",
                 header = TRUE)
GOl5<-read.delim("D:/3enrichment/GOlevel5.txt", 
                 sep = "\t", colClasses = "character",
                 header = TRUE)


GenesKOGpair.1v1<-read.delim(file = "eggnog-online/Gsoja.KOG.1v1.txt",
                             header = TRUE)
GenesKEGGpair.1v1<-read.delim(file = "eggnog-online/Gsoja.KEGG.1v1.txt",
                              header = TRUE)
GenesGOpair.1v1<-read.delim(file = "eggnog-online/Gsoja.GO.1v1.txt",
                            header = TRUE)

setwd("run/Gsoja/selected_miRNA/")
DEG<-read.delim("Gsoja_cor.tar.mrna2mirna.txt", header = TRUE)

DEG<-read.delim("run/Gmax/selected_miRNA/Gmax_cor.tar.mrna2mirna.txt", header = TRUE)
setwd("run/Gmax/selected_miRNA/")
names(DEG)[1]="Genes"
names(DEG)[2]="Groups"



##mir166m arm switching
DEG<-DEG[which(DEG$Groups != "gma-miR166m"),c(1:2)]

arm<-read.delim("target_mir166m-3p.txt", header = TRUE,skip = 1)
arm<-separate(arm, "Target_Acc.", c("UTR", "mRNA"),
                           sep = "\\.", remove = FALSE, convert = TRUE)
#sep = ";Parent=", remove = FALSE, convert = TRUE)
arm$Genes<-gsub(".\\d.Wm82.a2.v1", "", arm$mRNA)
arm$Groups<-gsub("gma-MIR166m", "gma-miR166m", arm$miRNA_Acc.)

arm<-subset(arm, select = c("Genes", "Groups"))
arm<-unique(arm)
DEG<-rbind(DEG, arm)

arm<-read.delim("target_mir166m-5p.txt", header = TRUE,skip = 1)
arm<-separate(arm, "Target_Acc.", c("UTR", "mRNA"),
              sep = "\\.", remove = FALSE, convert = TRUE)
#sep = ";Parent=", remove = FALSE, convert = TRUE)
arm$Genes<-gsub(".\\d.Wm82.a2.v1", "", arm$mRNA)
arm$Groups<-gsub("gma-MIR166m", "gma-miR166m", arm$miRNA_Acc.)

arm<-subset(arm, select = c("Genes", "Groups"))
arm<-unique(arm)
DEG<-rbind(DEG, arm)

arm<-read.delim("target_mir166m-cds.txt", header = TRUE,skip = 1)
arm<-separate(arm, "Target_Acc.", c("UTR", "mRNA"),
              sep = "\\.", remove = FALSE, convert = TRUE)
#sep = ";Parent=", remove = FALSE, convert = TRUE)
arm$Genes<-gsub(".\\d.Wm82.a2.v1", "", arm$mRNA)
arm$Groups<-gsub("gma-MIR166m", "gma-miR166m", arm$miRNA_Acc.)

arm<-subset(arm, select = c("Genes", "Groups"))
arm<-unique(arm)
DEG<-rbind(DEG, arm)

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
  KOG.summary<-as.data.frame(xtabs(~Description, plotdata))
  names(KOG.summary)[2]="Number of miRNA"
  write.table(KOG.summary, file = "KOGenrich.summary.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  list(plotdata$Group)
  a<-ggplot(plotdata, aes(x = Cluster, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "KOG Enrichment", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
    #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
    guides(size = guide_legend(order = 1))+
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
  #ggsave("KOGp1.bybackground.tiff", width = 20, height = 8, units = "in", dpi = 300)
  ggsave("KOGp1.bybackground.png", width = 15, height = 8, units = "in", dpi = 300)
  
  a<-ggplot(plotdata, aes(x = Cluster, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "KOG Enrichment", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
    #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
    guides(size = guide_legend(order = 1))+
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
  #ggsave("KOGp1.bygroup.tiff", width = 20, height = 8, units = "in", dpi = 300)
  ggsave("KOGp1.bygroup.png", width = 15, height = 8, units = "in", dpi = 300)
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
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE)
  names(kegg2ont)[2]="ID"
  plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = "KEGG.enrich-clean.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  KEGG.summary<-as.data.frame(xtabs(~Description, plotdata))
  plotdata$Description<-paste(plotdata$Description,"(", plotdata$BGnumerator, ")")
  
  write.table(plotdata, file = "KEGG.enrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  names(KEGG.summary)[2]="Number of miRNA"
  write.table(KEGG.summary, file = "KEGGenrich.summary.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  plotdata<-subset(plotdata, is.na(Description) == FALSE & p.adjust <= 1)
  plotdata$ONTOLOGY<-gsub(" ", "\n", plotdata$ONTOLOGY)
  list(plotdata$Group)
  
  a<-ggplot(subset(plotdata, ratio1 >= 0 & p.adjust < 1 & ONTOLOGY != "Human\nDiseases"), 
            aes(x = Cluster, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "KEGG Enrichment", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
    #                 labels = addline_format(c("Knot\n(117)", "Oidia\n(19)\n \nup", "Scl\n(286)", "Knot\n(151)", "Oidia\n(211)\n \ndown", "Scl\n(113)")))+
    guides(size = guide_legend(order = 1))+
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
  #ggsave("KEGGp1.bybackground.tiff", width = 20, height = 30, units = "in", dpi = 300)
  ggsave("KEGGp1.bybackground.png", width = 15, height = 25, units = "in", dpi = 300)
  
  a<-ggplot(subset(plotdata, ratio2 >= 0 & p.adjust < 1 & ONTOLOGY != "Human\nDiseases"),
            aes(x = Cluster, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "KEGG Enrichment", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
    #                 labels = addline_format(c("Knot\n(117)", "Oidia\n(19)\n \nup", "Scl\n(286)", "Knot\n(151)", "Oidia\n(211)\n \ndown", "Scl\n(113)")))+
    guides(size = guide_legend(order = 1))+
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
  #ggsave("KEGGp1.bygroup.r05.tiff", width = 20, height = 30, units = "in", dpi = 300)
  ggsave("KEGGp1.bygroup.r05.png", width = 15, height = 25, units = "in", dpi = 300)
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
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[3]<-"goClass"
  go2name$ONTOLOGY<-gsub("biological_process", "Biological\n Process", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("molecular_function", "Molecular\n Function", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("cellular_component", "Cellular\n Component", go2name$ONTOLOGY)
  
  plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  plotdata$Description<-paste(plotdata$Description,"(", plotdata$BGnumerator, ")")
  
  write.table(plotdata, file = "GOenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  GO.summary<-as.data.frame(xtabs(~goName, plotdata))
  names(GO.summary)[2]="Number of miRNA"
  write.table(GO.summary, file = "GOenrich.summary.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  #plotdata<-read.delim("GOenrich.txt", header = TRUE)
  ##all GO
  { 
    test<-plotdata[-which(plotdata$goClass %in% GOl3$goClass),]
    a<-ggplot(subset(plotdata, ratio1 >=0 & p.adjust < 1 & is.na(Description)==FALSE), 
              aes(x = Cluster, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      guides(size = guide_legend(order = 1))+
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
    #ggsave("GOp05.bybackground.r10p05.tiff", width = 20, height = 65, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("GOp05.bybackground.L4.png", width = 10, height = 15, units = "in", dpi = 300, limitsize = FALSE)
    
    a<-ggplot(subset(test, ratio2 >=0 & p.adjust < 1 & is.na(Description)==FALSE), 
              aes(x = Cluster, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      guides(size = guide_legend(order = 1))+
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
    #ggsave("GO.bygroup.r10p20.tiff", width = 30, height = 50, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("GO.bygroup.L4.png", width = 10, height = 15, units = "in", dpi = 300, limitsize = FALSE)
  }
  
  ##BP
  {
    a<-ggplot(subset(plotdata, ratio1 >=0.10 & p.adjust < 0.20 & ONTOLOGY == "Biological\n Process" ), 
              aes(x = Cluster, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment: Biological Process", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.20))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      guides(size = guide_legend(order = 1))+
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
    #ggsave("GO.BP.bybackground.r20p05.tiff", width = 12, height = 25, units = "in", dpi = 300)
    ggsave("GO.BP.bybackground.r10p20.png", width = 12, height = 25, units = "in", dpi = 300)
    
    a<-ggplot(subset(plotdata, ratio2 >= 0.1 & p.adjust < 0.2 & ONTOLOGY == "Biological\n Process" ), 
              aes(x = Cluster, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment: Biological Process", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.20))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      guides(size = guide_legend(order = 1))+
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
    #ggsave("GO.BP.bygroup.r10p20.tiff", width = 10, height = 30, units = "in", dpi = 300)
    ggsave("GO.BP.bygroup.r10p20.png", width = 12, height = 25, units = "in", dpi = 300, limitsize = FALSE)
  }
  
  ##MF
  {
    plotdata$Description[c(1223,1224)]<-"oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of one atom of oxygen ( 4 )"
                              
    a<-ggplot(subset(plotdata, ratio1 >=0.1 & p.adjust < 0.2 & ONTOLOGY == "Molecular\n Function" ), 
              aes(x = Cluster, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment: Molecular Function", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.05))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      guides(size = guide_legend(order = 1))+
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
    #ggsave("GO.MF.bybackground.r20p05.tiff", width = 18, height = 15, units = "in", dpi = 300)
    ggsave("GO.MF.bybackground.r10p20.png", width = 15, height = 12, units = "in", dpi = 300)
    
    a<-ggplot(subset(plotdata, ratio2 >= 0.1 & p.adjust < 0.20 & ONTOLOGY == "Molecular\n Function" ), 
              aes(x = Cluster, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment: Molecular Function", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.20))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      guides(size = guide_legend(order = 1))+
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
    #ggsave("GO.MF.bygroup.r10p20.tiff", width = 10, height = 8, units = "in", dpi = 300)
    ggsave("GO.MF.bygroup.r10p20.png", width = 18, height = 12, units = "in", dpi = 300)
  }
  
  ##CC
  {
    a<-ggplot(subset(plotdata, ratio1 > 0 & p.adjust < 0.2 & ONTOLOGY == "Cellular\n Component" ), 
              aes(x = Cluster, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment: Cellular Component", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      guides(size = guide_legend(order = 1))+
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
    #ggsave("GO.CC.bybackground.r10p05.tiff", width = 8, height = 15, units = "in", dpi = 300)
    ggsave("GO.CC.bybackground.p20.png", width = 12, height = 18, units = "in", dpi = 300)
    
    a<-ggplot(subset(plotdata, ratio2 >= 0 & p.adjust < 0.2 & ONTOLOGY == "Cellular\n Component" ), 
              aes(x = Cluster, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment: Cellular Component", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.20))+
      #scale_x_discrete(limits = c("Knot.up", "Oidia.up", "Scl.up", "Knot.down", "Oidia.down", "Scl.down"), 
      #                 labels = addline_format(c("Knot\n(515)", "Oidia\n(69)\n \nup", "Scl\n(1110)", "Knot\n(435)", "Oidia\n(131)\n \ndown", "Scl\n(433)")))+
      guides(size = guide_legend(order = 1))+
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
    #ggsave("GO.CC.bygroup.r10.tiff", width = 10, height = 15, units = "in", dpi = 300)
    ggsave("GO.CC.bygroup.p20.png", width = 12, height = 18, units = "in", dpi = 300, limitsize = FALSE)
  }
}

