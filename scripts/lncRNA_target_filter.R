library(tidyr)
library(stringr)
##Target prediction result
{
tar.lnc2mirna<-read.delim("PCC/psRNATarget_Gmax_miRNA_lncRNA.txt", header = TRUE)
tar.lnc2mirna$mID<-str_sub(tar.lnc2mirna$Target_Acc.,1,nchar(tar.lnc2mirna$Target_Acc.)-2)
tar.lnc2mirna$miRNA<-tar.lnc2mirna$miRNA_Acc.
ID<-read.table("PCC/Gmax_selected_id4Nong.txt")
names(ID)[1]="lncRNA"
names(ID)[2]="mID"
ID.1v1<-matrix(NA, nrow = 1, ncol = 2)
ID.1v1<-as.data.frame(ID.1v1)
names(ID.1v1)[1]="lncRNA"
names(ID.1v1)[2]="mID"

for (i in 1:nrow(ID)) {
  subtable<-ID[i,]
  rcnames<-list(c(strsplit(subtable$lncRNA[1], ',')[[1]]),c(strsplit(subtable$mID[1], ',')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$lncRNA<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(mID, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  ID.1v1<-rbind(ID.1v1, pairtable.new)
}
ID.1v1<-subset(ID.1v1, 
               is.na(ID.1v1$lncRNA)==FALSE,
               select = c("mID", "lncRNA"))
ID.1v1<-unique(ID.1v1)
rm(pairtable, pairtable.new, rcnames, subtable, ID)

ID.1v1$mID<-gsub("C08_MSTRG", "Gmax_MSTRG", ID.1v1$mID)
tar.lnc2mirna<-merge(tar.lnc2mirna, ID.1v1, by = "mID", all.x = TRUE)

##lncrna and mirna
cor.lncrna2mirna<-read.table("run/Gmax/selected_miRNA/PCC_Gmax_miRNA_lncRNA_corr.all.txt", header = TRUE)
cor.lncrna2mirna$lncRNA<-cor.lncrna2mirna$r.signif.lncRNA
cor.lncrna2mirna$miRNA<-cor.lncrna2mirna$r.signif.miRNA

cor.tar.lncrna2mirna<-merge(cor.lncrna2mirna, tar.lnc2mirna,
                            by = c("lncRNA", "miRNA"), all.x = TRUE)
cor.tar.lncrna2mirna<-subset(cor.tar.lncrna2mirna, is.na(Target_Acc.) ==FALSE)

write.table(cor.tar.lncrna2mirna, 
            file = "run/Gmax/selected_miRNA/Gmax_cor.tar.lncrna2mirna_all.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)

cor.lncrna2mirna<-read.table("run/Gmax/PCC_Gmax_miRNA_lncRNA_corr.all.txt", header = TRUE)
cor.lncrna2mirna$lncRNA<-cor.lncrna2mirna$r.signif.lncRNA
cor.lncrna2mirna$miRNA<-cor.lncrna2mirna$r.signif.miRNA

cor.tar.lncrna2mirna<-merge(cor.lncrna2mirna, tar.lnc2mirna,
                            by = c("lncRNA", "miRNA"), all.x = TRUE)
cor.tar.lncrna2mirna<-subset(cor.tar.lncrna2mirna, is.na(Target_Acc.) ==FALSE)

write.table(cor.tar.lncrna2mirna, 
            file = "run/Gmax/Gmax_cor.tar.lncrna2mirna_all.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)
}

##Target prediction result
{
  tar.lnc2mirna<-read.delim("PCC/psRNATarget_Gsoja_miRNA_lncRNA.txt", header = TRUE, skip = 1)
  tar.lnc2mirna$mID<-str_sub(tar.lnc2mirna$Target_Acc.,1,nchar(tar.lnc2mirna$Target_Acc.)-2)
  tar.lnc2mirna$mID<-gsub(".$", "", tar.lnc2mirna$mID)
  tar.lnc2mirna$miRNA<-tar.lnc2mirna$miRNA_Acc.
  ID<-read.table("PCC/Gsoja_selected_id4Nong.txt")
  names(ID)[1]="lncRNA"
  names(ID)[2]="mID"
  ID.1v1<-matrix(NA, nrow = 1, ncol = 2)
  ID.1v1<-as.data.frame(ID.1v1)
  names(ID.1v1)[1]="lncRNA"
  names(ID.1v1)[2]="mID"
  
  for (i in 1:nrow(ID)) {
    subtable<-ID[i,]
    rcnames<-list(c(strsplit(subtable$lncRNA[1], ',')[[1]]),c(strsplit(subtable$mID[1], ',')[[1]]))
    pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
    pairtable<-as.data.frame(pairtable)
    pairtable$lncRNA<-rownames(pairtable)
    rownames(pairtable)<-1:nrow(pairtable)
    pairtable<-as.data.frame(pairtable)
    pairtable.new<-pairtable %>% gather(mID, pair, c(1:ncol(pairtable)-1))
    pairtable.new<-pairtable.new[,c(1:2)]
    ID.1v1<-rbind(ID.1v1, pairtable.new)
  }
  ID.1v1<-subset(ID.1v1, 
                 is.na(ID.1v1$lncRNA)==FALSE,
                 select = c("mID", "lncRNA"))
  ID.1v1<-unique(ID.1v1)
  rm(pairtable, pairtable.new, rcnames, subtable, ID)
  
  #ID.1v1$mID<-gsub("C08_MSRG", "Gmax_MSTRG", ID.1v1$mID)
  tar.lnc2mirna<-merge(tar.lnc2mirna, ID.1v1, by = "mID", all.x = TRUE)
  
  ##lncrna and mirna
  cor.lncrna2mirna<-read.table("run/Gsoja/selected_miRNA/PCC_Gsoja_miRNA_lncRNA_corr.all.txt", header = TRUE)
  cor.lncrna2mirna$lncRNA<-cor.lncrna2mirna$r.signif.lncRNA
  cor.lncrna2mirna$miRNA<-cor.lncrna2mirna$r.signif.miRNA
  
  cor.tar.lncrna2mirna<-merge(cor.lncrna2mirna, tar.lnc2mirna,
                              by = c("lncRNA", "miRNA"), all.x = TRUE)
  cor.tar.lncrna2mirna<-subset(cor.tar.lncrna2mirna, is.na(Target_Acc.) ==FALSE)
  
  write.table(cor.tar.lncrna2mirna, 
              file = "run/Gsoja/selected_miRNA/Gsoja_cor.tar.lncrna2mirna_all.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  cor.lncrna2mirna<-read.table("run/Gsoja/PCC_Gsoja_miRNA_lncRNA_corr.all.txt", header = TRUE)
  cor.lncrna2mirna$lncRNA<-cor.lncrna2mirna$r.signif.lncRNA
  cor.lncrna2mirna$miRNA<-cor.lncrna2mirna$r.signif.miRNA
  
  cor.tar.lncrna2mirna<-merge(cor.lncrna2mirna, tar.lnc2mirna,
                              by = c("lncRNA", "miRNA"), all.x = TRUE)
  cor.tar.lncrna2mirna<-subset(cor.tar.lncrna2mirna, is.na(Target_Acc.) ==FALSE)
  
  write.table(cor.tar.lncrna2mirna, 
              file = "run/Gsoja/Gsoja_cor.tar.lncrna2mirna_all.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
}