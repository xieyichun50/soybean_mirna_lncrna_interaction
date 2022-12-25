library(psych)
library(tidyr)
library(dplyr)
exp.mRNA<-read.csv(file = "run/Gmax_mRNA_filter_gene_count_matrix.csv", 
                   header = TRUE)
DEG.mrna<-read.table("run/DEGlist.Gmax.txt", header = TRUE)
exp.mRNA<-exp.mRNA[which(exp.mRNA$gene_id %in% DEG.mrna$gene_id),]
row.names(exp.mRNA)=exp.mRNA$gene_id
exp.mRNA<-exp.mRNA[,2:ncol(exp.mRNA)]
exp.mRNA$max = apply(exp.mRNA, 1, max)
exp.mRNA$mean = apply(exp.mRNA, 1, mean)
exp.mRNA<-exp.mRNA[which(exp.mRNA$mean > 10 & exp.mRNA$max >= 50), 1:(ncol(exp.mRNA)-2)]
exp.mRNA<-log2(exp.mRNA+1)
exp.mRNA<-t(exp.mRNA)
  
exp.miRNA<-read.csv(file = "run/Gmax_miRNA_soybean_salt_smallRNA.miRBase_22.1_mature.gma.uniq_seq.norm.csv",
                    header = TRUE)
DEG.miRNA<-read.table("run/mirnadeglist.txt", header = TRUE)
exp.miRNA<-exp.miRNA[which(exp.miRNA$gene_id %in% DEG.miRNA[,1]),]
row.names(exp.miRNA)=exp.miRNA$gene_id
exp.miRNA<-exp.miRNA[,2:ncol(exp.miRNA)]
exp.miRNA$max = apply(exp.miRNA, 1, max)
exp.miRNA$mean = apply(exp.miRNA, 1, mean)
exp.miRNA<-exp.miRNA[which(exp.miRNA$mean > 1 & exp.miRNA$max >= 50), 1:(ncol(exp.miRNA)-2)]
exp.miRNA<-log2(exp.miRNA+1)
exp.miRNA<-t(exp.miRNA)

exp.lncRNA<-read.csv(file = "run/Gmax_lncRNA_NPC_Cufflinks_FPKM4R.csv",
                     header = TRUE)
DEG.lncRNA<-read.table("run/Gmax_NPC_TS_Cuffdiff_time_series_EI_SR_R_LI_SR_L_stranded_FPKM4R.txt", header = TRUE)
exp.lncRNA<-exp.lncRNA[which(exp.lncRNA$gene_id %in% DEG.lncRNA[,1]),]
row.names(exp.lncRNA)=exp.lncRNA$gene_id
exp.lncRNA<-exp.lncRNA[,2:ncol(exp.lncRNA)]
exp.lncRNA$max = apply(exp.lncRNA, 1, max)
exp.lncRNA$mean = apply(exp.lncRNA, 1, mean)
exp.lncRNA<-exp.lncRNA[which(exp.lncRNA$max > 1), 1:(ncol(exp.lncRNA)-2)]
exp.lncRNA<-log2(exp.lncRNA+1)
exp.lncRNA<-t(exp.lncRNA)

##Pearson correlation coefficient
##miRNA_lncRNA
cor.value.mi2lnc<-corr.test(exp.lncRNA, exp.miRNA, method = "pearson", adjust = "BH", alpha = 0.05)
ind.signif=which(cor.value.mi2lnc$p<0.05 & cor.value.mi2lnc$p.adj<0.05,arr.ind = T)
r.signif<-cor.value.mi2lnc$r[ind.signif]
p.signif<-cor.value.mi2lnc$p[ind.signif]
p.adj.signif<-cor.value.mi2lnc$p.adj[ind.signif]
r.signif.lncRNA<-dimnames(cor.value.mi2lnc$r)[[1]][ind.signif[,1]]
r.signif.miRNA<-dimnames(cor.value.mi2lnc$r)[[2]][ind.signif[,2]]
corr.filter<-cbind(r.signif.lncRNA, r.signif.miRNA, r.signif, p.signif, p.adj.signif)
corr.filter.merge<-subset(corr.filter, abs(r.signif) > 0.8)

ind.signif=which(cor.value.mi2lnc$p<=1 & cor.value.mi2lnc$p.adj<=1,arr.ind = T)
r.signif<-cor.value.mi2lnc$r[ind.signif]
p.signif<-cor.value.mi2lnc$p[ind.signif]
p.adj.signif<-cor.value.mi2lnc$p.adj[ind.signif]
r.signif.lncRNA<-dimnames(cor.value.mi2lnc$r)[[1]][ind.signif[,1]]
r.signif.miRNA<-dimnames(cor.value.mi2lnc$r)[[2]][ind.signif[,2]]
corr.filter<-cbind(r.signif.lncRNA, r.signif.miRNA, r.signif, p.signif, p.adj.signif)
corr.filter.merge<-subset(corr.filter, abs(r.signif) > 0)

write.table(corr.filter.merge, 
            file = "PCC/PCC_Gmax_miRNA_lncRNA_corr_all.txt", 
            sep = '\t', row.names = FALSE,
            quote = FALSE)

rm(cor.value.mi2lnc, ind.signif, 
   r.signif, p.signif, p.adj.signif, r.signif.miRNA, r.signif.lncRNA,
   corr.filter, corr.filter.merge)

##lncRNA_mRNA
cor.value.lnc2mrna<-corr.test(exp.mRNA, exp.lncRNA, method = "pearson", adjust = "BH", alpha = 0.05)
ind.signif=which(cor.value.lnc2mrna$p<0.05 & cor.value.lnc2mrna$p.adj<0.05,arr.ind = T)
r.signif<-cor.value.lnc2mrna$r[ind.signif]
p.signif<-cor.value.lnc2mrna$p[ind.signif]
p.adj.signif<-cor.value.lnc2mrna$p.adj[ind.signif]
r.signif.mRNA<-dimnames(cor.value.lnc2mrna$r)[[1]][ind.signif[,1]]
r.signif.lncRNA<-dimnames(cor.value.lnc2mrna$r)[[2]][ind.signif[,2]]
corr.filter<-cbind(r.signif.mRNA, r.signif.lncRNA, r.signif, p.signif, p.adj.signif)
corr.filter.merge<-subset(corr.filter, abs(r.signif) > 0.8)

write.table(corr.filter.merge, 
            file = "PCC_Gmax_mRNA_lncRNA_corr.txt", 
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(cor.value.lnc2mrna, ind.signif, 
   r.signif, p.signif, p.adj.signif, r.signif.mRNA, r.signif.lncRNA,
   corr.filter, corr.filter.merge)

##miRNA_mRNA
cor.value.mrna2mirna<-corr.test(exp.mRNA, exp.miRNA, method = "pearson", adjust = "BH", alpha = 0.05)
ind.signif=which(cor.value.mrna2mirna$p<0.05 & cor.value.mrna2mirna$p.adj<0.05,arr.ind = T)
r.signif<-cor.value.mrna2mirna$r[ind.signif]
p.signif<-cor.value.mrna2mirna$p[ind.signif]
p.adj.signif<-cor.value.mrna2mirna$p.adj[ind.signif]
r.signif.mRNA<-dimnames(cor.value.mrna2mirna$r)[[1]][ind.signif[,1]]
r.signif.miRNA<-dimnames(cor.value.mrna2mirna$r)[[2]][ind.signif[,2]]
corr.filter<-cbind(r.signif.mRNA, r.signif.miRNA,r.signif,p.signif, p.adj.signif)
corr.filter.merge<-subset(corr.filter, abs(r.signif) > 0.8)

write.table(corr.filter.merge,
            file = "PCC_Gmax_mRNA_miRNA_corr.txt", 
            sep = '\t', row.names = FALSE,
            quote = FALSE)

rm(cor.value.mrna2mirna, ind.signif, 
   r.signif, p.signif, p.adj.signif, r.signif.mRNA, r.signif.miRNA,
   corr.filter, corr.filter.merge)