setwd("C:/Users/Xinyi/Desktop/omics/human kidney gluconeogenesis")

# get KEGG data----
library(clusterProfiler)
human_pathway <- search_kegg_organism('hsa', by='kegg_code')

#load DEG list----
library(tidyverse)
library(rhdf5)
library(edgeR)

DEGList <- read_tsv("C:/Users/Xinyi/Desktop/omics/human kidney/PCT diabete DEG list.txt")

DEG_upper <- DEGList[(which(DEGList$avg_logFC > 0)),]

#get entrezid 
library(org.Hs.eg.db)
library(AnnotationDbi)
hs <- org.Hs.eg.db

#upregulated gene list----
my.symbols_upper <- DEG_upper$...1
entrezid_upper<- AnnotationDbi::select(hs, 
                                       keys = my.symbols_upper,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")

#gp analysis
goanalysis_upper <- enrichGO(entrezid_upper$ENTREZID, 
                             OrgDb = org.Hs.eg.db, ont="BP", pvalueCutoff=0.05)


goresult <- goanalysis_upper@result

#gluconeogenesis gene go analysis----
gluconeogenesis_gene <- goresult[which(goresult$Description == "gluconeogenesis"),]

(gluconeogenesis_gene$geneID)
gluconeogenesis_up_gene_entrezid <- c(2203,1407,5105,5771,8604,23516,8473,467,10165,57223,7874,
                                      10447,10776,2805,5230,1642,57486,8864,3953,5236,10891,
                                      2023,83667,57001,2597,229)
gluconeogenesis_gene_id <- entrezid_upper[(which(entrezid_upper$ENTREZID %in% gluconeogenesis_up_gene_entrezid)),]

gluconeogenesis_up <- DEG_upper[which(DEG_upper$...1 %in% gluconeogenesis_gene_id$SYMBOL),]

write_tsv(gluconeogenesis_up, "gluconeogenesis go analysis gene.txt")

# adrenergic gene go analysis (adrenergic receptor signaling pathway)----
adrenergic_gene <- goresult[which(goresult$Description == "adrenergic receptor signaling pathway"),]

(adrenergic_gene$geneID)
adrenergic_up_gene_entrezid <- c(5144,340481,57561,5142,493,92255,11214)
adrenergic_gene_id <- entrezid_upper[(which(entrezid_upper$ENTREZID %in% adrenergic_up_gene_entrezid)),]

adrenergic_up <- DEG_upper[which(DEG_upper$...1 %in% adrenergic_gene_id$SYMBOL),]

write_tsv(adrenergic_up, "adrenergic go analysis gene.txt")


#go analysis of all DEG
my.symbols <- DEGList$...1
entrezid<- AnnotationDbi::select(hs, 
                                       keys = my.symbols,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")

#gp analysis
goanalysis_all<- enrichGO(entrezid$ENTREZID, 
                             OrgDb = org.Hs.eg.db, ont="BP", pvalueCutoff=0.05)


goresult_all <- goanalysis_all@result

#gluconeogenesis gene go analysis----
#gluconeogenesis_gene_all <- goresult_all[which(goresult$Description == "gluconeogenesis"),]

#(gluconeogenesis_gene_all$geneID)
#gluconeogenesis_all_gene_entrezid <- c(27072,143187,4734,950,55187,57154,7879,54832,
#                                       30849,79158,9648,10490,79720,6653,8546,9765,
#                                      137492,55048,9570,23163,3920,23230,8943,1191,
#                                     8031,374666,200312,6272)
#gluconeogenesis_gene_all <- entrezid_upper[(which(entrezid$ENTREZID %in% gluconeogenesis_all_gene_entrezid)),]

#gluconeogenesis_all <- DEGList[which(DEGList$...1 %in% gluconeogenesis_gene_all$SYMBOL),]

#write_tsv(gluconeogenesis_all, "gluconeogenesis DEG.txt")


