setwd("C:/Users/Xinyi/Desktop/omics/human kidney gluconeogenesis")

# get KEGG data----
library(clusterProfiler)
human_pathway <- search_kegg_organism('hsa', by='kegg_code')

#load DEG list----
library(tidyverse)
library(rhdf5)
library(edgeR)

DEGList <- read_tsv("C:/Users/Xinyi/Desktop/omics/human kidney/dct&ct.txt")

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
adrenergic_up_gene_entrezid <- c(493,9693,5144,340481,11214)
adrenergic_gene_id <- entrezid_upper[(which(entrezid_upper$ENTREZID %in% adrenergic_up_gene_entrezid)),]

adrenergic_up <- DEG_upper[which(DEG_upper$...1 %in% adrenergic_gene_id$SYMBOL),]

write_tsv(adrenergic_up, "adrenergic go analysis gene collecting duct.txt")





