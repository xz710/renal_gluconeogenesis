setwd("C:/Users/Xinyi/Desktop/omics/human kidney gluconeogenesis")

# get KEGG data----
library(clusterProfiler)
human_pathway <- search_kegg_organism('hsa', by='kegg_code')

#load DEG list----
library(tidyverse)
library(rhdf5)
library(edgeR)

DEGList <- read_tsv("C:/Users/Xinyi/Desktop/omics/human kidney/PCT diabete DEG list.txt")

DEG_lower <- DEGList[(which(DEGList$avg_logFC < 0)),]

#get entrezid 
library(org.Hs.eg.db)
library(AnnotationDbi)
hs <- org.Hs.eg.db

#lowregulated gene list----
my.symbols_lower <- DEG_lower$...1
entrezid_lower<- AnnotationDbi::select(hs, 
                                       keys = my.symbols_lower,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")

#gp analysis
goanalysis_lower <- enrichGO(entrezid_lower$ENTREZID, 
                             OrgDb = org.Hs.eg.db, ont="BP", pvalueCutoff=0.05)


goresult_low <- goanalysis_lower@result

#gluconeogenesis gene go analysis----
gluconeogenesis_gene <- goresult_low[which(goresult$Description == "gluconeogenesis"),]

(gluconeogenesis_gene$geneID)
gluconeogenesis_low_gene_entrezid <- c(5903,2821,2033,84912,23411,5465,468,4190,8850,1196,
                                       2819,2806,10296,51547,2542,2875,7167,5164,2648,5091,
                                       5106,4485)
gluconeogenesis_gene_id <- entrezid_lower[(which(entrezid_lower$ENTREZID %in% gluconeogenesis_low_gene_entrezid)),]

gluconeogenesis_low <- DEG_lower[which(DEG_lower$...1 %in% gluconeogenesis_gene_id$SYMBOL),]

write_tsv(gluconeogenesis_low, "gluconeogenesis go analysis gene downregulated.txt")

# adrenergic gene go analysis (adrenergic receptor signaling pathway)----
adrenergic_gene <- goresult[which(goresult$Description == "adrenergic receptor signaling pathway"),]

(adrenergic_gene$geneID)
adrenergic_low_gene_entrezid <- c(115,9693,147)
adrenergic_gene_id <- entrezid_lower[(which(entrezid_lower$ENTREZID %in% adrenergic_low_gene_entrezid)),]

adrenergic_low <- DEG_lower[which(DEG_lower$...1 %in% adrenergic_gene_id$SYMBOL),]







