#GSM3823939	Control Sample 1 (snRNAseq)
#GSM3823940	Control Sample 2 (snRNAseq)
#GSM3823941	Control Sample 3 (snRNAseq)
#GSM3823942	Diabetes Sample 1 (snRNAseq)
#GSM3823943	Diabetes Sample 2 (snRNAseq)
#GSM3823944	Diabetes Sample 3 (snRNAseq)

setwd("C:/Users/Xinyi/Desktop/omics")

library(tidyverse)
library(rhdf5)
library(edgeR)

archs4.human <- "human_matrix_v11.h5" # if you placed the hdf5 file in your working directory, just use "human_matrix_v8.h5" as the path
all.samples.human <- h5read(archs4.human, name="meta/samples/geo_accession")
mySamples <- c("GSM3823939","GSM3823940","GSM3823941","GSM3823942","GSM3823943","GSM3823944") 

my.sample.locations <- which(all.samples.human %in% mySamples)
genes <- h5read(archs4.human, "meta/genes/gene_symbol")
expression <- h5read(archs4.human, "data/expression", 
                     index=list(my.sample.locations, 1:length(genes)))

expression <- t(expression)
rownames(expression) <- genes
colnames(expression) <- all.samples.human[my.sample.locations]
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)

keepers <- rowSums(archs4.cpm>1)>=3
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")
archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

sample_source_name <- h5read(archs4.human, "meta/samples/source_name_ch1")
sample_title <- h5read(archs4.human, name="meta/samples/title")
sample_characteristics<- h5read(archs4.human, name="meta/samples/characteristics_ch1")

studyDesign <- tibble(Sample_title = sample_title[my.sample.locations], 
                      Sample_source = sample_source_name[my.sample.locations],
                      Sample_characteristics = sample_characteristics[my.sample.locations])

StudyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      patient_condition = c("control", "control", "control", "diabetes", "diabetes"))

pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
pca.res.df <- as_tibble(pca.res$x)

pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color=StudyDesign$patient_condition,) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()


gluconeogenesis_DEG_name <- c(gluconeogenesis_low[1:9,]$...1,gluconeogenesis_up[1:25,]$...1)

gluconeogenesis_DEG_expression <- archs4.filtered.norm.log2.cpm[(which(rownames(archs4.filtered.norm.log2.cpm)%in% gluconeogenesis_DEG_name)),]
expression_data <- data.frame(gluconeogenesis_DEG_expression)
