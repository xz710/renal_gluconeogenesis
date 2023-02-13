#GSEA----
library(fgsea)

#wrangling data

which(DEG_upper$...1 != entrezid_upper_shaped$SYMBOL)

entrezid_upper_shaped <- entrezid_upper[-1326,]
entrezid_upper_shaped <- entrezid_upper_shaped[-3240,]

DEG_upper_withid <- cbind(DEG_upper, entrezid_upper_shaped[,2])

gseaDat <- filter(DEG_upper_withid, !is.na(DEG_upper_withid$`entrezid_upper_shaped[, 2]`))
dim(gseaDat)

# creat ranks
ranks <- gseaDat$avg_logFC
names(ranks) <- gseaDat$`entrezid_upper_shaped[, 2]`
head(ranks)
barplot(sort(ranks, decreasing = T))


# get pathways
pathway <- read_tsv("C:/Users/Xinyi/Downloads/c2.cp.kegg.v7.5.entrez.gmt")
pathways <- pathway%>%
  unite("gene id", "79868":"4245", remove=F, sep="/")
pathways <- pathways[,1:3]
colnames(pathways) <- c("pathway","description","Gene ID")



fgseaRes <- fgsea(pathways, ranks, minSize=1, maxSize = 500, scoreType = "pos")
