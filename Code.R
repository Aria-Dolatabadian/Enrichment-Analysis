#https://rpubs.com/jrgonzalezISGlobal/enrichment


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("airway")



if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")



if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")



if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("enrichplot")


library(airway)
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)



library(airway)
data(airway)
table(airway$dex)


library(edgeR)
dge <- DGEList(assay(airway), group = airway$dex)
dge <- calcNormFactors(dge)
mm <- model.matrix( ~ group, data=dge$samples)

keep.exprs <- filterByExpr(dge)
dge.filt <- dge[keep.exprs,]
dim(dge.filt)


mm <- model.matrix( ~ group, data=dge.filt$samples)
v <- voom(dge.filt, design = mm, plot = TRUE)


fit <- lmFit(v, mm)
fit <- eBayes(fit)
topTable(fit)

tt <- topTable(fit, n=Inf)
mask <- tt$adj.P.Val < 0.05 &
        abs(tt$logFC) > log2(2)
deGenes <- rownames(tt[mask, ])
head(deGenes)

length(deGenes)

geneUniverse <- rownames(tt)
length(geneUniverse)


library(org.Hs.eg.db)
deGenes <- unlist(mget(deGenes, envir=org.Hs.egENSEMBL2EG,
                       ifnotfound = NA))

geneUniverse <- unlist(mget(geneUniverse, envir=org.Hs.egENSEMBL2EG,
                       ifnotfound = NA))



library(clusterProfiler)
ans.go <- enrichGO(gene = deGenes, ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   universe = geneUniverse,
                   readable=TRUE,
                   pvalueCutoff = 0.05)
tab.go <- as.data.frame(ans.go)
tab.go<- subset(tab.go, Count>5)
tab.go[1:5, 1:6]



ans.kegg <- enrichKEGG(gene = deGenes,
                       organism = 'hsa',
                       universe = geneUniverse,
                       pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(ans.kegg)
tab.kegg<- subset(tab.kegg, Count>5)
tab.kegg[1:5, 1:6]


gda <- read.delim("curated_gene_disease_associations.tsv.gz")
disease2gene <- gda[, c("diseaseId", "geneId")]
disease2name <- gda[, c("diseaseId", "diseaseName")]


ans.dis <- enricher(deGenes, TERM2GENE=disease2gene,
                    TERM2NAME=disease2name)
tab.dis <- as.data.frame(ans.dis)
tab.dis<- subset(tab.dis, Count>5)
tab.dis[,1:6]


c3.tf <- read.gmt("c3.tft.v6.2.entrez.gmt")

ans.tf <- enricher(deGenes, TERM2GENE=c3.tf)
tab.tf <- as.data.frame(ans.tf)
tab.tf<- subset(tab.tf, Count>5)
tab.tf[1:5,1:5]


library(enrichplot)
p1 <- barplot(ans.dis, showCategory=10)
p1


p2 <- dotplot(ans.kegg, showCategory=20) + ggtitle("KEGG")
p3 <- dotplot(ans.dis, showCategory=20) + ggtitle("Disease")
plot_grid(p2, p3, nrow=2)


p4 <- upsetplot(ans.dis)
p4

p5 <- emapplot(ans.kegg)   #Error
p5


cowplot::plot_grid(p1, p3, p5, ncol=2, labels=LETTERS[1:3])

sessionInfo()

