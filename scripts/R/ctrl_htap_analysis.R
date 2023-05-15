## CTRL + HTAP analysis

# library
library(Seurat)
# library(dplyr)
# library(pheatmap)
# library(velocyto.R)
# library(tidyverse)
library(ggplot2)
# library(RColorBrewer)
# #library(DoubletFinder)
# library(gplots)
# #library(fgsea)
# #library(biomaRt)
# library(knitr)
# library(xtable)
# #library(ggrepel)
# #library(matchSCore2)
# library(cowplot)
# library(scater)
# library(scran)
# library(SingleCellExperiment)
# #library(DropletUtils)
# library(batchelor)
# library(BiocNeighbors)
# library(SeuratWrappers)
# library(edgeR)
# library(limma)
# library(EnhancedVolcano)
# library(DOSE)
# library(pathview)
# library(clusterProfiler)
# library(enrichplot)
# library(org.Hs.eg.db)
# library(Nebulosa)
# library(monocle3)


setwd("/data/data_mbouamboua/10x_tu/R_scripts")
#source("00.import.R")
ctrlhtap <- readRDS("/data/data_mbouamboua/10x_tu/Rmd/ctrlhtap.rds")

# Delete sexual sex‐determining gene 

ctrlhtap[["sex"]]=""
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "C686518"),] <- "M"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "C681725"),] <- "F"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "C681086"),] <- "F"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "C685116"),] <- "F"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "C676136"),] <- "F"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "H247"),] <- "M"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "H248"),] <- "M"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "H230"),] <- "F"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "H231"),] <- "F"
ctrlhtap[['sex']][which(ctrlhtap@meta.data$SAMPLE == "H250"),] <- "F"



sex <- subset(ctrlhtap, sex == "M") 
DefaultAssay(sex) <- "RNA"
sex <- NormalizeData(sex)
FeaturePlot(sex, features = c("DDX3Y", "XIST", "RPS4Y1", "EIF1AY"))

VlnPlot(sex, features = c("DDX3Y", "XIST", "RPS4Y1", "EIF1AY"), pt.size = 0, group.by  = "PATIENT")


sex <- subset(ctrlhtap, sex == "F") 
DefaultAssay(sex) <- "RNA"
sex <- NormalizeData(sex)
# FeaturePlot(sex, features = c("DDX3Y", "XIST", "RPS4Y1", "LARS2", "SUN2", "PCDH1", "PPP1R37", "KDM5D", "ZFY", "EIF2S3Y", "EIF1AY"))
# VlnPlot(sex, features = c("DDX3Y", "XIST", "RPS4Y1", "LARS2", "SUN2", "PCDH1", "PPP1R37", "KDM5D", "ZFY", "EIF2S3Y", "EIF1AY"), pt.size = 0, group.by  = "PATIENT")
FeaturePlot(sex, features = c("DDX3Y", "XIST", "RPS4Y1", "EIF1AY"))
VlnPlot(sex, features = c("DDX3Y", "XIST", "RPS4Y1",  "EIF1AY"), pt.size = 0, group.by  = "PATIENT")

rm(sex)

# Delete sexual sex‐determining gene 
counts <- GetAssayData(ctrlhtap, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c("DDX3Y", "XIST", "RPS4Y1",  "EIF1AY"))),]
ctrlhtap <- subset(ctrlhtap, features = rownames(counts))

rm(counts)


# DimPlot integrated
# Find cluster
DefaultAssay(ctrlhtap) <- "integrated"
ElbowPlot(ctrlhtap)
ctrlhtap <- FindNeighbors(object = ctrlhtap, do.plot=TRUE, dims = 1:20)
ctrlhtap <- FindClusters(object = ctrlhtap, resolution=0.2)
ctrlhtap[['seurat_clusters']][which(ctrlhtap@active.ident == 6),] <- 1
ctrlhtap[['seurat_clusters']][which(ctrlhtap@active.ident == 7),] <- 1
ctrlhtap <- SetIdent(ctrlhtap, value = "seurat_clusters")

ctrlhtap$pheno <- paste0(ctrlhtap$seurat_clusters, "_", ctrlhtap$PHENOTYPE)
ctrlhtap$pat <- paste0(ctrlhtap$seurat_clusters, "_", ctrlhtap$PATIENT)

DimPlot(object = ctrlhtap, label=TRUE) + labs(title = "Clustering", colour = "Cell clusters")
DimPlot(object = ctrlhtap, group.by = "PATIENT") + labs(title = "Integration 5 CTRL + 5 HTAP", colour = "Patient")



# DEA cluster 0 - cluster 1 contrast

ctrlhtap <- SetIdent(ctrlhtap, value = "seurat_clusters")
DefaultAssay(ctrlhtap) <- "RNA"
markers <- FindMarkers(object = ctrlhtap, ident.1=0, ident.2=1, only.pos = FALSE)
#print(kable(markers[1:20,]))

EnhancedVolcano(markers, lab = rownames(markers), x = 'avg_logFC', y = 'p_val_adj', title = "Volcano plot - DEA cluster 0 vs 1", subtitle = "", drawConnectors = TRUE,
                widthConnectors = 0.2)

# clusterProfiler
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

original_gene_list <- markers$avg_logFC
names(original_gene_list) <- rownames(markers)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
keytypes(org.Hs.eg.db)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
#emapplot(gse, showCategory = 10)
#cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
#ridgeplot(gse) + labs(x = "enrichment distribution")

#gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
#terms <- gse$Description[1:3]
#pmcplot(terms, 2010:2020, proportion=FALSE)

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
df <- markers[rownames(markers) %in% dedup_ids$SYMBOL,]

kegg_gene_list <- df$avg_logFC
names(kegg_gene_list) <- dedup_ids$ENTREZID
gene_list<-na.omit(kegg_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 20, title = "Enriched Pathways in cluster 0" , split=".sign") + facet_grid(.~.sign)
# emapplot(kk2)
# cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
# ridgeplot(kk2) + labs(x = "enrichment distribution")

