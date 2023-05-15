---
title: "R Notebook"
output: html_notebook
---

"CA4" = "Capillary"
"GJA5" = "Artery"
"LYVE1" = "Lymphatic"
"ACKR1" = "Veine"
"Bro-1-EC" = "MGP"
"Bro-2-EC" = "ARID5A"
"Cap EC" = "CAV1"
"CapA EC" = "FCN3"
"Cap-i1 EC" = "EDNRB"
"Cap-i2 EC" = "SPTBN1"

new.cluster.ids <- c("CA4", "GJA5", "LYVE1", "ACKR1", "FCN3", "EDNRB", "SPTBN1", "MGP", "ARID5A", "CAV1")



FeaturePlot(seurat.ctrl, reduction = "umap", label.size = 2, cols = c("grey", "red"), features = c("NRP1", "SELE", "PROX1", "LYVE1", "CA4", "IGFBP4", "IFI27", "MKI67", "CCL20"), min.cutoff = "q1")


DotPlot(seurat.ctrl, features =  c("PROX1", "LYVE1", "VWF", "ACKR3","ACKR3", "PDPN", "BGN", "CD93", "CA4",	"BMX", "GJA5", "IGFBP4", "IFI27", "MKI67", "CCL20","NRP1", "SELE", "PCDH1", "PLLP", "SGK1", "NOTCH4", "EFNB2","ID1", "EDN1", "KDR", "SOX18", "EPAS1", "RGCC", "HEY1"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis() + labs(x = "Markers", y = "Clusters")





# DotPlot(seurat.ctrl, features =  c("BTG2",
# "GM12216",
# "HILPDA",
# "IER3",
# "JUNB",
# "JUND",
# "PTMA",
# "SEMA3C",
# "SGK1",
# "TMSB4X",
# "UQCR11",
# "CSDE1",
# "GUK1",
# "STMN1",
# "TMSB10"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis() + labs(x = "Markers", y = "Clusters")
# 
# DotPlot(seurat.ctrl, features =c("PCDH1",
# "PHLDA3",
# "PLLP",
# "PMP22",
# "PTP4A3",
# "RGS12",
# "SEPT4",
# "STMN2",
# "TBX3",
# "TMCC2",
# "TMEM176A",
# "TPPP3",
# "TSPAN8",
# "HSPB1",
# "IGFBP7"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis() + labs(x = "Markers", y = "Clusters")


#my_levels <- features.EC
#aggr@active.ident <- factor(x = aggr@active.ident, levels = my_levels)
#saveRDS(aggr, "./output/htap.labels.rds")






# Cluster marker genes identification
library("scCATCH")

DefaultAssay(seurat.htap.ctrl) <- "RNA"
seurat.htap.ctrl <- NormalizeData(seurat.htap.ctrl)

# Step 1: prepare a Seurat object containing log1p normalized single-cell transcriptomic data matrix and the information of cell clusters.
# Note: please define the species for revising gene symbols. Human or Mouse. The default is to find potential marker genes for all clusters with the percentage of expressed cells (≥25%), using WRS test (P<0.05) and a log1p fold change of ≥0.25. These parameters are adjustable for users.

clu_markers <- findmarkergenes(seurat.htap.ctrl,
                               species = 'Human',
                               cluster = 'All',
                               match_CellMatch = TRUE,
                               cancer = NULL,
                               tissue = c("Adventitia", "Antecubital vein", "Artery", "Blood vessel", "Umbilical vein", "Airway epithelium", "Alveolus", "Bronchoalveolar system", "Lung"),
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)



# Step 2: evidence-based scoring and annotaion for identified potential marker genes of each cluster generated from findmarkergenes function.

clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Human',
                   cancer = NULL,
                   tissue = c("Adventitia", "Antecubital vein", "Artery", "Blood vessel", "Umbilical vein", "Airway epithelium", "Alveolus", "Bronchoalveolar system", "Lung"))












---
  title: "HTAP Project (Ly Tu) - Lib CTRL"
author: '[Yvon Mbouamboua](mailto:mbouamboua@ipmc.cnrs.fr) IPMC CNRS - Computational Biology'
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output:
  html_document:
  df_print: paged
toc: yes
toc_depth: '3'
html_notebook:
  code_folding: none
theme: journal
toc: yes
toc_depth: 3
toc_float: yes
pdf_document:
  toc: yes
toc_depth: '3'
---
  
  <style type="text/css">
  
  body, td {
    font-size: 15px;
  }
code.r{
  font-size: 15px;
}
pre {
  font-size: 15px
}
</style>
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(
  cache = FALSE,
  cache.lazy = FALSE,
  tidy = TRUE
)

setwd("/data/data_mbouamboua/10x_tu")
#annotations <- read.csv("data/annotation.csv")
patient_colors <- c("#b2df8a", "#20b2aa",  "#4292c6", "#bd0026", "#6a51a3")
#clust_colors <- c("#33a02c", "#ff7f00", "#08519c", "#c51bc7", "#e31a1c", "#b2df8a", "mediumblue", "#807dba")
# clust_colors <- c("#e31a1c",
#                   "#762a83", 
#                     "#08519c", 
#                     "mediumblue", 
#                     "#ff7f00", 
#                     "#33a02c", 
#                     "#b2df8a", 
#                     "#20b2aa", 
#                     "#c51b7d")


clust_colors <- c("#e31a1c",
                  "#bd0026", 
                  "#762a83",
                  "#6a51a3",
                  "#807dba",
                  "#08519c",
                  "#2171b5",
                  "#4292c6",
                  "#6baed6",
                  "#9ecae1",
                  "#c6dbef",
                  "#deebf7",
                  "mediumblue",
                  "#ff7f00", 
                  "#33a02c", 
                  "#b2df8a", 
                  "#20b2aa", 
                  "#c51b7d") 
source("00.import.R")

library(readxl)
EC <- read_excel("data/EC.xlsx")
EC_lambresch <- read_csv("data/Long_List_EC.csv")
#View(EC_lambresch)
```


# Loading data


# Loading data 10 samples

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}

A1 <- readRDS("/data/10x_data/10x_tu/output/A1.rds")
B1 <- readRDS("/data/10x_data/10x_tu/output/B1.rds")
C1 <- readRDS("/data/10x_data/10x_tu/output/C1.rds")
D1 <- readRDS("/data/10x_data/10x_tu/output/D1.rds")
E1 <- readRDS("/data/data_mbouamboua/10x_tu/output/E1.rds")

ctrl.list <- c(A1, B1, C1, D1, E1)
```

# Integration 5 CTRL 


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}
options(future.globals.maxSize = 80000 * 1024^2)
features <- SelectIntegrationFeatures(object.list = ctrl.list, nfeatures = 3000)
ctrl.list <- PrepSCTIntegration(object.list = ctrl.list, anchor.features = features)
ctrl.list <- lapply(X = ctrl.list, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = ctrl.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca", verbose = FALSE)
seurat.ctrl <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
seurat.ctrl <- RunPCA(seurat.ctrl, verbose = FALSE)
#seurat.ctrl <- RunFastMNN(object.list = SplitObject(seurat.ctrl, split.by = "SAMPLE"))
seurat.ctrl <- RunUMAP(seurat.ctrl, reduction = "pca", dims = 1:30)
seurat.ctrl <- RunTSNE(seurat.ctrl, reduction = "pca", dims = 1:30)
```


# Clustering, frequency

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=8}

seurat.ctrl <- FindNeighbors(seurat.ctrl, reduction = "pca", dims = 1:30)
seurat.ctrl <- FindClusters(seurat.ctrl, resolution = 0.2)
seurat.ctrl <- ScaleData(object = seurat.ctrl)

#saveRDS(seurat.ctrl, "./output/seurat.ctrl.rds")

# Visualize QC metrics as a violin plot 
VlnPlot(seurat.ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size = 0, ncol = 3, split.by = "PATIENT")

#DimPlot(object = seurat.ctrl, label = FALSE, group.by = "PHENOTYPE") + labs(colour = "Phenotype")
DimPlot(object = seurat.ctrl,  pt.size = 0.1, label = T, label.size = 6) + labs(colour = "Cell clusters")


DimPlot(object = seurat.ctrl, label = FALSE, group.by = "PATIENT", pt.size = 0.2, cols = patient_colors) + labs(colour = "Patients")



freq_table <- prop.table(x = table(seurat.ctrl@active.ident, seurat.ctrl@meta.data[, c("PATIENT")]), margin = 2)*100
freq_table <- round(freq_table, digits = 2)
#colnames(freq_table) <- c("CTRL Lib E", "CTRL Lib C", "CTRL Lib B", "CTRL Lib D", "CTRL Lib A" )

barplot(t(freq_table),
        las=1,
        horiz = F,
        ylim = c(0, max(freq_table)*1.1),
        #col = "#619CFF",
        ylab = "Cell numbers (%)",
        xlab = "Cell clusters",
        main = "Number of cells of control patients per cluster",
        legend = T,
        #col=c("#619CFF","#FF7A7A", "#55820f", "#CB760B", "#b2b210"),
        col= patient_colors,
        beside = T)
#box()
# B <- as.matrix(freq_table)
# text(B+2, A, labels = paste(as.character(B), "%", sep = ""), cex = 0.8, pos = 2)

```
# Cluster colrs

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=6}

# Define cluster colours
# Cluster 0
seurat.ctrl <- SetIdent(seurat.ctrl, value = "seurat_clusters")
seurat.ctrl[['cluster0']] <- " "
seurat.ctrl[['cluster0']][which(seurat.ctrl@active.ident == 0),] <- 0
seurat.ctrl <- SetIdent(seurat.ctrl, value = "cluster0")
DimPlot(seurat.ctrl, group.by = "cluster0", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)

seurat.ctrl <- SetIdent(seurat.ctrl, value = "seurat_clusters")
seurat.ctrl[['cluster1']] <- " "
seurat.ctrl[['cluster1']][which(seurat.ctrl@active.ident == 1),] <- 1
seurat.ctrl <- SetIdent(seurat.ctrl, value = "cluster1")
DimPlot(seurat.ctrl, group.by = "cluster1", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)

seurat.ctrl <- SetIdent(seurat.ctrl, value = "seurat_clusters")
seurat.ctrl[['cluster2']] <- " "
seurat.ctrl[['cluster2']][which(seurat.ctrl@active.ident == 2),] <- 2
seurat.ctrl <- SetIdent(seurat.ctrl, value = "cluster2")
DimPlot(seurat.ctrl, group.by = "cluster2", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)


seurat.ctrl <- SetIdent(seurat.ctrl, value = "seurat_clusters")
seurat.ctrl[['cluster3']] <- " "
seurat.ctrl[['cluster3']][which(seurat.ctrl@active.ident == 3),] <- 3
seurat.ctrl <- SetIdent(seurat.ctrl, value = "cluster3")
DimPlot(seurat.ctrl, group.by = "cluster3", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)


seurat.ctrl <- SetIdent(seurat.ctrl, value = "seurat_clusters")
seurat.ctrl[['cluster4']] <- " "
seurat.ctrl[['cluster4']][which(seurat.ctrl@active.ident == 4),] <- 4
seurat.ctrl <- SetIdent(seurat.ctrl, value = "cluster4")
DimPlot(seurat.ctrl, group.by = "cluster4", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)
```



# Clusters correlation heatmap

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}
seurat.ctrl <- SetIdent(seurat.ctrl, value = "seurat_clusters")
mnmat <- c()
uniq <- unique(seurat.ctrl@active.ident)
seurat.ctrl@meta.data$cluster <- seurat.ctrl@active.ident
for(i in 1:length(uniq)){
  mnmat <- cbind(mnmat, apply(as.matrix(seurat.ctrl@assays$RNA@data[, seurat.ctrl@meta.data$cluster==uniq[i]]), 1, mean))
}

colnames(mnmat) <- as.vector(unique(seurat.ctrl@active.ident))
ct=cor(mnmat)
pheatmap(ct, angle_col = 0, main = "Clusters correlation heatmap")
```


# Clusters gene markers heatmap

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=16, fig.width=8}
# Finding differentially expressed features (cluster biomarkers)
ctrl.markers <- FindAllMarkers(object = seurat.ctrl, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)

#View(markers)
top20 <- ctrl.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
DoHeatmap(seurat.ctrl, features = top20$gene, size=5, angle = 0) + labs(title = "CTRL dataset - Expression Heatmap of top 20 markers", colour = "Cell clusters")

write.table(markers, file="./output/seurat.ctrl.csv", sep="\t")


```




```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=12}
top <- markers %>% group_by(cluster) %>% top_n(nrow(markers), avg_logFC)
#View(top)
topFDR <- subset(top,
                 p_val_adj <= 0.05 & pct.1 >= 0.5 & cluster == 4)
topFDR <- topFDR[order(topFDR$pct.1, decreasing = T), ]
#View(topFDR)

top.genes <- topFDR$gene[1:40]

DotPlot(seurat.ctrl, features = top.genes, cols = c("green", "red"), dot.scale = 6) + RotatedAxis() + labs(title = "DGE of top 40 markers from CTRL dataset", x = "Markers", y = "Clusters")

#DotPlot(seurat.ctrl, features = c("ABCA4"), cols = c("blue", "red"), dot.scale = 6) + RotatedAxis()


```


```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=20, fig.width=12}
FeaturePlot(seurat.ctrl, reduction = "umap", cols = c("grey", "red"), features = top.genes, min.cutoff = "q1")
```




```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=8}
# top.markers <- c("TFPI", "MMNR1", "SSR3", "HOMER3", "PPARG", "HMGB2", "IFI27", "TNFAIP2")
# DotPlot(seurat.ctrl, features = top.markers, cols = c("blue", "red"), dot.scale = 6) + RotatedAxis()
# 
# cell.markers <- RenameIdents(seurat.ctrl, `0` = "TFPI", `1` = "SSR3", `2` = "HOMER3", `3` = "PPARG", `4` = "HMGB2", `5` = "IFI27", `6` = "TNFAIP2")
# DimPlot(object = cell.markers, label = TRUE, label.size = 5, pt.size = 0.5)+labs(colour = "Endothelial markers")
# 
# 
# 
# top.markers.clust <- c("FABP4", "FABP5",  "IGFBP7", "ALDH1A1", "MMRN1", "LMNA", "TXNRD1", "KRT7", "FTL", "MT2A", "S100A11", "S100A6", "GAPDH", "SRGN", "NQO1", "TKT")
# DotPlot(seurat.ctrl, features = top.markers.clust, cols = c("blue", "red"), dot.scale = 6) + RotatedAxis()

```




```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=20}
top <- ctrl.markers %>% group_by(cluster) %>% top_n(nrow(ctrl.markers), avg_logFC)

topFDR <- subset(top,
                 p_val_adj <= 0.05 & pct.1 >= 0.5)

top.CTRL <- topFDR %>% group_by(cluster) %>% top_n(10, avg_logFC)


DotPlot(seurat.ctrl, features = top.CTRL$gene, cols = c("green", "red"), dot.scale = 6, group.by = "cluster_ctrl") + RotatedAxis() + labs(title = "DGE of top 10 markers from CTRL dataset", x = "Markers", y = "Clusters")

```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}
christoph.markers <- c("LYVE1","VWF","PROX1", "PDPN", "NRP2", "EPHB2", "EPHB4", "ACKR3",  "BGN", "FLT4", "NRP1", "NR2F2",  "RGCC", "EPAS1", "SOX18", "NOTCH4", "KDR", "EDN1", "EFNB2", "ID1", "BMX", "HEY1", "IGFBP7")

DotPlot(seurat.ctrl, features = christoph.markers, cols = c("green", "red"), dot.scale = 6, group.by = "cluster_ctrl") + RotatedAxis() + labs(title = "", x = "Markers", y = "Clusters")


DotPlot(seurat.ctrl, features = christoph.markers, cols = c("green", "red"), dot.scale = 6) + RotatedAxis() + labs(title = "", x = "Markers", y = "Clusters")



```







```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=10}

# top1 <- markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
# genes.CTRL <- top1$gene 
# new.cluster.ids <- genes.CTRL
# names(x = new.cluster.ids) <- levels(x = seurat.ctrl)
# aggr <- RenameIdents(object = seurat.ctrl, new.cluster.ids)
# DimPlot(object = aggr, label = TRUE, label.size = 5, pt.size = 0.5)+labs(colour = "Endothelial markers")
# #DimPlot(object = seurat.ctrl, label = TRUE, label.size = 8, pt.size = 0.5) + labs(colour = "Cell clusters")
# DotPlot(seurat.ctrl, features =  genes.CTRL, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis() + labs(x = "Markers", y = "Clusters")
# FeaturePlot(seurat.ctrl, reduction = "umap", label.size = 2, cols = c("grey", "red"), features = genes.CTRL, min.cutoff = "q1")


# cell.types <- RenameIdents(seurat.ctrl, `0` = "BMP4", `1` = "NAV3", `2` = "REEP1", `3` = "MRPS6", `4` = "CD36", `5` = "BST2", `6` = "CTHRC1", `7` = "TUBB6", `8` = "TUBB6", `9` = "TNFAIP2", `10` = "FAM11B")
# DimPlot(object = cell.types, label = TRUE, label.size = 5, pt.size = 0.5)+labs(colour = "Endothelial markers")


```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=8}

# christoph.markers <- c("LYVE1","VWF","PROX1", "PDPN", "NRP2", "EPHB2", "EPHB4", "ACKR3",  "BGN", "FLT4", "NRP1", "NR2F2",  "RGCC", "EPAS1", "SOX18", "NOTCH4", "KDR", "EDN1", "EFNB2", "ID1", "BMX", "HEY1", "IGFBP7")
# 
# FeaturePlot(seurat.ctrl, reduction = "umap", label.size = 2, cols = c("grey", "red"), features = christoph.markers, min.cutoff = "q1")
# #FeaturePlot(seurat.ctrl, reduction = "umap", label.size = 2, cols = c("grey", "red"), features = "VWF", min.cutoff = "q1")
# #DimPlot(object = seurat.ctrl, label = TRUE, label.size = 5)+labs(colour = "Clusters")
# 
# christoph.ec.cell.types <- RenameIdents(seurat.ctrl, `0` = "EDN1", `1` = "NRP1", `2` = "PROX1", `3` = "LYVE1", `4` = "RGCC", `5` = "VWF", `6` = "NRP2", `7` = "EFNB2", `8` = "?", `9` = "HEY1", `10` = "NRP1")
# DimPlot(object = christoph.ec.cell.types, label = TRUE, label.size = 5, pt.size = 0.5)+labs(colour = "Endothelial markers")
# 
# christoph.ec <- RenameIdents(seurat.ctrl, `0` = "Capillary", `1` = "Artery", `2` = "Lymphatic", `3` = "Lymphatic", `4` = "Capillary", `5` = "Vein", `6` = "Lymphatic", `7` = "Capillary", `8` = "?", `9` = "Vein (ACKR3) ? Artery (HEY1) ?", `10` = "Artery")
# DimPlot(object = christoph.ec, label = F, pt.size = 0.5, label.size = 5)+labs(colour = "Endothelial cell types")
# DotPlot(seurat.ctrl, features =  christoph.markers, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis() + labs(x = "Markers", y = "Clusters")


cluster.ids <- christoph.markers
names(x = cluster.ids) <- levels(x = seurat.ctrl)
aggr <- RenameIdents(object = seurat.ctrl, cluster.ids)
DimPlot(object = aggr, label = TRUE, label.size = 5, pt.size = 0.5)+labs(colour = "Endothelial markers")



cell.types <- RenameIdents(seurat.ctrl, `0` = "Lymphatic", `1` = "?", `2` = "Cycling cells", `3` = "Vein")
DimPlot(object = cell.types, label = TRUE, label.size = 5, cols = clust_colors)+labs(colour = "Endothelial cells")

```


## Find the spécific markers per cluster

```{r}
# # Find markers of C1 (compared to a background of all other cells)
# c0.markers=FindMarkers(seurat.ctrl,1,6, only.pos = TRUE, test.use = "roc")
# 
# FeaturePlot(seurat.ctrl, reduction = "umap", cols = c("grey", "red"), features = "CDKN1C", min.cutoff = "q1")
# DotPlot(seurat.ctrl, features =  c("CDKN1C"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis() + labs(x = "Markers", y = "Clusters")
# VlnPlot(seurat.ctrl, features = c("CDKN1C"), pt.size = 0, ncol = 3)

```



# Session Info

```{r sessinf}
sessionInfo()
```






# HTAP.BMP9 vs HTAP


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=12}

A4 <- readRDS("/data/10x_data/10x_tu/output/A4.rds")
A3 <- readRDS("/data/10x_data/10x_tu/output/A3.rds")
B4 <- readRDS("/data/10x_data/10x_tu/output/B4.rds")
B3 <- readRDS("/data/10x_data/10x_tu/output/B3.rds")
C4 <- readRDS("/data/10x_data/10x_tu/output/C4.rds")
C3 <- readRDS("/data/10x_data/10x_tu/output/C3.rds")
D4 <- readRDS("/data/10x_data/10x_tu/output/D4.rds")
D3 <- readRDS("/data/10x_data/10x_tu/output/D3.rds")
E4 <- readRDS("/data/data_mbouamboua/10x_tu/output/E4.rds")
E3 <- readRDS("/data/data_mbouamboua/10x_tu/output/E3.rds")


htap.list <- c(A4, A3, B4, B3, C4, C3, D4, D3, E4, E3)

options(future.globals.maxSize = 80000 * 1024^2)
features <- SelectIntegrationFeatures(object.list = htap.list, nfeatures = 3000)
htap.list <- PrepSCTIntegration(object.list = htap.list, anchor.features = features)
htap.list <- lapply(X = htap.list, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = htap.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca", verbose = FALSE)
seurat.htap.bmp <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
seurat.htap.bmp <- RunPCA(seurat.htap.bmp, verbose = FALSE)
seurat.htap.bmp <- RunFastMNN(object.list = SplitObject(seurat.htap.bmp, split.by = "SAMPLE"))
seurat.htap.bmp <- RunUMAP(seurat.htap.bmp, reduction = "mnn", dims = 1:30)
seurat.htap.bmp <- RunTSNE(seurat.htap.bmp, reduction = "mnn", dims = 1:30)

metadata <- seurat.htap.bmp@meta.data
metadata$condition <- NA
metadata$condition[which(str_detect(metadata$TAG, "A3"))] <- "HTAP"
metadata$condition[which(str_detect(metadata$TAG, "B3"))] <- "HTAP"
metadata$condition[which(str_detect(metadata$TAG, "C3"))] <- "HTAP"
metadata$condition[which(str_detect(metadata$TAG, "D3"))] <- "HTAP"
metadata$condition[which(str_detect(metadata$TAG, "E3"))] <- "HTAP"

metadata$condition[which(str_detect(metadata$TAG, "A4"))] <- "HTAP.BMP9"
metadata$condition[which(str_detect(metadata$TAG, "B4"))] <- "HTAP.BMP9"
metadata$condition[which(str_detect(metadata$TAG, "C4"))] <- "HTAP.BMP9"
metadata$condition[which(str_detect(metadata$TAG, "D4"))] <- "HTAP.BMP9"
metadata$condition[which(str_detect(metadata$TAG, "E4"))] <- "HTAP.BMP9"

# Add metadata back to Seurat object
seurat.htap.bmp@meta.data <- metadata


# Clustering, frequency

DefaultAssay(seurat.htap.bmp) <- "integrated"
seurat.htap.bmp <- FindNeighbors(seurat.htap.bmp, reduction = "mnn", dims = 1:20)
seurat.htap.bmp <- FindClusters(seurat.htap.bmp, resolution = 0.15)

# Visualize QC metrics as a violin plot 
VlnPlot(seurat.htap.bmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size = 0, ncol = 3, split.by = "BMP9", c("#33a02c", "#ff7f00"))

#VlnPlot(seurat.htap.ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), split.by = "PHENOTYPE", cols = c("#33a02c", "#e31a1c"), pt.size = 0.0001, ncol = 3, label = TRUE)


DimPlot(object = seurat.htap.bmp, label = FALSE, group.by = "condition", cols = c("#33a02c", "#ff7f00")) + labs(colour = "Phenotype")
DimPlot(object = seurat.htap.bmp, label = TRUE) + labs(colour = "Cell clusters")
DimPlot(object = seurat.htap.bmp, label = FALSE, split.by = "condition") + labs(colour = "Cell clusters")

freq_table <- prop.table(x = table(seurat.htap.bmp@active.ident, seurat.htap.bmp@meta.data[, "condition"]), margin = 2)*100
freq_table <- round(freq_table, digits = 2)

barplot(t(freq_table),
        las=1,
        horiz = F,
        ylim = c(0, max(freq_table)*2.0),
        #col = "#619CFF",
        ylab = "Cell numbers (%)",
        xlab = "Cell clusters",
        main = "Number of cells by cluster",
        legend = T,
        col=c("#33a02c", "#ff7f00"),
        beside = F)
#B <- as.matrix(freq_table)
# text(B+2, A, labels = paste(as.character(B), "%", sep = ""), cex = 0.8, pos = 2)


DefaultAssay(seurat.htap.bmp) <- "RNA"
seurat.htap.bmp <- NormalizeData(seurat.htap.bmp)
sce.htap.bmp <- as.SingleCellExperiment(seurat.htap.bmp)
# DGE using limma
agg <- aggregateAcrossCells(sce.htap.bmp, ids = colData(sce.htap.bmp)[,c("seurat_clusters", "condition")])
colData(agg)[,c("seurat_clusters", "condition",  "ncells")]

group <- factor(agg$condition)
dge <- DGEList(assay(agg), group = group)

dge$samples

# Desing matrix HTAP and CTRL conditions
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

cont.matrix <- makeContrasts(HTAP.BMP.vs.HTAP = HTAP.BMP9 - HTAP,
                             #CTRL.vs.HTAP = CTRL - HTAP,
                             levels=colnames(design))

v <- voom(dge, design, plot=F)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
efit <- eBayes(vfit, robust = TRUE)
topTable <- limma::topTable(efit, number = "Inf")
#genes.HTAP.CTRL <- rownames(topTable)
#topTable$cluster <- rep("cluster8", nrow(topTable))

#top.markers.htap.vs.ctrl <- topTable[abs(topTable$logFC) >= 1.5 & topTable$adj.P.Val <= 0.05, ]
top.markers.htap.bmp.vs.htap <- topTable[abs(topTable$logFC) >= 1.5 & topTable$P.Value <= 10e-05, ]


#View(top.markers.htap.bmp.vs.htap)
top.genes.clusters <- rownames(top.markers.htap.bmp.vs.htap)




# Volcanoplot

library("EnhancedVolcano")


EnhancedVolcano(topTable,
                lab = rownames(topTable),
                x = 'logFC',
                y = 'P.Value',
                #xlim = c(-8, 8),
                ylim = c(0, -log10(10e-12)),
                title = "Differential Gene Expression \n
    HTAP.BMP vs HTAP",
                subtitle = "(linear model by limma package, 
    robust = true, FCcutoff = 1.5, pCutoff = 10e-05)" ,
                #selectLab = top.c8.genes,
                #xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.5,
                #labCol = 'black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                #drawConnectors = TRUE,
                #widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()
```








```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}

A1 <- readRDS("/data/10x_data/10x_tu/output/A1.rds")
A2 <- readRDS("/data/10x_data/10x_tu/output/A2.rds")
B1 <- readRDS("/data/10x_data/10x_tu/output/B1.rds")
B2 <- readRDS("/data/10x_data/10x_tu/output/B2.rds")
C1 <- readRDS("/data/10x_data/10x_tu/output/C1.rds")
C2 <- readRDS("/data/10x_data/10x_tu/output/C2.rds")
D1 <- readRDS("/data/10x_data/10x_tu/output/D1.rds")
D2 <- readRDS("/data/10x_data/10x_tu/output/D2.rds")
E1 <- readRDS("/data/data_mbouamboua/10x_tu/output/E1.rds")
E2 <- readRDS("/data/data_mbouamboua/10x_tu/output/E2.rds")


htap.list <- c(A1, A2, B1, B2, C1, C2, D1, D2, E1, E2)

options(future.globals.maxSize = 80000 * 1024^2)
features <- SelectIntegrationFeatures(object.list = htap.list, nfeatures = 3000)
htap.list <- PrepSCTIntegration(object.list = htap.list, anchor.features = features)
htap.list <- lapply(X = htap.list, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = htap.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca", verbose = FALSE)
seurat.ctrl.bmp <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
seurat.ctrl.bmp <- RunPCA(seurat.ctrl.bmp, verbose = FALSE)
seurat.ctrl.bmp <- RunFastMNN(object.list = SplitObject(seurat.ctrl.bmp, split.by = "SAMPLE"))
seurat.ctrl.bmp <- RunUMAP(seurat.ctrl.bmp, reduction = "mnn", dims = 1:30)
#seurat.ctrl.bmp <- RunTSNE(seurat.ctrl.bmp, reduction = "mnn", dims = 1:30)

metadata <- seurat.ctrl.bmp@meta.data
metadata$condition <- NA
metadata$condition[which(str_detect(metadata$TAG, "A1"))] <- "CTRL"
metadata$condition[which(str_detect(metadata$TAG, "B1"))] <- "CTRL"
metadata$condition[which(str_detect(metadata$TAG, "C1"))] <- "CTRL"
metadata$condition[which(str_detect(metadata$TAG, "D1"))] <- "CTRL"
metadata$condition[which(str_detect(metadata$TAG, "E1"))] <- "CTRL"

metadata$condition[which(str_detect(metadata$TAG, "A2"))] <- "CTRL.BMP9"
metadata$condition[which(str_detect(metadata$TAG, "B2"))] <- "CTRL.BMP9"
metadata$condition[which(str_detect(metadata$TAG, "C2"))] <- "CTRL.BMP9"
metadata$condition[which(str_detect(metadata$TAG, "D2"))] <- "CTRL.BMP9"
metadata$condition[which(str_detect(metadata$TAG, "E2"))] <- "CTRL.BMP9"

# Add metadata back to Seurat object
seurat.ctrl.bmp@meta.data <- metadata


# Clustering, frequency

#DefaultAssay(seurat.ctrl.bmp) <- "integrated"
seurat.ctrl.bmp <- FindNeighbors(seurat.ctrl.bmp, reduction = "mnn", dims = 1:20)
seurat.ctrl.bmp <- FindClusters(seurat.ctrl.bmp, resolution = 0.15)

# Visualize QC metrics as a violin plot 
VlnPlot(seurat.ctrl.bmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size = 0, ncol = 3, split.by = "BMP9", c("#762a83",  "#b2df8a"))

#VlnPlot(seurat.htap.ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), split.by = "PHENOTYPE", cols = c("#33a02c", "#e31a1c"), pt.size = 0.0001, ncol = 3, label = TRUE)


DimPlot(object = seurat.ctrl.bmp, label = FALSE, group.by = "condition", cols = c("#762a83",  "#b2df8a")) + labs(colour = "Phenotype", pt.size = 0.5)
DimPlot(object = seurat.ctrl.bmp, label = TRUE) + labs(colour = "Cell clusters")
DimPlot(object = seurat.ctrl.bmp, label = FALSE, split.by = "condition") + labs(colour = "Cell clusters")

freq_table <- prop.table(x = table(seurat.ctrl.bmp@active.ident, seurat.ctrl.bmp@meta.data[, "condition"]), margin = 2)*100
freq_table <- round(freq_table, digits = 2)

barplot(t(freq_table),
        las=1,
        horiz = F,
        ylim = c(0, max(freq_table)*2.1),
        #col = "#619CFF",
        ylab = "Cell numbers (%)",
        xlab = "Cell clusters",
        main = "Number of cells by cluster",
        legend = T,
        col=c("#762a83",  "#b2df8a"),
        beside = F)
#B <- as.matrix(freq_table)
# text(B+2, A, labels = paste(as.character(B), "%", sep = ""), cex = 0.8, pos = 2)


DefaultAssay(seurat.ctrl.bmp) <- "RNA"
seurat.ctrl.bmp <- NormalizeData(seurat.ctrl.bmp)
sce.htap.bmp <- as.SingleCellExperiment(seurat.ctrl.bmp)
# DGE using limma
agg <- aggregateAcrossCells(sce.htap.bmp, ids = colData(sce.htap.bmp)[,c("seurat_clusters", "condition")])
colData(agg)[,c("seurat_clusters", "condition",  "ncells")]

group <- factor(agg$condition)
dge <- DGEList(assay(agg), group = group)

dge$samples

# Desing matrix HTAP and CTRL conditions
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

cont.matrix <- makeContrasts(CTRL.BMP.vs.CTRL = CTRL.BMP9 - CTRL,
                             #CTRL.vs.HTAP = CTRL - HTAP,
                             levels=colnames(design))

v <- voom(dge, design, plot=F)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
efit <- eBayes(vfit, robust = TRUE)
topTable <- limma::topTable(efit, number = "Inf")
#genes.HTAP.CTRL <- rownames(topTable)
#topTable$cluster <- rep("cluster8", nrow(topTable))

#top.markers.htap.vs.ctrl <- topTable[abs(topTable$logFC) >= 1.5 & topTable$adj.P.Val <= 0.05, ]
top.markers.ctrl.bmp.vs.htap <- topTable[abs(topTable$logFC) >= 1.5 & topTable$P.Value <= 10e-05, ]


#View(top.markers.htap.bmp.vs.htap)
top.genes.clusters <- rownames(top.markers.ctrl.bmp.vs.htap)




# Volcanoplot

library("EnhancedVolcano")


EnhancedVolcano(topTable,
                lab = rownames(topTable),
                x = 'logFC',
                y = 'P.Value',
                #xlim = c(-8, 8),
                ylim = c(0, -log10(10e-12)),
                title = "Differential Gene Expression \n
    CTRL.BMP vs CTRL",
                subtitle = "(linear model by limma package, 
    robust = true, FCcutoff = 1.5, pCutoff = 10e-05)" ,
                #selectLab = top.c8.genes,
                #xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.5,
                #labCol = 'black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                #drawConnectors = TRUE,
                #widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()
```




# Cluster relabeling

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}

genes.HTAP.CTRL <- c("CDKN1C", "IFI27", "CDKN1C", "CCL20", "CA4")


FeaturePlot(seurat.htap.ctrl, reduction = "umap", label.size = 2, cols = c("grey", "red"), features = c("NRP1", "SELE", "PROX1"), min.cutoff = "q1")

DotPlot(seurat.htap.ctrl, features =   c("NRP1", "SELE", "PROX1"), cols = c("blue", "red"), dot.scale = 8, split.by = "PHENOTYPE") + RotatedAxis()

df <- data.frame("Symbols" = c("NRP1", "PROX1", "SELE"),
                 "Cell types" = c("Artery", "Lymphatic", "Vein"))
kable(df, caption = "Endothelial markers and cell types")

new.cluster.ids <- genes.HTAP.CTRL
names(x = new.cluster.ids) <- levels(x = seurat.htap.ctrl)
aggr <- RenameIdents(object = seurat.htap.ctrl, new.cluster.ids)
DimPlot(object = aggr, label = TRUE, label.size = 4, pt.size = 1)+labs(colour = "Endothelial markers")
DimPlot(object = seurat.htap.ctrl, label = TRUE, label.size = 8, pt.size = 1) + labs(colour = "Cell clusters")
DotPlot(seurat.htap.ctrl, features =  genes.HTAP.CTRL, cols = c("blue", "red"), dot.scale = 8, split.by = "PHENOTYPE") + RotatedAxis()


FeaturePlot(seurat.htap.ctrl, reduction = "umap", label.size = 2, pt.size = 1, cols = c("grey", "red"), features = genes.HTAP.CTRL, min.cutoff = "q1")


#my_levels <- features.EC

#aggr@active.ident <- factor(x = aggr@active.ident, levels = my_levels)

#saveRDS(aggr, "./output/htap.labels.rds")


cell.types <- RenameIdents(seurat.htap.ctrl, `0` = "Lymphatic", `1` = "Artery", `2` = "?", `3` = "Veine")
DimPlot(object = cell.types, label = TRUE, label.size = 5, pt.size = 1)+labs(colour = "Endothelial cells")


seurat.htap.ctrl.cell.types <- RenameIdents(seurat.htap.ctrl, `0` = "CDKN1C", `1` = "IFI27", `2` = "MKI67", 
                                            `3` = "CCL20")
DimPlot(object = seurat.htap.ctrl.cell.types, label = TRUE, label.size = 4)+labs(colour = "Markers")


DimPlot(object = seurat.htap.ctrl, label = TRUE, label.size = 4)+labs(colour = "Endothelial markers")

DimPlot(object = aggr, label = TRUE, label.size = 4)+labs(colour = "Endothelial markers")
DimPlot(seurat.htap.ctrl.cell.types, label = TRUE) + labs(colour = "Endothelial cells")


```



col <- c("#e31a1c",
         "#bd0026", 
         "#762a83",
         "#6a51a3",
         "#807dba",
         "#08519c",
         "#2171b5",
         "#4292c6",
         "#6baed6",
         "#9ecae1",
         "#c6dbef",
         "#deebf7",
         "mediumblue",
         "#ff7f00", 
         "#33a02c", 
         "#b2df8a", 
         "#20b2aa", 
         "#c51b7d") 



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=10, fig.width=25}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
top <- topTable[abs(topTable$logFC) >= 1.2 & topTable$P.Value <= 10e-03,]
DotPlot(seurat.htap.ctrl.female,features=rownames(head(top[order(top$AveExpr, decreasing = T),],30)), idents = c('0', '1', '2', '3', '4'), dot.scale = 8, group.by = "patient_htap.ctrl", cols = c("blue", "red")) + RotatedAxis() + labs(title = "HTAP vs CTRL in female dataset (Linear Model by limma, robust = true, FCcutoff = 1.5, Pvalue ≤ 0.05)", x = "Markers", y = "Clusters")

VlnPlot(seurat.htap.ctrl.female, features = c("ITGA4", "PSG4", "CRYBB2"), pt.size = 0, group.by = "patient_htap.ctrl")

DefaultAssay(seurat.htap.ctrl.male) <- "RNA"
seurat.htap.ctrl <- NormalizeData(seurat.htap.ctrl.male)
top <- topTable[abs(topTable$logFC) >= 1.2 & topTable$P.Value <= 10e-03,]
DotPlot(seurat.htap.ctrl.male,features=rownames(head(top[order(top$AveExpr, decreasing = T),],30)), idents = c('0', '1', '2', '3', '4'), dot.scale = 8, group.by = "patient_htap.ctrl", cols = c("blue", "red")) + RotatedAxis() + labs(title = "HTAP vs CTRL in male dataset (Linear Model by limma, robust = true, FCcutoff = 1.5, Pvalue ≤ 0.05)", x = "Markers", y = "Clusters")

```





## GSEA using gene sets from KEGG pathways
set.seed(123456)
gseaKEGG <- gseKEGG(geneList = genes.xx, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 100, # default number permutations
                    #minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    #pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
View(gseaKEGG_results)

## Output results from GO analysis to a table
enrichGO_summary_htap.vs.ctrl.female <- data.frame(ego)
#View(enrichGO_summary_htap.vs.ctrl.female)
write.table(enrichGO_summary_htap.vs.ctrl.female, file="./output/enrichGO_summary_htap.vs.ctrl.female.csv", sep=";")

# Visualization

#barplot(ego, showCategory=50) + ggtitle("GO enrichment chart (HTAP vs CTRL)")
dotplot(ego, showCategory=50) + ggtitle("GO enrichment chart (HTAP vs CTRL)")


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 50, pie_scale=0.8,layout="kk", line_scale = 0.1) + ggtitle("Enrichmap gathering the 50 most significant GO terms (by padj) \n to visualize the relations between the terms (HTAP vs CTRL)")


edo <- enrichDGN(names(geneList))

# Gene-Concept Network
## convert gene ID to Symbol

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
head(edox)
# Heatmap-like functional classification
#heatplot(edox, foldChange=geneList)


## Cnetplot details the genes associated with one or more terms 

cnetplot(ego,
         categorySize="pvalue",
         showCategory = 50,
         foldChange=geneList,
         vertex.label.font=6)

# Generate the plotting object

#cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
#cnetplot(edox, node_label="all") 

# Heatmap-like functional classification
#heatplot(edox, foldChange=geneList)











# Goplot

# david <- enrichGO_summary_htap.vs.ctrl.female[, c( "ONTOLOGY","ID","Description","geneID" , "p.adjust")]
# colnames(david) <- c("Category", "ID","Term", "Genes", "adj_pval")
# rownames(david) <- NULL
# # Generate the plotting object
# topTable.xx <- topTable.x[, c("gene", "logFC" ,"AveExpr","t" ,"P.Value","adj.P.Val" ,"B")]
# colnames(topTable.xx) <- c( "ID","logFC", "AveExpr","t",  "P.Value", "adj.P.Val","B")
# rownames(topTable.xx) <- NULL
# circ <- circle_dat(david, topTable.xx)
# #View(circ)
# #GOBar(subset(circ, category == 'BP'))
# 
# # Facet the barplot according to the categories of the terms 
# 
# GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))
# 
# #Generate a circular visualization of the results of gene- annotation enrichment analysis
# GOCircle(circ)
# 
# 
# #chord <- chord_dat(circ, topTable.xx[, c("ID", "logFC")], david$Term)
# 
# 
# # Generate the bubble plot with a label threshold of 3
# #GOBubble(circ)




View(gseaKEGG_results)
## Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa03040')

## Write GSEA results to file
View(gseaKEGG_results)
write.csv(gseaKEGG_results, file="./output/gseaKEGG_results.csv", sep=";", quote=F)

## Output images for a single significant KEGG pathway
# pathview(gene.data = geneList,
#               pathway.id = "hsa04657",
#               species = "hsa",
#               limit = list(gene = 2, # value gives the max/min limit for foldchanges
#               cpd = 1))


## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
  pathview(gene.data = geneList, pathway.id = gseaKEGG_results$ID[x], species = "hsa", 
           limit = list(gene = 2, cpd = 1))
}

purrr::map(1:length(gseaKEGG_results$ID), get_kegg_plots)


# GSEA using gene sets associated with ALL Gene Ontology terms
gseaGO <- gseGO(geneList = genes.xx, 
                OrgDb = org.Hs.eg.db, 
                ont = 'ALL', 
                nPerm = 1000, 
                #minGSSize = 20, 
                pvalueCutoff = 0.05,
                verbose = FALSE) 

View(gseaGO@result)
write.table(gseaGO_results_ADA2_NOMO1, file="./output/gseaGO_results_ADA2_NOMO1.csv", sep=";")

gseaplot(gseaGO, geneSetID = 'GO:0030246')

# Group GO
ggo <- groupGO(gene     = names(geneList),
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

View(data.frame(ggo))





# Marking differences between clusters 4 and 0

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=3.5, fig.width=16}
markers_4_0 <- FindMarkers(object = seurat.htap.ctrl, ident.1=4, ident.2=0, only.pos = FALSE)
markers_4 <- head(markers_4_0,50)
DotPlot(seurat.htap.ctrl, features=rownames(head(markers_4[order(markers_4$avg_logFC),],50)), idents = c('4','0'),  group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis()  + labs(title = "Specific markers in cluster 4 (HTAP + CTRL)", x = "Markers", y = "Clusters")


markers_0_4 <- FindMarkers(object = seurat.htap.ctrl, ident.1=0, ident.2=4, only.pos = FALSE)
markers_0 <- head(markers_0_4,50)
DotPlot(seurat.htap.ctrl, features=rownames(head(markers_0[order(markers_0$avg_logFC),],50)), idents = c('4','0'),  group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis()  + labs(title = "Specific markers in cluster 0 (HTAP + CTRL)", x = "Markers", y = "Clusters")
```



# Marking differences between clusters 2 and 0 and 1

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=16}
markers_2_0 <- FindMarkers(object = seurat.htap.ctrl, ident.1=2, ident.2=0, only.pos = FALSE)
markers_2.0 <- head(markers_2_0,50)
DotPlot(seurat.htap.ctrl, features=rownames(head(markers_2.0[order(markers_2.0$avg_logFC),],50)), idents = c('2','0'),  group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis()  + labs(title = "Cluster 2 vs 0 (HTAP + CTRL)", x = "Markers", y = "Clusters")


markers_2_1 <- FindMarkers(object = seurat.htap.ctrl, ident.1=2, ident.2=1, only.pos = FALSE)
markers_2.1 <- head(markers_2_1,50)
DotPlot(seurat.htap.ctrl, features=rownames(head(markers_2.1[order(markers_2.1$avg_logFC),],50)), idents = c('2','1'),  group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis()  + labs(title = "Cluster 2 vs 1 (HTAP + CTRL)", x = "Markers", y = "Clusters")
```






# Marking differences between clusters 0 and 1

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=3.8, fig.width=16}

markers_0_1 <- FindMarkers(object = seurat.htap.ctrl.female, ident.1=0, ident.2=1, only.pos = TRUE)
markers_0 <- head(markers_0_1,50)
DotPlot(seurat.htap.ctrl.female,features=rownames(head(markers_0[order(markers_0$avg_logFC),],50)), idents = c('0','1'), group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis() + labs(title = "Specific markers in cluster 0 (HTAP + CTRL)", x = "Markers", y = "Clusters")

```

# Marking differences between clusters 1 and 3

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=16}

markers_1.vs.0.3 <- FindMarkers(object = seurat.htap.ctrl, ident.1=1, ident.2=c(0,3), only.pos = FALSE)
markers_1.vs <- head(markers_1.vs.0.3,50)
DotPlot(seurat.htap.ctrl, features=rownames(head(markers_1.vs[order(markers_1.vs$avg_logFC),],50)), idents = c('1','0', '3'), group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis()  + labs(title = "Specific markers in cluster 1 (HTAP + CTRL)", x = "Markers", y = "Clusters")




markers_1_3 <- FindMarkers(object = seurat.htap.ctrl, ident.1=1, ident.2=3, only.pos = FALSE)
markers_1 <- head(markers_1_3,50)
DotPlot(seurat.htap.ctrl, features=rownames(head(markers_1[order(markers_1$avg_logFC),],50)), idents = c('1','3'), group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis()  + labs(title = "Specific markers in cluster 1 (HTAP + CTRL)", x = "Markers", y = "Clusters")

```




# # Select the samples by sex
# DefaultAssay(seurat.htap.ctrl) <- "RNA"
# seurat.htap.ctrl <- NormalizeData(seurat.htap.ctrl)
# VlnPlot(seurat.htap.ctrl, features = c("XIST", "DDX3Y",  "EIF1AY", "RPS4Y1"), group.by = "PATIENT", pt.size = 0, ncol = 3)

#DimPlot(object = seurat.htap.ctrl.male, reduction = "umap", label = TRUE) + labs(title = "Male", colour = "Cell clusters")











































@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  
  
  
  
  # Clustering, Frequency 
  
  ```{r}
DefaultAssay(seurat.htap.ctrl) <- "integrated"
seurat.htap.ctrl <- FindNeighbors(seurat.htap.ctrl, reduction = "mnn", dims = 1:20)
seurat.htap.ctrl <- FindClusters(seurat.htap.ctrl, resolution = 0.15)
seurat.htap.ctrl[['seurat_clusters']][which(seurat.htap.ctrl@active.ident == 4),] <- 0
seurat.htap.ctrl <- SetIdent(seurat.htap.ctrl, value = "seurat_clusters")
seurat.htap.ctrl <- RenameIdents(seurat.htap.ctrl, `0` = "0", `1` = "1", `2` = "Cell cycling", `3` = "3")
seurat.htap.ctrl$pheno <- paste0(seurat.htap.ctrl$seurat_clusters, "_", seurat.htap.ctrl$PHENOTYPE)
seurat.htap.ctrl$pat <- paste0(seurat.htap.ctrl$seurat_clusters, "_", seurat.htap.ctrl$PATIENT)

seurat.htap.ctrl[["sex"]]=""
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C686518"),] <- "M"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C681725"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C681086"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C685116"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C676136"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H247"),] <- "M"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H248"),] <- "M"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H230"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H231"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H250"),] <- "F"
#seurat.htap.ctrl.male <- subset(seurat.htap.ctrl, sex == "M") 
seurat.htap.ctrl.female <- subset(seurat.htap.ctrl, sex == "F") 

DimPlot(object = seurat.htap.ctrl.female, label = T, label.size = 5)+labs(color = "Cell clusters")
DimPlot(object = seurat.htap.ctrl.female, split.by = "PHENOTYPE", label = F, label.size = 5)+labs(color = "Cell clusters")

# # Define cluster colors
# # Cluster 0
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
seurat.htap.ctrl.female[['cluster0']] <- " "
seurat.htap.ctrl.female[['cluster0']][which(seurat.htap.ctrl.female@active.ident == 0),] <- 0
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster0")
DimPlot(seurat.htap.ctrl.female, group.by = "cluster0", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)

seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
seurat.htap.ctrl.female[['cluster1']] <- " "
seurat.htap.ctrl.female[['cluster1']][which(seurat.htap.ctrl.female@active.ident == 1),] <- 1
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster1")
DimPlot(seurat.htap.ctrl.female, group.by = "cluster1", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)

seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
seurat.htap.ctrl.female[['cluster2']] <- " "
seurat.htap.ctrl.female[['cluster2']][which(seurat.htap.ctrl.female@active.ident == 2),] <- 2
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster2")
DimPlot(seurat.htap.ctrl.female, group.by = "cluster2", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)


seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
seurat.htap.ctrl.female[['cluster3']] <- " "
seurat.htap.ctrl.female[['cluster3']][which(seurat.htap.ctrl.female@active.ident == 3),] <- 3
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster3")
DimPlot(seurat.htap.ctrl.female, group.by = "cluster3", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)


# seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
# seurat.htap.ctrl.female[['cluster4']] <- " "
# seurat.htap.ctrl.female[['cluster4']][which(seurat.htap.ctrl.female@active.ident == 4),] <- 4
# seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster3")
# DimPlot(seurat.htap.ctrl.female, group.by = "cluster4", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)

```







```{r}
VlnPlot(seurat.htap.ctrl.female, features = c("XIST", "DDX3Y",  "EIF1AY", "RPS4Y1"), group.by = "PATIENT", pt.size = 0, ncol = 2)

VlnPlot(seurat.htap.ctrl.male, features = c("XIST", "DDX3Y",  "EIF1AY", "RPS4Y1"), group.by = "PATIENT", pt.size = 0, ncol = 2)

FeaturePlot(object = seurat.htap.ctrl.female, features = c("XIST", "DDX3Y",  "EIF1AY", "RPS4Y1"), ncol=2) 

FeaturePlot(object = seurat.htap.ctrl.male, features = c("XIST", "DDX3Y",  "EIF1AY", "RPS4Y1"), ncol=2)


# Visualize QC metrics as a violin plot 
VlnPlot(seurat.htap.ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size = 0, ncol = 3)

VlnPlot(seurat.htap.ctrl.female, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), split.by = "PHENOTYPE", cols = c("#33a02c", "#e31a1c"), pt.size = 0, ncol = 3)

FeaturePlot(object = seurat.htap.ctrl, features = c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol=3) 

DimPlot(object = seurat.htap.ctrl, label = FALSE, group.by = "PHENOTYPE", cols = c("#33a02c", "#e31a1c")) + labs(colour = "Phenotype")

DimPlot(object = seurat.htap.ctrl, label = TRUE) + labs(colour = "Cell clusters")

DimPlot(object = seurat.htap.ctrl, label = FALSE, split.by = "PHENOTYPE") + labs(colour = "Cell clusters")
DimPlot(object = seurat.htap.ctrl, label = FALSE, group.by = "PATIENT") + labs(colour = "Patients")



seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
freq_table <- prop.table(x = table(seurat.htap.ctrl.female@active.ident, seurat.htap.ctrl.female@meta.data[, "PHENOTYPE"]), margin = 2)*100
freq_table <- round(freq_table, digits = 2)

barplot(t(freq_table),
        las=1,
        horiz = F,
        ylim = c(0, max(freq_table)*1.2),
        #col = "#619CFF",
        ylab = "Cell numbers (%)",
        xlab = "Cell clusters",
        main = "Number of cells by cluster",
        legend = T,
        args.legend=list(title="Patients"),
        col=c(rep("#33a02c", 1), rep("#e31a1c",1)),
        beside = T)
#B <- as.matrix(freq_table)
# text(B+2, A, labels = paste(as.character(B), "%", sep = ""), cex = 0.8, pos = 2)


freq_table <- prop.table(x = table(seurat.htap.ctrl.female@active.ident, seurat.htap.ctrl.female@meta.data[, "PATIENT"]), margin = 2)*100
freq_table <- round(freq_table, digits = 2)

barplot(t(freq_table),
        las=1,
        horiz = F,
        ylim = c(0, max(freq_table)*1.2),
        #col = "#619CFF",
        ylab = "Cell numbers (%)",
        xlab = "Cell clusters",
        main = "Number of cells by patient",
        legend = T,
        args.legend=list(title="Patients"),
        col=c(rep("#33a02c", 4), rep("#e31a1c",3)),
        beside = T)
#B <- as.matrix(freq_table)
# text(B+2, A, labels = paste(as.character(B), "%", sep = ""), cex = 0.8, pos = 2)




# freq_table <- prop.table(x = table(seurat.htap.ctrl.male@active.ident, seurat.htap.ctrl.male@meta.data[, "PATIENT"]), margin = 2)*100
# freq_table <- round(freq_table, digits = 2)
# 
# barplot(t(freq_table),
#         las=1,
#         horiz = F,
#         ylim = c(0, max(freq_table)*1.2),
#         #col = "#619CFF",
#         ylab = "Cell numbers (%)",
#         xlab = "Cell clusters",
#         main = "Number of cells by cluster in male dataset",
#         legend = T,
#         args.legend=list(title="Patients"),
#         col=c("#33a02c", rep("#e31a1c",2)),
#         beside = T)

```


# Clusters gene markers heatmap

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=10, fig.width=10}
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
seurat.htap.ctrl.female <- ScaleData(seurat.htap.ctrl.female)
markers <- FindAllMarkers(object = seurat.htap.ctrl.female, only.pos = TRUE)
top20 <- markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

DoHeatmap(seurat.htap.ctrl.female, features=top20$gene, size=5.5, label = F, angle = 0) + labs(title = " CTRL + HTAP dataset - Expression Heatmap", colour = "Cell clusters")

write.table(markers, file="./output/seurat.htap.ctrl.female.csv", sep="\t")

```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=14, fig.width=8}
DefaultAssay(seurat.htap.ctrl.male) <- "RNA"
seurat.htap.ctrl.male <- NormalizeData(seurat.htap.ctrl.male)
seurat.htap.ctrl.male <- ScaleData(seurat.htap.ctrl.male)
markers <- FindAllMarkers(object = seurat.htap.ctrl.male, only.pos = TRUE)
top20 <- markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

DoHeatmap(seurat.htap.ctrl.male, features=top20$gene, size=5.5, label = T, angle = 0) + labs(title = "Male (CTRL + HTAP) ", colour = "Cell clusters")

write.table(markers, file="./output/seurat.htap.ctrl.male.csv", sep="\t")


male <-  markers
DoHeatmap(seurat.htap.ctrl.male, features=top20$gene, size=5.5, label = T, angle = 0) + labs(title = "Male (CTRL + HTAP) ", colour = "Cell clusters")


```




# Clusters correlation heatmap

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}

mnmat <- c()
uniq <- unique(seurat.htap.ctrl.female@active.ident)
seurat.htap.ctrl.female@meta.data$cluster <- seurat.htap.ctrl.female@active.ident
for(i in 1:length(uniq)){
  mnmat <- cbind(mnmat, apply(as.matrix(seurat.htap.ctrl.female@assays$RNA@data[, seurat.htap.ctrl.female@meta.data$cluster==uniq[i]]), 1, mean))
}

colnames(mnmat) <- as.vector(unique(seurat.htap.ctrl.female@active.ident))
ct=cor(mnmat)
pheatmap(ct, angle_col = 0, main = "Clusters correlation heatmap")


mnmat <- c()
uniq <- unique(seurat.htap.ctrl.male@active.ident)
seurat.htap.ctrl.male@meta.data$cluster <- seurat.htap.ctrl.male@active.ident
for(i in 1:length(uniq)){
  mnmat <- cbind(mnmat, apply(as.matrix(seurat.htap.ctrl.male@assays$RNA@data[, seurat.htap.ctrl.male@meta.data$cluster==uniq[i]]), 1, mean))
}

colnames(mnmat) <- as.vector(unique(seurat.htap.ctrl.male@active.ident))
ct=cor(mnmat)
pheatmap(ct, angle_col = 0, main = "Clusters correlation heatmap")
```



# FeaturePlot

```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=10, fig.width=12}
FeaturePlot(seurat.htap.ctrl.female, features = "IFI27", label = T)

DimPlot(object = seurat.htap.ctrl, label = TRUE) + labs(colour = "Cell clusters")
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=20, fig.width=25}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
FeaturePlot(seurat.htap.ctrl.female, features = c("LYVE1","PROX1", "PDPN", "NRP2", "FLT4", "SOX18", "IGFBP7"), split.by = "SAMPLE")


```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=40, fig.width=20}

christoph.markers <- c("LYVE1","VWF","PROX1", "PDPN", "NRP2", "EPHB2", "EPHB4", "ACKR3",  "BGN", "FLT4", "NRP1", "NR2F2",  "RGCC", "EPAS1", "SOX18", "NOTCH4", "KDR", "EDN1", "EFNB2", "ID1", "BMX", "HEY1", "IGFBP7")

# Female

DimPlot(object = seurat.htap.ctrl, label = TRUE) + labs(colour = "Cell clusters")
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
#FeaturePlot(seurat.htap.ctrl.female, features = christoph.markers)
FeaturePlot(seurat.htap.ctrl.female, features = christoph.markers, split.by = "SAMPLE")


# Male
DefaultAssay(seurat.htap.ctrl.male) <- "integrated"
seurat.htap.ctrl <- FindNeighbors(seurat.htap.ctrl.male, reduction = "mnn", dims = 1:20)
seurat.htap.ctrl <- FindClusters(seurat.htap.ctrl.male, resolution = 0.8)
DimPlot(object = seurat.htap.ctrl, label = TRUE) + labs(colour = "Cell clusters")
DefaultAssay(seurat.htap.ctrl.male) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.male)
# FeaturePlot(seurat.htap.ctrl.male, features = christoph.markers, split.by = "SAMPLE")
FeaturePlot(seurat.htap.ctrl.male, features = christoph.markers)
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=20, fig.width=10}
DefaultAssay(seurat.htap.ctrl.male) <- "RNA"
seurat.htap.ctrl.male <- NormalizeData(seurat.htap.ctrl.male)

FeaturePlot(seurat.htap.ctrl.male, features = c("LYVE1","PROX1", "PDPN", "NRP2", "FLT4", "SOX18", "IGFBP7"), split.by = "SAMPLE")
```





```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}
# Enrichment analysis
markers_0_1$SYMBOL <- rownames(markers_0_1)
gene.df <- bitr(rownames(markers_0_1), fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

d <- merge(gene.df, markers_0_1, by="SYMBOL")

geneList <- d[,"avg_logFC"]
names(geneList) <- as.character(d[,"ENTREZID"])
geneList <- sort(geneList, decreasing = TRUE)

edo <- enrichDGN(names(geneList))
barplot(edo, showCategory=50) + ggtitle("Disease terms founded in HTAP vs CTRL DGE analysis")
#dotplot(edo, showCategory=50) + ggtitle("Disease terms founded in HTAP vs CTRL DGE analysis")

# Gene-Concept Network
## convert gene ID to Symbol

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

cnetplot(edox, node_label="all") 

# Heatmap-like functional classification
heatplot(edox, foldChange=geneList)

emapplot(edo, pie_scale=0.8,layout="kk", line_scale = 0.1) 



# GSEA analysis (Gene Set Enrichment Analysis GSEA)

## Create ranks
#barplot(sort(geneList, decreasing = T))


# KEGG enrichment analysis
sigGenes <- d$ENTREZID
sigGenes <- na.exclude(sigGenes)
kk <- enrichKEGG(gene = sigGenes, organism = 'hsa')
head(kk, n=10)

#browseKEGG(kk, 'hsa04512')


pathview(gene.data = geneList, 
         pathway.id = "hsa04151", 
         species = "hsa", 
         limit = list(gene=5, cpd=1))




# markers_1_0 <- FindMarkers(object = seurat.htap.ctrl.female, ident.1=1, ident.2=0, only.pos = FALSE)
# markers_1 <- head(markers_1_0,50)
# DotPlot(seurat.htap.ctrl.female,features=rownames(head(markers_1[order(markers_1$avg_logFC),],50)), idents = c('0','1'), group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis() + labs(title = "Specific markers in cluster 1 (HTAP + CTRL)", x = "Markers", y = "Clusters")
# 
# 
# 
# markers_0_1 <- FindMarkers(object = seurat.htap.ctrl.male, ident.1=0, ident.2=1, only.pos = FALSE)
# markers_0 <- head(markers_0_1,50)
# DotPlot(seurat.htap.ctrl.male,features=rownames(head(markers_0[order(markers_0$avg_logFC),],50)), idents = c('0','1'), group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis() + labs(title = "Specific markers in cluster 0 (HTAP + CTRL)", x = "Markers", y = "Clusters")
# 
# 
# markers_1_0 <- FindMarkers(object = seurat.htap.ctrl.male, ident.1=1, ident.2=0, only.pos = FALSE)
# markers_1 <- head(markers_1_0,50)
# DotPlot(seurat.htap.ctrl.male,features=rownames(head(markers_1[order(markers_1$avg_logFC),],50)), idents = c('0','1'), group.by = "cluster_htap.ctrl", cols = c("blue", "red")) + RotatedAxis() + labs(title = "Specific markers in cluster 1 (HTAP + CTRL)", x = "Markers", y = "Clusters")

```







```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=10, fig.width=12}

markers_0.vs.1.4$SYMBOL <- rownames(markers_0.vs.1.4)

gene.df <- bitr(markers_0.vs.1.4$SYMBOL, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

d <- merge(gene.df, markers_0.vs.1.4, by="SYMBOL")

geneList <- d$avg_logFC
names(geneList) <- as.character(d$ENTREZID)
geneList <- sort(geneList, decreasing = TRUE)

edo <- enrichDGN(names(geneList))
barplot(edo, showCategory=50) + ggtitle("Disease terms founded in HTAP vs CTRL DGE analysis")
dotplot(edo, showCategory=50) + ggtitle("Disease terms founded in HTAP vs CTRL DGE analysis")

# Gene-Concept Network
## convert gene ID to Symbol

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

cnetplot(edox, node_label="all") 

# Heatmap-like functional classification
heatplot(edox, foldChange=geneList)

emapplot(edo, pie_scale=0.8,layout="kk", line_scale = 0.1) 
```








## Marking differences between clusters 0 and 1

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=15}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
markers_0.vs.1 <- FindMarkers(object = seurat.htap.ctrl.female, ident.1 = 0, ident.2 = 1, only.pos = TRUE)

DotPlot(seurat.htap.ctrl.female, features = rownames(head(markers_0.vs.1[order(markers_0.vs.1$avg_logFC, decreasing = TRUE),],50)), idents = c('0', '1'), group.by = "pheno") + RotatedAxis()  + labs(title = "Top 50 specific markers in cluster 0 vs 1 (HTAP + CTRL)", x = "Markers", y = "Clusters")

# Data
gprofiler_0 <- gprofiler2::gconvert(query = rownames(markers_0.vs.1), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
markers_0.vs.1$SYMBOL <- rownames(markers_0.vs.1)
genes.c0 <- bitr(gprofiler_0$target, fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)
## Remove any Entrez duplicates
genes.c0 <- genes.c0[which(duplicated(genes.c0$ENTREZID) == F), ]
genes.c0.List <- markers_0.vs.1$avg_logFC
names(genes.c0.List) <- as.character(genes.c0$ENTREZID)
## Remove any NA values
genes.c0.List <- genes.c0.List[!is.na(names(genes.c0.List))]
#genes.c0.List <- genes.c0.List[abs(genes.c0.List) >= 1]
genes.c0.List <- sort(genes.c0.List, decreasing = TRUE)
genes.c0.List
```






#### Enrich GO terms


```{r}
#data(geneList, package="DOSE")
enrichGO0 <- enrichGO(gene = names(genes.c0.List),
                      #universe      = names(geneList),
                      keyType = "ENTREZID",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE,
                      pool = TRUE)
enrichGO_0 <- setReadable(enrichGO0, 'org.Hs.eg.db')


#View(data.frame(enrichGO0))


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN0 <- enrichDGN(names(genes.c0.List))
enrichDGN_0 <- setReadable(enrichDGN0, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO0 <- enrichDO(gene          = names(genes.c0.List),
                      ont           = "DO",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      #universe      = names(genes.c0.List),
                      minGSSize     = 5,
                      maxGSSize     = 500,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)
enrichDO_0 <- setReadable(enrichDO0, 'org.Hs.eg.db')

```


#### Barplot

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=12, fig.width=14}

barplot(enrichGO0, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 0 (HTAP + CTRL)", x = "GO terms", y = "Genes related to GO term ")

```

#### Emaplot
```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=12, fig.width=14}
## Plotting terms of interest
cnetplot(enrichGO_0, 
         categorySize="pvalue", 
         foldChange=genes.c0.List, 
         showCategory = 10, 
         vertex.label.font=6,
         colorEdge = TRUE,
         circular = FALSE)

```


#### Heatmap plot

```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=18}
# Heatmap-like functional classification

heatplot(enrichDO_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015) in cluster 0 (HTAP + CTRL)", x = "Genes", y = "Disease Ontology")


heatplot(enrichDGN_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Gene-disease association (DisGeNET, Janet et al. 2015)   in cluster 0 (HTAP + CTRL)", y = "Diseases", x = "Genes")
```

```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=15}
# Heatmap-like functional classification
heatplot(enrichGO_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Top 50 enriched Gene Ontology terms in cluster 0 (HTAP + CTRL)", x = "Genes", y = "GO terms")

```

#### GSEA using gene sets from KEGG pathways

```{r}
c0 <- genes.c0.List[genes.c0.List >= 1]
## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = c0, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    #nPerm = 1000, # default number permutations
                    minGSSize = 5, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    #pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
```


#### UMAP

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=8}

DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)

FeaturePlot(seurat.htap.ctrl.female, features = c("ALDH1A1", "FABP4", "TFPI", "BGN", "TIMP3", "MMRN1"), label = F)

DotPlot(seurat.htap.ctrl.female, features = c("ALDH1A1", "FABP4", "TFPI", "BGN", "TIMP3", "MMRN1"), group.by = "pat") + RotatedAxis() + labs(title = "Specific markers in cluster 0", x = "Gene", y = "Patients cluster")


```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=15, fig.width=25}
VlnPlot(seurat.htap.ctrl.female, features =  c("ALDH1A1", "FABP4", "TFPI", "BGN", "TIMP3", "MMRN1"), group.by = "pat", pt.size = 0, ncol = 2)
```





#### KEGG over-representation test

```{r}
enrichKEGG_0 <- enrichKEGG(gene         = names(genes.c0.List),
                           organism     = 'hsa',
                           pvalueCutoff = 0.05)
enrichKEGG_0 <- data.frame(enrichKEGG_0)
enrichKEGG_0$cluster <- rep("cluster0", nrow(enrichKEGG_0))
View(data.frame(enrichKEGG_0))

## Output images for a single significant KEGG pathway
pathview(gene.data = genes.c0.List,
         pathway.id = "hsa05165",
         species = "hsa",
         limit = list(gene = max(abs(genes.c0.List)), # value gives the max/min limit for foldchanges
                      cpd = 1))


```



## Marking differences between clusters 2 and 1

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=15}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)

markers_2.vs.1 <- FindMarkers(seurat.htap.ctrl.female, ident.1 = 2, ident.2 = 1,  only.pos = TRUE)

DotPlot(seurat.htap.ctrl.female, features=rownames(head(markers_2.vs.1[order(markers_2.vs.1$avg_logFC, decreasing = TRUE),],50)), idents = c('1', '2'), group.by = "pheno") + RotatedAxis()  + labs(title = "Top 50 specific markers in cluster 2 vs 1 (HTAP + CTRL)", x = "Markers", y = "Clusters")


# Data
gprofiler_2 <- gprofiler2::gconvert(query = rownames(markers_2.vs.1), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
gprofiler_2 <- gprofiler_2[which(duplicated(gprofiler_2$target) == F), ]
markers_2.vs.1$SYMBOL <- rownames(markers_2.vs.1)
genes.c2 <- bitr(gprofiler_2$target, fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)

## Remove any Entrez duplicates
genes.c2 <- genes.c2[which(duplicated(genes.c2$SYMBOL) == F), ]

genes.c2.List <- markers_2.vs.1$avg_logFC
names(genes.c2.List) <- as.character(genes.c2$ENTREZID)
## Remove any NA values
genes.c2.List <- genes.c2.List[!is.na(names(genes.c2.List))]
#genes.c2.List <- genes.c2.List[genes.c2.List >= 0.5]
genes.c2.List <- sort(genes.c2.List, decreasing = TRUE)
genes.c2.List
```

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=50, fig.width=20}
FeaturePlot(seurat.htap.ctrl.female, features = rownames(head(markers_2.vs.1[order(markers_2.vs.1$avg_logFC, decreasing = TRUE),],50)))
```


#### Enrich GO terms


```{r}
#data(geneList, package="DOSE")
enrichGO2 <- enrichGO(gene = names(genes.c2.List),
                      #universe      = names(geneList),
                      keyType = "ENTREZID",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE,
                      pool = TRUE)
enrichGO_2 <- setReadable(enrichGO2, 'org.Hs.eg.db')
#View(data.frame(enrichGO_1))

# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN2 <- enrichDGN(names(genes.c2.List))
enrichDGN_2 <- setReadable(enrichDGN2, 'org.Hs.eg.db')

# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO2 <- enrichDO(gene          = names(genes.c2.List),
                      ont           = "DO",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      #universe      = names(genes.c0.List),
                      minGSSize     = 5,
                      maxGSSize     = 500,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)
enrichDO_2 <- setReadable(enrichDO2, 'org.Hs.eg.db')


```


#### Barplot

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=12, fig.width=14}

barplot(enrichGO2, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 1 (HTAP + CTRL)", x = "GO terms", y = "Genes related to GO term ")

```



#### Heatmap plot

```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=18}
# Heatmap-like functional classification
heatplot(enrichGO_2, foldChange=genes.c2.List, showCategory=50) + labs(title = "Top 50 enriched Gene Ontology terms \n in cluster 2 (HTAP + CTRL)", x = "Genes", y = "GO terms")


heatplot(enrichDO_2, foldChange=genes.c2.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015) in cluster 2 (HTAP + CTRL)", x = "Genes", y = "Disease Ontology")


heatplot(enrichDGN_2, foldChange=genes.c2.List, showCategory=50) + labs(title = "Gene-disease association (DisGeNET, Janet et al. 2015) \n  in cluster 2 (HTAP + CTRL)", y = "Diseases", x = "Genes")
```

#### UMAP

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=8}

DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)

FeaturePlot(seurat.htap.ctrl.female, features = c("IFI27", "RGS5", "BST2", "FN1", "MMP1", "SLCTA11", "VWF"), label = F)

DotPlot(seurat.htap.ctrl.female, features = c("IFI27", "RGS5", "BST2", "FN1", "MMP1", "SLCTA11", "VWF"), group.by = "pat") + RotatedAxis() + labs(title = "Specific markers in cluster 2", x = "Gene", y = "Patients cluster")


```




#### KEGG over-representation test

```{r}

enrichKEGG_2 <- enrichKEGG(gene         = names(genes.c2.List),
                           organism     = 'hsa',
                           pvalueCutoff = 0.05)
enrichKEGG_2 <- data.frame(enrichKEGG_2)
enrichKEGG_2$cluster <- rep("cluster2", nrow(enrichKEGG_2))
View(data.frame(enrichKEGG_2))

## Output images for a single significant KEGG pathway
pathview(gene.data = genes.c1.List,
         pathway.id = "hsa05205",
         species = "hsa",
         limit = list(gene = max(abs(genes.c2.List)), # value gives the max/min limit for foldchanges
                      cpd = 1))


## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
  pathview(gene.data = genes.c1.List, pathway.id = enrichKEGG_1$ID[x], species = "hsa", 
           limit = list(gene = max(abs(genes.c1.List)), cpd = 1))
}

purrr::map(1:length(enrichKEGG_1$ID), get_kegg_plots)

```


## Marking differences between clusters 1 vs 0 and 2

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=15}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")

markers_1.vs.0 <- FindMarkers(seurat.htap.ctrl.female, ident.1 = 1, ident.2 = 0,   min.pct = 0.25, only.pos = F)

markers_1.vs.2 <- FindMarkers(seurat.htap.ctrl.female, ident.1 = 1, ident.2 = 2,   min.pct = 0.25, only.pos = F)

rbind_1 <- rbind(markers_1.vs.0, markers_1.vs.2)

markers.c1.0.2 <- intersect(rownames(markers_1.vs.0), rownames(markers_1.vs.2))

DotPlot(seurat.htap.ctrl.female, features = markers.c1.0.2, idents = c('0', '1', '2')) + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 0 and 2 (HTAP + CTRL)", x = "Markers", y = "Clusters")

markers.c1 <- c("SSR3", "SLC12A2", "FDFT1", "SQLE")

markers_1 <- rbind_1[markers.c1, ]



# genedot1=markers_1.vs.0.2 %>% top_n(-10, avg_logFC)
# genedot2=markers_1.vs.0.2 %>% top_n(10, avg_logFC)
# genedot=c(rownames(genedot1),rownames(genedot2))
# DotPlot(seurat.htap.ctrl.female, features = genedot, idents = c('2', '1', '0')) + RotatedAxis()


# Data
gprofiler_1 <- gprofiler2::gconvert(query = rownames(markers_1), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
gprofiler_1 <- gprofiler_1[which(duplicated(gprofiler_1$target) == F), ]
markers_1$SYMBOL <- rownames(markers_1)
genes.c1 <- bitr(gprofiler_1$target, fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)

## Remove any Entrez duplicates
genes.c1 <- genes.c1[which(duplicated(genes.c1$SYMBOL) == F), ]
genes.c1.List <- markers_1$avg_logFC
names(genes.c1.List) <- as.character(genes.c1$ENTREZID)
## Remove any NA values
genes.c1.List <- genes.c1.List[!is.na(names(genes.c1.List))]
#genes.c1.List <- genes.c1.List[abs(genes.c1.List) >= 1]
genes.c1.List <- sort(genes.c1.List, decreasing = TRUE)
genes.c1.List
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}
FeaturePlot(seurat.htap.ctrl.female, features = markers.c1)
#VlnPlot(seurat.htap.ctrl.female, features = markers.c1, pt.size = 0)
DotPlot(seurat.htap.ctrl.female, features = markers.c1, group.by = "pat") + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 0 and 2 (HTAP + CTRL)", x = "Markers", y = "Clusters")
```

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=10}

EnhancedVolcano(markers_1,
                lab = rownames(markers_1),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 1 vs 0 and 2", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()


EnhancedVolcano(markers_1.vs.0,
                lab = rownames(markers_1.vs.0),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 1 vs 0", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```


#### Enrich GO terms


```{r}
#data(geneList, package="DOSE")
enrichGO1 <- enrichGO(gene = names(genes.c1.List),
                      #universe      = names(geneList),
                      keyType = "ENTREZID",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE,
                      pool = TRUE)
enrichGO_1 <- setReadable(enrichGO1, 'org.Hs.eg.db')


#View(data.frame(enrichGO0))


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN1 <- enrichDGN(names(genes.c1.List))
enrichDGN_1 <- setReadable(enrichDGN1, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO1 <- enrichDO(gene          = names(genes.c1.List),
                      ont           = "DO",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      #universe      = names(genes.c0.List),
                      minGSSize     = 5,
                      maxGSSize     = 500,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)
enrichDO_1 <- setReadable(enrichDO1, 'org.Hs.eg.db')

```


#### Barplot

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=12, fig.width=20}

barplot(enrichGO1, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 1 vs 0 and 2 (HTAP + CTRL)", x = "GO terms", y = "Genes related to GO term ")

```


#### Heatmap plot

```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=15}
# Heatmap-like functional classification
heatplot(enrichGO_1, foldChange=genes.c1.List, showCategory=50) + labs(title = "Top 50 enriched Gene Ontology terms \n in cluster 1 vs 0 and 2 (HTAP + CTRL)", x = "Genes", y = "GO terms")

heatplot(enrichDO_1, foldChange=genes.c1.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015) in cluster 1 (HTAP + CTRL)", x = "Genes", y = "Disease Ontology")


```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=8}

heatplot(enrichDGN_1, foldChange=genes.c1.List, showCategory=50) + labs(title = "Gene-disease association  (DisGeNET, Janet et al. 2015) \n   in cluster 1 (HTAP + CTRL)", y = "Diseases", x = "Genes")
```


```{r}
enrichKEGG_1 <- enrichKEGG(gene         = names(genes.c1.List),
                           organism     = 'hsa',
                           pvalueCutoff = 0.05)
enrichKEGG_1 <- data.frame(enrichKEGG_1)
enrichKEGG_1$cluster <- rep("cluster2", nrow(enrichKEGG_1))
View(data.frame(enrichKEGG_1))

## Output images for a single significant KEGG pathway
pathview(gene.data = genes.c1.List,
         pathway.id = "hsa00100",
         species = "hsa",
         limit = list(gene = max(abs(genes.c1.List)), # value gives the max/min limit for foldchanges
                      cpd = 1))
```




#### Cluster 1 vs 0

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")

markers_1.vs.0 <- FindMarkers(seurat.htap.ctrl.female, ident.1 = 1, ident.2 = 0,   min.pct = 0.25, only.pos = T)

DotPlot(seurat.htap.ctrl.female, features=rownames(head(markers_1.vs.0[order(markers_1.vs.0$avg_logFC, decreasing = TRUE),],50)), idents = c('0', '1')) + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 0 (HTAP + CTRL)", x = "Markers", y = "Clusters")



gprofiler_1.0 <- gprofiler2::gconvert(query = rownames(markers_1.vs.0), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
gprofiler_1.0 <- gprofiler_1.0[which(duplicated(gprofiler_1.0$target) == F), ]
markers_1.vs.0$SYMBOL <- rownames(markers_1.vs.0)
genes.c1.0 <- bitr(gprofiler_1.0$target, fromType = "ENSEMBL",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

## Remove any Entrez duplicates
genes.c1.0 <- genes.c1.0[which(duplicated(genes.c1.0$SYMBOL) == F), ]
genes.c1.0.List <- markers_1.vs.0$avg_logFC
names(genes.c1.0.List) <- as.character(genes.c1.0$ENTREZID)
## Remove any NA values
genes.c1.0.List <- genes.c1.0.List[!is.na(names(genes.c1.0.List))]
genes.c1.0.List <- genes.c1.0.List[abs(genes.c1.0.List) >= 1]
genes.c1.0.List <- sort(genes.c1.0.List, decreasing = TRUE)
genes.c1.0.List



#data(geneList, package="DOSE")
enrichGO.1.0 <- enrichGO(gene = names(genes.c1.0.List),
                         #universe      = names(geneList),
                         keyType = "ENTREZID",
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE,
                         pool = TRUE)
enrichGO_1_0 <- setReadable(enrichGO.1.0, 'org.Hs.eg.db')


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN.1.0 <- enrichDGN(names(genes.c1.0.List))
enrichDGN_1_0 <- setReadable(enrichDGN.1.0, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO.1.0 <- enrichDO(gene          = names(genes.c1.0.List),
                         ont           = "DO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         minGSSize     = 5,
                         maxGSSize     = 500,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
enrichDO_1_0 <- setReadable(enrichDO.1.0, 'org.Hs.eg.db')




```




```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}

barplot(enrichGO.1.0, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 1 vs 0 (HTAP + CTRL)", x = "GO terms", y = "Genes related to GO term ")

# Heatmap-like functional classification
heatplot(enrichGO_1_0, foldChange=genes.c1.0.List, showCategory=50) + labs(title = "Enriched GO terms \n in cluster 1 vs 0 (HTAP + CTRL)", x = "Genes", y = "GO terms")

heatplot(enrichDO_1_0, foldChange=genes.c1.0.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015)\n in cluster 1 vs 0 (HTAP + CTRL)", x = "Genes", y = "Disease Ontology")


heatplot(enrichDGN_1_0, foldChange=genes.c1.0.List, showCategory=50) + labs(title = "Gene-disease association \n  (DisGeNET, Janet et al. 2015) \n   in cluster 1 vs 0 (HTAP + CTRL)", y = "Diseases", x = "Genes")
```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(markers_1.vs.0,
                lab = rownames(markers_1.vs.0),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 1 vs 0", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=8}
FeaturePlot(seurat.htap.ctrl.female, features = c("KRT7", "MMP1", "MT2A", "NQO1"), label = F)

DotPlot(seurat.htap.ctrl.female, features = c("KRT7", "MMP1", "MT2A", "NQO1"), group.by = "pat") + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 0", x = "Gene", y = "Patients cluster")


```



#### Cluster 1 vs 2

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=12}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")

markers_1.vs.2 <- FindMarkers(seurat.htap.ctrl.female, ident.1 = 1, ident.2 = 2,   min.pct = 0.25, only.pos = T)

DotPlot(seurat.htap.ctrl.female, features=rownames(head(markers_1.vs.2[order(markers_1.vs.2$avg_logFC, decreasing = TRUE),],50)), idents = c('0', '2')) + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 2 (HTAP + CTRL)", x = "Markers", y = "Clusters")



gprofiler_1.2 <- gprofiler2::gconvert(query = rownames(markers_1.vs.2), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
gprofiler_1.2 <- gprofiler_1.2[which(duplicated(gprofiler_1.2$target) == F), ]
markers_1.vs.2$SYMBOL <- rownames(markers_1.vs.2)
genes.c1.2 <- bitr(gprofiler_1.2$target, fromType = "ENSEMBL",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

## Remove any Entrez duplicates
genes.c1.2 <- genes.c1.2[which(duplicated(genes.c1.2$SYMBOL) == F), ]
genes.c1.2.List <- markers_1.vs.2$avg_logFC
names(genes.c1.2.List) <- as.character(genes.c1.2$ENTREZID)
## Remove any NA values
genes.c1.2.List <- genes.c1.2.List[!is.na(names(genes.c1.2.List))]
#genes.c1.2.List <- genes.c1.2.List[abs(genes.c1.2.List) >= 1]
genes.c1.2.List <- sort(genes.c1.2.List, decreasing = TRUE)
genes.c1.2.List



#data(geneList, package="DOSE")
enrichGO.1.2 <- enrichGO(gene = names(genes.c1.2.List),
                         #universe      = names(geneList),
                         keyType = "ENTREZID",
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE,
                         pool = TRUE)
enrichGO_1_2 <- setReadable(enrichGO.1.2, 'org.Hs.eg.db')


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN.1.2 <- enrichDGN(names(genes.c1.2.List))
enrichDGN_1_2 <- setReadable(enrichDGN.1.2, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO.1.2 <- enrichDO(gene          = names(genes.c1.2.List),
                         ont           = "DO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         minGSSize     = 5,
                         maxGSSize     = 500,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
enrichDO_1_2 <- setReadable(enrichDO.1.2, 'org.Hs.eg.db')




```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=15}

barplot(enrichGO.1.2, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 1 vs 2 (HTAP + CTRL)", x = "GO terms", y = "Genes related to GO term ")

# Heatmap-like functional classification
heatplot(enrichGO_1_2, foldChange=genes.c1.2.List, showCategory=50) + labs(title = "Enriched GO terms \n in cluster 1 vs 2  (HTAP + CTRL)", x = "Genes", y = "GO terms")

heatplot(enrichDO_1_2, foldChange=genes.c1.2.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015)\n in cluster 1 vs 2 (HTAP + CTRL)", x = "Genes", y = "Disease Ontology")


heatplot(enrichDGN_1_2, foldChange=genes.c1.2.List, showCategory=50) + labs(title = "Gene-disease association \n  (DisGeNET, Janet et al. 2015) \n   in cluster 1 vs 2 (HTAP + CTRL)", y = "Diseases", x = "Genes")
```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(markers_1.vs.2,
                lab = rownames(markers_1.vs.2),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 1 vs 2", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=8}
FeaturePlot(seurat.htap.ctrl.female, features = rownames(markers_1.vs.2[markers_1.vs.2$avg_logFC>=0.7,]), label = F)

FeaturePlot(seurat.htap.ctrl.female, features = c("CCND2", "PDPN", "PPARG", "PROX1", "TFF3", "CCL2"), label = F)

DotPlot(seurat.htap.ctrl.female, features = c("FABP4", "FN1", "IFI27", "RGS5", "RNASE1", "TFF3"), group.by = "pat") + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 2", x = "Gene", y = "Patients cluster")


```



























## Tables

```{r}
enrichGO_0 <- data.frame(enrichGO_0)
enrichGO_0$cluster <- rep("cluster0", nrow(enrichGO_0))
enrichDGN_0 <- data.frame(enrichDGN_0)
enrichDGN_0$cluster <- rep("cluster0", nrow(enrichDGN_0))
enrichDO_0 <- data.frame(enrichDO_0)
enrichDO_0$cluster <- rep("cluster0", nrow(enrichDO_0))


enrichGO_1 <- data.frame(enrichGO_1)
enrichGO_1$cluster <- rep("cluster1", nrow(enrichGO_1))
enrichDGN_1 <- data.frame(enrichDGN_1)
enrichDGN_1$cluster <- rep("cluster1", nrow(enrichDGN_1))
enrichDO_1 <- data.frame(enrichDO_1)
enrichDO_1$cluster <- rep("cluster1", nrow(enrichDO_1))

enrichGO_3 <- data.frame(enrichGO_3)
enrichGO_3$cluster <- rep("cluster3", nrow(enrichGO_3))
enrichDGN_3 <- data.frame(enrichDGN_3)
enrichDGN_3$cluster <- rep("cluster3", nrow(enrichDGN_3))
enrichDO_3 <- data.frame(enrichDO_3)
enrichDO_3$cluster <- rep("cluster3", nrow(enrichDO_3))

enrichGO_4 <- data.frame(enrichGO_4)
enrichGO_4$cluster <- rep("cluster4", nrow(enrichGO_4))
enrichDGN_4 <- data.frame(enrichDGN_4)
enrichDGN_4$cluster <- rep("cluster4", nrow(enrichDGN_4))
enrichDO_4 <- data.frame(enrichDO_4)
enrichDO_4$cluster <- rep("cluster4", nrow(enrichDO_4))

enrichGO <- rbind(enrichGO_0, enrichGO_1, enrichGO_3, enrichGO_4)
write.table(enrichGO, file="./output/enrichGO_HTAP_CTRL.csv", sep=";")

enrichDO <- rbind(enrichDO_0, enrichDO_1, enrichDO_3, enrichDO_4)
enrichDO$Database <- rep("Disease Ontology (Yu et al. 2015)", nrow(enrichDO))
enrichDGN <- rbind(enrichDGN_0, enrichDGN_1, enrichDGN_3, enrichDGN_4)
enrichDGN$Database <- rep("DisGeNET (Janet et al. 2015)", nrow(enrichDGN))
enrichDiseases <- rbind(enrichDO, enrichDGN)
write.table(enrichDO, file="./output/enrichDiseases_HTAP_CTRL.csv", sep=";")

enrichKEGG <- rbind(enrichKEGG_0, enrichKEGG_1, enrichKEGG_3, enrichKEGG_4)
write.table(enrichKEGG, file="./output/enrichKEGG_HTAP_CTRL.csv", sep=";")
```







# DGE using limma package


```{r}

DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)
sce.htap.ctrl.female <- as.SingleCellExperiment(seurat.htap.ctrl.female)

agg.clust <- aggregateAcrossCells(sce.htap.ctrl.female, ids = colData(sce.htap.ctrl.female)[,c("seurat_clusters")])
agg.clust <- assay(agg.clust)
agg.clust <- round(agg.clust/rowSums(agg.clust)*100, 2)
agg.clust[is.na(agg.clust)] = 0.00
agg.clust <- data.frame(agg.clust)
colnames(agg.clust) <- c("cluster0", "cluster1", "cluster2", "cluster3", "cluster4")
agg.clust$gene <- rownames(agg.clust)
head(agg.clust)


agg.patient <- aggregateAcrossCells(sce.htap.ctrl.female, ids = colData(sce.htap.ctrl.female)[,c("PATIENT")])
agg.patient <- assay(agg.patient)
# agg.patient <- round(agg.patient/rowSums(agg.patient)*100, 2)
# agg.patient[is.na(agg.patient)] = 0.00
agg.patient <- data.frame(agg.patient)
agg.patient$gene <- rownames(agg.patient)
#colnames(agg.patient) <- c("cluster0", "cluster1", "cluster2", "cluster3", "cluster4")
head(agg.patient)

agg.clust.patient <- merge(agg.patient, agg.clust, by = "gene")
#View(agg.clust.patient)

agg <- aggregateAcrossCells(sce.htap.ctrl.female, ids = colData(sce.htap.ctrl.female)[,c("seurat_clusters", "PHENOTYPE")])
#colData(agg)[,c("seurat_clusters", "PHENOTYPE",  "ncells")]

group <- factor(agg$PHENOTYPE)
dge <- DGEList(assay(agg), group = group)

# keep <- filterByExpr(dge, group=group)
# dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# plotMD(cpm(dge, log=TRUE), column=1)
# abline(h=0, col="red", lty=2, lwd=2)

# Desing matrix HTAP and CTRL conditions
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
con <- makeContrasts(HTAP.vs.CTRL = HTAP - CTRL,
                     #CTRL.vs.HTAP = CTRL - HTAP,
                     levels=colnames(design))

v <- voom(dge, design, plot=F)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=con)
efit <- eBayes(vfit, robust = TRUE)
topTable.x <- limma::topTable(efit, number = "Inf")
topTable.x$gene <- rownames(topTable.x)
plotMD(efit, main = "HTAP vs CTRL")

limma.topTable.HTAP.vs.CTRL.female <- merge(agg.clust.patient, topTable.x, by = "gene")

#View(limma.topTable.HTAP.vs.CTRL.female)


# Maplot
# maplot <- topTable.x[, c("AveExpr", "logFC", "adj.P.Val")]
# colnames(maplot) <- c("baseMean", "log2FoldChange", "padj")
# maplot <- maplot[!apply(is.na(maplot) | maplot == "", 1, all),]
# 
# #library("ggpubr")
# ggpubr::ggmaplot(maplot,
#                  main = "MA plot of HTAP vs CTRL (FDR 5%)",
#                  fdr = 0.05,
#                  fc = 1.5,
#                  size = 1.5,
#                  palette = c("#B31B21", "#1465AC", "darkgray"),
#                  genenames = as.vector(rownames(maplot)),
#                  legend = "right",
#                  top = 15,
#                  alpha = 1,
#                  xlab = "Log2 mean expression",
#                  ylab = "Log2 fold change",
#                  select.top.method = c("fc"),
#                  font.label = c("bold", 6),
#                  label.rectangle = TRUE,
#                  font.legend = "bold",
#                  #font.main = "bold",
#                  ggtheme = ggplot2::theme_classic())






# Male

# DefaultAssay(seurat.htap.ctrl.male) <- "RNA"
# seurat.htap.ctrl <- NormalizeData(seurat.htap.ctrl.male)
# sce.htap.ctrl.male <- as.SingleCellExperiment(seurat.htap.ctrl.male)
# 
# agg.clust <- aggregateAcrossCells(sce.htap.ctrl.male, ids = colData(sce.htap.ctrl.male)[,c("seurat_clusters")])
# agg.clust <- assay(agg.clust)
# agg.clust <- round(agg.clust/rowSums(agg.clust)*100, 2)
# agg.clust[is.na(agg.clust)] = 0.00
# agg.clust <- data.frame(agg.clust)
# colnames(agg.clust) <- c("cluster0", "cluster1", "cluster2", "cluster3", "cluster4")
# agg.clust$gene <- rownames(agg.clust)
# head(agg.clust)
# 
# 
# agg.patient <- aggregateAcrossCells(sce.htap.ctrl.male, ids = colData(sce.htap.ctrl.male)[,c("PATIENT")])
# agg.patient <- assay(agg.patient)
# #agg.patient <- round(agg.patient/rowSums(agg.patient)*100, 3)
# #agg.patient[is.na(agg.patient)] = 0.000
# agg.patient <- data.frame(agg.patient)
# agg.patient$gene <- rownames(agg.patient)
# #colnames(agg.patient) <- c("cluster0", "cluster1", "cluster2", "cluster3", "cluster4")
# head(agg.patient)
# 
# agg.clust.patient <- merge(agg.patient, agg.clust, by = "gene")
# 
# 
# agg <- aggregateAcrossCells(sce.htap.ctrl.male, ids = colData(sce.htap.ctrl.male)[,c("seurat_clusters", "PHENOTYPE")])
# #colData(agg)[,c("seurat_clusters", "PHENOTYPE",  "ncells")]
# 
# group <- factor(agg$PHENOTYPE)
# dge <- DGEList(assay(agg), group = group)
# 
# # keep <- filterByExpr(dge, group=group)
# # dge <- dge[keep, , keep.lib.sizes=FALSE]
# dge <- calcNormFactors(dge)
# 
# # plotMD(cpm(dge, log=TRUE), column=1)
# # abline(h=0, col="red", lty=2, lwd=2)
# 
# # Desing matrix HTAP and CTRL conditions
# design <- model.matrix(~0+group)
# colnames(design) <- gsub("group", "", colnames(design))
# con <- makeContrasts(HTAP.vs.CTRL = HTAP - CTRL,
#                              #CTRL.vs.HTAP = CTRL - HTAP,
#                              levels=colnames(design))
# 
# v <- voom(dge, design, plot=F)
# vfit <- lmFit(v, design)
# vfit <- contrasts.fit(vfit, contrasts=con)
# efit <- eBayes(vfit, robust = TRUE)
# topTable.y <- limma::topTable(efit, number = "Inf")
# topTable.y$gene <- rownames(topTable.y)
# plotMD(efit, main = "HTAP vs CTRL")
# 
# 
# limma.topTable.HTAP.vs.CTRL.male <- merge(agg.clust.patient, topTable.y, by = "gene")
# 
# #View(limma.topTable.HTAP.vs.CTRL.male)
# 
# 
# 
# # Merge female and male dateset
# limma.topTable.HTAP.vs.CTRL <- merge(limma.topTable.HTAP.vs.CTRL.female , limma.topTable.HTAP.vs.CTRL.male, by = "gene" )
# 
# View(limma.topTable.HTAP.vs.CTRL)
# 
# write.table(limma.topTable.HTAP.vs.CTRL, file="./output/limma.topTable.HTAP.vs.CTRL.csv", sep=";")


# Volcanoplot









#counts <- seurat.htap.ctrl@assays$RNA@counts


# Molecular function of the interest genes in HTAP vs CTRL
# library("gprofiler2")
# 
# gprofiler2.htap.vs.ctrl <- gprofiler2::gconvert(query = rownames(topLimma), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
# View(gprofiler2.htap.vs.ctrl)
# gprofiler2.htap.vs.ctrl <- gprofiler2.htap.vs.ctrl[, -c(1:3)]
# gprofiler2.htap.vs.ctrl <- gprofiler2.htap.vs.ctrl[, c("name","target", "description")]
#View(gprofiler2.htap.vs.ctrl)


# Gene list functional enrichment analysis with gost

# gostres <- gost(query = top.genes.clusters, 
#                 organism = "hsapiens", ordered_query = FALSE, 
#                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
#                 measure_underrepresentation = FALSE, evcodes = TRUE, 
#                 user_threshold = 0.05, correction_method = "g_SCS", 
#                 domain_scope = "annotated", custom_bg = NULL, 
#                 numeric_ns = "", sources = NULL, as_short_link = FALSE)
# 
# 
# names(gostres)
# #View(gostres$result)
# 
# p <- gostplot(gostres, capped = TRUE, interactive =  FALSE)
# 
# publish_gostplot(p, highlight_terms = gostres$result$term_id, 
#                        width = 8, height = 20, filename = NULL )


```






```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}
View(topTable.x[abs(topTable.x$logFC) >= 1.5 & topTable.x$P.Value <= 10e-05, ])

top.markers <- rownames(topTable.x[topTable.x$adj.P.Val <= 0.05,])

library("EnhancedVolcano")
EnhancedVolcano(topTable.x,
                lab = rownames(topTable.x),
                x = 'logFC',
                y = 'P.Value',
                #xlim = c(-8, 8),
                ylim = c(0, -log10(10e-12)),
                title = "Differential Gene Expression \n
    HTAP vs CTRL in female dataset",
                subtitle = "(linear model by limma package, 
    robust = true, FCcutoff = 1.5, pCutoff = 10e-05)" ,
                #selectLab = top.c8.genes,
                #xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                #selectLab = top.markers,
                #shade = top.markers,
                #shadeAlpha = 1/2,
                # shadeFill = 'purple',
                # borderWidth = 1.0,
                # borderColour = 'black',
                #boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()


# EnhancedVolcano(topTable.y,
#     lab = rownames(topTable.y),
#     x = 'logFC',
#     y = 'P.Value',
#      #xlim = c(-8, 8),
#     ylim = c(0, -log10(10e-12)),
#     title = "Differential Gene Expression \n
#     HTAP vs CTRL in male dataset",
#     subtitle = "(linear model by limma package, 
#     robust = true, FCcutoff = 1.2, pCutoff = 10e-03)" ,
#     #selectLab = top.c8.genes,
#     #xlab = bquote(~Log[2]~ 'fold change'),
#     pCutoff = 10e-05,
#     FCcutoff = 1.5,
#     pointSize = 1.0,
#      labSize = 3.5,
#     #labCol = 'black',
#     labFace = 'bold',
#     #boxedLabels = TRUE,
#     colAlpha = 1,
#     legendPosition = 'right',
#     legendLabSize = 8,
#    legendIconSize = 2.0,
#     drawConnectors = TRUE,
#     widthConnectors = 0.2,
#     #colConnectors = 'black'
#    ) + ggplot2::theme_linedraw()
```

# Gene enrichment

```{r}
## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
```






# ADA2 and NOMO1

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)

# Enrichment analysis of ADA2 and NOMO1
sig <- subset(topTable.x,
              gene == c("NOMO1", "ADA2"))
colnames(sig) <- c("logFC","AveExpr","t" , "P.Value" ,"adj.P.Val", "B", "SYMBOL")
sig.df <- bitr(sig$SYMBOL, fromType = "SYMBOL",
               toType = c("ENSEMBL", "ENTREZID"),
               OrgDb = org.Hs.eg.db)
d <- merge(sig.df, sig, by="SYMBOL")
## Remove any Entrez duplicates
d <- d[which(duplicated(d$ENTREZID)==F), ]
ada2_nomo1 <- d$logFC
names(ada2_nomo1) <- as.character(d$ENTREZID)
ada2_nomo1 <- sort(ada2_nomo1, decreasing = TRUE)

## GSEA using gene sets from KEGG pathways

## Output images for a single significant KEGG pathway

#data(paths.hsa)
#View(paths.hsa)


enrichKEGG_htap.ctrl <- enrichKEGG(gene         = names(ada2_nomo1),
                                   organism     = 'hsa',
                                   pvalueCutoff = 0.05)
enrichKEGG_HTAP.vs.CTRL <- data.frame(enrichKEGG_htap.ctrl)
enrichKEGG_HTAP.vs.CTRL$Gene_expr_analyse <- rep("HTAP.vs.CTRL", nrow(enrichKEGG_HTAP.vs.CTRL))
View(data.frame(enrichKEGG_HTAP.vs.CTRL))


## Output images for a single significant KEGG pathway
pathview(gene.data = ada2_nomo1,
         pathway.id = "hsa00230",
         species = "hsa",
         limit = list(gene = max(abs(ada2_nomo1)), # value gives the max/min limit for foldchanges
                      cpd = 1))




## Run GO enrichment analysis 
#ada2_nomo1_DOSE <- data(ada2_nomo1, package="DOSE")

enrichGO_htap.ctrl <- enrichGO(gene          = unique(sig.df$ENTREZID),
                               #universe      = names(ada2_nomo1_DOSE),
                               #keyType = "ENSEMBL",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE,
                               pool = TRUE)
enrichGO_HTAP.vs.CTRL <- setReadable(enrichGO_htap.ctrl, 'org.Hs.eg.db')
heatplot(enrichGO_HTAP.vs.CTRL, foldChange=ada2_nomo1) + ggtitle("Gene-disease association")
emapplot(enrichGO_HTAP.vs.CTRL, showCategory=50, pie_scale=0.8,layout="kk", line_scale = 0.2) + ggtitle("Relationship between the top 50 most significantly \n enriched GO terms in HTAP vs CTRL")
enrichGO_HTAP.vs.CTRL <- data.frame(enrichGO_HTAP.vs.CTRL)
enrichGO_HTAP.vs.CTRL$Gene_expr_analyse <- rep("HTAP.vs.CTRL", nrow(enrichGO_HTAP.vs.CTRL))
View(enrichGO_HTAP.vs.CTRL)


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN_HTAP.vs.CTRL <- enrichDGN(names(ada2_nomo1))
enrichDGN_HTAP.vs.CTRL <- setReadable(enrichDGN_HTAP.vs.CTRL, 'org.Hs.eg.db')
enrichDGN_HTAP.vs.CTRL <- data.frame(enrichDGN_HTAP.vs.CTRL)
enrichDGN_HTAP.vs.CTRL$Database <- rep("DisGeNET (Janet et al. 2015)", nrow(enrichDGN_HTAP.vs.CTRL))
# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO_HTAP.vs.CTRL <- enrichDO(gene          = names(ada2_nomo1),
                                  ont           = "DO",
                                  pvalueCutoff  = 0.05,
                                  pAdjustMethod = "BH",
                                  #universe      = names(genes.c0.List),
                                  minGSSize     = 5,
                                  maxGSSize     = 500,
                                  qvalueCutoff  = 0.05,
                                  readable      = FALSE)
enrichDO_HTAP.vs.CTRL <- setReadable(enrichDO_HTAP.vs.CTRL, 'org.Hs.eg.db')
enrichDO_HTAP.vs.CTRL <- data.frame(enrichDO_HTAP.vs.CTRL)
enrichDO_HTAP.vs.CTRL$Database <- rep("Disease Ontology (Yu et al. 2015)", nrow(enrichDO_HTAP.vs.CTRL))



#GOplot

david <- data.frame(enrichGO_htap.ctrl)[, c( "ONTOLOGY","ID","Description","geneID" , "p.adjust")]
colnames(david) <- c("Category", "ID","Term", "Genes", "adj_pval")
rownames(david) <- NULL
# Generate the plotting object
Sig_ada2_nomo1 <- sig[, c( "SYMBOL", "logFC" ,"AveExpr","t" ,"P.Value","adj.P.Val" ,"B")]
colnames(Sig_ada2_nomo1) <- c( "ID","logFC", "AveExpr","t",  "P.Value", "adj.P.Val","B")
rownames(Sig_ada2_nomo1) <- NULL
circ <- circle_dat(david, Sig_ada2_nomo1)
#View(circ)
#GOBar(subset(circ, category == 'BP'))

# Facet the barplot according to the categories of the terms

GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))

# Generate the bubble plot with a label threshold of 3
GOBubble(circ, label = 1)


# Display of the relationship between genes and terms (GOChord)

genes <- circ[, c("genes","logFC")]
colnames(genes) <- c("ID", "logFC")
genes <- genes[!duplicated(genes),]
process <- circ$term

# Now it is time to generate the binary matrix
chord <- chord_dat(circ, genes, process)
head(chord)

# Generate the matrix with a list of selected genes
chord <- chord_dat(data = circ, genes = genes)
# Generate the matrix with selected processes
chord <- chord_dat(data = circ, process = process)

# Create the plot
chord <- chord_dat(data = circ, genes = genes, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)


barplot(ego, showCategory=50) + ggtitle("GO enrichment chart of ADA2 and NOMO1 markers (HTAP vs CTRL)")
dotplot(ego, showCategory=100) + ggtitle("Top 50 GO enrichment dotplot \n of ADA2 and NOMO1 in HTAP vs CTRL")

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 50, pie_scale=0.8,layout="kk", line_scale = 0.1) + ggtitle("Enrichmap gathering the 50 most significant GO terms (by padj) \n to visualize the relations between the terms")

edo <- enrichDGN(names(ada2_nomo1))
barplot(edo, showCategory=30) + ggtitle("Gene-disease association")




# Heatmap-like functional and disease classifications
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
heatplot(edox, foldChange=ada2_nomo1) + ggtitle("Gene-disease association")



## Cnetplot details the genes associated with one or more terms 

cnetplot(ego,
         categorySize="pvalue",
         showCategory = 50,
         foldChange=ada2_nomo1,
         vertex.label.font=6)

cnetplot(ego, foldChange=ada2_nomo1, circular = TRUE, colorEdge = TRUE)



#gene_ada2_nomo1 <- names(ada2_nomo1)[abs(ada2_nomo1) > 1]
kk <- enrichKEGG(gene         = names(ada2_nomo1),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
View(data.frame(kk))

# Pathview
pathview(gene.data  = ada2_nomo1,
         pathway.id = "hsa00230",
         species    = "hsa",
         limit      = list(gene=max(abs(ada2_nomo1)), cpd=1))




```






```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=15, fig.width=28}
#GOChord(chord, space = 0.2, gene.order = 'logFC', gene.space = 0.2, gene.size = 10)

# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord[,-46])
dim(chord)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

GOBubble(circ, labels = 1)


```








```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6.5, fig.width=8}
DotPlot(seurat.htap.ctrl.female, features = c("ADA2", "NOMO1"), group.by = "patient_htap.ctrl") + RotatedAxis() + labs(title = "DGE in HTAP vs CTRL (LogFC ≥ 1.5, FDR ≤ 0.05)", x = "Gene", y = "Patients cluster")

```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=12}
FeaturePlot(seurat.htap.ctrl.female, features = c("NOMO1", "ADA2"), split.by = "PHENOTYPE",  pt.size = 1, min.cutoff = "q1", cols = c("lightgrey", "blue"))
```










```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=4}
#top <- topTable[order(topTable$AveExpr, decreasing = T),]
top.x <- topTable.x %>% mutate(Color = ifelse(logFC < 0, "blue", "red"))

top.x <- head(top.x[order(top.x$AveExpr, decreasing = T),],30)
barplot(top.x$logFC, 
        main="HTAP vs CTRL",
        names.arg = top.x$gene,
        horiz = T,
        las=1, 
        cex.names = 0.5,
        #xlim = c(-10, 10), 
        col = top.x$Color,
        #legend = c("Up", "Down"),
        beside = F, 
        #xlab = "Genes",
        xlab = "LogFC")
# text(p, srt = 45, xpd = TRUE,
#      labels = paste(top$SYMBOL), cex = 0.65)

```



# Functional gene annotations


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}
# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(fgsea)
library(GOplot)


#top.x <- topTable.x[abs(topTable.x$logFC) > 1.0 & topTable.x$P.Value <= 0.05,]
topTable.xx <- topTable.x
#colnames(topTable.xx) <- c("logFC","AveExpr","t" , "P.Value" ,"adj.P.Val", "B", "SYMBOL")
topTable.xx.df <- bitr(topTable.xx$SYMBOL, fromType = "SYMBOL",
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = org.Hs.eg.db)
v <- merge(topTable.xx.df, topTable.xx, by="SYMBOL")
## Remove any Entrez duplicates
v <- v[which(duplicated(v$ENTREZID)==F), ]
genes.xx <- v$logFC
names(genes.xx) <- as.character(v$ENTREZID)
genes.xx <- sort(genes.xx, decreasing = TRUE)
genes.xx <- na.exclude(genes.xx)


enrichGO <- enrichGO(gene          = topTable.xx.df$ENTREZID,
                     universe      = names(geneList_DOSE),
                     #keyType = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE,
                     pool = TRUE)

head(enrichGO, n=10)
View(data.frame(enrichGO))



## GSEA using gene sets from KEGG pathways
genes <- names(genes.xx)[abs(genes.xx) > 1]
genes <- sort(genes, decreasing = TRUE)
kk <- enrichKEGG(gene         = genes,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)


gseaKEGG <- gseKEGG(geneList = genes.xx, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
View(gseaKEGG_results)
```





# HTAP vs CTRL

```{r}
# Visualization
p1 <- DimPlot(seurat.htap.ctrl.female, reduction = "umap", group.by = "PHENOTYPE")
p2 <- DimPlot(seurat.htap.ctrl.female, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(seurat.htap.ctrl.female, reduction = "umap", split.by = "PHENOTYPE")
```



```{r}
Idents(seurat.htap.ctrl.female) <- "PHENOTYPE"
avg.cells <- log1p(AverageExpression(seurat.htap.ctrl.female, verbose = FALSE)$RNA)
avg.cells$gene <- rownames(avg.cells)

#genes.to.label <- "TPM2"
p1 <- ggplot(avg.cells, aes(HTAP, CTRL)) + geom_point() + ggtitle("Average Expression (HTAP vs CTRL)")
#p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=6}
Idents(seurat.htap.ctrl.female) <- "pheno"

pheno_0 <- FindMarkers(seurat.htap.ctrl.female, ident.1 = "0_HTAP", ident.2 = "0_CTRL", verbose = FALSE, only.pos = F)
DotPlot(seurat.htap.ctrl.female, features=rownames(head(pheno_0[order(pheno_0$avg_logFC, decreasing = TRUE),],50)), group.by = "pheno", idents = c('0_HTAP', '0_CTRL')) + RotatedAxis()  +labs(title = "HTAP vs CTRL in cluster 0", x = "Genes", y = "Clusters")

EnhancedVolcano(pheno_0,
                lab = rownames(pheno_0),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "HTAP vs CTRL in cluster 0", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()



pheno_1 <- FindMarkers(seurat.htap.ctrl.female, ident.1 = "1_HTAP", ident.2 = "1_CTRL", verbose = FALSE, only.pos = F)
DotPlot(seurat.htap.ctrl.female, features=rownames(head(pheno_1[order(pheno_1$avg_logFC, decreasing = TRUE),],50)), group.by = "pheno", idents = c('1_HTAP', '1_CTRL')) + RotatedAxis()  +labs(title = "HTAP vs CTRL in cluster 1", x = "Genes", y = "Clusters")

EnhancedVolcano(pheno_1,
                lab = rownames(pheno_1),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "HTAP vs CTRL in cluster 1", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()




pheno_3 <- FindMarkers(seurat.htap.ctrl.female, ident.1 = "3_HTAP", ident.2 = "3_CTRL", verbose = FALSE, only.pos = F)
head(pheno_3, n = 15)
DotPlot(seurat.htap.ctrl.female, features=rownames(head(pheno_3[order(pheno_3$avg_logFC, decreasing = TRUE),],50)), group.by = "pheno", idents = c('3_HTAP', '3_CTRL')) + RotatedAxis() +labs(title = "HTAP vs CTRL in cluster 3", x = "Genes", y = "Clusters")


EnhancedVolcano(pheno_3,
                lab = rownames(pheno_3),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "HTAP vs CTRL in cluster 3", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```




```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=12, fig.width=5}

FeaturePlot(seurat.htap.ctrl.female, features=rownames(head(pheno_0[order(pheno_0$avg_logFC, decreasing = TRUE),],50)), split.by = "PHENOTYPE")


plots <- VlnPlot(seurat.htap.ctrl.female, features = rownames(head(pheno_0[order(pheno_0$avg_logFC, decreasing = TRUE),],50)), split.by = "PHENOTYPE", group.by = "pheno", 
                 pt.size = 0, combine = FALSE)
library(patchwork)
wrap_plots(plots = plots, ncol = 1)
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=20, fig.width=6}
DefaultAssay(seurat.htap.ctrl.female) <- "RNA"
seurat.htap.ctrl.female <- NormalizeData(seurat.htap.ctrl.female)

FeaturePlot(seurat.htap.ctrl.female, features=rownames(head(pheno_1[order(pheno_1$avg_logFC, decreasing = TRUE),],50)), split.by = "PHENOTYPE")


plots <- VlnPlot(seurat.htap.ctrl.female, features = rownames(head(pheno_1[order(pheno_1$avg_logFC, decreasing = TRUE),],50)), split.by = "PHENOTYPE", group.by = "pheno", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)
```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=20, fig.width=6}

FeaturePlot(seurat.htap.ctrl.female, features=rownames(head(pheno_3[order(pheno_3$avg_logFC, decreasing = TRUE),],50)), split.by = "PHENOTYPE")

plots <- VlnPlot(seurat.htap.ctrl.female, features = rownames(head(pheno_3[order(pheno_3$avg_logFC, decreasing = TRUE),],50)), split.by = "PHENOTYPE", group.by = "pheno", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
```





```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(pheno_0,
                lab = rownames(pheno_0),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "HTAP vs CTRL in cluster 2", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```



# HTAP.BMP9 vs HTAP

```{r}
seurat.htap.ctrl <- readRDS("/data/data_mbouamboua/10x_tu/output/seurat.htap.ctrl.rds")

DefaultAssay(seurat.htap.ctrl) <- "integrated"
seurat.htap.ctrl <- FindNeighbors(seurat.htap.ctrl, reduction = "mnn", dims = 1:20)
seurat.htap.ctrl <- FindClusters(seurat.htap.ctrl, resolution = 0.15)
DimPlot(object = seurat.htap.ctrl, reduction = "umap", label = TRUE) + labs(colour = "Cell clusters")
seurat.htap.ctrl$cluster_htap.ctrl <- paste0(seurat.htap.ctrl$seurat_clusters, "_", seurat.htap.ctrl$PHENOTYPE)
seurat.htap.ctrl$patient_htap.ctrl <- paste0(seurat.htap.ctrl$seurat_clusters, "_", seurat.htap.ctrl$PATIENT)


seurat.htap.ctrl[["sex"]]=""
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C686518"),] <- "M"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C681725"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C681086"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C685116"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "C676136"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H247"),] <- "M"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H248"),] <- "M"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H230"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H231"),] <- "F"
seurat.htap.ctrl[['sex']][which(seurat.htap.ctrl@meta.data$SAMPLE == "H250"),] <- "F"

#seurat.htap.ctrl.male <- subset(seurat.htap.ctrl, sex == "M") 
seurat.htap.ctrl.F <- subset(seurat.htap.ctrl, sex == "F") 
seurat.htap.ctrl.female <- FindNeighbors(seurat.htap.ctrl.F, reduction = "mnn", dims = 1:20)
seurat.htap.ctrl.female <- FindClusters(seurat.htap.ctrl.female, resolution = 0.10)
DimPlot(object = seurat.htap.ctrl.female, reduction = "umap", label = TRUE)
seurat.htap.ctrl.female$pheno <- paste0(seurat.htap.ctrl.female$seurat_clusters, "_", seurat.htap.ctrl.female$PHENOTYPE)
seurat.htap.ctrl.female$pat <- paste0(seurat.htap.ctrl.female$seurat_clusters, "_", seurat.htap.ctrl.female$PATIENT)
seurat.htap.ctrl.female <- RenameIdents(seurat.htap.ctrl.female, `0` = "0", `1` = "1", `2` = "Cell cycling", `3` = "3")
DimPlot(object = seurat.htap.ctrl.female, label = F, label.size = 5)



# # Define cluster colours
# # Cluster 0
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
seurat.htap.ctrl.female[['cluster0']] <- " "
seurat.htap.ctrl.female[['cluster0']][which(seurat.htap.ctrl.female@active.ident == 0),] <- 0
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster0")
DimPlot(seurat.htap.ctrl.female, group.by = "cluster0", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)

seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
seurat.htap.ctrl.female[['cluster1']] <- " "
seurat.htap.ctrl.female[['cluster1']][which(seurat.htap.ctrl.female@active.ident == 1),] <- 1
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster1")
DimPlot(seurat.htap.ctrl.female, group.by = "cluster1", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)

seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
seurat.htap.ctrl.female[['cluster2']] <- " "
seurat.htap.ctrl.female[['cluster2']][which(seurat.htap.ctrl.female@active.ident == 3),] <- 3
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster2")
DimPlot(seurat.htap.ctrl.female, group.by = "cluster2", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)


seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "seurat_clusters")
seurat.htap.ctrl.female[['cluster3']] <- " "
seurat.htap.ctrl.female[['cluster3']][which(seurat.htap.ctrl.female@active.ident == 2),] <- 2
seurat.htap.ctrl.female <- SetIdent(seurat.htap.ctrl.female, value = "cluster3")
DimPlot(seurat.htap.ctrl.female, group.by = "cluster3", cols = c("#C0C0C0", "#e31a1c"), reduction = "umap", label = T, label.size = 6)
```



19/10/2020

# Marking differences between clusters 0 and 1

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=16}
markers_0.vs.1 <- FindMarkers(object = ctrl, ident.1=0, ident.2=1, only.pos = FALSE)
# kable(head(markers_0.vs.1,10))
# markers_0 <- head(markers_0.vs.1,50)
DotPlot(ctrl,features=rownames(head(markers_0.vs.1[order(markers_0.vs.1$avg_logFC, decreasing = T),],50)), idents = c('0','1')) + RotatedAxis() + labs(title = "Clusters 0 vs 1 (CTRL dataset)", x = "Markers", y = "Clusters")
```



```{r}
gprofiler_0 <- gprofiler2::gconvert(query = rownames(markers_0.vs.1), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
markers_0.vs.1$SYMBOL <- rownames(markers_0.vs.1)
genes.c0 <- bitr(gprofiler_0$target, fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)
## Remove any Entrez duplicates
genes.c0 <- genes.c0[which(duplicated(genes.c0$ENTREZID) == F), ]
genes.c0.List <- markers_0.vs.1$avg_logFC
names(genes.c0.List) <- as.character(genes.c0$ENTREZID)
## Remove any NA values
genes.c0.List <- genes.c0.List[!is.na(names(genes.c0.List))]
#genes.c0.List <- genes.c0.List[abs(genes.c0.List) >= 1]
genes.c0.List <- sort(genes.c0.List, decreasing = TRUE)
genes.c0.List


#data(geneList, package="DOSE")
enrichGO0 <- enrichGO(gene = names(genes.c0.List),
                      #universe      = names(geneList),
                      keyType = "ENTREZID",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE,
                      pool = TRUE)
enrichGO_0 <- setReadable(enrichGO0, 'org.Hs.eg.db')


#View(data.frame(enrichGO0))


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN0 <- enrichDGN(names(genes.c0.List))
enrichDGN_0 <- setReadable(enrichDGN0, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO0 <- enrichDO(gene          = names(genes.c0.List),
                      ont           = "DO",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      #universe      = names(genes.c0.List),
                      minGSSize     = 5,
                      maxGSSize     = 500,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)
enrichDO_0 <- setReadable(enrichDO0, 'org.Hs.eg.db')
```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(markers_0.vs.1,
                lab = rownames(markers_0.vs.1),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 0 vs 1", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```


#### Barplot

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=12, fig.width=14}

barplot(enrichGO0, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 0 (HTAP + CTRL)", x = "GO terms", y = "Genes related to GO term ")

```



#### Emaplot

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}
## Plotting terms of interest
cnetplot(enrichGO_0, 
         categorySize="pvalue", 
         foldChange=genes.c0.List, 
         showCategory = 50, 
         vertex.label.font=6,
         colorEdge = TRUE,
         circular = TRUE)

```


#### Heatmap plot


```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=25}
# Heatmap-like functional classification
heatplot(enrichGO_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Top 50 enriched  Gene Ontology terms  in cluster 0", x = "Genes", y = "GO terms")

heatplot(enrichDO_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015) in cluster 0 ", x = "Genes", y = "Disease Ontology")

heatplot(enrichDGN_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Gene-disease association (DisGeNET, Janet et al. 2015)   in cluster 0 ", y = "Diseases", x = "Genes")

```

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=6}
heatplot(enrichGO_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Top 50 enriched Gene Ontology terms \n in cluster 0 (CTRL dataset)", x = "Genes", y = "GO terms")

```

```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=6}
# Heatmap-like functional classification

heatplot(enrichDO_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Gene-disease association \n (Disease Ontology, Yu et al. 2015) \n in cluster 0 (CTRL dataset)", x = "Genes", y = "Disease Ontology")


heatplot(enrichDGN_0, foldChange=genes.c0.List, showCategory=50) + labs(title = "Gene-disease association \n (DisGeNET, Janet et al. 2015) \n  in cluster 0 (CTRL dataset)", y = "Diseases", x = "Genes")
```

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=7}

FeaturePlot(ctrl, features = c("ALDH1A1", "FABP4"), label = F)

DotPlot(ctrl, features = c("ALDH1A1", "FABP4"), group.by = "pat", idents = c('0', '1', '3')) + RotatedAxis() + labs(title = "Specific markers in cluster 0 ", x = "Gene", y = "Patients cluster")


```
## KEGG over-representation

```{r}

enrichKEGG_0 <- enrichKEGG(gene         = names(genes.c0.List),
                           organism     = 'hsa',
                           pvalueCutoff = 0.05)
enrichKEGG_0 <- data.frame(enrichKEGG_0)
enrichKEGG_0$cluster <- rep("cluster0", nrow(enrichKEGG_0))
View(data.frame(enrichKEGG_0))

## Output images for a single significant KEGG pathway
pathview(gene.data = genes.c0.List,
         pathway.id = "hsa05205",
         species = "hsa",
         limit = list(gene = max(abs(genes.c0.List)), # value gives the max/min limit for foldchanges
                      cpd = 1))


## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
  pathview(gene.data = genes.c0.List, pathway.id = enrichKEGG_0$ID[x], species = "hsa", 
           limit = list(gene = max(abs(genes.c0.List)), cpd = 1))
}

purrr::map(1:length(enrichKEGG_0$ID), get_kegg_plots)

```

# Making differences between cluster 1 and all


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=12}
DefaultAssay(ctrl) <- "RNA"
ctrl <- NormalizeData(ctrl)
ctrl <- SetIdent(ctrl, value = "seurat_clusters")

markers_1 <- FindMarkers(ctrl, ident.1 = 1, ident.2 = c(0, 3),   min.pct = 0.25, only.pos = F)

DotPlot(ctrl,features=rownames(head(markers_1[order(markers_1$avg_logFC, decreasing = T),],50)), idents = c('0','1', '3')) + RotatedAxis() + labs(title = "Clusters 1 vs 0 and 3 (CTRL dataset)", x = "Markers", y = "Clusters")

# Data
gprofiler_1 <- gprofiler2::gconvert(query = rownames(markers_1), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
gprofiler_1 <- gprofiler_1[which(duplicated(gprofiler_1$target) == F), ]
markers_1$SYMBOL <- rownames(markers_1)
genes.c1 <- bitr(gprofiler_1$target, fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)

## Remove any Entrez duplicates
genes.c1 <- genes.c1[which(duplicated(genes.c1$SYMBOL) == F), ]
genes.c1.List <- markers_1$avg_logFC
names(genes.c1.List) <- as.character(genes.c1$ENTREZID)
## Remove any NA values
genes.c1.List <- genes.c1.List[!is.na(names(genes.c1.List))]
#genes.c1.List <- genes.c1.List[abs(genes.c1.List) >= 1]
genes.c1.List <- sort(genes.c1.List, decreasing = TRUE)
genes.c1.List


#data(geneList, package="DOSE")
enrichGO1 <- enrichGO(gene = names(genes.c1.List),
                      #universe      = names(geneList),
                      keyType = "ENTREZID",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE,
                      pool = TRUE)
enrichGO_1 <- setReadable(enrichGO1, 'org.Hs.eg.db')


#View(data.frame(enrichGO0))


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN1 <- enrichDGN(names(genes.c1.List))
enrichDGN_1 <- setReadable(enrichDGN1, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO1 <- enrichDO(gene          = names(genes.c1.List),
                      ont           = "DO",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      #universe      = names(genes.c0.List),
                      minGSSize     = 5,
                      maxGSSize     = 500,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)
enrichDO_1 <- setReadable(enrichDO1, 'org.Hs.eg.db')

```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(markers_1,
                lab = rownames(markers_1),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 1 vs 0 and 3", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```



```{r  message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=25}
# Heatmap-like functional classification
heatplot(enrichGO_1, foldChange=genes.c1.List, showCategory=50) + labs(title = "Top 50 enriched  Gene Ontology terms  in cluster 1", x = "Genes", y = "GO terms")

heatplot(enrichDO_1, foldChange=genes.c1.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015) in cluster 1 ", x = "Genes", y = "Disease Ontology")

heatplot(enrichDGN_1, foldChange=genes.c1.List, showCategory=50) + labs(title = "Gene-disease association (DisGeNET, Janet et al. 2015)   in cluster 1 ", y = "Diseases", x = "Genes")

```




```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=6}

FeaturePlot(ctrl, features = c("KRT7"), label = F)

DotPlot(ctrl, features = c("KRT7"), group.by = "pat", idents = c('0', '1', '3')) + RotatedAxis() + labs(title = "Specific marker in cluster 1", x = "Gene", y = "Patients cluster")


```

```{r}

enrichKEGG_1 <- enrichKEGG(gene         = names(genes.c1.List),
                           organism     = 'hsa',
                           pvalueCutoff = 0.05)
enrichKEGG_1 <- data.frame(enrichKEGG_1)
enrichKEGG_1$cluster <- rep("cluster1.vs.0.3", nrow(enrichKEGG_1))
View(data.frame(enrichKEGG_1))

## Output images for a single significant KEGG pathway
pathview(gene.data = genes.c1.List,
         pathway.id = "hsa05205",
         species = "hsa",
         limit = list(gene = max(abs(genes.c1.List)), # value gives the max/min limit for foldchanges
                      cpd = 1))


## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
  pathview(gene.data = genes.c1.List, pathway.id = enrichKEGG_1$ID[x], species = "hsa", 
           limit = list(gene = max(abs(genes.c1.List)), cpd = 1))
}

purrr::map(1:length(enrichKEGG_1$ID), get_kegg_plots)

```



# Making differences between clusters 1 and 0

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}
DefaultAssay(ctrl) <- "RNA"
ctrl <- NormalizeData(ctrl)
ctrl <- SetIdent(ctrl, value = "seurat_clusters")

markers_1.vs.0 <- FindMarkers(ctrl, ident.1 = 1, ident.2 = 0,   min.pct = 0.25, only.pos = F)


gprofiler_1.0 <- gprofiler2::gconvert(query = rownames(markers_1.vs.0), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
gprofiler_1.0 <- gprofiler_1.0[which(duplicated(gprofiler_1.0$target) == F), ]
markers_1.vs.0$SYMBOL <- rownames(markers_1.vs.0)
genes.c1.0 <- bitr(gprofiler_1.0$target, fromType = "ENSEMBL",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

## Remove any Entrez duplicates
genes.c1.0 <- genes.c1.0[which(duplicated(genes.c1.0$SYMBOL) == F), ]
genes.c1.0.List <- markers_1.vs.0$avg_logFC
names(genes.c1.0.List) <- as.character(genes.c1.0$ENTREZID)
## Remove any NA values
genes.c1.0.List <- genes.c1.0.List[!is.na(names(genes.c1.0.List))]
#genes.c1.0.List <- genes.c1.0.List[abs(genes.c1.0.List) >= 1]
genes.c1.0.List <- sort(genes.c1.0.List, decreasing = TRUE)
genes.c1.0.List



#data(geneList, package="DOSE")
enrichGO.1.0 <- enrichGO(gene = names(genes.c1.0.List),
                         #universe      = names(geneList),
                         keyType = "ENTREZID",
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE,
                         pool = TRUE)
# View(data.frame(enrichGO.1.0))
enrichGO_1_0 <- setReadable(enrichGO.1.0, 'org.Hs.eg.db')


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN.1.0 <- enrichDGN(names(genes.c1.0.List))
enrichDGN_1_0 <- setReadable(enrichDGN.1.0, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO.1.0 <- enrichDO(gene          = names(genes.c1.0.List),
                         ont           = "DO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         minGSSize     = 5,
                         maxGSSize     = 500,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
enrichDO_1_0 <- setReadable(enrichDO.1.0, 'org.Hs.eg.db')




```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=12}
DotPlot(ctrl, features=rownames(head(markers_1.vs.0[order(markers_1.vs.0$avg_logFC, decreasing = TRUE),],50)), idents = c('0', '1')) + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 0 (CTRL dataset)", x = "Markers", y = "Clusters")
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=30}

#barplot(enrichGO.1.0, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 1 vs 0 (HTAP + CTRL)", x = "GO terms", y = "Genes related to GO term ")

# Heatmap-like functional classification
heatplot(enrichGO_1_0, foldChange=genes.c1.0.List, showCategory=50) + labs(title = "Enriched GO terms in cluster 1 vs 0 (CTRL dataset)", x = "Genes", y = "GO terms")

heatplot(enrichDO_1_0, foldChange=genes.c1.0.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015) in cluster 1 vs 0 (CTRL dataset)", x = "Genes", y = "Disease Ontology")


heatplot(enrichDGN_1_0, foldChange=genes.c1.0.List, showCategory=50) + labs(title = "Gene-disease association   (DisGeNET, Janet et al. 2015) in cluster 1 vs 0 (CTRL dataset)", y = "Diseases", x = "Genes")
```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(markers_1.vs.0,
                lab = rownames(markers_1.vs.0),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 1 vs 0", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=8}
FeaturePlot(ctrl, features = c("KRT7", "MMP1", "MT2A", "NQO1"), label = F)

DotPlot(ctrl, features = c("KRT7", "MMP1", "MT2A", "NQO1"), group.by = "pat") + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 0", x = "Gene", y = "Patients cluster")


```


# Marking differences between clusters 1 and 3

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=8, fig.width=10}
DefaultAssay(ctrl) <- "RNA"
ctrl <- NormalizeData(ctrl)
ctrl <- SetIdent(ctrl, value = "seurat_clusters")

markers_1.vs.3 <- FindMarkers(ctrl, ident.1 = 1, ident.2 = 3,   min.pct = 0.25, only.pos = F)


gprofiler_1.3 <- gprofiler2::gconvert(query = rownames(markers_1.vs.3), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
gprofiler_1.3 <- gprofiler_1.3[which(duplicated(gprofiler_1.3$target) == F), ]
markers_1.vs.3$SYMBOL <- rownames(markers_1.vs.3)
genes.c1.3 <- bitr(gprofiler_1.3$target, fromType = "ENSEMBL",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

## Remove any Entrez duplicates
genes.c1.3 <- genes.c1.3[which(duplicated(genes.c1.3$SYMBOL) == F), ]
genes.c1.3.List <- markers_1.vs.3$avg_logFC
names(genes.c1.3.List) <- as.character(genes.c1.3$ENTREZID)
## Remove any NA values
genes.c1.3.List <- genes.c1.3.List[!is.na(names(genes.c1.3.List))]
#genes.c1.3.List <- genes.c1.3.List[abs(genes.c1.3.List) >= 1]
genes.c1.3.List <- sort(genes.c1.3.List, decreasing = TRUE)
genes.c1.3.List



#data(geneList, package="DOSE")
enrichGO.1.3 <- enrichGO(gene = names(genes.c1.3.List),
                         #universe      = names(geneList),
                         keyType = "ENTREZID",
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE,
                         pool = TRUE)
# View(data.frame(enrichGO.1.0))
enrichGO_1_3 <- setReadable(enrichGO.1.3, 'org.Hs.eg.db')


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN.1.3 <- enrichDGN(names(genes.c1.3.List))
enrichDGN_1_3 <- setReadable(enrichDGN.1.3, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO.1.3 <- enrichDO(gene          = names(genes.c1.3.List),
                         ont           = "DO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         minGSSize     = 5,
                         maxGSSize     = 500,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
enrichDO_1_3 <- setReadable(enrichDO.1.3, 'org.Hs.eg.db')

```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=12}
DotPlot(ctrl, features=rownames(head(markers_1.vs.3[order(markers_1.vs.3$avg_logFC, decreasing = TRUE),],50)), idents = c('1', '3')) + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 3 (CTRL dataset)", x = "Markers", y = "Clusters")
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(markers_1.vs.3,
                lab = rownames(markers_1.vs.3),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 1 vs 3 (CTRL dataset)", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=6}

#barplot(enrichGO.1.3, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 1 vs 3 (CTRL dataset)", x = "GO terms", y = "Genes related to GO term ")

# Heatmap-like functional classification
heatplot(enrichGO_1_3, foldChange=genes.c1.3.List, showCategory=50) + labs(title = "Enriched GO terms \n in cluster 1 vs 3 (CTRL dataset)", x = "Genes", y = "GO terms")

#heatplot(enrichDO_1_3, foldChange=genes.c1.3.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015) in cluster 1 vs 3 (CTRL dataset)", x = "Genes", y = "Disease Ontology")


heatplot(enrichDGN_1_3, foldChange=genes.c1.3.List, showCategory=50) + labs(title = "Gene-disease association   \n (DisGeNET, Janet et al. 2015) \n in cluster 1 vs 3 (CTRL dataset)", y = "Diseases", x = "Genes")
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=3, fig.width=6}

DefaultAssay(ctrl) <- "RNA"
ctrl <- NormalizeData(ctrl)
ctrl <- SetIdent(ctrl, value = "seurat_clusters")

FeaturePlot(ctrl, features = c("PTX3", "TFF3"), label = F)

DotPlot(ctrl, features = c("PTX3", "TFF3"), idents = c('1', '0', '3'), group.by = "pat") + RotatedAxis() + labs(title = "Specific markers in cluster 1 vs 3", x = "Gene", y = "Patients by cluster")


```
#### KEGG over-representation test

```{r}

enrichKEGG_1_3 <- enrichKEGG(gene         = names(genes.c1.3.List),
                             organism     = 'hsa',
                             pvalueCutoff = 0.05)
enrichKEGG_1_3 <- data.frame(enrichKEGG_1_3)
enrichKEGG_1_3$cluster <- rep("cluster1.vs.3", nrow(enrichKEGG_1_3))
View(data.frame(enrichKEGG_1_3))

## Output images for a single significant KEGG pathway
pathview(gene.data = genes.c1.3.List,
         pathway.id = "hsa05205",
         species = "hsa",
         limit = list(gene = max(abs(genes.c1.3.List)), # value gives the max/min limit for foldchanges
                      cpd = 1))


## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
  pathview(gene.data = genes.c1.3.List, pathway.id = enrichKEGG_1_3$ID[x], species = "hsa", 
           limit = list(gene = max(abs(genes.c1.3.List)), cpd = 1))
}

purrr::map(1:length(enrichKEGG_1_3$ID), get_kegg_plots)

```


# Marking differences between clusters 3 and 1

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=12}
DefaultAssay(ctrl) <- "RNA"
ctrl <- NormalizeData(ctrl)
ctrl <- SetIdent(ctrl, value = "seurat_clusters")

markers_3.vs.1 <- FindMarkers(ctrl, ident.1 = 3, ident.2 = 1,   min.pct = 0.25, only.pos = F)

DotPlot(ctrl, features=rownames(head(markers_3.vs.1[order(markers_3.vs.1$avg_logFC, decreasing = TRUE),],50)), idents = c('1', '3')) + RotatedAxis() + labs(title = "Specific markers in cluster 3 vs 1 (CTRL dataset)", x = "Markers", y = "Clusters")



gprofiler_3.1 <- gprofiler2::gconvert(query = rownames(markers_3.vs.1), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
gprofiler_3.1 <- gprofiler_3.1[which(duplicated(gprofiler_3.1$target) == F), ]
markers_3.vs.1$SYMBOL <- rownames(markers_3.vs.1)
genes.c3.1 <- bitr(gprofiler_3.1$target, fromType = "ENSEMBL",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

## Remove any Entrez duplicates
genes.c3.1 <- genes.c3.1[which(duplicated(genes.c3.1$SYMBOL) == F), ]
genes.c3.1.List <- markers_3.vs.1$avg_logFC
names(genes.c3.1.List) <- as.character(genes.c3.1$ENTREZID)
## Remove any NA values
genes.c3.1.List <- genes.c3.1.List[!is.na(names(genes.c3.1.List))]
genes.c3.1.List <- genes.c3.1.List[abs(genes.c3.1.List) >= 1]
genes.c3.1.List <- sort(genes.c3.1.List, decreasing = TRUE)
genes.c3.1.List



#data(geneList, package="DOSE")
enrichGO.3.1 <- enrichGO(gene = names(genes.c3.1.List),
                         #universe      = names(geneList),
                         keyType = "ENTREZID",
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE,
                         pool = TRUE)
enrichGO_3_1 <- setReadable(enrichGO.3.1, 'org.Hs.eg.db')


# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN.3.1 <- enrichDGN(names(genes.c3.1.List))
enrichDGN_3_1 <- setReadable(enrichDGN.3.1, 'org.Hs.eg.db')


# DOSE (Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis.

enrichDO.3.1 <- enrichDO(gene          = names(genes.c3.1.List),
                         ont           = "DO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         minGSSize     = 5,
                         maxGSSize     = 500,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
enrichDO_3_1 <- setReadable(enrichDO.3.1, 'org.Hs.eg.db')




```

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(markers_3.vs.1,
                lab = rownames(markers_3.vs.1),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 3 vs 1 (CTRL dataset)", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=6}

#barplot(enrichGO.3.1, showCategory=50) + labs(title = "Number of genes associated with the first 50 GO terms \n in cluster 3 vs 1 (CTRL dataset)", x = "GO terms", y = "Genes related to GO term ")

# Heatmap-like functional classification
heatplot(enrichGO_3_1, foldChange=genes.c3.1.List, showCategory=50) + labs(title = "Enriched GO terms \n in in cluster 3 vs 1 (CTRL dataset)", x = "Genes", y = "GO terms")

#heatplot(enrichDO_3_1, foldChange=genes.c3.1.List, showCategory=50) + labs(title = "Gene-disease association (Disease Ontology, Yu et al. 2015)\n in cluster 3 vs 1 (CTRL dataset)", x = "Genes", y = "Disease Ontology")


heatplot(enrichDGN_3_1, foldChange=genes.c3.1.List, showCategory=50) + labs(title = "Gene-disease association \n  (DisGeNET, Janet et al. 2015) \n   in cluster 3 vs 1 (CTRL dataset)", y = "Diseases", x = "Genes")
```

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=4, fig.width=6}

FeaturePlot(ctrl, features = c("IFI27", "RGS5", "RNASE1"), label = F)

DotPlot(ctrl, features = c("IFI27", "RGS5", "RNASE1"), group.by = "pat", idents = c('0', '1', '3')) + RotatedAxis() + labs(title = "Specific markers in cluster 3 vs 1", x = "Gene", y = "Patients cluster")


```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=5, fig.width=10}
EnhancedVolcano(markers_3.vs.1,
                lab = rownames(markers_3.vs.1),
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Differential Gene Expression",
                subtitle = "Cluster 3 vs 1 in CTRL dataset", 
                pointSize = 2.0,
                labSize = 3.5,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()

```


#### KEGG over-representation test

```{r}

enrichKEGG_3_1 <- enrichKEGG(gene         = names(genes.c3.1.List),
                             organism     = 'hsa',
                             pvalueCutoff = 0.05)
enrichKEGG_3_1 <- data.frame(enrichKEGG_3_1)
enrichKEGG_3_1$cluster <- rep("cluster3.vs.1", nrow(enrichKEGG_3_1))
View(data.frame(enrichKEGG_3_1))

## Output images for a single significant KEGG pathway
# pathview(gene.data = genes.c1.3.List,
#               pathway.id = "hsa05144",
#               species = "hsa",
#               limit = list(gene = max(abs(genes.c3.1.List)), # value gives the max/min limit for foldchanges
#               cpd = 1))


## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
  pathview(gene.data = genes.c1.3.List, pathway.id = enrichKEGG_3_1$ID[x], species = "hsa", 
           limit = list(gene = max(abs(genes.c3.1.List)), cpd = 1))
}

purrr::map(1:length(enrichKEGG_1_3$ID), get_kegg_plots)

```








## DGE HTAP vs CTRL using Limma packages



### Limma



```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=6, fig.width=8}

DefaultAssay(htap.ctrl) <- "RNA"
htap.ctrl <- NormalizeData(htap.ctrl)
Idents(htap.ctrl) <- "cluster.ids"
sce.htap.ctrl <- as.SingleCellExperiment(htap.ctrl)

agg.clust <- aggregateAcrossCells(sce.htap.ctrl, ids = colData(sce.htap.ctrl)[,c("cluster.ids")])
agg.clust <- assay(agg.clust)
agg.clust <- round(agg.clust/rowSums(agg.clust)*100, 3)
agg.clust[is.na(agg.clust)] = 0.000
agg.clust <- data.frame(agg.clust)
colnames(agg.clust) <- c("cluster 0", "cluster 1", "cycling cells", "cluster 3")
agg.clust$gene <- rownames(agg.clust)
head(agg.clust)


agg.patient <- aggregateAcrossCells(sce.htap.ctrl, ids = colData(sce.htap.ctrl)[,c("PATIENT")])
agg.patient <- assay(agg.patient)
agg.patient <- round(agg.patient/rowSums(agg.patient)*100, 3)
agg.patient[is.na(agg.patient)] = 0.000
agg.patient <- data.frame(agg.patient)
agg.patient$gene <- rownames(agg.patient)
#colnames(agg.patient) <- c("cluster0", "cluster1", "cluster2", "cluster3", "cluster4")
head(agg.patient)

agg.clust.patient <- merge(agg.clust, agg.patient, by = "gene")




agg <- aggregateAcrossCells(sce.htap.ctrl, ids = colData(sce.htap.ctrl)[,c("cluster.ids", "PHENOTYPE")])
#colData(agg)[,c("seurat_clusters", "PHENOTYPE",  "ncells")]

group <- factor(agg$PHENOTYPE)
dge <- DGEList(assay(agg), group = group)

# keep <- filterByExpr(dge, group=group)
# dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# plotMD(cpm(dge, log=TRUE), column=1)
# abline(h=0, col="red", lty=2, lwd=2)

# Desing matrix HTAP and CTRL conditions
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
con <- makeContrasts(HTAP.vs.CTRL = HTAP - CTRL,
                     #CTRL.vs.HTAP = CTRL - HTAP,
                     levels=colnames(design))

v <- voom(dge, design, plot=F)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=con)
efit <- eBayes(vfit, robust = TRUE)
topTable <- limma::topTable(efit, number = "Inf")
topTable$gene <- rownames(topTable)
View(topTable)
plotMD(efit, main = "HTAP vs CTRL")


limma.topTable.HTAP.vs.CTRL <- merge(agg.clust.patient, topTable, by = "gene")

#View(limma.topTable.HTAP.vs.CTRL)

write.table(limma.topTable.HTAP.vs.CTRL, file="./limma.HTAP.vs.CTRL.csv", sep="\t")


# Volcanoplot

library("EnhancedVolcano")


EnhancedVolcano(topTable,
                lab = rownames(topTable),
                x = 'logFC',
                y = 'P.Value',
                #xlim = c(-8, 8),
                ylim = c(0, -log10(10e-12)),
                title = "Differential Gene Expression \n
    HTAP vs CTRL",
                subtitle = "(linear model by limma package, 
    robust = true, FCcutoff = 1.5, pCutoff = 10e-05)" ,
                #selectLab = top.c8.genes,
                #xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.5,
                #labCol = 'black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                #colConnectors = 'black'
) + ggplot2::theme_linedraw()


# Maplot
# maplot <- topTable[, c("AveExpr", "logFC", "adj.P.Val")]
# colnames(maplot) <- c("baseMean", "log2FoldChange", "padj")
# maplot <- maplot[!apply(is.na(maplot) | maplot == "", 1, all),]
# 
# #library("ggpubr")
# ggpubr::ggmaplot(maplot,
#                  main = "MA plot of HTAP vs CTRL (FDR 5%)",
#                  fdr = 0.05,
#                  fc = 1.5,
#                  size = 1.5,
#                  palette = c("#B31B21", "#1465AC", "darkgray"),
#                  genenames = as.vector(rownames(maplot)),
#                  legend = "right",
#                  top = 15,
#                  alpha = 1,
#                  xlab = "Log2 mean expression",
#                  ylab = "Log2 fold change",
#                  select.top.method = c("fc"),
#                  font.label = c("bold", 6),
#                  label.rectangle = TRUE,
#                  font.legend = "bold",
#                  #font.main = "bold",
#                  ggtheme = ggplot2::theme_classic())





#counts <- seurat.htap.ctrl@assays$RNA@counts


# Molecular function of the interest genes in HTAP vs CTRL
library("gprofiler2")

gprofiler2.htap.vs.ctrl <- gprofiler2::gconvert(query = rownames(topLimma), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
View(gprofiler2.htap.vs.ctrl)
# gprofiler2.htap.vs.ctrl <- gprofiler2.htap.vs.ctrl[, -c(1:3)]
# gprofiler2.htap.vs.ctrl <- gprofiler2.htap.vs.ctrl[, c("name","target", "description")]
#View(gprofiler2.htap.vs.ctrl)


# Gene list functional enrichment analysis with gost

# gostres <- gost(query = top.genes.clusters, 
#                 organism = "hsapiens", ordered_query = FALSE, 
#                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
#                 measure_underrepresentation = FALSE, evcodes = TRUE, 
#                 user_threshold = 0.05, correction_method = "g_SCS", 
#                 domain_scope = "annotated", custom_bg = NULL, 
#                 numeric_ns = "", sources = NULL, as_short_link = FALSE)
# 
# 
# names(gostres)
# #View(gostres$result)
# 
# p <- gostplot(gostres, capped = TRUE, interactive =  FALSE)
# 
# publish_gostplot(p, highlight_terms = gostres$result$term_id, 
#                        width = 8, height = 20, filename = NULL )


```

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=12, fig.width=8}
top <- topTable[abs(topTable$AveExpr) >= 1.5 & topTable$P.Value <= 10e-05,]
#dim(top)
View(top)
#top <- head(top,50)
DotPlot(htap.ctrl,features=rownames(head(top[order(top$AveExpr, decreasing = T),],50)), dot.scale = 8, group.by = "pat") + RotatedAxis() + labs(title = "HTAP vs CTRL", x = "Markers", y = "Clusters")

```

```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=10, fig.width=10}
DefaultAssay(htap.ctrl) <- "RNA"
htap.ctrl <- NormalizeData(htap.ctrl)
Idents(htap.ctrl) <- "cluster.ids"
FeaturePlot(htap.ctrl, features = rownames(top))
```


```{r message=FALSE, warning=FALSE, echo=FALSE, eval=T, fig.height=10, fig.width=12}
top.xy <- topTable[abs(topTable$logFC) > 1.0 & topTable$P.Value <= 0.05,]
top.xy$SYMBOL <- rownames(top.xy)
#colnames(topTable.xx) <- c("logFC","AveExpr","t" , "P.Value" ,"adj.P.Val", "B", "SYMBOL")
topTable.xy.df <- bitr(top.xy$SYMBOL, fromType = "SYMBOL",
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = org.Hs.eg.db)

v <- merge(topTable.xy.df, top.xy, by="SYMBOL")
## Remove any Entrez duplicates
v <- v[which(duplicated(v$ENTREZID)==F), ]
genes.xx <- v$logFC
names(genes.xx) <- as.character(v$ENTREZID)
genes.xx <- sort(genes.xx, decreasing = TRUE)
genes.xx <- na.exclude(genes.xx)



# Enrichment analysis of ADA2 and NOMO1
sig <- subset(topTable,
              gene == c("NOMO1"))
colnames(sig) <- c("logFC","AveExpr","t" , "P.Value" ,"adj.P.Val", "B", "SYMBOL")
sig.df <- bitr(sig$SYMBOL, fromType = "SYMBOL",
               toType = c("ENSEMBL", "ENTREZID"),
               OrgDb = org.Hs.eg.db)
d <- merge(sig.df, sig, by="SYMBOL")
## Remove any Entrez duplicates
d <- d[which(duplicated(d$ENTREZID)==F), ]
nomo1 <- d$logFC
names(nomo1) <- as.character(d$ENTREZID)
nomo1 <- sort(nomo1, decreasing = TRUE)

enrichGO_htap.ctrl <- enrichGO(gene          = unique(sig.df$ENTREZID),
                               #universe      = names(ada2_nomo1_DOSE),
                               #keyType = "ENSEMBL",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE,
                               pool = TRUE)
enrichGO_HTAP.vs.CTRL <- setReadable(enrichGO_htap.ctrl, 'org.Hs.eg.db')
heatplot(enrichGO_HTAP.vs.CTRL, foldChange=nomo1) + ggtitle("NOMO1 Gene Ontology")


emapplot(enrichGO_HTAP.vs.CTRL, showCategory=50, pie_scale=0.8,layout="kk", line_scale = 0.2) + ggtitle("Relationship between the top 50 most significantly \n enriched GO terms in HTAP vs CTRL")
enrichGO_HTAP.vs.CTRL <- data.frame(enrichGO_HTAP.vs.CTRL)
enrichGO_HTAP.vs.CTRL$Gene_expr_analyse <- rep("HTAP.vs.CTRL", nrow(enrichGO_HTAP.vs.CTRL))
View(enrichGO_HTAP.vs.CTRL)



# Disease analysis
# DisGeNET(Janet et al. 2015)
enrichDGN_HTAP.vs.CTRL <- enrichDGN(names(nomo1))
enrichDGN_HTAP.vs.CTRL <- setReadable(enrichDGN_HTAP.vs.CTRL, 'org.Hs.eg.db')
enrichDGN_HTAP.vs.CTRL <- data.frame(enrichDGN_HTAP.vs.CTRL)
enrichDGN_HTAP.vs.CTRL$Database <- rep("DisGeNET (Janet et al. 2015)", nrow(enrichDGN_HTAP.vs.CTRL))
heatplot(enrichDGN_HTAP.vs.CTRL, foldChange=nomo1) + ggtitle("Disease-gene association (DisGeNET, Janet et al. 2015)")















enrichGO <- enrichGO(gene          = names(nomo1),
                     #universe      = names(geneList_DOSE),
                     #keyType = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE,
                     pool = TRUE)

head(enrichGO, n=10)
View(data.frame(enrichGO_HTAP.vs.CTRL))

enrichGO_HTAP.vs.CTRL <- setReadable(enrichGO, 'org.Hs.eg.db')
heatplot(enrichGO_HTAP.vs.CTRL, foldChange=genes.xx) + ggtitle("GO in HTAP vs CTRL")


## GSEA using gene sets from KEGG pathways
genes <- names(genes.xx)[abs(genes.xx) >= 1]
genes <- sort(genes, decreasing = TRUE)
kk <- enrichKEGG(gene         = genes,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)




# gseaKEGG <- gseKEGG(geneList = genes.xx, # ordered named vector of fold changes (Entrez IDs are the associated names)
#                     organism = "hsa", # supported organisms listed below
#                     nPerm = 1000, # default number permutations
#                     minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
#                     pvalueCutoff = 0.05, # padj cutoff value
#                     verbose = FALSE)
# 
# ## Extract the GSEA results
# gseaKEGG_results <- gseaKEGG@result
# View(gseaKEGG_results)
```
