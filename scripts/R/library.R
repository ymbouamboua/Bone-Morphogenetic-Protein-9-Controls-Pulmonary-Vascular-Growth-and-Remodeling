
library(Seurat)
library(SeuratWrappers)
library(dittoSeq)
library(dplyr)
library(plyr)
library(pheatmap)
library(velocyto.R)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(grid)
#library(DoubletFinder)
library(gplots)
#library(fgsea)
#library(biomaRt)
library(knitr)
library(xtable)
#library(ggrepel)
#library(matchSCore2)
library(cowplot)
library(scater)
library(scran)
library(SingleCellExperiment)
#library(DropletUtils)
library(batchelor)
library(BiocNeighbors)
library(SeuratWrappers)
library(edgeR)
library(limma)
library(EnhancedVolcano)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(Nebulosa)
library(monocle3)
library(ReactomeGSA)

# custom.pal <- c("grey","#FF9B00","#EC380B")
# 
# cc.genes <- readLines(con = "/data/10x_data/cell_cycle_vignette_files/regev_lab_cell_cycle_genes_mouse.txt")
# #We can segregate this list into markers of G2/M phase and markers of S phase
# s.genes <- cc.genes[1:43]
# g2m.genes <- cc.genes[44:97]
# 
# ggplotColours <- function(n = 6, h = c(0, 360) + 15){
#   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
#   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
# }
# col=ggplotColours(n = 10)
# 
# cell_type_color <- c('Cycling Basal' = '#F3766E',
#                      'Basal' = '#A31E22',
#                      'Suprabasal' = '#FCCC0A',
#                      'Secretory like' = '#a9c653',
#                      'Secretory' = '#53c653',
#                      'Pneumo precursor' = '#a7d6a7',
#                      'Goblet' = '#006600',
#                      'Deuterosomal' = '#2da9d2',
#                      'Multiciliated' = '#466cb9',
#                      'Ionocytes' = '#b34d6d',
#                      'Neuro-endocrine' = '#904fd1',
#                      'Glandulaire' = '#3e937f',
#                      'Fibroblast' = '#aea9ce',
#                      'Smooth muscle' = '#7c5959',
#                      'Endothelial' = '#281c1c',
#                      'Myofibroblast' = '#a98989',
#                      'Chondrocyte' = '#9b9b98',
#                      'Immune' = '#e2bae2',
#                      'Macrophage' = '#c997e7',
#                      'LB / Plasmocytes' = '#7e6384',
#                      'Mastocytes' = '#e366ff',
#                      'Lymphocytes (NK/NKT + LT)' = '#ff1aff',
#                      'Pneumocytes' = '#fffbb7')
