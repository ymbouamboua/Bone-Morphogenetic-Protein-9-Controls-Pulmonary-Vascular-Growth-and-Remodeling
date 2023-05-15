# Violin plot for muliple genes with comparison test : wilcox

violin.plot.wilcox <- function(object, group.by, gene.signature, file.name, test.sign){
  
  require(Seurat)
  require(ggplot2)
  
  if ("Seurat" != class(object)[1]) {
    stop("object should be of class Seurat")
  }
  
  if (!missing(x = group.by)) {
    object <- SetIdent(object = object, value = group.by)
  }
  
  plot <- function(signature, y.max = NULL){
    p <- VlnPlot(object = object, 
                 group.by = group.by, 
                 features = signature,
                 cols = c("#88CCEE","#117733"),
                 pt.size = 0, 
                 y.max = y.max)
    
    p <- p + theme( plot.title = element_text(hjust = 0.5, size = 8, face="bold.italic"),
                    text =  element_text(size = 5, family = "Helvetica"), 
                    legend.position = "none",
                    panel.grid.major = ggplot2::element_blank(), 
                    panel.grid.minor = ggplot2::element_blank(),
                    panel.background = ggplot2::element_blank(),
                    axis.line =   element_line(colour = "black"), 
                    axis.ticks =  element_line(colour = "black", size = 5/12), 
                    axis.title =  element_text(face = "plain"), 
                    axis.title.x = element_blank(),
                    axis.text =   element_text(size = 5), 
                    axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) 
    
    p <- p + stat_compare_means(comparisons = test.sign, label = "p.signif", size = 3)
    
    #labels <- c(expression("ALK1"^italic("Low")), expression("ALK1"^italic("High")))
    
    #p + scale_x_discrete(labels = labels)
    
  }
  
  plot.list <- list()
  y.max.list <- list()
  
  for (gene in gene.signature) {
    plot.list[[gene]] <- plot(gene)
    y.max.list[[gene]] <- max(plot.list[[gene]]$data[[gene]])
    plot.list[[gene]] <- plot(gene, y.max = (y.max.list[[gene]] + 1) )
  }
  
  cowplot::plot_grid(plotlist = plot.list, ncol = 5)
  file.name <- paste0(file.name, "_rplot.pdf")
  ggsave(file.name, width = 7, height = 4)
}

# comparisons <- list(c("celltype1", "celltype2"))
# 
# violin.plot.wilcox(object = sobj,
#                    group.by = "celltype",
#                    gene.signature = gene.signature, 
#                    file.name = "seurat",
#                    test.sign = comparisons)
