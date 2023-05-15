asDataframeSeurat <- function(x, genes = Seurat::VariableFeatures(x), fix_names = TRUE, ...) {
  # TODO possibly also a warning if it is not a variable gene and an error
  # if it does not exist also an argument to force though the error ...
  tmp <- Seurat::FetchData(x, vars = genes, ...)
  tmp <-  as.data.frame(tmp, stringsAsFactors = FALSE)
  tmp$ident <- Seurat::Idents(x)[rownames(tmp)]
  
  if (fix_names) {
    colnames(tmp) <- make.names(colnames(tmp))
  }
  return(tmp)
}
