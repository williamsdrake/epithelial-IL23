downsampleSeurat <- function(seuratObj,
                             target_cell = 10000) 
{
  # randomly select 'target_cell' number of cells from the passed seurat object
  seuratObj <- seuratObj[, sample(colnames(seuratObj), size = target_cell, replace = F)]

  return(seuratObj)
  
  
}