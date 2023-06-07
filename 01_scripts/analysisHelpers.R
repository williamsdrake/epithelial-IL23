downsampleSeurat <- function(seuratObj,
                             target_cell = 10000) 
{
  # randomly select 'target_cell' number of cells from the passed seurat object
  seuratObj <- seuratObj[, sample(colnames(seuratObj), size = target_cell, replace = F)]

  return(seuratObj)
  
  
}

convertEnsemblRownames <- function(seu) {
  # obtain the genes in the current object
  genes <- as.data.frame(rownames(seu))
  # convert ENSEMBL to SYMBOL
  master <- mapIds(org.Hs.eg.db, 
                   keys = genes$`rownames(seu)`, 
                   keytype = "ENSEMBL", 
                   column="SYMBOL")
  master <- as.data.frame(master)
  # match ENSEMBL to SYMBOL
  genes$geneID <- master$master
  # remove ENSEMBL genes that do not map to SYMBOLs
  seu <- seu[!is.na(genes$geneID),]
  genes <- genes[!is.na(genes$geneID),]
  # remove duplicate SYMBOLs
  dis <- duplicated(genes$geneID)
  seu <- seu[!dis,]
  genes <- genes[!dis,]
  # rename rows in appropriate slots (note: scale.data slot is empty since we have not done any scaling so do not rename them here)
  rownames(seu@assays$RNA@counts) <- genes$geneID
  rownames(seu@assays$RNA@data) <- genes$geneID
  return(seu)
}

getBulkExpression <- function(seu,
                              cellType,
                              donor,
                              disease,
                              diseaseLevels){
  bulk <- AggregateExpression(seu, 
                              return.seurat = T, 
                              slot = "counts", 
                              assays = "RNA", 
                              group.by = c(cellType,
                                           donor,
                                           disease))
  
  bulk$celltype <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
  bulk$donor <- sapply(strsplit(Cells(bulk), split = "_"), "[", 2)
  bulk$disease <- sapply(strsplit(Cells(bulk), split = "_"), "[", 3)
  bulk$disease <- factor(x = bulk$disease, levels = diseaseLevels)
  return(bulk)
}

# Data came in the format of an individual RDS
# Function to load all data in the folder into a list of seurat objects
read_rds_files <- function(folder) {
  # Get a list of all RDS files in the folder
  rds_files <- list.files(folder, pattern = "\\.rds$")
  
  # Add the path to each RDS file
  rds_files <- file.path(folder, rds_files)
  
  # Read each RDS file and create a Seurat object
  seurat_objects <- lapply(rds_files, readRDS)
  
  # Return a list of Seurat objects, where the name of each object is the name of the file
  return(setNames(seurat_objects, rds_files))
}

