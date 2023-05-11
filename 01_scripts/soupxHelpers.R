desoupSingle <- function(dir_raw,
                     dir_filt,
                     varFeat = 4000,
                     num_pc = 50,
                     pc_dim = 1:50) {
  # follow desouping vignette
  tod <- Seurat::Read10X(data.dir = dir_raw)
  toc <- Seurat::Read10X(data.dir = dir_filt)
  sc <- SoupX::SoupChannel(tod, toc)
  srat <- Seurat::CreateSeuratObject(sc$toc)
  srat <- Seurat::NormalizeData(srat)
  srat <- Seurat::FindVariableFeatures(srat, nfeatures = varFeat)
  # Cell-cycle scoring and regression
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  srat <- Seurat::CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes,
                           set.ident = TRUE)
  srat <- Seurat::ScaleData(srat, vars.to.regress = c("S.Score", "G2M.Score"), verbose = TRUE)
  srat <- Seurat::RunPCA(srat, npcs = num_pc, verbose = TRUE)
  srat <- Seurat::FindNeighbors(srat, reduction = "pca", dims = pc_dim)
  srat <- Seurat::FindClusters(srat, resolution = 1)
  sc <- SoupX::setClusters(sc,srat@active.ident)
  sc <- SoupX::autoEstCont(sc)
  return(SoupX::adjustCounts(sc))
}

desoupGroup <- function(dir_vector,
                        varFeat = 4000,
                        num_pc = 50,
                        pc_dim = 1:50) {
  desoupedList <- list()
  raw_vector <- dir_vector[grepl("raw_feature_bc_matrix|raw_gene_bc_matrices", dir_vector)]
  filt_vector <- dir_vector[grepl("filtered_gene_bc_matrices|filtered_feature_bc_matrix", dir_vector)]
  for(i in seq_along(raw_vector)){
    print(paste0("Desouping sample ", i, " of ", length(raw_vector), "."))
    desoupedList[i] <- desoupSingle(raw_vector[i], filt_vector[i])
    print(paste0("Sample ", i, " of ", length(raw_vector), " desouped."))
  }
  print("Desoup complete.")
  return(desoupedList)
}

dcgToSeurat <- function(matrixList,
                        dir_vector) 
  {
  seuratList <- list()
  # extract patient ID from file path
  dir_vector <- substring(dir_vector,57)
  print(dir_vector)
  dir_vector <- unique(sub("/.*", "", dir_vector))
  print(dir_vector)
  for(i in seq_along(matrixList))
    {
    seuratList[i] <- CreateSeuratObject(unlist(matrixList[[i]]), project = dir_vector[i])
    }
  return(seuratList)
}

