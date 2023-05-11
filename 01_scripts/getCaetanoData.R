library(GEOquery)
library(scCustomize)
library(dplyr)
caetano2021 <- "GSE152042"
myDir <- "/home/williamsdrw/Immunodeficient/02_raw_data/"

# get GEO information
gse <- getGEO(caetano2021)

# get GEO files
getGEOSuppFiles(GEO = caetano2021, baseDir = myDir)
myDir <- paste0(myDir,caetano2021,"/")
files <- list.files(myDir) %>%
  .[grepl("tar", .)]

#untar files
untar(paste0(myDir,files),
      exdir = paste0(myDir))

# make a list of files that have RNA counts
files <- list.files(myDir) %>%
  .[!grepl("GSE", .)]


# get/tidy metadata from GEO object to match other sample metadata
metadata <- as.data.frame(pData(gse[[1]])) %>%
  .[c("title", "extract_protocol_ch1.1", "disease state:ch1")]
colnames(metadata) <- c("orig.ident", "tech", "status")
statusChange <- unique(metadata$status)
metadata[metadata == statusChange[1]] <- "H"
metadata[metadata == statusChange[2]] <- "MP"
metadata[metadata == statusChange[3]] <- "P"
metadata$sex <- c("M", "M", "F", "M")
metadata$tissue <- "gingiva"
metadata$source <- caetano2021
metadata$ref <- "Caetano 2021"
metadata$age <- "41-65"
write.csv(metadata, file = paste0("/home/williamsdrw/Immunodeficient/03_meta_data/", Sys.Date(), "_GSE152042_metadata.csv"))

seurat_list <- list()

for(i in seq_along(files)){
  temp <- Read10X_h5(paste0(myDir, files[i]))
  seurat_list[i] <- CreateSeuratObject(temp,
                                       project = metadata$orig.ident[i])
}

caetano_srat <- Merge_Seurat_List(list_seurat = seurat_list)
caetano_srat <- Add_Sample_Meta(caetano_srat, meta_data = metadata, join_by_meta = "orig.ident", join_by_seurat = "orig.ident")
caetano_srat@meta.data

# Create a label for this tissue only
saveRDS(caetano_srat, file = paste0("/home/williamsdrw/Immunodeficient/04_data_objects/01_SoupX_objects/", Sys.Date(), "_GSE152042_withMeta.RDS"))
