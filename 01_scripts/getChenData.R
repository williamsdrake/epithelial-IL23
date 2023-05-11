library(GEOquery)
library(scCustomize)
library(dplyr)
chen2022 <- "GSE171213"
projDir <- "/home/williamsdrw/Immunodeficient/02_raw_data/"

# get GEO information
gse <- getGEO(chen2022)

# get GEO files
getGEOSuppFiles(GEO = chen2022, baseDir = myDir)
myDir <- paste0(myDir,chen2022,"/")
files <- list.files(myDir) %>%
  .[grepl("tar", .)]

#untar files
untar(paste0(myDir,files),
      exdir = paste0(myDir))
      
# make a list of files that have RNA counts
files <- list.files(myDir) %>%
  .[!grepl("cellname|GSE171213", .)]


# get/tidy metadata from GEO object to match other sample metadata
metadata <- as.data.frame(pData(gse[[1]])) %>%
  .[c("title", "extract_protocol_ch1.1", "age:ch1", "gender:ch1", "group:ch1", "tissue:ch1")]
colnames(metadata) <- c("orig.ident", "tech", "age", "sex", "status", "tissue")
statusChange <- unique(metadata$status)
metadata[metadata == statusChange[1]] <- "H"
metadata[metadata == statusChange[2]] <- "P"
metadata[metadata == statusChange[3]] <- "P_treated"
metadata$sex <- substr(metadata$sex, 0, 1)
metadata$tissue <- "gingiva"
metadata$tech <- "BD Rhapsody"
metadata$source <- chen2022
metadata$ref <- "Chen 2022"
write.csv(metadata, file = paste0("/home/williamsdrw/Immunodeficient/03_meta_data/", Sys.Date(), "_GSE171213_metadata.csv"))

seurat_list <- list()

for(i in seq_along(files)){
  tab <- read.table(paste0(myDir, files[i]), 
                    header = T, 
                    sep = "\t", 
                    row.names = 1, 
                    check.names = T)
  seurat_list[i] <- CreateSeuratObject(tab,
                                       project = metadata$orig.ident[i])
}

chen_srat <- Merge_Seurat_List(list_seurat = seurat_list)
# Create a label for this tissue only
saveRDS(chen_srat, file = paste0("/home/williamsdrw/Immunodeficient/04_data_objects/01_SoupX_objects/", Sys.Date(), "_GSE171213_withMeta.RDS"))
