seuratQC <- function(seuratObj,
                     feat_min = 200, 
                     feat_max = 5000,
                     pct_mt = 15,
                     study_group = "Group_A") 
{
  
  # Store a dataframe with the number of cells/sample before filtering
  preQC_data <- as.data.frame(table(seuratObj$orig.ident))
  
  # Generate mitochrondial gene metadata
  seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj,
                                                    pattern = "^MT-")
  
  # Make a set of violin plots that show the number of genes (features), the number of raw RNA counts,
  # and the percentage of mitochondrial genes per cell (pre QC)
  preQC_vln <- VlnPlot(seuratObj,
                       features = c("nFeature_RNA",
                                    "nCount_RNA",
                                    "percent.mt"),
                       ncol = 3,
                       group.by = "orig.ident",
                       pt.size = 0) +
    theme(axis.title.x = element_blank())
  
  # Filter out bad cells (fewer than e.g. 200 genes expressed [dying cell or RBC], greater than e.g. 5000 genes [may be a doublet],
  # and more than e.g. 15% mitochondrial gene expression [probably a dying cell])
  seuratObj <- subset(seuratObj, 
                      subset = nFeature_RNA > feat_min & 
                        nFeature_RNA < feat_max & 
                        percent.mt < pct_mt)
  
  # Store a dataframe of the number of cells/sample postQC
  postQC_data <- as.data.frame(table(seuratObj$orig.ident))
  
  # Make a set of violin plots that show the number of genes (features), the number of raw RNA counts,
  # and the percentage of mitochondrial genes per cell (post QC)
  postQC_vln <- VlnPlot(seuratObj,
                        features = c("nFeature_RNA",
                                     "nCount_RNA",
                                     "percent.mt"),
                        ncol = 3,
                        group.by = "orig.ident",
                        pt.size = 0) +
    theme(axis.title.x = element_blank())
  
  filename = paste0("~/updatedAtlas-hu/05_plots/01_QC_plots/", 
                    Sys.Date(),
                    "_",
                    study_group,
                    "_featureSummary.pdf")
  
  # Save the QC feature plots
  CairoPDF(file = filename)
  print(preQC_vln/postQC_vln)
  dev.off()
  
  
  # Print message after file save
  print(preQC_vln/postQC_vln)
  print(paste0("Feature summary plot is located at ", 
               filename))
  
  # Merge pre and postQC dataframes, format for graphing
  qcSum <- merge(preQC_data, postQC_data, by = "Var1")
  names(qcSum)[1] <- "Patient"
  names(qcSum)[2] <- "preQC"
  names(qcSum)[3] <- "postQC"
  qcSum <- reshape2::melt(qcSum[,c('Patient','preQC','postQC')],id.vars = 1)
  names(qcSum)[3] <- "Cells"
  
  # Generate plot to visualize how many cells were removed by QC
  cellNumSummary <- ggplot(qcSum, aes(x = Patient, y = Cells)) + 
    geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
    theme_cowplot() +
    theme(legend.title = element_blank(), 
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_blank()) + 
    scale_fill_manual(values = c("grey", "black"))
  filename = paste0("~/updatedAtlas-hu/05_plots/01_QC_plots/", 
                    Sys.Date(),
                    "_",
                    study_group,
                    "_cellNumSummary.pdf")
  
  CairoPDF(file = filename)
  print(cellNumSummary)
  dev.off()
  
  # Print message after file save
  print(cellNumSummary)
  print(paste0("Plot showing number of cells removed by QC is located at ", 
               filename))
  
  # Print message to signal QC completion
  print(paste0("QC for ", study_group, " is complete."))
  
  return(seuratObj)
  
  
}