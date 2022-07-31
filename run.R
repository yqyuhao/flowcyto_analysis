library(flowCore)
library(Seurat)
library(readxl)
library(cowplot)
library(ggplot2)
library(dplyr)
setwd("/data/analysis")

# define number
cell_id=0

# loading the panel information
panel_filename <- "config/Righton_CyTOF_panel.xlsx"

# Select Biomarker for using normalized expression from fcs files
panel <- read_excel(panel_filename,sheet="Panel2")

### sample list
samples <- read_excel("metadata/patients_metadata.xlsx", range = cell_cols("A:A")) %>% .$sample_id

### import fcs files from different data sets
for (i in seq_along(samples)){
  
  # read fcs data
  assign(paste0("fcs_data_", i), read.flowSet(paste0("data/", samples[i], ".fcs"),transformation = FALSE, truncate_max_range = FALSE))

  # read the panel file which provides information on the markers
  assign(paste0("panel_fcs_data_", i), pData(parameters(eval(parse(text = paste0("fcs_data_", i)))[[1]])))
  
  # Using Arcsinh transformation to normalize the fcs files
  assign(paste0("fcs_data_normalized_", i),fsApply(eval(parse(text = paste0("fcs_data_", i))),
                                                     function(x, cofactor = 6)
                                                     {
                                                       colnames(x) <-  eval(parse(text = paste0("panel_fcs_data_", i)))$name
                                                       expr <- exprs(x)
                                                       expr <- asinh(expr[,c(eval(parse(text = paste0("panel_fcs_data_", i)))$name)] / cofactor)
                                                       exprs(x) <- expr
                                                     }
  )[,c(panel$Biomarker)])
  #Create Seurat Object using normalized expression from fcs files
  assign(paste0("fcs_data_normalized_matrix_", i) , as.matrix(t(eval(parse(text = paste0("fcs_data_normalized_", i))))))
  
  Cells <- c()
  for (f in 1:ncol(eval(parse(text = paste0("fcs_data_normalized_matrix_", i)))))
  {
    cell_id=cell_id+1
    a <- paste("cell",cell_id, sep="")
    Cells <- c(Cells,a)
  }
  matrix_ID=eval(parse(text = paste0("fcs_data_normalized_matrix_", i)))
  colnames(matrix_ID) <- Cells
  assign(paste0("Obj_fcs_data_", i), CreateSeuratObject(matrix_ID, project = samples[i]))
}
remove(matrix_ID)
remove(a)
remove(cell_id)
remove(Cells)
remove(f)

library("future")
plan("multiprocess", workers = 60)
options(future.globals.maxSize = 100000 * 1024^60)
for (i in seq_along(samples)){
  
  # Variable feature selection
  Obj_fcs_data = eval(parse(text = paste0("Obj_fcs_data_", i)))
  VariableFeatures(Obj_fcs_data) <- rownames(Obj_fcs_data)
  Obj_fcs_data <- ScaleData(Obj_fcs_data)
  
  # Add metadata
  metatable <- read_excel("metadata/patients_metadata.xlsx")
  metadata <- FetchData(Obj_fcs_data, "orig.ident")
  metadata$cell_id <- rownames(metadata)
  metadata$sample_id <- metadata$orig.ident
  metadata <- left_join(x = metadata, y = metatable, by = "sample_id")
  rownames(metadata) <- metadata$cell_id
  Obj_fcs_data <- AddMetaData(Obj_fcs_data, metadata = metadata)
  
  # Dimensionality reduction
  Obj_fcs_data <- RunPCA(Obj_fcs_data, verbose = TRUE)
  Obj_fcs_data <- RunUMAP(Obj_fcs_data, dims = 1:10, verbose = T)
  Obj_fcs_data <- FindNeighbors(Obj_fcs_data, dims = 1:10)
  assign(paste0("Obj_fcs_data_dimensionality_reduction_", i) , RunTSNE(Obj_fcs_data, dims = 1:10, verbose = TRUE))
  
  # save seurat_objects
  saveRDS(eval(parse(text = paste0("Obj_fcs_data_dimensionality_reduction_", i))), file = paste0("objects/mutiplex_",samples[i],"_single.RDS"))
}

# Integration of Obj fcs data dimensionality reduction to uncover similarities and differences among cell types from each sample
Integration.anchors <- FindIntegrationAnchors(object.list = list(Obj_fcs_data_dimensionality_reduction_1, Obj_fcs_data_dimensionality_reduction_2,Obj_fcs_data_dimensionality_reduction_3,Obj_fcs_data_dimensionality_reduction_4,Obj_fcs_data_dimensionality_reduction_5,Obj_fcs_data_dimensionality_reduction_6,Obj_fcs_data_dimensionality_reduction_7,Obj_fcs_data_dimensionality_reduction_8,Obj_fcs_data_dimensionality_reduction_9,Obj_fcs_data_dimensionality_reduction_10, 
                                                                 Obj_fcs_data_dimensionality_reduction_11,Obj_fcs_data_dimensionality_reduction_12,Obj_fcs_data_dimensionality_reduction_13,Obj_fcs_data_dimensionality_reduction_14,Obj_fcs_data_dimensionality_reduction_15,Obj_fcs_data_dimensionality_reduction_16,Obj_fcs_data_dimensionality_reduction_17,Obj_fcs_data_dimensionality_reduction_18,Obj_fcs_data_dimensionality_reduction_19,Obj_fcs_data_dimensionality_reduction_20,
                                                                 Obj_fcs_data_dimensionality_reduction_21,Obj_fcs_data_dimensionality_reduction_22,Obj_fcs_data_dimensionality_reduction_23,Obj_fcs_data_dimensionality_reduction_24,Obj_fcs_data_dimensionality_reduction_25,Obj_fcs_data_dimensionality_reduction_26,Obj_fcs_data_dimensionality_reduction_27,Obj_fcs_data_dimensionality_reduction_28,Obj_fcs_data_dimensionality_reduction_29,Obj_fcs_data_dimensionality_reduction_30,
                                                                 Obj_fcs_data_dimensionality_reduction_31, Obj_fcs_data_dimensionality_reduction_32,Obj_fcs_data_dimensionality_reduction_33,Obj_fcs_data_dimensionality_reduction_34,Obj_fcs_data_dimensionality_reduction_35,Obj_fcs_data_dimensionality_reduction_36,Obj_fcs_data_dimensionality_reduction_37,Obj_fcs_data_dimensionality_reduction_38,Obj_fcs_data_dimensionality_reduction_39,Obj_fcs_data_dimensionality_reduction_40,
                                                                 Obj_fcs_data_dimensionality_reduction_41, Obj_fcs_data_dimensionality_reduction_42,Obj_fcs_data_dimensionality_reduction_43,Obj_fcs_data_dimensionality_reduction_44,Obj_fcs_data_dimensionality_reduction_45,Obj_fcs_data_dimensionality_reduction_46,Obj_fcs_data_dimensionality_reduction_47,Obj_fcs_data_dimensionality_reduction_48,Obj_fcs_data_dimensionality_reduction_49,Obj_fcs_data_dimensionality_reduction_50,
                                                                 Obj_fcs_data_dimensionality_reduction_51, Obj_fcs_data_dimensionality_reduction_52,Obj_fcs_data_dimensionality_reduction_53,Obj_fcs_data_dimensionality_reduction_54,Obj_fcs_data_dimensionality_reduction_55,Obj_fcs_data_dimensionality_reduction_56,Obj_fcs_data_dimensionality_reduction_57,Obj_fcs_data_dimensionality_reduction_58,Obj_fcs_data_dimensionality_reduction_59,Obj_fcs_data_dimensionality_reduction_60,
                                                                 Obj_fcs_data_dimensionality_reduction_61, Obj_fcs_data_dimensionality_reduction_62,Obj_fcs_data_dimensionality_reduction_63,Obj_fcs_data_dimensionality_reduction_64,Obj_fcs_data_dimensionality_reduction_65,Obj_fcs_data_dimensionality_reduction_66,Obj_fcs_data_dimensionality_reduction_67,Obj_fcs_data_dimensionality_reduction_68,Obj_fcs_data_dimensionality_reduction_69,Obj_fcs_data_dimensionality_reduction_70,
                                                                 Obj_fcs_data_dimensionality_reduction_71
                                                                 ), dims = 1:10)
Integration.combined <- IntegrateData(anchorset = Integration.anchors, dims = 1:10)
DefaultAssay(Integration.combined) <- "integrated"
Obj_fcs_data <- ScaleData(Integration.combined)
Obj_fcs_data <- RunPCA(Obj_fcs_data, verbose = TRUE)
ElbowPlot(Obj_fcs_data, ndims = 10)
ggsave2(paste0("Righton_flow_cytometry_Analysis.Integration.combined.ElbowPlot.pdf"), path = "results/Dimensionality_reduction", width = 30, height = 30, units = "cm")
Obj_fcs_data <- RunUMAP(Obj_fcs_data, reduction = "pca", dims = 1:10)
Obj_fcs_data <- FindNeighbors(Obj_fcs_data, reduction = "pca", dims = 1:10)

# Resolution Analysis
for (j in c(0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,1,2)) {
  Obj_fcs_data <- FindClusters(Obj_fcs_data, resolution = j)
  print(DimPlot(Obj_fcs_data, reduction = "umap") + labs(title = paste0("resolution: ", j)))
  DimPlot(Obj_fcs_data, reduction = "umap",label=TRUE) + labs(title = paste0("resolution: ", j))
  ggsave2(paste0("Righton_flow_cytometry_Analysis.Integration.combined.resolution.",j,".DimPlot.pdf"), path = "results/Dimensionality_reduction", width = 20, height = 20, units = "cm")
}

# clustree Graph
clustree(Obj_fcs_data)
ggsave2("Righton_flow_cytometry_Analysis.Integration.combined.clustree.pdf", path = "results/Dimensionality_reduction", width = 20, height = 20, units = "cm")

# Feature Plot
BioMarker = as.character(panel$Biomarker)
Idents(Obj_fcs_data) <- Obj_fcs_data@meta.data$integrated_snn_res.0.5
for (x in seq_along(BioMarker)) {
  print(FeaturePlot(Obj_fcs_data, features = BioMarker[x], coord.fixed = T, sort.cell = T, label=TRUE, label.size = 2, cols=rev(rainbow(10))))
  FeaturePlot(Obj_fcs_data, features = BioMarker[x], coord.fixed = T, sort.cell = T, label=TRUE, label.size = 2, cols=rev(rainbow(10)))
  ggsave2(paste0("Righton_flow_cytometry_Analysis.Integration.combined.",BioMarker[x],".FeaturePlot.png"), path = "results/Dimensionality_reduction", width = 30, height = 20, units = "cm")
}
saveRDS(Obj_fcs_data, file = "objects/mutiplex_integrated.RDS")
