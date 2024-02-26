# Lung adenocarcinoma
# from Bischoff P, et al. Oncogene 40:6748-6748 (2021).
# The following files were downloaded for analysis: "cellranger", "curated_annotation", and "metadata"

# Processing of raw data (filtering, normalization, dimensionality reduction, scaled, and cell type annotation) was performed using code provided in citation above ("code/lung_all_final.R")
# Analysis was performed using Seurat package version 4.9.9.9060

# Processed data saved as "all_SCTransform_cluster_id.RDS"
# Processed data from tumor samples saved as "tumor_tissue.RDS"


# Load libraries
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)

# Load data files
setwd("/Users/deannaedwards/Downloads/Lung_adenocarcinoma/data")

  #All samples (both tumor and normal tissues)
  seu_obj <- readRDS("all_SCTransform_cluster_id.RDS")
  
  #Tumor samples only
  tumor <- readRDS("tumor_tissue.RDS")

# Check to endothelial cluster
  Idents(seu_obj) <- seu_obj$SCT_snn_res.0.2
  DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)
  DotPlot(seu_obj, assay = 'RNA', features = c('PECAM1','CDH5','TIE1','ROBO4','CLDN5','VWF'))

  Idents(tumor) <- tumor@meta.data$SCT_snn_res.0.2
  DimPlot(tumor, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)
  DotPlot(tumor, assay = 'RNA', features = c('PECAM1','CDH5','TIE1','ROBO4','CLDN5','VWF'))


# Generate Figures
  
  # CLSTN1 dot plot in pre-identified cell populations
    levels(tumor) <- c('Endothelial', 'T Cells', 'B Cells', 'NK Cells', 'Fibroblasts', 'Myeloid', 'Normal Epithelial', 'Cancer Epithelial')
    DotPlot(tumor, assay = 'RNA', features = c('CLSTN1'))


  # Run ssGSEA gene set enrichment for selected gene sets
    gene.sets <- list(mTORC1_mediated_pathway = c("AKT1S1","EEF2K","EIF4B","EIF4E","EIF4EBP1","EIF4G1","LAMTOR1","LAMTOR2","LAMTOR3","LAMTOR4","LAMTOR5","MLST8","MTOR","RHEB","RPS6","RPS6KB1","RPTOR","RRAGA","RRAGB","RRAGD","SLC38A9","YWHAB"),
                      RAB_Regulation_of_Trafficking = c("AKT1","AKT2","AKT3","ALS2","ALS2CL","ANKRD27","ARF6","CCZ1","CCZ1B","CHM","CHML","DENND1A","DENND1B","DENND1C","DENND2A","DENND2B","DENND2C","DENND2D","DENND3","DENND4A","DENND4B","DENND4C","DENND5A","DENND5B","DENND6A","DENND6B","GABARAP","GABARAPL2","GAPVD1","GDI1","GDI2","GGA1","GGA2","GGA3","HPS1","HPS4","MADD","MAP1LC3B","MON1A","MON1B","OPTN","RAB10","RAB11A","RAB11B","RAB12","RAB13","RAB14","RAB18","RAB1A","RAB1B","RAB21","RAB27A","RAB27B","RAB31","RAB32","RAB33A","RAB33B","RAB35","RAB38","RAB39A","RAB39B","RAB3A","RAB3GAP1","RAB3GAP2","RAB3IL1","RAB3IP","RAB4A","RAB5A","RAB5B","RAB5C","RAB6A","RAB6B","RAB7A","RAB7B","RAB8A","RAB8B","RAB9A","RABEP1","RABGAP1","RABGEF1","RGP1","RIC1","RIN1","RIN2","RIN3","RINL","SBF1","SBF2","SYTL1","TBC1D10A","TBC1D10B","TBC1D10C","TBC1D13","TBC1D14","TBC1D15","TBC1D16","TBC1D17","TBC1D2","TBC1D20","TBC1D24","TBC1D25","TBC1D3","TBC1D7","TRAPPC1","TRAPPC10","TRAPPC11","TRAPPC12","TRAPPC13","TRAPPC2","TRAPPC2L","TRAPPC3","TRAPPC4","TRAPPC5","TRAPPC6A","TRAPPC6B","TRAPPC8","TRAPPC9","TSC1","TSC2","ULK1","YWHAE"))

    library(escape)
    ES_data <- enrichIt(obj = seu_obj,
                        gene.sets = gene.sets,
                        min.size = NULL)
    ES_data <- Seurat::AddMetaData(seu_obj, ES_data)
    write.csv(ES_data@meta.data, file = "LC_Patient_ssGSEA_enrichment_data.csv")
    # Gene Signatures for endothelial cells were plotted in Prism
  

  #Identify mTORC1 signaling low and high groups
    #stratify by mTORC1_pathway median score for endothelial cells
    ES_data_mTORC1_low <- subset(ES_data, subset = mTORC1_mediated_pathway < 5691.124076)
    ES_data_mTORC1_high <- subset(ES_data, subset = mTORC1_mediated_pathway >= 5691.124076)
    ES_data@meta.data$mTORC1_level <- ifelse(rownames(ES_data@meta.data) %in% colnames(ES_data_mTORC1_low),
                                            "mTORC1_low",
                                            "mTORC1_high")

    #Extract endothelial cell subset
    Idents(ES_data) <- ES_data@meta.data$SCT_snn_res.0.2
    endothelial <- subset(ES_data, idents = "Endothelial")


    #Dot Plot of CLSTN1 in endothelial cells based on mTORC1_pathway stratification
    #stratified by median mTORC1_pathway score for endothelial cells
    DotPlot(endothelial, assay = 'RNA', 
            features = c('CLSTN1'), 
            split.by = "tissue_type",
            group.by = "mTORC1_level", 
            cluster.idents = TRUE,
            cols = 'PRGn')










