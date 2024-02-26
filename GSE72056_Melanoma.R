# Melanoma dataset GSE72056
# downloaded GSE72056_melanoma_single_cell_revised_v2.txt.gz
# Code evaluates gene expression and gene signatures from pre-identified cell populations, as defined in Tirosh I, et al. Science. 352(6282): 189-196 (2016).
# Analysis was performed using Seurat package version 4.9.9.9060 

#load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(conflicted)


# Read in data file
setwd(".../GSE72056_Melanoma/data")
data <- read.delim("GSE72056_melanoma_single_cell_revised_v2.txt", header = T, stringsAsFactors = F)

# Configure data for analysis
  data_genes <- data[-1:-3,]

  gene_list <- data_genes %>% 
    pull("Cell") %>% 
    make.unique(sep = ".")

  rownames(data_genes) <- gene_list
  data_genes <- data_genes[, -1]
  data_meta <- data[1:3,]
  rownames(data_meta) <- data_meta[, 1]
  data_meta <- data_meta[, -1]
  data_meta_transpose <- data.frame(t(data_meta))

  data_2 <- CreateSeuratObject(counts = data_genes, assay = "RNA", meta.data = data_meta_transpose)
  View(data_2@meta.data)
  # Note: malignant cells and non-malignant cell id headers changed
  # malignant column: malignant(1=no,2=yes,0=unresolved)
  # Non-malignant cells defined by authors: 1= T cells, 2= B cells, 3= Macrophages, 4= Endothelial cells, 5=CAF, 6=NK



# set identity of clusters and validate cluster 4 as endothelial cells
  Idents(data_2)
  Idents(data_2) <- "cluster_id"
  Idents(data_2)

  levels(data_2) <- c('0', '1', '2', '3', '4', '5', '6')
  DotPlot(data_2, assay = 'RNA', features = c('PECAM1','CDH5','TIE1','ROBO4','CLDN5','VWF'))


  
# Generate Figures

  # CLSTN1 dot plot in cell populations 
    DotPlot(data_2, assay = 'RNA', features = c('CLSTN1'))
  
  

  #Run ssGSEA gene set enrichment for selected gene sets 
    gene.sets <- list(mTORC1_mediated_pathway = c("AKT1S1","EEF2K","EIF4B","EIF4E","EIF4EBP1","EIF4G1","LAMTOR1","LAMTOR2","LAMTOR3","LAMTOR4","LAMTOR5","MLST8","MTOR","RHEB","RPS6","RPS6KB1","RPTOR","RRAGA","RRAGB","RRAGD","SLC38A9","YWHAB"),
                      RAB_Regulation_of_Trafficking = c("AKT1","AKT2","AKT3","ALS2","ALS2CL","ANKRD27","ARF6","CCZ1","CCZ1B","CHM","CHML","DENND1A","DENND1B","DENND1C","DENND2A","DENND2B","DENND2C","DENND2D","DENND3","DENND4A","DENND4B","DENND4C","DENND5A","DENND5B","DENND6A","DENND6B","GABARAP","GABARAPL2","GAPVD1","GDI1","GDI2","GGA1","GGA2","GGA3","HPS1","HPS4","MADD","MAP1LC3B","MON1A","MON1B","OPTN","RAB10","RAB11A","RAB11B","RAB12","RAB13","RAB14","RAB18","RAB1A","RAB1B","RAB21","RAB27A","RAB27B","RAB31","RAB32","RAB33A","RAB33B","RAB35","RAB38","RAB39A","RAB39B","RAB3A","RAB3GAP1","RAB3GAP2","RAB3IL1","RAB3IP","RAB4A","RAB5A","RAB5B","RAB5C","RAB6A","RAB6B","RAB7A","RAB7B","RAB8A","RAB8B","RAB9A","RABEP1","RABGAP1","RABGEF1","RGP1","RIC1","RIN1","RIN2","RIN3","RINL","SBF1","SBF2","SYTL1","TBC1D10A","TBC1D10B","TBC1D10C","TBC1D13","TBC1D14","TBC1D15","TBC1D16","TBC1D17","TBC1D2","TBC1D20","TBC1D24","TBC1D25","TBC1D3","TBC1D7","TRAPPC1","TRAPPC10","TRAPPC11","TRAPPC12","TRAPPC13","TRAPPC2","TRAPPC2L","TRAPPC3","TRAPPC4","TRAPPC5","TRAPPC6A","TRAPPC6B","TRAPPC8","TRAPPC9","TSC1","TSC2","ULK1","YWHAE"))

    library(escape)
    ES_data <- enrichIt(obj = data_genes,
                          gene.sets = gene.sets,
                          min.size = NULL)

    ES_data <- Seurat::AddMetaData(data_2, ES_data)
    View(ES_data@meta.data)
    write.csv(ES_data@meta.data, file = "Melanoma_Patient_ssGSEA_enrichment_data.csv")
    # Gene Signatures for endothelial cells were plotted in Prism



  # Identify mTORC1 signaling low and high groups to examine CLSTN1 expression
    # Stratify by mTORC1_pathway median score for endothelial cells
    ES_data_mTORC1_low <- subset(ES_data, subset = mTORC1_mediated_pathway < 3966.35797)
    ES_data_mTORC1_high <- subset(ES_data, subset = mTORC1_mediated_pathway >= 3966.35797)
    ES_data@meta.data$mTORC1_level <- ifelse(rownames(ES_data@meta.data) %in% colnames(ES_data_mTORC1_low),
                                              "mTORC1_low",
                                              "mTORC1_high")

    # Extract endothelial cell subset
    Idents(ES_data) <- ES_data@meta.data$cluster_id
    endothelial <- subset(ES_data, idents = "4")
    
    

    # Dot Plot of CLSTN1 in endothelial cells based on mTORC1_pathway stratification
    DotPlot(endothelial, assay = 'RNA', 
            features = c('CLSTN1'), 
            group.by = "mTORC1_level", 
            cluster.idents = TRUE,
            cols = 'PRGn')
