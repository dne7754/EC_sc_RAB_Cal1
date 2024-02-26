# Breast cancer dataset GSE176078 (available at NCBI GEO)
# single cell seq breast cancer patients = 26

# downloaded "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
# pre-processed data has been filtered, normalized, and scaled
# cell clusters previously identified in Wu SZ, et al. Nat Genet. 53(9):1334-1347 (2021).
# Analysis was performed using Seurat package version 4.9.9.9060


# Load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)

# Configure data for analysis
setwd(".../GSE176078_Breast_Cancer/data")

  #If using Seurat v5
  options(Seurat.object.assay.version = "v3")

  expression_matrix <- ReadMtx(
    mtx = "count_matrix_sparse.mtx",
    cells = "count_matrix_barcodes.tsv",
    features = "count_matrix_genes.tsv",
    cell.column = 1,
    feature.column = 1)

  data_2 <- CreateSeuratObject(counts = expression_matrix, assay = "RNA")
  View(data_2@meta.data)

  data_2$sample <- rownames(data_2@meta.data)
  View(data_2@meta.data)

  data_2@meta.data <- separate(data_2@meta.data, col = 'sample', into = c('Patient', 'Barcode'), 
                              sep = '_')
  View(data_2@meta.data)
  unique(data_2@meta.data$Patient)


# Quality control-check to ensure only cells with low mtDNA content were used
  data_2$mitoPercent <- PercentageFeatureSet(data_2, pattern='^MT-')
  View(data_2@meta.data)

  VlnPlot(data_2, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
  FeatureScatter(data_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = 'lm')


# add Sample metadata
  metadata_file <- read.csv("metadata.csv")
  data_2[['subtype']] <-metadata_file$subtype[match(rownames(data_2@meta.data), metadata_file$X)] 
  data_2[['celltype_subset']] <-metadata_file$celltype_subset[match(rownames(data_2@meta.data), metadata_file$X)] 
  data_2[['celltype_minor']] <-metadata_file$celltype_minor[match(rownames(data_2@meta.data), metadata_file$X)] 
  data_2[['celltype_major']] <-metadata_file$celltype_major[match(rownames(data_2@meta.data), metadata_file$X)] 
  View(data_2@meta.data)

  Idents(data_2)
  Idents(data_2) <- "celltype_major"
  Idents(data_2)


# Scale counts to TP10K (transcripts per 10k)
  tp10k = NormalizeData(obj = data_2, scale.factor = 1e4, normalization.method="RC")@assays$RNA@data
  tp10k_data <- CreateSeuratObject(counts = tp10k)

  tp10k_data[['subtype']] <-metadata_file$subtype[match(rownames(tp10k_data@meta.data), metadata_file$X)]  
  tp10k_data[['celltype_subset']] <-metadata_file$celltype_subset[match(rownames(tp10k_data@meta.data), metadata_file$X)] 
  tp10k_data[['celltype_minor']] <-metadata_file$celltype_minor[match(rownames(tp10k_data@meta.data), metadata_file$X)] 
  tp10k_data[['celltype_major']] <-metadata_file$celltype_major[match(rownames(tp10k_data@meta.data), metadata_file$X)] 
  View(tp10k_data@meta.data)


# Validate pre-identified endothelial cells by generating a dot plot of endothelial cell markers 
  Idents(tp10k_data)
  Idents(tp10k_data) <- "celltype_major"
  Idents(tp10k_data)

  DotPlot(tp10k_data, assay = 'RNA', features = c('PECAM1','CDH5','TIE1','ROBO4','CLDN5','VWF'))


# Generate Figures
  
  # CLSTN1 dot plot in pre-identified cell populations
    DotPlot(tp10k_data, assay = 'RNA', features = c('CLSTN1'))



  # Run ssGSEA gene set enrichment for selected gene sets
    # get expression data of endothelial cells
    cells.use <- WhichCells(object = tp10k_data, ident = 'Endothelial')
    expr <- LayerData(object = tp10k_data, assay = "RNA", layer = "data")[, cells.use]
    expr <- as(Class = 'matrix', object = expr)
    write.csv(x = expr, file = "expresssion_endothelial.csv", quote = FALSE)

    gene.sets <- list(mTORC1_mediated_pathway = c("AKT1S1","EEF2K","EIF4B","EIF4E","EIF4EBP1","EIF4G1","LAMTOR1","LAMTOR2","LAMTOR3","LAMTOR4","LAMTOR5","MLST8","MTOR","RHEB","RPS6","RPS6KB1","RPTOR","RRAGA","RRAGB","RRAGD","SLC38A9","YWHAB"),
                      RAB_Regulation_of_Trafficking = c("AKT1","AKT2","AKT3","ALS2","ALS2CL","ANKRD27","ARF6","CCZ1","CCZ1B","CHM","CHML","DENND1A","DENND1B","DENND1C","DENND2A","DENND2B","DENND2C","DENND2D","DENND3","DENND4A","DENND4B","DENND4C","DENND5A","DENND5B","DENND6A","DENND6B","GABARAP","GABARAPL2","GAPVD1","GDI1","GDI2","GGA1","GGA2","GGA3","HPS1","HPS4","MADD","MAP1LC3B","MON1A","MON1B","OPTN","RAB10","RAB11A","RAB11B","RAB12","RAB13","RAB14","RAB18","RAB1A","RAB1B","RAB21","RAB27A","RAB27B","RAB31","RAB32","RAB33A","RAB33B","RAB35","RAB38","RAB39A","RAB39B","RAB3A","RAB3GAP1","RAB3GAP2","RAB3IL1","RAB3IP","RAB4A","RAB5A","RAB5B","RAB5C","RAB6A","RAB6B","RAB7A","RAB7B","RAB8A","RAB8B","RAB9A","RABEP1","RABGAP1","RABGEF1","RGP1","RIC1","RIN1","RIN2","RIN3","RINL","SBF1","SBF2","SYTL1","TBC1D10A","TBC1D10B","TBC1D10C","TBC1D13","TBC1D14","TBC1D15","TBC1D16","TBC1D17","TBC1D2","TBC1D20","TBC1D24","TBC1D25","TBC1D3","TBC1D7","TRAPPC1","TRAPPC10","TRAPPC11","TRAPPC12","TRAPPC13","TRAPPC2","TRAPPC2L","TRAPPC3","TRAPPC4","TRAPPC5","TRAPPC6A","TRAPPC6B","TRAPPC8","TRAPPC9","TSC1","TSC2","ULK1","YWHAE"))

    library(escape)
    ES_data <- enrichIt(obj = tp10k_data,
                        gene.sets = gene.sets,
                        min.size = NULL)
    ES_data <- Seurat::AddMetaData(tp10k_data, ES_data)
    write.csv(ES_data@meta.data, file = "BC_Patient_ssGSEA_enrichment_data.csv")
    # Gene Signatures for endothelial cells were plotted in Prism
  


  # Identify mTORC1 signaling low and high groups
    # stratify by mTORC1_pathway median score for endothelial cells
    ES_data_mTORC1_low <- subset(ES_data, subset = mTORC1_mediated_pathway < 7671.53297)
    ES_data_mTORC1_high <- subset(ES_data, subset = mTORC1_mediated_pathway >= 7671.53297)
    ES_data@meta.data$mTORC1_level <- ifelse(rownames(ES_data@meta.data) %in% colnames(ES_data_mTORC1_low),
                                            "mTORC1_low",
                                            "mTORC1_high")

    #Extract endothelial cell subset
    Idents(ES_data) <- ES_data@meta.data$celltype_major
    endothelial <- subset(ES_data, idents = "Endothelial")


    #Dot Plot of CLSTN1 in endothelial cells based on mTORC1_pathway stratification
    DotPlot(endothelial, assay = 'RNA', 
            features = c('CLSTN1'), 
            group.by = "mTORC1_level", 
            cluster.idents = TRUE,
            cols = 'PRGn')

















