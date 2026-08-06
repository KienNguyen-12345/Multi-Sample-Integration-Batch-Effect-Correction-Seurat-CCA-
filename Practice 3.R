# Study Design:
# Identification of distinct tumor cell populations and key genetic mehcanism in hepatoblastoma (HB) through single cell RNA sequencing (scRNA-seq) - GSE180665
#Design: samples from human tumor, background liver and patient derived xenografts were colected to demostrate gene expression patterns within tumor and to identify intratumor cell subtype heterogenety to define differeing roles in pahotgenesis based on intracellular signlaiing in pediatric HB.
#Goal: To integrate data from different patients and correct for batch effects.

# What is integration?
#Integration: In single-cell RNA sequencing (scRNA-seq), integration is a computational process that 
# aligns and harmonizes multiple datasets to correct for batch effects—unwanted technical variations caused by differences 
# in experimental protocols, sequencing platforms, or donors.  The primary goal is to match shared cell types and states across
# these datasets so that cells with similar biological features cluster together, while preserving true biological differences 
# between experimental conditions. 

# When to integrate?
# Integrating >1 scRNA-seq datasets (treated/untrated, KO/WT, different samples) -- Can be from same cells or different cells
# Cell label transfer - transfer cell type classifications from a refrence to query dataset.
# Integration of multimodal single cell data (e.g. scRNA-Seq and ScATAC-Seq) - integrate, into a single0cell multi-omic dataset, signals collected from separate assays.
# Integration of scRNA-Seq and spatial expressiond ata - integrate topological arragnement of cells in tissues with gene expression data.


# Types of integration:
  # Horizontal integration
# * Same modality from independent cells
# * E.g. scRNA-seq from same tissue but in different patients/donors/sequencing technologies
# * Assays are anchored by common gene set

  #Vertical Integration:
# Multiple modalities profiled simutaneously from same cells
# E.g. scRNA-seq and scATAC-seq from same cells
# Assays are anchored by cells

  #Diagonal integration:
# Different modalities from different cells
# scRNA-seq and scATAC-seq performed on seqparate group of cells.

# Batch correction methods:
# 7 ways: MNN, Seurat v3, LIGER, Harmony, BBKNN, scVI, Conos, Scmap, Scanorama, scAlign
# --> This time we will use seurat V3
setwd("/Users/nguyenkien/Desktop/Practice R/Lesson 3/Practice 3")
installed.packages("gridExtra")
#Load libraries:
library(gridExtra)
library(Seurat)
library(ggplot2)
library(tidyverse)
#Increase memory (from 16Gb to 32Gb):
install.packages("usethis") # Run if you don't have this installed
usethis::edit_r_environ()
mem.maxVSize()
#Download data:
# +) After loading down, the file names end with .tar.gz.
# +) Remember to open these tar files to get access to the files
# Get data location:
dirs <- list.dirs(path = "/Users/nguyenkien/Desktop/Practice R/Lesson 3/Practice 3/Data", 
                  recursive = FALSE, 
                  full.names = FALSE
                  )
for(x in dirs){
  # 1. Clean the folder name to use as an R variable name
  name <- gsub("_filtered_feature_bc_matrix", '', x)
  
  # 2. Build the full path to the sample folder
  sample_path <- paste0('/Users/nguyenkien/Desktop/Practice R/Lesson 3/Practice 3/Data/', x, '/')
  
  # 3. Read the 10x matrix files
  cts <- ReadMtx(
    mtx      = paste0(sample_path, 'matrix.mtx.gz'),
    features = paste0(sample_path, 'features.tsv.gz'),
    cells    = paste0(sample_path, 'barcodes.tsv.gz')
  )
  
  #Explain: when we open each unzipped file, there will be 3 files. One is called "features.tsv.gz"
  # and one is called "barcodes.tsv.gz", we will assign feature with "feature" file and cells with "barcodes" files
  #Explain: paste0() joins multiple pieces of text together into a single string with zero spaces between them.
  
  # 4. Create Seurat Object and assign it to the clean variable name
  assign(name, CreateSeuratObject(counts = cts, project = name))
}


#Merge datasets:
merged_surat <- merge(HB17_background, y= c(HB17_PDX,HB17_tumor, HB30_PDX, HB30_tumor, HB53_background, HB53_tumor),
      add.cell.ids = ls () [3:9],
      project = "HB")
#This function merges HB17_background with all of other seurat files then assign the newly created version with variable called "merged_seurat"

#QC and filtering:
view(merged_surat@meta.data)

#Create sample column:
merged_surat$sample <- rownames(merged_surat@meta.data)
#What we did was just to replicate the row names into a new column called "sample"

#Split sample column:
merged_surat@meta.data <- separate(merged_surat@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'),
         sep = '_')
# What we did: in the merged_surat meta data, we use the separate function to separate the sample column into 3 smaller columns
# named "Patient", "Type", and "Barcode". The way the code understands where it separate is by using the "_" as a sign that the next phrase belongs to a different column

#Calculate mitochondrial percentage:
merged_surat@meta.data$percent.mt <- PercentageFeatureSet(merged_surat, pattern= "^MT-")

#Explore QC:
VlnPlot(merged_surat, features = c("nCount_RNA","nFeature_RNA", "percent.mt"), ncol=3)
#Draws a violin plot of single cell data (gene expression, metrics, PC scores, etc.)

FeatureScatter(merged_surat, feature1 = "nCount_RNA", feature2 = 'nFeature_RNA') + geom_smooth(method= "lm")
              
#Creates a scatter plot of two features (typically feature expression), across a set of single cells. Cells are colored by their identity class. Pearson correlation between the two features is displayed above the plot.
#Geom_smooth: Aids the eye in seeing patterns in the presence of overplotting. geom_smooth() and stat_smooth() are effectively aliases: they both use the same arguments. Use stat_smooth() if you want to display the results with a non-standard geom.

#Filtering:

merged_surat_filtered <- subset(merged_surat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   percent.mt < 10)
#Perform standard worflow steps to figure out if we see any batch effects:

#1. Normalize data:

merged_surat_filtered <- NormalizeData(merged_surat_filtered)

#2. Identify highly variable features:

merged_surat_filtered <- FindVariableFeatures(merged_surat_filtered, selection.method= "vst")

# Identify 10 most highly variable genes:
top10 <- head(VariableFeatures(merged_surat_filtered), 10)

#Plot variables:
plot1 <- VariableFeaturePlot(merged_surat_filtered)
LabelPoints(plot = plot1, points = top10, repel = T)

#3. Scaling:
all.gene <- rownames(merged_surat_filtered)
merged_surat_filtered <- ScaleData(merged_surat_filtered)

#4. Perform linear dimensionality reduction:
merged_surat_filtered <- RunPCA(merged_surat_filtered, features = VariableFeatures(merged_surat_filtered))
DimHeatmap(merged_surat_filtered, dims= 1, cells = 500, balanced= T)

#Determine dimensionality of the data:
ElbowPlot(merged_surat_filtered)

#5. Clustering:
merged_surat_filtered <- FindNeighbors(merged_surat_filtered, dims = 1: 20)

#Understanding resolution:
merged_surat_filtered <- FindClusters(merged_surat_filtered)
view(merged_surat_filtered@meta.data)


#UMAP (For non-dimensional reduction plot):
merged_surat_filtered <- RunUMAP(merged_surat_filtered, dims = 1:20)

#Plot:
p1 <- DimPlot(merged_surat_filtered, reduction= "umap", group.by = "Patient")

p2 <- DimPlot(merged_surat_filtered, reduction = "umap",group.by = "Type", cols = c("red", "green","blue"))

grid.arrange(p1, p2, nrow= 2)

#Comment:in this graph, we see 2 plots of dimsplot showing clustering groups based on clustering using tissue type or using original patients.
# We can see that the graphs of "patient", the clusterings are due to different in patients and not due to different in biological differences --> masking
# --> must correct for the batch effects

#Batch effect adjustment method:
#Perform integration steps:
obj.list <- SplitObject(merged_surat_filtered, split.by = 'Patient')

for(i in 1:length(obj.list)){
  obj.list[[i]][["RNA"]] <- as(object = obj.list[[i]][["RNA"]], Class = "Assay")
}

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# Select integration features:
features <- SelectIntegrationFeatures(object.list = obj.list)
# Find integration anchors (CCA):
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features,
                                  reduction = "cca")

saveRDS(anchors, file = "cca_anchors_HB.rds")

#Integrate data:
seurat.integrated <- IntegrateData(anchorset = anchors)

# Ccale the data:
seurat.integrated <- ScaleData(seurat.integrated)

# 2. Run PCA on the integrated assay
seurat.integrated <- RunPCA(seurat.integrated, assay = "integrated")

# 3. Run UMAP using the first 50 principal components (or choose based on ElbowPlot)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:50)

#4 Visualize:
p3 <- DimPlot(seurat.integrated, reduction = "umap", group.by =  "Patient")
p4 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Type")

grid.arrange(p3,p4, ncol=2)

