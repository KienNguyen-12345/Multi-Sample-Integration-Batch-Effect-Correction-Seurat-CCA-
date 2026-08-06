# Single-Cell RNA Sequencing Integration and Batch Correction in Hepatoblastoma (GSE180665)

An end-to-end single-cell RNA sequencing (scRNA-Seq) data integration and batch correction pipeline implemented in R using the **Seurat v3** framework. This workflow processes single-cell transcriptomic profiles derived from primary human tumor tissue, adjacent background liver, and patient-derived xenografts (PDX) from pediatric hepatoblastoma (HB) patients (GSE180665). By leveraging Canonical Correlation Analysis (CCA), the pipeline harmonizes data across patient cohorts to correct for patient-specific batch effects while preserving true biological heterogeneity.

---

## 📋 Overview of Workflow Steps

1. **Data Loading & Preprocessing:** Iteratively load 10x Genomics raw feature-barcode matrices (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`) across sample directories using `ReadMtx()` and construct individual `SeuratObject` containers.


2. **Dataset Merging & Metadata Parsing:** Merge individual sample objects into a combined dataset and parse sample identities using `tidyr::separate()` to extract `Patient`, `Type`, and `Barcode` attributes.


3. **Quality Control (QC) & Filtering:** Compute mitochondrial read percentages (`percent.mt`) using `PercentageFeatureSet()` and filter low-quality cells, empty droplets, and damaged cells based on expression thresholds (`nCount_RNA > 800`, `nFeature_RNA > 500`, and `percent.mt < 10%`).


4. **Unintegrated Pipeline Evaluation:** Execute standard single-cell processing steps (`NormalizeData`, `FindVariableFeatures`, `ScaleData`, `RunPCA`, and `RunUMAP`) on the uncorrected data to assess patient-specific clustering masks and technical variation.


5. **Multi-Dataset Splitting:** Split the merged object by patient ID (`SplitObject`) into individual sub-objects to prepare for dataset harmonization.


6. **Integration Feature Selection:** Normalize individual sub-objects independently and identify shared highly variable integration features using `SelectIntegrationFeatures()`.


7. **Anchor Identification & CCA Integration:** Compute integration anchors across patient datasets using Canonical Correlation Analysis (`FindIntegrationAnchors`) and run `IntegrateData()` to harmonize gene expression profiles across batches.


8. **Downstream Scaling & Dimensionality Reduction:** Perform data scaling (`ScaleData`), Principal Component Analysis (`RunPCA`), and non-linear dimensionality reduction (`RunUMAP`) on the integrated assay to visualize batch-corrected cell populations.

---

## Workflow Architecture


```

```
                            [ Raw 10x Genomics Matrices ]
                                          │
                                          ▼
                             [ Data Import & Merging ]
                                          │
                                          ▼
                           [ QC & Filtering (MT% / Counts) ]
                                          │
                                          ▼
                           [ Unintegrated Pipeline ] 
                           (Normalization -> PCA -> UMAP)
                                          │
                                          ▼
                         [ Batch Effect Identification ]
                                          │
                                          ▼
                         [ Seurat v3 Anchor Integration ]
                         (Canonical Correlation Analysis)
                                          │
                                          ▼
                            [ Post-Integration UMAP ]

```

```

---

## Installation & Requirements

### R Version & Core Packages
Ensure you have **R 4.0** or higher installed. Install the required R dependencies using the code below:

```R
# Install core CRAN packages
install.packages(c("tidyverse", "ggplot2", "gridExtra", "usethis"))

# Install Seurat
install.packages("Seurat")

```

### Memory Allocation

Because single-cell datasets require significant RAM, it is recommended to expand memory capacity in your `.Renviron` file:

```R
library(usethis)
usethis::edit_r_environ()
# Add: R_MAX_VSIZE=32Gb

```

---

## Pipeline Breakdown

### 1. Data Import and Formatting

Iterates through sample directories, imports 10x Genomics matrix files via `ReadMtx()`, and constructs individual `SeuratObject` containers for each sample condition.

### 2. Merging & Metadata Extraction

Merges all individual datasets into a single object and parses sample names using `tidyr::separate()` to break down barcodes into distinct metadata columns: `Patient`, `Type`, and `Barcode`.

### 3. Quality Control (QC) & Filtering

Filtering thresholds applied to remove low-quality cells, empty droplets, and dying cells:

* **UMI Count (`nCount_RNA`):** > 800
* **Gene Count (`nFeature_RNA`):** > 500
* **Mitochondrial Percentage (`percent.mt`):** < 10%

### 4. Unintegrated Workflow Evaluation

Standard Seurat workflow (`NormalizeData`, `FindVariableFeatures`, `ScaleData`, `RunPCA`, `RunUMAP`) was run to visualize uncorrected batch effects across patient samples.

### 5. Seurat v3 CCA Integration

To correct for patient-specific clustering masks:

* Split data by patient (`SplitObject`).
* Normalized and identified variable features independently for each sub-object.
* Selected integration features using `SelectIntegrationFeatures()`.
* Identified integration anchors via Canonical Correlation Analysis (`FindIntegrationAnchors()`).
* Integrated datasets with `IntegrateData()`.

---

## Results & Visualization

### Quality Control Inspection

```R
# Visualize distributions
VlnPlot(merged_surat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
FeatureScatter(merged_surat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")

```

### Unintegrated vs. Integrated Comparison

Prior to integration, dimension reduction plots demonstrate strong clustering driven primarily by patient identity rather than cell lineage:

```R
# Pre-integration visualization
p1 <- DimPlot(merged_surat_filtered, reduction = "umap", group.by = "Patient")
p2 <- DimPlot(merged_surat_filtered, reduction = "umap", group.by = "Type", cols = c("red", "green", "blue"))
grid.arrange(p1, p2, nrow = 2)

```

Post-integration using CCA anchors aligns shared cell types across patients while preserving distinct tissue-type differences (Tumor vs. Background vs. PDX):

```R
# Post-integration visualization
p3 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Patient")
p4 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Type")
grid.arrange(p3, p4, ncol = 2)

```

---

## Key Takeaways

* **Batch Effect Elimination:** Unintegrated scRNA-seq analysis leads to patient-driven clustering, which masks underlying biological insights.
* **Anchor Identification:** Seurat v3's CCA integration successfully aligns shared cellular states across diverse donor backgrounds.
* **Biological Context Retention:** Tumor-specific signatures and background liver traits are cleanly separated post-integration without patient-level technical bias.

```

```
