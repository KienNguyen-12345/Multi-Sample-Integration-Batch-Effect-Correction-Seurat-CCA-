```markdown
# Single-Cell RNA Sequencing Integration and Batch Correction in Hepatoblastoma (GSE180665)

## Table of Contents
- [Study Overview](#study-overview)
- [Dataset Information](#dataset-information)
- [Workflow Architecture](#workflow-architecture)
- [Installation & Requirements](#installation--requirements)
- [Pipeline Breakdown](#pipeline-breakdown)
- [Results & Visualization](#results--visualization)
- [Key Takeaways](#key-takeaways)

---

## Study Overview

Hepatoblastoma (HB) is the most common primary liver tumor in pediatric populations. Understanding intratumoral heterogeneity and distinct tumor cell populations is critical to elucidating key genetic mechanisms and cellular signaling pathways involved in HB pathogenesis.

This project focuses on the computational integration of scRNA-seq datasets derived from primary human tumor tissue, adjacent background liver, and patient-derived xenografts (PDX). Using standard batch-effect correction strategies in **Seurat v3**, this workflow aligns biological states across patient cohorts while mitigating patient-specific technical variation.

---

## Dataset Information

* **Accession:** GSE180665
* **Organism:** *Homo sapiens*
* **Sample Types:** 
  * Primary Hepatoblastoma Tumor
  * Adjacent Background Liver
  * Patient-Derived Xenografts (PDX)
* **Patients Included:** HB17, HB30, HB53
* **Format:** 10x Genomics raw feature-barcode matrices (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`)

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
