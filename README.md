

---

```markdown
# Single-Cell RNA Sequencing (scRNA-Seq) Data Integration & Batch Effect Correction in Hepatoblastoma

A single-cell bioinformatics pipeline designed to load, quality-control (QC), process, and integrate scRNA-seq datasets from human hepatoblastoma (HB) samples, background liver tissue, and Patient-Derived Xenografts (PDX) to correct for batch effects using **Seurat v5 CCA Integration**.

---

## 📌 Table of Contents
- [Study Overview](#-study-overview)
- [Dataset Information](#-dataset-information)
- [Workflow Architecture](#-workflow-architecture)
- [Installation & Requirements](#-installation--requirements)
- [Pipeline Breakdown](#-pipeline-breakdown)
- [Results & Visualization](#-results--visualization)
- [Key Takeaways](#-key-takeaways)

---

## 🔬 Study Overview

* **Publication Context:** Based on the single-cell RNA-seq analysis of pediatric Hepatoblastoma (HB) tissue samples.
* **GEO Accession:** `GSE180665` / `GSE180666`
* **Primary Objective:** Identify intratumoral cell subtype heterogeneity and driver populations across multiple patients while harmonizing non-biological batch effects introduced by independent sequencing runs and biological donors.

---

## 📊 Dataset Information

The pipeline analyzes scRNA-seq data collected across **3 Patient IDs** and **3 Tissue Types**:

| Metadata Column | Values / Description |
| :--- | :--- |
| **Patient ID** | `HB17`, `HB30`, `HB53` |
| **Sample Type** | `background` (Normal liver), `tumor` (Primary HB tumor), `PDX` (Patient-derived xenograft) |

---

## 🛠 Workflow Architecture


```

10X Matrix Files (.mtx.gz, .tsv.gz)
│
▼
Data Import Loop & Seurat Object Creation
│
▼
Dataset Merging & Metadata Extraction (Patient / Type)
│
▼
Quality Control & Filtering (nCount_RNA, nFeature_RNA, percent.mt)
│
▼
Unintegrated Workflow (Normalize ➔ Variable Features ➔ Scale ➔ PCA ➔ UMAP)
│
├─► [QC Plot]: Verify presence of Patient Batch Effects
│
▼
Batch Effect Adjustment (SplitObject ➔ CCA Integration Anchors ➔ IntegrateData)
│
▼
Post-Integration Dimensionality Reduction (Scale ➔ PCA ➔ UMAP)
│
▼
[Final Plot]: Verify successful batch correction & biological preservation

```

---

## 💻 Installation & Requirements

### Required R Packages

Ensure you have R (>= 4.0) and the following libraries installed:

```R
install.packages(c("ggplot2", "tidyverse", "gridExtra", "usethis"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("Seurat")

```

---

## 🧬 Pipeline Breakdown

### 1. Environment & Data Import

* Configures memory management using `usethis::edit_r_environ()`.
* Automatically iterates over sample directories containing 10X Genomics output files (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`).
* Instantiates individual Seurat objects using dynamic variable assignment.

### 2. QC & Filtering

* Merges individual samples into a single dataset (`merged_surat`).
* Parses metadata string IDs into discrete columns (`Patient`, `Type`, `Barcode`).
* Calculates mitochondrial content percentage (`percent.mt`).
* Filters out low-quality cells based on strict thresholds:
* `nCount_RNA > 800`
* `nFeature_RNA > 500`
* `percent.mt < 10%`



### 3. Evaluating Pre-Integration Batch Effects

* Runs classical linear and non-linear dimensionality reduction (PCA, UMAP).
* Generates baseline plots (`p1`, `p2`) grouped by `Patient` and `Type`.
* **Diagnostic Rule:** If cells cluster primarily by `Patient` rather than cell identity, a **batch effect** is present and requires computational integration.

### 4. Integration & Batch Effect Correction (Seurat v3 CCA)

* Splits the merged dataset by `Patient` into an object list (`obj.list`).
* Applies a compatibility fix converting v5 assay structures to standard v3 `Assay` objects.
* Normalizes and selects highly variable features across each batch independently.
* Identifies integration anchors via **Canonical Correlation Analysis (CCA)** (`FindIntegrationAnchors`).
* Saves generated anchors (`cca_anchors_HB.rds`) for reproducibility.
* Projects all cells into a unified, batch-corrected expression assay (`seurat.integrated`).

### 5. Final Integrated Visualization

* Scales integrated data, runs PCA, and calculates a 50-dimension UMAP embedding.
* Generates final diagnostic plots (`p3`, `p4`) using `grid.arrange()`.

---

## 📈 Results & Visualization

| Analysis Stage | `Patient` Plot Observation | `Type` Plot Observation | Verdict |
| --- | --- | --- | --- |
| **Pre-Integration** | Cells form isolated, distinct clusters based on donor ID (`HB17`, `HB30`, `HB53`). | Sample types are fragmented across patient boundaries. | ❌ **Batch Effect Present** |
| **Post-Integration** | Cells from all 3 patients smoothly intermingle within shared clusters. | `background` tissue forms a distinct cluster, while `tumor` and `PDX` overlap. | ✅ **Batch Effect Corrected** |

---

## 🔑 Key Takeaways

1. **Successful Integration:** Overlapping patient distributions in the post-integration UMAP confirm that non-biological technical variation has been effectively removed.
2. **Biological Truth Preserved:** Healthy `background` liver tissue remains distinctly separate from malignant cells, confirming that integration preserves genuine biological differences.
3. **PDX Model Validation:** The strong spatial alignment between `tumor` and `PDX` cells indicates that patient-derived xenografts faithfully replicate primary human hepatoblastoma gene expression landscapes.

```

```
