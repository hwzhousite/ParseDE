# PARSE-DE Pipeline

This repository implements the **PARSE-DE** workflow for differential expression (DE) analysis in single-cell data.  
The project is organized into modular scripts for reproducibility and flexibility.

---

## 📂 File Structure and Descriptions

| File | Description |
|------|--------------|
| **`0_data_loading.R`** | Handles **data input and preprocessing**. <br>Update the data file path here to analyze your own dataset. |
**`1_parse_fit.R`** | Performs **PARSE model fitting** on the input data. <br>Includes **parameter initialization**, **model optimization**, and **fitted object output**. <br><br>🔧 **Important:** The tuning parameters `lambda` and `pca_adjust` must be carefully selected. Start with a **small `lambda`** and a **`pca_adjust` value close to 1** to ensure informative genes are properly captured. |
| **`2_de_gene_procedure.R`** | Conducts **differential expression (DE) analysis** using the fitted PARSE model.<br> In the Seurat part, resolution shall be adjusted. <br>Generates DE gene lists, statistical summaries, and visualization-ready results. |
| **`utility.R`** | Contains **helper functions** used across the workflow (e.g., normalization, plotting, matrix operations, and evaluation metrics). |

---

## 🚀 Usage

1. **Set up your environment**
   ```r
   install.packages(c("PARSE", "Seurat", "mvtnorm", "ggplot2", "Matrix","knockoff"))
