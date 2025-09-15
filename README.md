<div align="center">

# 🧬 Generalized Additive Modeling of _TRPM4_-Ribo Transcriptional Space in Prostate Cancer
### Modeling of _TRPM4_-Ribo relationship with generalized additive model (GAM), with focus on modeling optimization, interpretability, mapping for top-ranked explanatory power (TREP) cells, and identification of their biology.

![Project Status](https://img.shields.io/badge/status-active-brightgreen)
![GitHub last commit](https://img.shields.io/github/last-commit/kahkengwong/GAM_PCa_Project)
![GitHub languages](https://img.shields.io/github/languages/count/kahkengwong/GAM_PCa_Project)
![Contributors](https://img.shields.io/github/contributors/kahkengwong/GAM_PCa_Project)

**Author:** **Kah** **Keng** **Wong**  

</div>

---

## 📋Overview
Analysis of an scRNA-seq dataset of prostate cancer (PCa) and benign prostate (non-Ca) samples using the **Seurat** package, with subsequent modeling via generalized additive models (GAMs) through the **mgcv** package. The GAM is implemented with thin-plate regression splines (TPRS), regularized by penalized residual sum of squares (PRSS) based on the smoothing parameter $\lambda$ derived from restricted maximum likelihood (REML). 

## 🎯Aims of the Project
The following are the project's objectives:
1) To optimize and interpret the modeling with mgcv package including *k*, $\lambda$, and $\gamma$ optimization process, PRSS and REML convergence, and visualization of splines' formation.
2) To uncover how much does _TRPM4_ explain variability in Ribo (average expression of seven ribosomal genes) expression via deviance explained, GAM's performance metrics. 
3) To identify which cells hold strong _TRPM4_-Ribo relationship, and to uncover their biology via Gene Ontology (GO) enrichment analysis. 

This processed dataset is based on publicly available data from: 
Wong HY, Sheng Q, Hesterberg AB, Croessmann S et al. Single cell analysis of cribriform prostate cancer reveals cell intrinsic and tumor microenvironmental pathways of aggressive disease. Nat Commun 2022;13(1):6036. https://doi.org/10.1038/s41467-022-33780-1

## 💻 Requirements
### Software Requirements
- **R:** ≥ 4.0.0 (tested with R 4.4.2 "Pile of Leaves")
- **RStudio:** Recommended (tested with RStudio 2025.05.0 Build 496)

### Hardware Requirements
- **RAM:** 16GB (32GB recommended)
- **CPU:** 8 cores (12+ cores recommended)
- **Storage:** ~12GB free disk space for data and results

### Performance Notes
- All analyses are CPU-based; GPU acceleration not required
- Higher specs are primarily for scRNA-seq preprocessing (e.g., doublet removal, Kendall's τ computation)
- GAM modeling with [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html) is computationally efficient and does not require high-end specs
- Understanding the algorithm/modeling is key; minimum specs are sufficient for GAM analysis using the processed Seurat object

### Development Environment
- **Developed and tested on:** Intel Core i9-14900KF, 64GB DDR5 RAM, Windows 11 Pro, RStudio 2025.05.0 Build 496, R 4.4.2

---

## 🔀Workflow of the Project
![Workflow](https://raw.githubusercontent.com/kahkengwong/GAM_PCa_Project/main/Project_Workflow.jpg)

Key methdologies of this study include **decomposition of model's DE into cell-level**: If a cell highly influences the shape of the fitted spline, then it should have high **explanatory power (EP)**. These cells are termed as **top-ranked EP (TREP) cells**. The decomposition steps are detailed in part (V) of the figure above. 

Another key method is to binarize the transcriptional space by **inflection point (IP)** into **pre-IP and post-IP regions**. These regions exhibit distinct distribution pattern of TREP cells, producing **quadrants of four subpopulation of cells with different biology** - These are detailed in the next section.  

---

## 📊Project Key Findings
![Project Key Findings](https://raw.githubusercontent.com/kahkengwong/GAM_PCa_Project/main/Project_Key_Findings.jpg)


**Hence, the CEP-IP framework based on pairwise genes produces quadrants of cell subpopulations, enabling distinct biology identification that may guide therapy.**


---

## 📊scRNA-seq Analysis and GAM-PRSS-REML Modeling Scripts
The processed Seurat object `GSE185344_Seurat_processed.RData` (9.52 GB) is available on HuggingFace: 
https://huggingface.co/datasets/kahkengwong/GAM_PCa_Project

The results of the GAM modeling in this study can be replicated by analyzing the processed Seurat object `GSE185344_Seurat_processed.RData` by following the code block `Part_3.01_Mean_Expression_Justifications.r` until `Part_3.15_Monocle3_Pre-IP_vs_Post-IP_TREP.r`

For the complete workflow, the scripts should be used in the following sequence, corresponding to the flow of the main project/manuscript:

## 📜Descriptions of the Scripts
| No | Script File                                      | Description                                                                                     |
|------|--------------------------------------------------|-------------------------------------------------------------------------------------------------|
| 1    | `Part_1_scRNAseq_preprocessing_` <br> `and_UMAP_clusters.r` | Preprocesses scRNA-seq data (PCa and non-Ca samples). Includes QC steps: removing low-quality cells, regressing out cell cycle phase effects, and correcting batch effects. Performs UMAP clustering to identify cell populations. |
| 2    | `Part_2_UMAP_Heatmap_Spearman-` <br> `Kendall's-matrix.r` | Visualizes _TRPM4_ expression in UMAP clusters for PCa and non-Ca samples. Calculates Spearman’s and Kendall’s correlation matrices and creates a heatmap to compare TRPM4 with other significant gene sets. |
| 3    | `Part_3.01_Mean_Expression_` <br> `Justifications.r`                | Analyzes reliability of Ribo and AR gene sets in tumor and non-cancerous (BP) samples. Computes Cronbach’s α, McDonald’s ω, and KMO scores, with summary statistics and gene set reliability metrics. |
| 4    | `Part_3.02_Family_Distribution_` <br> `Analysis.r`             | Performs GAM diagnostics for PCa and non-Ca samples, testing multiple distribution families. Generates diagnostic plots and exports AIC/BIC metrics and best family results. |
| 5    | `Part_3.03_GAM_REML_PRSS_` <br> `Setup.r`    | GAM-REML-PRSS modeling with `mgcv` for PCa and non-Ca samples for multiple gene sets. Extracts basis functions, model details, and PRSS metrics, analyzing REML convergence. |
| 6    | `Part_3.04_GAM_REML_PRSS_` <br> `Analysis.r` | Conducts GAM-REML-PRSS analysis for PCa and non-Ca samples across gene sets. Extracts model parameters, calculates PRSS, and convergence metrics. |
| 7    | `Part_3.05_REML_Extraction_` <br> `and_Convergence.r` | Extracts detailed REML convergence metrics for GAM models in PCa and non-Ca samples. Processes convergence details, including scores and iterations, and exports results to Excel for analysis. |
| 8    | `Part_3.06_REML_PRSS_Plots_` <br> `and_EDF_Analysis.r` | Creates visualizations of REML and PRSS iterations for PCa and non-Ca samples. Generates detailed EDF reports, applies FDR correction, and non-linear relationship summaries. |
| 9    | `Part_3.07_Validation_of_` <br> `k_and_Lambda_Selection.r`  | Validates GAM models by refitting with varied k and lambda values for PCa samples. Analyzes PRSS and REML scores, generating plots to identify optimal parameters for model robustness. |
| 10    | `Part_3.08_Visualize_TPRS_` <br> `and_GAM_Components.r`  | Visualizes TPRS and GAM components for PCa samples, plotting basis functions, knot placement, and data distribution to illustrate spline formation and model contributions. |
| 11    | `Part_3.09_Extract_GAM's_` <br> `ND_MD_DE.r`  | Extracts and analyzes null deviance, model deviance, and deviance explained for GAM models in PCa and non-Ca samples. |
| 12    | `Part_3.10_Extract_TRPM4-` <br> `Ribo_EP.r`  | Extracts TRPM4-Ribo explanatory power for PCa samples at the cell level. Computes deviance contributions, sorts cells by explanatory power (EP), and exports detailed results. |
| 13    | `Part_3.11_MCCV_of_CEP_` <br> `Classification.r`  | Monte Carlo cross-validation (MCCV) for cell EP in PCa samples. Classifies deviance cells, generates leverage-based random control plots, and validates GAM model performance. |
| 14    | `Part_3.12_TREP_and_` <br> `non-TREP_GAM_Plots.r`  | Generates scatter plots for PCa samples, distinguishing TRPM4-Ribo explanatory power (TREP) cells in purple from non-TREP cells in gray, and exports cell-level data. |
| 15    | `Part_3.13_TREP_and_non-TREP_` <br> `Mosaic_and_Raincloud_Plots.r`  | Creates mosaic and raincloud plots for Ribo expression in PCa samples, comparing pre-IP and post-IP cells. Exports contingency data and distribution visualizations. |
| 16    | `Part_3.14_TREP_vs_` <br> `non-TREP_DEGs_Analysis.r`  | Identifies differentially expressed genes (DEGs) between TREP and non-TREP cells in pre-IP and post-IP groups for PCa samples. |
| 17    | `Part_3.15_Monocle3_` <br> `Pre-IP_vs_Post-IP_TREP.r`  | Generates Monocle3 trajectory visualizations for TREP and non-TREP cells in PCa samples. Performs quantitative clustering analysis, comparing UMAP1 distributions. |

## 🛠️Packages and Dependencies
The R packages and versions used in this study were saved in the `renv.lock` file included in the GitHub repository. This lockfile contains version information for all 37 packages used with their dependencies. Download the repository and run `renv::restore()` to install the identical package versions used in this study. Note that renv creates an isolated library and will not modify a system's existing R packages and setup. 

For manual or selected package installation, a simplified `r_packages_info.json` file is provided with the package names, versions, and sources (CRAN, Bioconductor, or GitHub).

## 🚀 Quick Start
1. **Get the repository:** 
   - **Option A:** Fork this repository on GitHub, then download to your local machine
   - **Option B:** Click the green "< > Code" button (top right of the file list), select "Download ZIP", and extract to a folder.
2. **Open the project:** Navigate to the extracted folder and double-click any `.r` file to launch RStudio
3. **Set working directory:** In RStudio, go to Session → Set Working Directory → To Source File Location (or use `setwd`)
4. **Install dependencies:** Type `renv::restore()` in the RStudio console and press Enter to install all required packages
5. **Download data:** Get the processed Seurat object from HuggingFace: [Processed_Seurat_Object](https://huggingface.co/datasets/kahkengwong/GAM_PCa_Project)
6. **Place data file:** Save the downloaded `.RData` file in the same folder as your scripts
7. **Run analysis:** Execute scripts Part_3.01 through Part_3.15 sequentially for main analysis with `mgcv`

---

## 🧾License
This project is licensed under the [MIT License](https://github.com/kahkengwong/GAM_PCa_Project/blob/main/LICENSE), an open-source license to encourage collaboration and reuse, while ensuring proper attribution to the original author(s). For the full details, please refer to the [LICENSE](https://github.com/kahkengwong/GAM_PCa_Project/blob/main/LICENSE) file in this repository.


---

## 🤝🏻Contributing
Contributions are welcome! Please open an issue or submit a pull request if you have suggestions or improvements.

---

## 📩Contact
All analyses, modeling, and interpretations were conducted by Kah Keng Wong [kahkeng@usm.my](mailto:kahkeng@usm.my), [kahkeng3@gmail.com](mailto:kahkeng3@gmail.com)

---

## ✨ Support the Project
If you find this project valuable, please consider **starring** ⭐ or **forking** 🍴 the repository. Your support helps others discover these analyses and insights 🙌
