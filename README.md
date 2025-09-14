<div align="center">

# Generalized Additive Modeling of _TRPM4_-Ribo Transcriptional Space in Prostate Cancer
### _Modeling of _TRPM4_-Ribo relationship with GAM, with focus on modeling optimization, interpretability, mapping for top-ranked explanatory power (TREP) cells, and identification of their biology._

![Project Status](https://img.shields.io/badge/status-active-brightgreen)
![GitHub last commit](https://img.shields.io/github/last-commit/kahkengwong/GAM_PCa_Project)
![GitHub repo size](https://img.shields.io/github/repo-size/kahkengwong/GAM_PCa_Project)
![GitHub languages](https://img.shields.io/github/languages/count/kahkengwong/GAM_PCa_Project)
![GitHub top language](https://img.shields.io/github/languages/top/kahkengwong/GAM_PCa_Project)
![Contributors](https://img.shields.io/github/contributors/kahkengwong/GAM_PCa_Project)

Author: Kah Keng Wong  

</div>

---

## Overview
Analysis of an scRNA-seq dataset of prostate cancer (PCa) and benign prostate (non-Ca) samples using the **Seurat** package, with subsequent modeling via generalized additive models (GAMs) through the **mgcv** package. The GAM is implemented with thin-plate regression splines (TPRS), regularized by penalized residual sum of squares (PRSS) based on the smoothing parameter $\lambda$ derived from restricted maximum likelihood (REML). 

## Aims of the Project
The following are the project's objectives:
1) To optimize and interpret the modeling with mgcv package including *k*, $\lambda$, and $\gamma$ optimization process, PRSS and REML convergence, and visualization of splines' formation.
2) To uncover how much does _TRPM4_ explain variability in Ribo (average expression of seven ribosomal genes) expression via deviance explained, GAM's performance metrics. 
3) To identify which cells hold strong _TRPM4_-Ribo relationship, and to uncover their biology via Gene Ontology (GO) enrichment analysis. 

Key findings include _TRPM4_-Ribo modeling with validated parameters and detailed GAM, PRSS, and REML mechanistic insights.

---

# scRNA-seq Analysis and GAM-PRSS-REML Modeling Scripts
The scripts should be used in the following sequence, corresponding to the flow of the main project/manuscript:

## Descriptions of the Scripts
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
| 9    | `Part_3.08_Visualize_TPRS_` <br> `and_GAM_Components.r`  | Visualizes TPRS and GAM components for PCa samples, plotting basis functions, knot placement, and data distribution to illustrate spline formation and model contributions. |
| 9    | `Part_3.09_Extract_GAM's_` <br> `ND_MD_DE.r`  | Extracts and analyzes null deviance, model deviance, and deviance explained for GAM models in PCa and non-Ca samples. |
| 9    | `Part_3.10_Extract_TRPM4-` <br> `Ribo_EP.r`  | Extracts TRPM4-Ribo explanatory power for PCa samples at the cell level. Computes deviance contributions, sorts cells by explanatory power (EP), and exports detailed results. |
| 9    | `Part_3.11_MCCV_of_CEP_` <br> `Classification.r`  | Monte Carlo cross-validation (MCCV) for cell EP in PCa samples. Classifies deviance cells, generates leverage-based random control plots, and validates GAM model performance. |
| 9    | `Part_3.12_TREP_and_` <br> `non-TREP_GAM_Plots.r`  | Generates scatter plots for PCa samples, distinguishing TRPM4-Ribo explanatory power (TREP) cells in purple from non-TREP cells in gray, and exports cell-level data. |
| 9    | `Part_3.13_TREP_and_` <br> `non-TREP_Mosaic_and_Raincloud_Plots.r`  | Creates mosaic and raincloud plots for Ribo expression in PCa samples, comparing pre-IP and post-IP cells. Exports contingency data and distribution visualizations. |
| 9    | `Part_3.14_TREP_vs_` <br> `non-TREP_DEGs_Analysis.r`  | Identifies differentially expressed genes (DEGs) between TREP and non-TREP cells in pre-IP and post-IP groups for PCa samples. |
| 9    | `Part_3.15_Monocle3_` <br> `Pre-IP_vs_Post-IP_TREP.r`  | Generates Monocle3 trajectory visualizations for TREP and non-TREP cells in PCa samples. Performs quantitative clustering analysis, comparing UMAP1 distributions. |

## Dependencies
This project requires the following R packages:
- **Seurat**: scRNA-seq preprocessing and clustering.
- **mgcv**: Generalized additive modeling with PRSS and REML.
- **purrr**, **dplyr**, **tidyr**, **ggplot2**: Data manipulation and visualization of results.
- **Matrix**, **SparseArray**: Sparse matrix handling.
- **openxlsx**, **writexl**: Excel file I/O.
- **parallel**, **pbapply**: Parallel processing.

Additional libraries (e.g., `circlize`, `ComplexHeatmap`, `monocle3`, `viper`) are used for specific analyses and visualizations. See script headers in `Part-1-scRNAseq-Preprocessing-and-UMAP-Clusters.r` for full list of packages used.

---

## License
This project is licensed under the [MIT License](https://github.com/kahkengwong/GAM_PRSS_REML_Project/blob/main/LICENSE), an open-source license to encourage collaboration and reuse, while ensuring proper attribution to the original author(s). For the full details, please refer to the [LICENSE](https://github.com/kahkengwong/GAM_PRSS_REML_Project/blob/main/LICENSE) file in this repository.


---

## Contributing
Contributions are welcome! Please open an issue or submit a pull request if you have suggestions or improvements.

---

## Contact
All analyses, modeling, and interpretations were conducted by KK Wong [kahkeng3@gmail.com](mailto:kahkeng3@gmail.com)

---
