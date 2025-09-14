<div align="center">

# Generalized Additive Modeling of _TRPM4_-Ribo Transcriptional Space in Prostate Cancer
### _Modeling of _TRPM4_-Ribo relationship with GAM, with focus on modeling optimization and interpretability._

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
1) To predict the functions of the oncogene _TRPM4_ in high-grade prostate cancer (scRNA-seq dataset) via GAM, PRSS, and REML modeling. 
2) To investigate GAM’s touted interpretability using scRNA-seq’s continuous expression profiles, navigating the limits of its interpretability such as mgcv’s (R package for GAM) lambda optimization process (e.g., gradients, Hessians, and transformed REML formula). 

The scripts include extraction and visualization of PRSS and REML components, and validating parameters (*k*, $\lambda$, $\gamma$) through manual recalculation, and visualizations. Key goals include:
- Preprocessing and clustering scRNA-seq data with UMAP.
- Analyzing gene expressions (e.g., *TRPM4*, *KLK4*) and their gene set associations.
- Modeling *TRPM4* with GAMs, optimizing *k*, $\lambda$, and $\gamma$ (manually set, affecting $\lambda$, and validated qualitatively and quantitatively).
- Ensuring interpretability through validation and visualization.

Key findings include *TRPM4* modeling with validated parameters and detailed GAM, PRSS, and REML mechanistic insights. A manuscript is in preparation for a Q1 journal (target: April 2025), with results shared here upon acceptance.

---

# scRNA-seq Analysis and GAM-PRSS-REML Modeling Scripts
The scripts should be used in the following sequence, corresponding to the flow of the main project/manuscript:

## Descriptions of the Scripts
| No | Script File                                      | Description                                                                                     |
|------|--------------------------------------------------|-------------------------------------------------------------------------------------------------|
| 1    | `Part-1-scRNAseq-Preprocessing-and-UMAP-Clusters.r` | Preprocesses scRNA-seq data (PCa and non-Ca samples). Includes QC steps: removing low-quality cells, regressing out cell cycle phase effects, and correcting batch effects. Performs UMAP clustering to identify cell populations. |
| 2    | `Part-2-UMAP-Heatmap-and-SK-Matrix.r`           | Visualizes *TRPM4* and *KLK4* expression across UMAP clusters in PCa and non-Ca samples. Computes Spearman’s and Kendall’s correlation matrices and generates a heatmap to compare *TRPM4* against other relevant gene sets. |
| 3    | `Part-3.1-GAM-PRSS-REML-Setup.r`                | Sets up the Generalized Additive Model (GAM) with Restricted Maximum Likelihood (REML) and assesses model convergence. |
| 4    | `Part-3.2-GAM-PRSS-REML-Analysis.r`             | Performs GAM analysis, extracts best-fitting models, and summarizes key statistics, including PRSS and REML convergence. |
| 5    | `Part-3.3-REML-Extraction-and-Convergence.r`    | Extracts and analyzes detailed optimization data from REML processes for **interpretability**. Recalculates REML components manually to verify `mgcv`-computed scores. |
| 6    | `Part-3.4-GAM-PRSS-REML-Plots-and-EDF-Analysis.r` | Generates visual plots to analyze GAM modeling of gene expression. Tracks PRSS and REML optimization to understand how $k$ and $\lambda$ are selected, emphasizing **interpretability**. |
| 7    | `Part-3.5-TPRS-Visualization-and-GAM-Components.r` | Visualizes TPRS and GAM components for **interpretability**. Explains spline basis construction around knots, penalization by $\lambda$, weighting by coefficients, and how regularized splines combine to form the full GAM fit. |
| 8    | `Part-3.6-Validation-of-k-and-Lambda-Selection.r` | Validates the selection of $k$ and $\lambda$ using the nested REML approach. Refits models with independent $k$ or $\lambda$ values (without wrapping REML in PRSS, as in the `analyze_sample` function) and performs 10-fold CV (RMSE and RSS) to ensure reproducibility. |
| 9    | `Part-3.7-Gamma-Consequences-on-GAM-Fitting.r`  | Assesses the impact of different $\gamma$ values on GAM fitting. Compares models with varying $\gamma$ to a reference model (default $\gamma$) using 10-fold CV (RMSE). |

## Additional Details of the Scripts
| No | Script File                                      | Code Lines | Main Functions                                                                                          | Main Outputs                                                                                                                    |
|----|--------------------------------------------------|------------|----------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| 1  | `Part-1-scRNAseq-Preprocessing-and-UMAP-Clusters.r` | 322        | `process_seurat()`, `filter_ribosomal()`, `filter_mitochondrial()`, `remove_doublets()`, `correct_batch_effects()`, `generate_elbow_plot()`, `downstream_analyses()` | Filtered Seurat objects, UMAP plots, elbow plots, cluster marker TSV files, and a saved workspace |
| 2  | `Part-2-UMAP-Heatmap-and-SK-Matrix.r`           | 197        | `plot_umap_gene_expression()`, `get_cluster_averages()`, `plot_correlation_matrix()`, `calculate_correlations()`, `scale_to_zscore()` | UMAP plots for _TRPM4_ and _KLK4_ expression, correlation matrices, TRPM4 correlation Excel files, and a heatmap |
| 3  | `Part-3.1-GAM-PRSS-REML-Setup.r`                | 1113       | `extract_basis_functions()`, `extract_model_details()`, `calculate_manual_prediction()`, `calculate_prss()`, `extract_reml_iterations_enhanced()`, `analyze_sample()`, `analyze_reml_convergence()` | Best GAM model, PRSS data, REML iterations, best model parameters, sample data with predictions, REML convergence analysis, and REML convergence statistics |
| 4  | `Part-3.2-GAM-PRSS-REML-Analysis.r`             | 1023       | `analyze_multiple_gene_sets()`, `format_equation()`, `create_summary_df()`, `extract_best_models_and_params()`, `export_detailed_results()`, `export_prss_data()`, `add_reml_convergence_columns()` | Summary Excel file, detailed PCa results, detailed non-Ca results, PRSS data for PCa, PRSS data for non-Ca, REML convergence analysis, and REML summary |
| 5  | `Part-3.3-REML-Extraction-and-Convergence.r`    | 1530       | `extract_reml_info()`, `extract_lambda_components()`, `process_reml_info()`, `extract_reml_components_simple()`, `process_reml_components_simple()`, `analyze_reml_discrepancies()`, `extract_reml_optimization_details()` | REML summary and components Excel file, REML information with legend, REML components with discrepancies, PCa REML information, non-Ca REML information, convergence details Excel file, and optimization details |
| 6  | `Part-3.4-GAM-PRSS-REML-Plots-and-EDF-Analysis.r` | 873        | `prepare_data()`, `reml_optimization_visualization()`, `generate_visualizations()`, `save_plot_to_files()`, `generate_and_save_visualizations()`, `plot_lambda_values()`, `create_gam_plot()` | PRSS and REML plots for PCa and non-Ca, GAM plots for PCa and non-Ca, EDF summary Excel file, detailed basis functions report, non-linear relationships summary, non-linear stats, and $\lambda$ value plots |
| 7  | `Part-3.5-TPRS-Visualization-and-GAM-Components.r` | 1134       | `visualize_tprs_improved()`, `plot_cumulative_components_revised()`, `plot_gam_with_examples()`, `visualize_phi1_emergence()`, `predict()`, `gam()`, `ggplot()` | TPRS basis functions plot, weighted basis functions plot, GAM components plot, variance contribution plots, cumulative variance plot, cumulative smooth components plot, GAM components with examples plot |
| 8  | `Part-3.6-Validation-of-k-and-Lambda-Selection.r` | 1929       | `analyze_prss_for_k_values()`, `analyze_eigenvalues_for_lambdas()`, `run_kfold_cv_analysis_for_k()`, `run_kfold_cv_analysis_improved()`, `calculate_prss()`, `extract_eigenvalues()`, `perform_kfold_cv_for_k()` | PRSS versus $k$ plot, eigenvalue bar plots, eigenvalue distribution plots, additional eigenvalue comparison plots, CV $k$ versus RMSE and deviance plots, PRSS versus CV $k$ comparison plot, CV $\lambda$ versus RMSE and deviance plots |
| 9  | `Part-3.7-Gamma-Consequences-on-GAM-Fitting.r`  | 472        | `visualize_gam_gammas_free()`, `evaluate_gamma_boxplot()`, `gam()`, `predict()`, `ggplot()`, `pivot_longer()` | GAM full fits plot, difference from gamma = 1 plot, selected $\lambda$ values bar plot, selected $\lambda$ values data, CV RMSE boxplot, CV summary statistics, CV RMSE raw data |

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
