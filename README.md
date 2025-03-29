# GAM_PRSS_REML_Project
GAM with PRSS and REML in scRNA-seq prostate cancer and benign cases dataset; Detailed GAM-PRSS-REML tracking for interpretability.

![Project Status](https://img.shields.io/badge/status-active-brightgreen)
![GitHub last commit](https://img.shields.io/github/last-commit/kahkengwong/GAM_PRSS_REML_Project)
![GitHub languages](https://img.shields.io/github/languages/count/kahkengwong/GAM_PRSS_REML_Project)
![GitHub top language](https://img.shields.io/github/languages/top/kahkengwong/GAM_PRSS_REML_Project)
![Contributors](https://img.shields.io/github/contributors/kahkengwong/GAM_PRSS_REML_Project)

Author: Kah Keng Wong  

---

## Overview
Comprehensive analysis of an scRNA-seq dataset of prostate cancer and benign prostate samples using the Seurat package, with subsequent modeling utilizing a generalized additive model (GAM) through the mgcv package. The GAM is implemented via thin-plate regression splines (TPRS) for smoothing, which are regularized by the penalized residual sum of squares (PRSS) based on the smoothing parameter $\lambda$ derived from restricted maximum likelihood (REML). The codes are designed not only to achieve the modeling but also to emphasize the **interpretability** of the modeling process through detailed extraction of components that constitute PRSS and REML, as well as validation of the parameters selected by PRSS ($k$) and REML ($\lambda$) through manual recalculation, 10-fold cross-validation, and visualization of relevant plots. *Code generation was assisted by Claude (Sonnet 3.5 and 3.7) from Anthropic, and refined by the author.*

The scripts are designed to be run sequentially, following the workflow of the main project/manuscript. Key goals include:
- Preprocessing and clustering scRNA-seq data using UMAP.
- Analyzing gene expressions (e.g., *TRPM4* and *KLK4*) and their associations with the relevant gene sets.
- Modeling of _TRPM4_ with the relevant gene sets by GAMs with PRSS and REML, optimizing parameters like $k$, $\lambda$, and $\gamma$.
- Ensuring interpretability by validating and visualizing the modeling process.

Key findings include robust modeling of _TRPM4_ expression with validated $k$ and $\lambda$ parameters, alongside detailed interpretations of GAM, PRSS, and REML mechanisms. A manuscript is in preparation for submission to a Q1 journal, with example results to be shared here upon acceptance.

---

# scRNA-seq Analysis and GAM-PRSS-REML Modeling Scripts
The scripts should be used in the following sequence, corresponding to the flow of the main project/manuscript:

## Descriptions of the Scripts
| No | Script File                                      | Description                                                                                     |
|------|--------------------------------------------------|-------------------------------------------------------------------------------------------------|
| 1    | `Part-1-scRNAseq-Preprocessing-and-UMAP-Clusters.r` | Preprocesses scRNA-seq data (prostate cancer and benign samples). Includes QC steps: removing low-quality cells, regressing out cell cycle phase effects, and correcting batch effects. Performs UMAP clustering to identify cell populations. |
| 2    | `Part-2-UMAP-Heatmap-and-SK-Matrix.r`           | Visualizes *TRPM4* and *KLK4* expression across UMAP clusters in prostate cancer and benign samples. Computes Spearman’s and Kendall’s correlation matrices and generates a heatmap to compare *TRPM4* against other relevant gene sets. |
| 3    | `Part-3.1-GAM-PRSS-REML-Setup.r`                | Sets up the Generalized Additive Model (GAM) with Restricted Maximum Likelihood (REML) and assesses model convergence. |
| 4    | `Part-3.2-GAM-PRSS-REML-Analysis.r`             | Performs GAM analysis, extracts best-fitting models, and summarizes key statistics, including PRSS and REML convergence. |
| 5    | `Part-3.3-REML-Extraction-and-Convergence.r`    | Extracts and analyzes detailed optimization data from REML processes for **interpretability**. Recalculates REML components manually to verify `mgcv`-computed scores. |
| 6    | `Part-3.4-GAM-PRSS-REML-Plots-and-EDF-Analysis.r` | Generates visual plots to analyze GAM modeling of gene expression. Tracks PRSS and REML optimization to understand how $k$ and $\lambda$ are selected, emphasizing **interpretability**. |
| 7    | `Part-3.5-TPRS-Visualization-and-GAM-Components.r` | Visualizes TPRS and GAM components for **interpretability**. Explains spline basis construction around knots, penalization by $\lambda$, weighting by coefficients, and how regularized splines combine to form the full GAM fit. |
| 8    | `Part-3.6-Validation-of-k-and-Lambda-Selection.r` | Validates the selection of $k$ and $\lambda$ using the nested REML approach. Refits models with independent $k$ or $\lambda$ values (without wrapping REML in PRSS, as in the `analyze_sample` function) and performs 10-fold cross-validation (RMSE and RSS) to ensure reproducibility. |
| 9    | `Part-3.7-Gamma-Consequences-on-GAM-Fitting.r`  | Assesses the impact of different $\gamma$ values on GAM fitting. Compares models with varying $\gamma$ to a reference model (default $\gamma$) using 10-fold cross-validation (RMSE). |

## Additional Details of the Scripts
| No | Script File                                      | Code Lines | Main Functions                                                                                          | Main Outputs                                                                                                                    |
|----|--------------------------------------------------|------------|----------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| 1  | `Part-1-scRNAseq-Preprocessing-and-UMAP-Clusters.r` | 322        | `process_seurat()`, `filter_ribosomal()`, `filter_mitochondrial()`, `remove_doublets()`, `correct_batch_effects()`, `generate_elbow_plot()`, `downstream_analyses()` | Filtered Seurat objects, UMAP plots, elbow plots, cluster marker TSV files, and a saved workspace |
| 2  | `Part-2-UMAP-Heatmap-and-SK-Matrix.r`           | 197        | `plot_umap_gene_expression()`, `get_cluster_averages()`, `plot_correlation_matrix()`, `calculate_correlations()`, `scale_to_zscore()` | UMAP plots for TRPM4 and KLK4 expression, correlation matrices, TRPM4 correlation Excel files, and a heatmap |
| 3  | `Part-3.1-GAM-PRSS-REML-Setup.r`                | 1113       | `extract_basis_functions()`, `extract_model_details()`, `calculate_manual_prediction()`, `calculate_prss()`, `extract_reml_iterations_enhanced()`, `analyze_sample()`, `analyze_reml_convergence()` | Best GAM model, PRSS data, REML iterations, best model parameters, sample data with predictions, REML convergence analysis, and REML convergence statistics |
| 4  | `Part-3.2-GAM-PRSS-REML-Analysis.r`             | 1023       | `analyze_multiple_gene_sets()`, `format_equation()`, `create_summary_df()`, `extract_best_models_and_params()`, `export_detailed_results()`, `export_prss_data()`, `add_reml_convergence_columns()` | Summary Excel file, detailed PCa results, detailed non-Ca results, PRSS data for PCa, PRSS data for non-Ca, REML convergence analysis, and REML summary |
| 5  | `Part-3.3-REML-Extraction-and-Convergence.r`    | 1530       | `extract_reml_info()`, `extract_lambda_components()`, `process_reml_info()`, `extract_reml_components_simple()`, `process_reml_components_simple()`, `analyze_reml_discrepancies()`, `extract_reml_optimization_details()` | REML summary and components Excel file, REML information with legend, REML components with discrepancies, PCa REML information, non-Ca REML information, convergence details Excel file, and optimization details |
| 6  | `Part-3.4-GAM-PRSS-REML-Plots-and-EDF-Analysis.r` | 873        | `prepare_data()`, `reml_optimization_visualization()`, `generate_visualizations()`, `save_plot_to_files()`, `generate_and_save_visualizations()`, `plot_lambda_values()`, `create_gam_plot()` | PRSS and REML plots for PCa and non-Ca, GAM plots for PCa and non-Ca, EDF summary Excel file, detailed basis functions report, non-linear relationships summary, non-linear stats, and lambda value plots |
| 7  | `Part-3.5-TPRS-Visualization-and-GAM-Components.r` | 1134       | `visualize_tprs_improved()`, `plot_cumulative_components_revised()`, `plot_gam_with_examples()`, `visualize_phi1_emergence()`, `predict()`, `gam()`, `ggplot()` | TPRS basis functions plot, weighted basis functions plot, GAM components plot, variance contribution plots, cumulative variance plot, cumulative smooth components plot, GAM components with examples plot |
| 8  | `Part-3.6-Validation-of-k-and-Lambda-Selection.r` | 1929       | `analyze_prss_for_k_values()`, `analyze_eigenvalues_for_lambdas()`, `run_kfold_cv_analysis_for_k()`, `run_kfold_cv_analysis_improved()`, `calculate_prss()`, `extract_eigenvalues()`, `perform_kfold_cv_for_k()` | PRSS versus k plot, eigenvalue bar plots, eigenvalue distribution plots, additional eigenvalue comparison plots, CV k versus RMSE and deviance plots, PRSS versus CV k comparison plot, CV lambda versus RMSE and deviance plots |
| 9  | `Part-3.7-Gamma-Consequences-on-GAM-Fitting.r`  | 472        | `visualize_gam_gammas_free()`, `evaluate_gamma_boxplot()`, `gam()`, `predict()`, `ggplot()`, `pivot_longer()` | GAM full fits plot, difference from gamma = 1 plot, selected lambda values bar plot, selected lambda values data, CV RMSE boxplot, CV summary statistics, CV RMSE raw data |

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
For further information or questions, please email [kahkeng@usm.my](mailto:kahkeng@usm.my)

---
