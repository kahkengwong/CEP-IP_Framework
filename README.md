# GAM_PRSS_REML_Project
GAM with PRSS and REML in scRNA-seq prostate cancer and benign cases dataset; Detailed GAM-PRSS-REML tracking for interpretability.

![GitHub last commit](https://img.shields.io/github/last-commit/kahkengwong/GAM_PRSS_REML_Project)
![GitHub languages](https://img.shields.io/github/languages/count/kahkengwong/GAM_PRSS_REML_Project)
![GitHub top language](https://img.shields.io/github/languages/top/kahkengwong/GAM_PRSS_REML_Project)
![Contributors](https://img.shields.io/github/contributors/kahkengwong/GAM_PRSS_REML_Project)

Author: Kah Keng Wong  

---

## Overview
Comprehensive analysis of an scRNA-seq dataset of prostate cancer and benign prostate samples using the Seurat package, with subsequent modeling utilizing a generalized additive model (GAM) through the mgcv package. The GAM is implemented via thin-plate regression splines (TPRS) for smoothing, which are regularized by the penalized residual sum of squares (PRSS) based on the smoothing parameter $\lambda$ derived from restricted maximum likelihood (REML). The codes are designed not only to achieve the modeling but also to emphasize the **interpretability** of the modeling process through detailed extraction of components that constitute PRSS and REML, as well as validation of the parameters selected by PRSS ($k$) and REML ($\lambda$) through manual recalculation, 10-fold cross-validation, and visualization of relevant plots.

The scripts are designed to be run sequentially, following the workflow of the main project/manuscript. Key goals include:
- Preprocessing and clustering scRNA-seq data using UMAP.
- Analyzing gene expressions (e.g., *TRPM4* and *KLK4*) and their associations with the relevant gene sets.
- Modeling of _TRPM4_ with the relevant gene sets by GAMs with PRSS and REML, optimizing parameters like $k$, $\lambda$, and $\gamma$.
- Ensuring interpretability by validating and visualizing the modeling process.

---

## scRNA-seq Analysis and GAM-PRSS-REML Modeling Scripts
The scripts should be used in the following sequence, corresponding to the flow of the main project/manuscript:

# Table format
| No | Script File                                      | Description                                                                                     | Code Lines | Functions                                                                                          | Outputs                                                                                                                    |
|----|--------------------------------------------------|-------------------------------------------------------------------------------------------------|------------|----------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| 1  | `Part-1-scRNAseq-Preprocessing-and-UMAP-Clusters.r` | Preprocesses scRNA-seq data (prostate cancer and benign samples). Includes QC steps: removing low-quality cells, regressing out cell cycle phase effects, and correcting batch effects. Performs UMAP clustering to identify cell populations. | 322        | `process_seurat()`, `filter_ribosomal()`, `filter_mitochondrial()`, `remove_doublets()`, `correct_batch_effects()`, `generate_elbow_plot()`, `downstream_analyses()` | Filtered Seurat objects (prostate_ca_seurat_integrated, non_cancerous_seurat_integrated), UMAP plots (pre- and post-integration), elbow plots (prostate_ca_ElbowPlot.pdf, non_cancerous_ElbowPlot.pdf), cluster marker TSV files (prostate_ca_top_markers_for_each_cluster_vRibo.tsv, non_cancerous_top_markers_for_each_cluster_vRibo.tsv), and a saved workspace (Dt2_scRNAseq_workspace_vRibo_v2.RData) |
| 2  | `Part-2-UMAP-Heatmap-and-SK-Matrix.r`           | Visualizes *TRPM4* and *KLK4* expression across UMAP clusters in prostate cancer and benign samples. Computes Spearman’s and Kendall’s correlation matrices and generates a heatmap to compare *TRPM4* against other relevant gene sets. | 197        | `plot_umap_gene_expression()`, `get_cluster_averages()`, `plot_correlation_matrix()`, `calculate_correlations()`, `scale_to_zscore()` | UMAP plots for TRPM4 and KLK4 expression, correlation matrices (cluster_correlation_matrix_spearman.csv, cluster_correlation_matrix_kendall.csv), TRPM4 correlation Excel files (TRPM4_Correlations_PCa_Clusters_Combined.xlsx, TRPM4_Correlations_PCa_Cluster_16.xlsx, TRPM4_Correlations_BPNonCa.xlsx), and a heatmap (TRPM4_heatmap.pdf) |
| 3  | `Part-3.1-GAM-PRSS-REML-Setup.r`                | Sets up the Generalized Additive Model (GAM) with Restricted Maximum Likelihood (REML) and assesses model convergence. | 1113       | `extract_basis_functions()`, `extract_model_details()`, `calculate_manual_prediction()`, `calculate_prss()`, `extract_reml_iterations_enhanced()`, `analyze_sample()`, `analyze_reml_convergence()` | Best GAM model (best_model), PRSS data (prss_data), REML iterations (all_reml_iterations), best model parameters (best_params_df), sample data with predictions (gam_data), REML convergence analysis (convergence_analysis), and REML convergence statistics (convergence_stats) |
| 4  | `Part-3.2-GAM-PRSS-REML-Analysis.r`             | Performs GAM analysis, extracts best-fitting models, and summarizes key statistics, including PRSS and REML convergence. | 1023       | `analyze_multiple_gene_sets()`, `format_equation()`, `create_summary_df()`, `extract_best_models_and_params()`, `export_detailed_results()`, `export_prss_data()`, `add_reml_convergence_columns()` | Summary Excel file (Ribo_AR_All_Analysis_Summary_g1.5-all-lin.xlsx), detailed PCa results (Ribo_AR_All_Detailed_PCa_g1.5-all-lin.xlsx), detailed non-Ca results (Ribo_AR_All_Detailed_NonCa_g1.5-all-lin.xlsx), PRSS data for PCa (Ribo_AR_All_vs_Iteration_PCa_g1.5-all-lin.xlsx), PRSS data for non-Ca (Ribo_AR_All_vs_Iteration_NonCa_g1.5-all-lin.xlsx), REML convergence analysis (REML_Convergence_Analysis), and REML summary (REML_Summary) |
| 5  | `Part-3.3-REML-Extraction-and-Convergence.r`    | Extracts and analyzes detailed optimization data from REML processes for **interpretability**. Recalculates REML components manually to verify `mgcv`-computed scores. | 1530       | `extract_reml_info()`, `extract_lambda_components()`, `process_reml_info()`, `extract_reml_components_simple()`, `process_reml_components_simple()`, `analyze_reml_discrepancies()`, `extract_reml_optimization_details()` | REML summary and components Excel file (Ribo_AR_REML_Results_g1.5-all-lin.xlsx), REML information with legend (reml_info_with_legend), REML components with discrepancies (all_reml_components), PCa REML information (pca_reml_info), non-Ca REML information (non_ca_reml_info), convergence details Excel file (Ribo_AR_Convergence_Details_g1.5-all-lin.xlsx), and optimization details (all_details) |
| 6  | `Part-3.4-GAM-PRSS-REML-Plots-and-EDF-Analysis.r` | Generates visual plots to analyze GAM modeling of gene expression. Tracks PRSS and REML optimization to understand how $k$ and $\lambda$ are selected, emphasizing **interpretability**. | 873        | `prepare_data()`, `reml_optimization_visualization()`, `generate_visualizations()`, `save_plot_to_files()`, `generate_and_save_visualizations()`, `plot_lambda_values()`, `create_gam_plot()` | PRSS and REML plots for PCa and non-Ca (gam_optimization_pca_*.pdf, gam_optimization_pca_*.jpg, gam_optimization_nonca_*.pdf, gam_optimization_nonca_*.jpg), GAM plots for PCa and non-Ca (Ribo_AR_All_Plot_PCa_g1.5-all_*.pdf, Ribo_AR_All_Plot_NonCa_g1.5-all_*.pdf), EDF summary Excel file (Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx), detailed basis functions report (all_detailed), non-linear relationships summary (nonlinear_summary), non-linear stats (nonlinear_stats), and lambda value plots (output_plots/*.pdf, output_plots/*.jpg) |
| 7  | `Part-3.5-TPRS-Visualization-and-GAM-Components.r` | Visualizes TPRS and GAM components for **interpretability**. Explains spline basis construction around knots, penalization by $\lambda$, weighting by coefficients, and how regularized splines combine to form the full GAM fit. | 1134       | `visualize_tprs_improved()`, `plot_cumulative_components_revised()`, `plot_gam_with_examples()`, `visualize_phi1_emergence()`, `predict()`, `gam()`, `ggplot()` | TPRS basis functions plot (tprs_viz$basis_functions_plot), weighted basis functions plot (tprs_viz$weighted_basis_plot), GAM components plot (tprs_viz$components_plot), variance contribution plots (tprs_viz$importance_page1), cumulative variance plot (tprs_viz$cumulative_variance_plot), cumulative smooth components plot (p_revised$plot), GAM components with examples plot (examples_viz$main_plot) |
| 8  | `Part-3.6-Validation-of-k-and-Lambda-Selection.r` | Validates the selection of $k$ and $\lambda$ using the nested REML approach. Refits models with independent $k$ or $\lambda$ values (without wrapping REML in PRSS, as in the `analyze_sample` function) and performs 10-fold cross-validation (RMSE and RSS) to ensure reproducibility. | 1929       | `analyze_prss_for_k_values()`, `analyze_eigenvalues_for_lambdas()`, `run_kfold_cv_analysis_for_k()`, `run_kfold_cv_analysis_improved()`, `calculate_prss()`, `extract_eigenvalues()`, `perform_kfold_cv_for_k()` | PRSS versus k plot (prss_results$plot), eigenvalue bar plots (create_eigenvalue_plots output), eigenvalue distribution plots (eigen_dist_results), additional eigenvalue comparison plots (additional_plots), CV k versus RMSE and deviance plots (cv_analysis_for_k$cv_plots$rmse_plot, cv_analysis_for_k$cv_plots$deviance_plot), PRSS versus CV k comparison plot (p_comparison), CV lambda versus RMSE and deviance plots (cv_analysis_improved$cv_plots$rmse_plot, cv_analysis_improved$cv_plots$deviance_plot) |
| 9  | `Part-3.7-Gamma-Consequences-on-GAM-Fitting.r`  | Assesses the impact of different $\gamma$ values on GAM fitting. Compares models with varying $\gamma$ to a reference model (default $\gamma$) using 10-fold cross-validation (RMSE). | 472        | `visualize_gam_gammas_free()`, `evaluate_gamma_boxplot()`, `gam()`, `predict()`, `ggplot()`, `pivot_longer()`, `set.seed()` | GAM full fits plot (result$main_plot), difference from gamma = 1 plot (result$difference_plot), selected lambda values bar plot (result$lambda_plot), selected lambda values data (result$selected_lambdas), CV RMSE boxplot (result$cv_plot), CV summary statistics (result$cv_summary), CV RMSE raw data (result$rmse_data) |

---

## License
This project is licensed under the [MIT License](https://github.com/kahkengwong/GAM_PRSS_REML_Project/blob/main/LICENSE), a permissive open-source license designed to encourage collaboration and reuse, while ensuring proper attribution to the original author(s). For the full details, please refer to the [LICENSE](https://github.com/kahkengwong/GAM_PRSS_REML_Project/blob/main/LICENSE) file in this repository.



---

## Contributing
Contributions are welcome! Please open an issue or submit a pull request if you have suggestions or improvements.

---

## Contact
For further information or questions, please email [kahkeng@usm.my](mailto:kahkeng@usm.my)

---
