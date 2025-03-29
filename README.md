# GAM_PRSS_REML_Project
GAM with PRSS and REML in scRNA-seq PCa dataset

# PDAC ML Paper Repository
![GitHub last commit](https://img.shields.io/github/last-commit/kahkengwong/GAM_PRSS_REML_Project)
![GitHub repo size](https://img.shields.io/github/repo-size/kahkengwong/GAM_PRSS_REML_Project)
![GitHub languages](https://img.shields.io/github/languages/count/kahkengwong/GAM_PRSS_REML_Project)
![GitHub top language](https://img.shields.io/github/languages/top/kahkengwong/GAM_PRSS_REML_Project)
![Contributors](https://img.shields.io/github/contributors/kahkengwong/GAM_PRSS_REML_Project)

Author: Kah Keng Wong  

---

## ML Scripts

| Step | Script File                                      | Description                                                                                     |
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




Scripts should be used in the following sequence, corresponding to the flow of the main project/manuscript:
1. `Part_1_scRNAseq_preprocessing_&_UMAP_clusters.r` - Scripts to conduct scRNA-seq dataset preprocessing (prostate cancer and benign samples; QC steps such as removing of low-quality cells, regressing out cell cycle phases effects, batch effects correction) and UMAP clustering.
   
2. `Part_2_UMAP_Heatmap_&_SK-matrix.r` - scRNA-seq analysis by visualizing _TRPM4_ and _KLK4_ across different clusters in prostate cancer and benign samples; Calculates Sperman's and Kendall's correlation matrices; Generates heatmap to compare _TRPM4_ vs other relevant gene sets.
   
3. `Part_3.1_GAM_PRSS_REML_Setup.r` - Setup for generalized additive model (GAM) and restricted maximum likelihood (REML), and to assess model convergence.
   
4. `Part_3.2_GAM_PRSS_REML_Analysis.r` - GAM analysis codes, and to extract best-fitting models and the summary of the key statistics including PRSS and REML convergence.
   
5. `Part_3.3_REML_Extraction_&_Convergence.r` - Extracts and analyzes detailed optimization data REML processes with the aim for **interpretability**. This is achieved by extracting the values of each component of REML formula, and to manually recalculate the relevant scores, ensuring that they tally with mgcv-computed scores.
    
6. `Part_3.4_GAM-PRSS-REML_Plots_&_EDF_Analysis.r` - Visual plots to analyze GAM modeling of gene expression in prostate cancer and benign samples. To track optimization processes i.e., PRSS and REML, with the main aim for **interpretability** by tracking how each PRSS and REML reaches their selected $k$ and $\lambda$ values, respectively.
    
7. `Part_3.5_TPRS_Visualization_and_GAM_Components.r` - Visualization of TPRS and GAM components. This is also aimed at **interpretability** by understanding how each spline basis is initially constructed surrounding the knots, and how each spline is penalized by $\lambda$ and weighted by their coefficient, and how each of these regularized splines are combined to produce the smooth terms that ultimately combine with linear terms to produce the full GAM fit.
    
8. `Part_3.6_Validation_of_k_&_Lambda_Selection.r` - Validates the selection of key parameters ($k$ and $\lambda$) by the nested REML approach. Validation is conducted by refitting using $k$ or $\lambda$ values independently (without wrapping REML within PRSS as conducted in `analyze_sample` function), and also by 10-fold cross-validation (RMSE and RSS), ensuring the initial approach adopted by `analyze_sample` function is reproducible. 

9. `Part_3.7_Gamma_Consequences_on_GAM_Fitting.r` - This block of codes aims to assess how different $\gamma$ values affect the fit of GAMs. The GAM fitting is conducted using different $\gamma$ values, compares them to the reference model that adopts the default $\gamma$ value, and to evaluate their performance via 10-fold cross-validation (RMSE). 

---

## Contact
For further information or questions, please email [kahkeng@usm.my](mailto:kahkeng@usm.my)

---
