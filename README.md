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

| Step | Script File                                      | Description                                                                                     |
|------|--------------------------------------------------|-------------------------------------------------------------------------------------------------|
| 1    | `Part_1_scRNAseq_preprocessing_&_UMAP_clusters.r` | Preprocesses scRNA-seq data (prostate cancer and benign samples). Includes QC steps: removing low-quality cells, regressing out cell cycle phase effects, and correcting batch effects. Performs UMAP clustering to identify cell populations. |
| 2    | `Part_2_UMAP_Heatmap_&_SK-matrix.r`           | Visualizes *TRPM4* and *KLK4* expression across UMAP clusters in prostate cancer and benign samples. Computes Spearman’s and Kendall’s correlation matrices and generates a heatmap to compare *TRPM4* against other relevant gene sets. |
| 3    | `Part_3.1_GAM_PRSS_REML_Setup.r`                | Sets up the Generalized Additive Model (GAM) with Restricted Maximum Likelihood (REML) and assesses model convergence. |
| 4    | `Part_3.2_GAM_PRSS_REML_Analysis.r`             | Performs GAM analysis, extracts best-fitting models, and summarizes key statistics, including PRSS and REML convergence. |
| 5    | `Part_3.3_REML_Extraction_&_Convergence.r`    | Extracts and analyzes detailed optimization data from REML processes for **interpretability**. Recalculates REML components manually to verify `mgcv`-computed scores. |
| 6    | `Part_3.4_GAM-PRSS-REML_Plots_&_EDF_Analysis.r` | Generates visual plots to analyze GAM modeling of gene expression. Tracks PRSS and REML optimization to understand how $k$ and $\lambda$ are selected, emphasizing **interpretability**. |
| 7    | `Part_3.5_TPRS_Visualization_and_GAM_Components.r` | Visualizes TPRS and GAM components for **interpretability**. Explains spline basis construction around knots, penalization by $\lambda$, weighting by coefficients, and how regularized splines combine to form the full GAM fit. |
| 8    | `Part_3.6_Validation_of_k_&_Lambda_Selection.r` | Validates the selection of $k$ and $\lambda$ using the nested REML approach. Refits models with independent $k$ or $\lambda$ values (without wrapping REML in PRSS, as in the initial `analyze_sample` function) and performs 10-fold cross-validation (RMSE and RSS) to ensure reproducibility. |
| 9    | `Part_3.7_Gamma_Consequences_on_GAM_Fitting.r`  | Assesses the impact of different $\gamma$ values on GAM fitting. Compares models with varying $\gamma$ to a reference model (default $\gamma$) using 10-fold cross-validation (RMSE). |

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
