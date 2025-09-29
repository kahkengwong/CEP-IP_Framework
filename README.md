<div align="center">

# 🧬 CEP-IP: An Explainable Framework for Cell Subpopulation Identification in Single-cell Transcriptomics
### ✨ Introducing the CEP-IP framework that separates cell subpopulations into quadrants, each with distinct biology.
### ✨ Modeling of _TRPM4_-Ribo relationship with generalized additive model (GAM), and subsequent stratification by the CEP-IP framework. 

![Project Status](https://img.shields.io/badge/status-active-brightgreen?logo=check&logoColor=white)
[![Project Page](https://img.shields.io/badge/Code-GitHub-4E81BE?logo=github&logoColor=white)](https://github.com/kahkengwong/GAM_PCa_Project)
![GitHub languages](https://img.shields.io/github/languages/count/kahkengwong/GAM_PCa_Project)
![GitHub top language](https://img.shields.io/github/languages/top/kahkengwong/GAM_PCa_Project?logo=R&logoColor=486FBA)
![GitHub last commit](https://img.shields.io/github/last-commit/kahkengwong/GAM_PCa_Project)
[![Dataset](https://img.shields.io/badge/Dataset-HuggingFace-B08C00?logo=huggingface&logoColor=F3D34E)](https://huggingface.co/datasets/kahkengwong/GAM_PCa_Project)
[![Paper](https://img.shields.io/badge/Paper-arXiv-red?logo=arxiv&logoColor=white)](https://arxiv.org/abs/2509.12073)

**Author:** **Kah** **Keng** **Wong**  

</div>

---

## 📋Overview
- Analysis of an scRNA-seq dataset of prostate cancer (PCa) and benign prostate (non-Ca) samples using the **Seurat** package, with subsequent modeling via generalized additive models (GAMs) through the **mgcv** package. 

- The GAM fitting is implemented with thin-plate regression splines (TPRS), regularized by penalized residual sum of squares (PRSS) based on the smoothing parameter $\lambda$ derived from restricted maximum likelihood (REML).

- Introducing the **CEP-IP framework** **[cell explanatory power (CEP) with inflection point (IP)]** that subsets the cell subpopulations into quadrants, and Gene Ontology (GO) analysis uncovers distinct biological pathways for each cell subpopulation. 

- **The methodology is generalizable** - **[Fork this repository](https://github.com/kahkengwong/CEP-IP_Framework/fork)** and adapt this framework for your own gene pairs and cell types.

## 🎯Aims of the Project
The following are the project's objectives:
1) To optimize and explain the modeling with mgcv package including *k*, $\lambda$, and $\gamma$ optimization process, PRSS and REML convergence, and visualization of splines' formation.
2) To determine how much does _TRPM4_ explain variability in Ribo (average expression of seven ribosomal genes) expression via deviance explained, GAM's performance metrics. 
3) To identify which cells hold strong _TRPM4_-Ribo relationship via the CEP-IP framework, and to uncover their potentially distinctive biology.

This processed dataset is based on publicly available data from: 
Wong HY, Sheng Q, Hesterberg AB, Croessmann S et al. Single cell analysis of cribriform prostate cancer reveals cell intrinsic and tumor microenvironmental pathways of aggressive disease. Nat Commun 2022;13(1):6036. https://doi.org/10.1038/s41467-022-33780-1

---

## 🔀Workflow of the Project
![Workflow](https://raw.githubusercontent.com/kahkengwong/CEP-IP_Framework/main/Project_Workflow.jpg)

- Key methodologies of this study include identification of cells most well-predicted by the model: If a cell is well-predicted, it should have high **explanatory power (EP)**. These cells are termed as **top-ranked EP (TREP) cells**. These crucial steps are detailed in part (V) of the figure above.

- Another key method is to binarize the transcriptional space by **inflection point (IP)** into **pre-IP and post-IP regions**. These regions exhibit distinct distribution pattern of TREP cells, producing **quadrants of four subpopulation of cells with different biology**.

- Collectively, this forms the **CEP-IP framework** detailed in the next section.  

---

## 📊Project Key Findings 
![Project Key Findings](https://raw.githubusercontent.com/kahkengwong/GAM_PCa_Project/main/Project_Key_Findings.jpg)


- **The CEP-IP framework transforms pairwise gene relationships into clinically actionable cell subpopulations, each with distinct biology and therapeutic potential.**

---

## 📊scRNA-seq Analysis and GAM Modeling Scripts
- The processed Seurat object `GSE185344_Seurat_processed.RData` (9.52 GB) is available on [HuggingFace](https://huggingface.co/datasets/kahkengwong/GAM_PCa_Project/tree/main). *Note: HuggingFace displays this as 9.74 GB due to platform metadata - this is the same file.*

- The results of the GAM modeling in this study can be replicated by analyzing the processed Seurat object `GSE185344_Seurat_processed.RData` by following the code block `Part_3.01_Mean_Expression_Justifications.r` until `Part_3.15_Monocle3_Pre-IP_vs_Post-IP_TREP.r`

- For the complete workflow, the scripts should be used in the following sequence, corresponding to the flow of the main project/manuscript:

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
- The R packages and versions used in this study were saved in the `renv.lock` file included in the GitHub repository. This lockfile contains version information for all 37 packages used with their dependencies.

- Download or fork the repository and run `renv::restore()` to install the identical package versions used in this study. Note that renv creates an isolated library and will not modify a system's existing R packages and setup. 

- For manual or selected package installation, a simplified `r_packages_info.json` file is provided with the package names, versions, and sources (CRAN, Bioconductor, or GitHub).

---

## 🚀 Quick Start
1. Download repository (ZIP or git clone)
2. Optional: Install dependencies with `renv::restore()` or selected packages according to `r_packages_info.json`
3. Download processed Seurat object from [HuggingFace](https://huggingface.co/datasets/kahkengwong/GAM_PCa_Project/tree/main) 
4. Run scripts Part_3.01 through Part_3.15


## 🔬 Detailed Setup Instructions (Optional)
<details>
<summary><strong>📜 Click to view all the setup descriptions (7 steps)</strong></summary>

1. **Get the Repository**:
   - **Option A (most users)**: Click the green "<> Code" button, select "Download ZIP", and extract to a folder (e.g., `C:/Users/YourName/GAM_PCa_Project`).
   - **Option B (advanced)**: Fork this repository on GitHub, then clone locally.

2. **Open the Project**:
   - Navigate to the extracted folder and double-click any `.r` file to launch RStudio.

3. **Set Working Directory**:
   - In RStudio, go to `Session > Set Working Directory > To Source File Location` or run `setwd("path/to/repository")`.

4. **Install Packages and Dependencies**:
   - Run `renv::restore()` in the RStudio console to install all 37 required R packages and dependencies as specified in `renv.lock`. This uses an isolated library and preserves your system R setup.
   - For manual installation, refer to `r_packages_info.json` for package names, versions, and sources.

5. **Download Data**:
   - Obtain the processed Seurat object (`GSE185344_Seurat_processed.RData`) from [HuggingFace](https://huggingface.co/datasets/kahkengwong/GAM_PCa_Project/tree/main).

6. **Place Data File**:
   - Save the downloaded `.RData` file in the repository root directory (same folder as the scripts).

7. **Run Analysis**:
   - Execute the scripts sequentially from `Part_3.01_Mean_Expression_Justifications.r` to `Part_3.15_Monocle3_Pre-IP_vs_Post-IP_TREP.r` for the main GAM analysis.
   - This pipeline implements the PRSS-REML optimization and generates CEP-IP quadrants.
   - Scripts have been stress-tested for consistent and reproducible results.
</details>


## 📋 Supplementary Tables & Results (Zenodo)
**Contents on Zenodo and aims:** 
- 13 Supplementary Tables (Excel format) available from [Zenodo](https://zenodo.org/records/17114394)
- Detailed tabular results from analyses conducted in the study
- Refer these tables to understand the results without running the full pipeline, or to validate your own results

**Contents organized by analysis stage:**
- **scRNA-seq QC & Clustering** (Supp. Tables 1-2): Cell filtering metrics, cluster statistics
- **Gene Selection & Enrichment** (Supp. Tables 3-5): Gene set reliability, correlation matrices, GO enrichment
- **GAM Modeling & Optimization** (Supp. Tables 6-10): Model parameters, REML convergence, deviance metrics, k-optimization, λ-optimization
- **Cell Classification & GO Enrichment** (Supp. Tables 11-13): TREP vs non-TREP comparisons, DEG analysis, pathway enrichment

---

## 💻 Requirements
### Software Requirements
- **R:** ≥ 4.0.0 (tested with R 4.4.2 "Pile of Leaves")
- **RStudio:** Recommended (tested with RStudio 2025.05.0 Build 496)

### Hardware Requirements
- **RAM:** 16GB (32GB recommended)
- **CPU:** 8 cores (12+ cores recommended)
- **Storage:** ~12GB free disk space for data and results

### Development Environment
- **Developed and tested on:** Intel Core i9-14900KF, 64GB DDR5 RAM, RTX 4090, Windows 11 Pro, RStudio 2025.05.0 Build 496, R 4.4.2

### Performance Notes
- All analyses are CPU-based; GPU acceleration not required
- Higher specs are primarily for scRNA-seq preprocessing (e.g., doublet removal, Kendall's τ computation)
- GAM modeling with [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html) is computationally efficient and does not require high-end specs
- Understanding the algorithms and mathematical principles is key; minimum specs are sufficient for GAM analysis using the processed Seurat object
> [!NOTE]
> 💡**To students and young researchers:** Limited computational resources should never prevent one from learning the key methodologies and mathematical principles in machine learning; understanding them matters more than hardware specs.

---

## 🔬 Adapting CEP-IP for Your Research
This framework is designed to be generalizable beyond prostate cancer. Consider forking if you want to:

- **Apply to different cancer types** - Adapt for breast, lung, colorectal, or other cancers
- **Explore different gene pairs** - Replace _TRPM4_-Ribo with your genes of interest
- **Integrate with your pipeline** - Use CEP-IP as a module in larger workflows
- **Benchmark against other methods** - Compare with your current cell classification approaches
- **Educational purposes** - Learn GAM, PRSS, and REML methodologies with working code

**[Fork this repository](https://github.com/kahkengwong/CEP-IP_Framework/fork)** to start customizing for your needs. 

---

## 🧾License
This project is licensed under the [MIT License](https://github.com/kahkengwong/GAM_PCa_Project/blob/main/LICENSE), an open-source license to encourage collaboration and reuse, while ensuring proper attribution to the original author(s). For the full details, please refer to the [LICENSE](https://github.com/kahkengwong/GAM_PCa_Project/blob/main/LICENSE) file in this repository.


---

## 🤝🏻 Contributing
Contributions, issues, and feature requests are welcome!

### 🐛 Found a Bug?
[Open an issue](https://github.com/kahkengwong/CEP-IP_Framework/issues) describing the problem, including your error message, R version, and operating system.

### 💡 Want to Contribute Improvements?
If you've made enhancements you'd like to share with others, you can contribute back via a **pull request** (i.e., a way to propose merging your improvements back into the main repository so others can benefit).

**How to submit:**
1. Make changes in your forked copy
2. Click "Contribute" → "Open pull request" on your fork's GitHub page
3. Describe your changes and submit

I'll review and may ask questions before merging. **Or simply:** [Open an issue](https://github.com/kahkengwong/CEP-IP_Framework/issues) describing your idea.

---

## 📚 Citation
If you found this study useful, please cite the preprint:

**BibTeX:**
```bibtex
@misc{wong2025cep-ip,
      title={CEP-IP: An Explainable Framework for Cell Subpopulation Identification in Single-cell Transcriptomics}, 
      author={Kah Keng Wong},
      year={2025},
      eprint={2509.12073},
      archivePrefix={arXiv},
      primaryClass={q-bio.GN},
      url={https://arxiv.org/abs/2509.12073}, 
}
```
Wong KK (2025). CEP-IP: An Explainable Framework for Cell Subpopulation Identification in Single-cell Transcriptomics. arXiv preprint arXiv:2509.12073. https://arxiv.org/abs/2509.12073

Please also cite the source dataset:
Wong HY, Sheng Q, Hesterberg AB, Croessmann S et al (2022). Single cell analysis of cribriform prostate cancer reveals cell intrinsic and tumor microenvironmental pathways of aggressive disease. Nat Commun 13(1):6036. https://doi.org/10.1038/s41467-022-33780-1

---

## 📩Contact
All analyses, modeling, and interpretations were conducted by Kah Keng Wong: [kahkeng@usm.my](mailto:kahkeng@usm.my)

---
