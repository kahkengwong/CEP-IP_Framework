<div align="center">

# 🧬 CEP-IP: An Explainable Framework for Cell Subpopulation Identification in Single-cell Transcriptomics
### ✨ **CEP-IP** is a novel explainable AI framework that identifies cell subpopulations harboring strong pairwise monotonic **gene-of-interest (GOI)** and **dual-filtered gene (DFG)** module relationships in scRNA-seq data.
### ✨ Strength of the **GOI-DFG module** is quantified by generalized additive model (GAM), and individual cells harboring strong GOI-DFG module are mapped by the **CEP-IP framework**.
### ✨ First demonstrated with the *TRPM4*-Ribo module in prostate cancer, then validated with the *CARM1P1*-DFG module in Allen brain MTG and *FOXM1*-DFG module in Neftel GBM datasets. 

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
- Analysis of multiple scRNA-seq datasets using the **Seurat** package, followed by modeling of pairwise monotonic **GOI-DFGs** relationships via **generalized additive models (GAMs)** implemented in the **mgcv** package.

- The strength of each **GOI-DFGs module** relationship is quantified by **deviance explained (DE)**. This overall DE is then mapped to individual cells via **cell explanatory power (CEP)** classification, identifying **top-ranked explanatory power (TREP)** cells that strongly harbor the GOI-DFGs signal.

- The **CEP-IP framework** (Cell Explanatory Power with Inflection Point) then stratifies the transcriptional space into biologically distinct subpopulations using inflection-point (IP) analysis. DEGs, GO, and Monocle3 analyses then uncover the biology of each subpopulation.

**Three GOI-DFGs module pairings were tested by CEP-IP in this study:**
1. **Prostate cancer (PCa) dataset** - *TRPM4*-Ribo module (7 dual-filtered ribosomal genes averaged as “Ribo”).
2. **Allen Human Middle Temporal Gyrus (MTG) dataset** - *CARM1P1*-DFG module (validation).
3. **Neftel glioblastoma multiforme (GBM) dataset** - *FOXM1*-DFG module (validation).

---

## 🎯Aims of the Project
The project objectives are:
1. To optimize and explain GAM modeling (including *k*, λ, and γ optimization, PRSS/REML convergence, and visualization of thin-plate regression splines).
2. To quantify how much a **GOI** explains variability in its monotonically co-expressed **DFGs** via deviance explained (DE).
3. To identify cells harboring the strongest **GOI-DFGs** relationship via the CEP-IP framework and uncover their distinctive biology through GO enrichment and Monocle3 analysis.

**CEP-IP was validated in two independent brain datasets**, each using a different **GOI-DFGs module pairing**:
- Allen MTG dataset (*CARM1P1*-DFG module)
- Neftel GBM dataset (*FOXM1*-DFG module)

---
 
### 📂 Datasets Used
 
The datasets used in this study are based on publicly available data from the following sources:
 
**1. PCa dataset (GSE185344)**
H.Y. Wong, Q. Sheng, A.B. Hesterberg, S. Croessmann, B.L. Rios, et al., Single cell analysis of cribriform prostate cancer reveals cell intrinsic and tumor microenvironmental pathways of aggressive disease, *Nat Commun*, 13 (2022) 6036. https://doi.org/10.1038/s41467-022-33780-1
 
**2. Allen Human Middle Temporal Gyrus (MTG) dataset**
R.D. Hodge, T.E. Bakken, J.A. Miller, K.A. Smith, E.R. Barkan, et al., Conserved cell types with divergent features in human versus mouse cortex, *Nature*, 573 (2019) 61-68. https://doi.org/10.1038/s41586-019-1506-7
 
**3. Neftel Glioblastoma Multiforme (GBM) dataset**
C. Neftel, J. Laffy, M.G. Filbin, T. Hara, M.E. Shore, et al., An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma, *Cell*, 178 (2019) 835-849 e821. https://doi.org/10.1016/j.cell.2019.06.024

---

## 🔀Workflow of the Project
![Workflow](https://raw.githubusercontent.com/kahkengwong/CEP-IP_Framework/main/Project_Workflow.jpg)

- The workflow begins with identification of the **GOI** and its monotonically co-expressed **DFGs** via Spearman-Kendall dual-filtering, establishing the GOI-DFG module.

- GAM then quantifies the strength of their relationship (DE), and PRSS-REML optimization ensures an appropriately fitted model.
 
- Key methodologies of this study include identification of cells most well-predicted by the model: If a cell is well-predicted, it should have high **explanatory power (EP)**. These cells are termed as **top-ranked EP (TREP) cells** — those most strongly harboring the GOI-DFG module (part (V) of the figure above)
 
- Another key method is to binarize the GOI-DFG transcriptional space by **inflection point (IP)** into **pre-IP and post-IP regions**. These regions exhibit distinct distribution patterns of TREP cells, producing **quadrants of four subpopulations of cells with different biology**, revealed through GO enrichment and Monocle3 analysis.
 
- Collectively, this forms the **CEP-IP framework** detailed in the next section.

---

## ⭐Project Key Findings 
![Project Key Findings](https://raw.githubusercontent.com/kahkengwong/GAM_PCa_Project/main/Project_Key_Findings.jpg)


- **The CEP-IP framework transforms pairwise GOI-DFG module relationships into clinically actionable cell subpopulations, each with distinct biology and trajectories.**
  
---

## 📊scRNA-seq Analysis and GAM Modeling Scripts
Three processed Seurat objects are required to reproduce all analyses and are available on [HuggingFace](https://huggingface.co/datasets/kahkengwong/CEP-IP_Framework/tree/main):
 
| Seurat Object | Size | Dataset | GOI-DFG Module |
|---|---|---|---|
| `GSE185344_Seurat_processed.RData` | 9.52 GB* | PCa (GSE185344) | _TRPM4_-Ribo |
| `AllenMTG_Seurat_processed.RData` | 8.22 GB | Allen Human MTG | _CARM1P1_-DFG |
| `NeftelGBM_SS2_AdultMalignant_Seurat_processed.RData` | 4.02 GB | Neftel GBM | _FOXM1_-DFG |
 
*Note: HuggingFace displays `GSE185344_Seurat_processed.RData` as 9.74 GB due to platform metadata — this is the same file.
 
- The results of the GAM modeling and CEP-IP framework application in the PCa dataset (_TRPM4_-Ribo GOI-DFG module) can be replicated by analyzing `GSE185344_Seurat_processed.RData` following the code blocks `Part_3.01_Mean_Expression_Justifications.r` through `Part_3.15_CEP-IP_in_Monocle3_Trajectory.r`.
 
- Validation of the CEP-IP framework with the _CARM1P1_-DFG module in cortical neurons can be reproduced by analyzing `AllenMTG_Seurat_processed.RData` using `Part_3.16_CEP-IP_Validation_Allen_MTG_dataset.r`.
 
- Validation with the _FOXM1_-DFG module in GBM cells can be reproduced by analyzing `NeftelGBM_SS2_AdultMalignant_Seurat_processed.RData` using `Part_3.17_CEP-IP_Validation_Neftel_GBM_dataset.r`.
 
- For the complete workflow, the scripts should be used in the following sequence, corresponding to the flow of the main project/manuscript:
 
## 📜 Descriptions of the Scripts
| No | Script File | Description |
|------|------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1 | `Part_1_scRNAseq_preprocessing_` <br> `and_UMAP_clusters.r` | Preprocesses scRNA-seq data (PCa, MTG, GBM). Includes QC (gene/ribosomal/mitochondrial filtering, cell-cycle regression, doublet removal, batch correction) and UMAP clustering. |
| 2 | `Part_2_UMAP_Heatmap_Spearman-` <br> `Kendall's-matrix.r` | Visualizes GOI expression in UMAPs and computes Spearman–Kendall dual-filter correlation matrices for GOI–DFG identification. |
| 3 | `Part_3.01_Mean_Expression_` <br> `Justifications.r` | Assesses internal reliability of gene sets (Cronbach's α, McDonald's ω, KMO) for downstream averaging into composite scores. |
| 4 | `Part_3.02_Family_Distribution_` <br> `Analysis.r` | Performs GAM family diagnostics and selects optimal distribution. |
| 5 | `Part_3.03_GAM_REML_PRSS_` <br> `Setup.r` | Sets up PRSS–REML optimization pipeline for GOI–DFG modeling. |
| 6 | `Part_3.04_GAM_REML_PRSS_` <br> `Analysis.r` | Runs full GAM-REML-PRSS optimization across samples and gene sets. |
| 7 | `Part_3.05_REML_Extraction_` <br> `and_Convergence.r` | Extracts detailed REML convergence metrics, gradients, and Hessian values. |
| 8 | `Part_3.06_REML_PRSS_Plots_` <br> `and_EDF_Analysis.r` | Generates PRSS/REML convergence plots and effective degrees-of-freedom reports (FDR-corrected). |
| 9 | `Part_3.07_Validation_of_` <br> `k_and_Lambda_Selection.r` | Validates optimal k and λ choices through refitting and robustness checks. |
| 10 | `Part_3.08_Visualize_TPRS_` <br> `and_GAM_Components.r` | Visualizes thin-plate regression spline basis functions, knots, and model component contributions. |
| 11 | `Part_3.09_Extract_GAM's_` <br> `ND_MD_DE.r` | Extracts null deviance, model deviance, and deviance explained for each GOI–DFG model. |
| 12 | `Part_3.10_Extract_TRPM4-` <br> `Ribo_EP.r` | Computes cell-level explanatory power (EP) for the _TRPM4_–Ribo module. |
| 13 | `Part_3.11_CEP-IP_MCCV_of_` <br> `CEP_Classification.r` | Monte Carlo cross-validation (MCCV) of CEP classification (TREP vs non-TREP) with random and leverage-based controls. |
| 14 | `Part_3.12_CEP-IP_GAM_Plots.r` | Generates GOI–DFG GAM scatter plots highlighting TREP (purple) vs non-TREP (gray) cells. |
| 15 | `Part_3.13_CEP-IP_Mosaic_and_` <br> `Raincloud_Plots.r` | Creates mosaic plots (TREP distribution above/below GAM curve) and raincloud plots for pre-IP vs post-IP regions. |
| 16 | `Part_3.14_CEP-IP_DEGs_Analysis.r` | Performs differential expression (TREP vs non-TREP) within pre-IP and post-IP regions. |
| 17 | `Part_3.15_CEP-IP_in_` <br> `Monocle3_Trajectory.r` | Runs Monocle3 trajectory analysis and quantitative UMAP1 distribution comparisons for the four CEP-IP subpopulations. |
| 18 | `Part_3.16_CEP-IP_Validation_` <br> `Allen_MTG_dataset.r` | Full CEP-IP pipeline on Allen MTG dataset (_CARM1P1_–DFG module). |
| 19 | `Part_3.17_CEP-IP_Validation_` <br> `Neftel_GBM_dataset.r` | Full CEP-IP pipeline on Neftel GBM dataset (_FOXM1_–DFG module). |


## 🛠️Packages and Dependencies
- The R packages and versions used in this study were saved in the `renv.lock` file included in the GitHub repository. This lockfile contains version information for all 37 packages used with their dependencies.
 
- Download or fork the repository and run `renv::restore()` to install the identical package versions used in this study. Note that renv creates an isolated library and will not modify a system's existing R packages and setup.
 
- For manual or selected package installation, a simplified `r_packages_info.json` file is provided with the package names, versions, and sources (CRAN, Bioconductor, or GitHub).

---

## 🚀 Quick Start
1. Download repository (ZIP or git clone)
2. Optional: Install dependencies with `renv::restore()` or selected packages according to `r_packages_info.json`
3. Download all three processed Seurat objects from [HuggingFace](https://huggingface.co/datasets/kahkengwong/CEP-IP_Framework/tree/main):
   - `GSE185344_Seurat_processed.RData` (9.52 GB) — for the main PCa analysis
   - `AllenMTG_Seurat_processed.RData` (8.22 GB) — for the Allen MTG validation
   - `NeftelGBM_SS2_AdultMalignant_Seurat_processed.RData` (4.02 GB) — for the Neftel GBM validation
4. Run scripts `Part_3.01` through `Part_3.15` using `GSE185344_Seurat_processed.RData` for the main PCa analysis (_TRPM4_-Ribo GOI-DFG module)
5. Run `Part_3.16` using `AllenMTG_Seurat_processed.RData` for CEP-IP validation with the _CARM1P1_-DFG module in the Allen MTG dataset
6. Run `Part_3.17` using `NeftelGBM_SS2_AdultMalignant_Seurat_processed.RData` for CEP-IP validation with the _FOXM1_-DFG module in the Neftel GBM dataset


## 🔬 Detailed Setup Instructions (Optional)
<details>
<summary><strong>📜 Click to view all the setup descriptions (7 steps)</strong></summary>

1. **Get the Repository**:
   - **Option A (most users)**: Click the green "<> Code" button, select "Download ZIP", and extract to a folder (e.g., `C:/Users/YourName/CEP-IP_Framework`).
   - **Option B (advanced)**: Fork this repository on GitHub, then clone locally.
 
2. **Open the Project**:
   - Navigate to the extracted folder and double-click any `.r` file to launch RStudio.
 
3. **Set Working Directory**:
   - In RStudio, go to `Session > Set Working Directory > To Source File Location` or run `setwd("path/to/repository")`.
 
4. **Install Packages and Dependencies**:
   - Run `renv::restore()` in the RStudio console to install all 37 required R packages and dependencies as specified in `renv.lock`. This uses an isolated library and preserves your system R setup.
   - For manual installation, refer to `r_packages_info.json` for package names, versions, and sources.
 
5. **Download Data**:
   - Obtain all three processed Seurat objects from [HuggingFace](https://huggingface.co/datasets/kahkengwong/CEP-IP_Framework/tree/main):
     - `GSE185344_Seurat_processed.RData` (9.52 GB) — PCa dataset for the _TRPM4_-Ribo GOI-DFG module analysis.
     - `AllenMTG_Seurat_processed.RData` (8.22 GB) — Allen Human MTG dataset for the _CARM1P1_-DFG module validation.
     - `NeftelGBM_SS2_AdultMalignant_Seurat_processed.RData` (4.02 GB) — Neftel GBM dataset for the _FOXM1_-DFG module validation.
 
6. **Place Data File**:
   - Save the downloaded `.RData` file in the repository root directory (same folder as the scripts).
 
7. **Run Analysis**:
   - Execute the scripts sequentially from `Part_3.01_Mean_Expression_Justifications.r` to `Part_3.15_CEP-IP_in_Monocle3_Trajectory.r` for the main PCa GAM analysis with the _TRPM4_-Ribo GOI-DFG module.
   - This pipeline implements the PRSS-REML optimization and generates CEP-IP quadrants.
   - Run `Part_3.16_CEP-IP_Validation_Allen_MTG_dataset.r` for CEP-IP validation with the _CARM1P1_-DFG module in human cortical neurons.
   - Run `Part_3.17_CEP-IP_Validation_Neftel_GBM_dataset.r` for CEP-IP validation with the _FOXM1_-DFG module in adult malignant GBM cells.
   - Scripts have been stress-tested for consistent and reproducible results.
</details>


## 📋 Supplementary Tables & Results (Zenodo)
**Contents on Zenodo and aims:**
- 17 Supplementary Tables (Excel format) available from [Zenodo](https://zenodo.org/records/17114394)
- Detailed tabular results from analyses conducted in the study
- Refer to these tables to understand the results without running the full pipeline, or to validate your own results
 
**Contents organized by analysis stage:**
- **QC & Clustering** (Supp. Tables 1–2): Cell filtering metrics, cluster statistics
- **GOI–DFG Module Selection** (Supp. Tables 3–5): DFG identification, composite reliability metrics, correlation matrices, GO enrichment
- **GAM Modeling & Optimization** (Supp. Tables 6–10): Model parameters, REML convergence, DE metrics, k- and λ-optimization
- **CEP-IP Classification & GO Enrichment** (Supp. Tables 11–14): MCCV results, HVG analysis, CEP-IP quadrant cell counts, DEG analysis and GO enrichment per subpopulation
- **IP Reliability & Validation** (Supp. Tables 15–17): IPRS scores for all datasets; GAM metrics and GO enrichment for Allen MTG (_CARM1P1_–DFG) and Neftel GBM (_FOXM1_–DFG) validations; GBM GOI screening and within-positive monotonicity results

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
K.K. Wong. CEP-IP: An Explainable Framework for Cell Subpopulation Identification in Single-cell Transcriptomics. arXiv preprint arXiv:2509.12073. (2025) https://arxiv.org/abs/2509.12073

Please also cite the source dataset:

H.Y. Wong, Q. Sheng, A.B. Hesterberg, S. Croessmann, B.L. Rios, et al., Single cell analysis of cribriform prostate cancer reveals cell intrinsic and tumor microenvironmental pathways of aggressive disease, Nat Commun, 13 (2022) 6036. https://doi.org/10.1038/s41467-022-33780-1
R.D. Hodge, T.E. Bakken, J.A. Miller, K.A. Smith, E.R. Barkan, et al., Conserved cell types with divergent features in human versus mouse cortex, Nature, 573 (2019) 61-68. https://doi.org/10.1038/s41586-019-1506-7
C. Neftel, J. Laffy, M.G. Filbin, T. Hara, M.E. Shore, et al., An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma, Cell, 178 (2019) 835-849 e821. https://doi.org/10.1016/j.cell.2019.06.024 

---

## 📩Contact
All analyses, modeling, and interpretations were conducted by Kah Keng Wong: [kahkeng@usm.my](mailto:kahkeng@usm.my)

---
