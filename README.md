# Temporal Drift Confound Correction in Multi-Cohort Microbiome Studies

A six-layer cross-validated machine learning framework for detecting and
correcting temporal drift confounds in harmonised multi-cohort 16S rRNA
microbiome datasets.

**Reference:**  
Desh SS, Jose N, Anudeep PP, Nayak S. Temporal drift confound correction in
multi-cohort microbiome studies: a cross-validated machine learning framework.
*Bioinformatics Advances*. 2026; (in press). <!-- to update after acceptance -->

---

## Data Access

This repository contains **code only**. No sequence data or OTU tables are
included.

- **Primary dataset:** Rodriguez et al. (2023), BioProject
  [PRJNA891951](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA891951)
  — 11 dietary-fibre intervention studies, 2,283 samples, 529 subjects,
  9,612 OTUs.
- **SCFA annotation reference:** Frolova et al. (2022), Table 1. Download the
  file (`Table 1.XLSX`) from the original publication and place it in the
  repository root. This file is not redistributable and is excluded from the
  repo via `.gitignore`.

## Requirements

### Python

Python ≥ 3.10. Install dependencies:

```bash
pip install -r requirements.txt
```

**PyTorch note:** The pinned version pulls CPU-only wheels. If you need a
different build:

```bash
pip install torch==2.2.0 --index-url https://download.pytorch.org/whl/cpu
```

### R

R ≥ 4.3. Install the following packages:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr",
  "lme4", "lmerTest", "pbapply",
  "tempted", "ggplot2", "pROC"
))
```

---

## Execution Order

The pipeline has a branching dependency structure. Steps marked ∥ can run in
parallel once their prerequisites are met.

| Step | Notebook / Script | Depends on |
|------|-------------------|------------|
| 1 | `01_data_prep.ipynb` | — |
| 2 | `02_preprocessing.ipynb` | 1 |
| 3 | `03_random_forest.ipynb` | 2 |
| 4 | `delta_clr_classification.ipynb` | 2 |
| 5 ∥ | `lme_analysis.R` | 4 |
| 6 ∥ | `tempted_analysis.R` | 4 |
| 7 | `07_tempted_fibertype_replot.ipynb` | 6 |
| 8 | `08_lme_shap_convergence.ipynb` | 4, 5 |
| 9 ∥ | `04_autoencoder_dim_reduction.ipynb` | 2 |
| 10 ∥ | `05_clustering.ipynb` | 9 |
| 11 | `06_nutrigenomic_pathway_layer.ipynb` | 4, 5, 10 |
| 12 | `09_responder_analysis.ipynb` | 2, 4 |

---

## Repository Structure

```
├── 01_data_prep.ipynb
├── 02_preprocessing.ipynb
├── 03_random_forest.ipynb
├── 04_autoencoder_dim_reduction.ipynb
├── 05_clustering.ipynb
├── 06_nutrigenomic_pathway_layer.ipynb
├── 07_tempted_fibertype_replot.ipynb
├── 08_lme_shap_convergence.ipynb
├── 09_responder_analysis.ipynb
├── delta_clr_classification.ipynb
├── lme_analysis.R
├── tempted_analysis.R
├── requirements.txt
├── .gitignore
├── LICENSE
└── README.md
```

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE).

## Contact

Sachin S Desh (corresponding author)  
ICMR-National Institute of Child Health and Development Research, New Delhi  
<!-- drsachinsdesh@gmail.com -->
