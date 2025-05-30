[![PyPI](https://img.shields.io/pypi/v/scMalignantFinder)](https://pypi.org/project/scMalignantFinder)
[![PyPI Downloads](https://img.shields.io/pepy/dt/scMalignantFinder?logo=pypi)](https://pepy.tech/project/scMalignantFinder)
[![Stars](https://img.shields.io/github/stars/Jonyyqn/scMalignantFinder?style=flat&logo=GitHub&color=blue)](https://github.com/Jonyyqn/scMalignantFinder/stargazers)

# scMalignantFinder: Distinguishing Malignant Cells in Single-Cell and Spatial Transcriptomics by Leveraging Cancer Signatures

**scMalignantFinder** is a Python package designed for analyzing cancer single-cell RNA-seq and spatial transcriptomics datasets to distinguish malignant cells from their normal counterparts. Trained on over 400,000 high-quality single-cell transcriptomes, scMalignantFinder uses curated pan-cancer gene signatures for training set calibration and selects features by taking the union of differentially expressed genes across each dataset.

![workflow](docs/workflow.png)

## Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->

<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Latest updates](#latest-updates)
- [Installation](#installation)
- [Data preparation](#data-preparation)
- [User guidance](#user-guidance)
  - [Identify malignant cells from scRNA-seq data](#identify-malignant-cells-from-scRNA-seq-data)
  - [Identify malignant regions from spatial transcriptomics](#identify-malignant-regions-from-spatial-transcriptomics)
  - [Analyze cancer cell states using curated gene sets](#Analyze-cancer-cell-states-using-curated-gene-sets)
  
<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Latest updates
### **Version 1.1.7 (2025-05-27)**
- Enhanced downstream analysis capabilities with a new module for cancer cell state analysis

### **Version 1.1.6 (2025-04-26)**

**Expanded input file format compatibility**

- Extended support for test data input to include `.txt`, `.tsv`, and `.csv` files, as well as their compressed `.gz` versions.

### **Version 1.1.5 (2025-04-20)**

**New spatial region identification feature**

- Add the ability to identify malignant regions from spatial transcriptomics data by applying a clustering-based method.

### **Version 1.0.5 (2025-01-11)**

**Enhanced Flexibility for Test Input**

- Test data can now be provided as a path to an .h5ad file or directly as an AnnData object.

### Version 1.0.0 (2024-12-24)

**New features**

- Introduced malignancy probability output

## Installation

We recommend using a conda environment to install scMalignantFinder.

1. Create and activate a conda environment

```bash
conda create -n scmalignant python=3.10.10
conda activate scmalignant
```

2. Install `scMalignantFinder` from PyPI:

```bash
pip install scMalignantFinder
```

Optional: scMalignantFinder includes a built-in pan-cancer cell type annotation tool, scATOMIC. If you want to perform basic cell type annotation before identifying malignant cells, follow the [scATOMIC official tutorial](https://github.com/abelson-lab/scATOMIC) to complete its installation in the same conda environment.

## Data preparation

A pretrained model and a list of ordered features are provided in the `model` directory:
- Pretrained model: `model/model.joblib`
- Ordered feature list: `model/ordered_feature.tsv`

Users can also download or use the training data for training the model. 
1. **Training data**: Download the training data used in the original study from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/combine_training.h5ad), or use your own dataset to train the model.
2. **Feature file**: The feature list file can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/combined_tumor_up_down_degs.txt).
3. **Example test data**: 
   - Cancer cell line data containing malignant cells can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/test_cancerCellLine.h5ad).
   - Healthy tissue data containing normal epithelial cells can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/test_TabulaSapiens.h5ad).

In addition, for malignant region identification in spatial transcriptomics data, the malignant cell gene signature file is provided as `model/sc_malignant_deg.gmt`.
Alternatively, users can supply their own gene signature file, provided it follows the same .gmt format.

## User guidance

### Identify malignant cells from scRNA-seq data

```python
### Load package
from scMalignantFinder import classifier

# Initialize model
model = classifier.scMalignantFinder(
    test_input="path/to/test_data.h5ad",  # Test data input. Supports:
                                          # - AnnData object
                                          # - Path to .h5ad file
                                          # - Path to tab-delimited .txt/.tsv file
                                          # - Path to comma-delimited .csv file
                                          # (.txt, .tsv, and .csv can also be gzip-compressed as
                                          # .txt.gz, .tsv.gz, or .csv.gz)
                                          # For text-based files (.txt, .tsv, .csv), 
                                          # rows must correspond to gene symbols
                                          # and columns to cell barcodes.
                                          
    pretrain_dir=None,                    # Directory containing pretrained model and feature list. 
                                          # If None, the model will be trained from scratch. Default: None.
                                           
    train_h5ad_path="/path/to/train_data.h5ad",  # Path to training dataset (.h5ad)
    feature_path="/path/to/features.txt",        # Path to ordered feature list (.txt)
    model_method="LogisticRegression",           # Classifier method: LogisticRegression, RandomForest, or XGBoost
    norm_type=True,                              # Whether to normalize test data. Default: True.
    n_thread=1                                   # Number of threads for parallel processing
)

# Load data
model.load()

# Predict malignancy
result_adata = model.predict()

# View results
print(result_adata.obs["scMalignantFinder_prediction"].head())

## Example output for scMalignantFinder_prediction:
Index
KUL01-T_AAACCTGGTCTTTCAT     Malignant
KUL01-T_AAACGGGTCGGTTAAC     Malignant
KUL01-T_AAAGATGGTATAGGGC     Normal
KUL01-T_AAAGATGGTGGCCCTA     Malignant
KUL01-T_AAAGCAAGTAAACACA     Malignant
Name: scMalignantFinder_prediction, dtype: category
Categories (2, object): ['Normal', 'Malignant']

print(result_adata.obs["malignancy_probability"].head())
## Example output for malignancy_probability:
Index
KUL01-T_AAACCTGGTCTTTCAT     0.98578
KUL01-T_AAACGGGTCGGTTAAC     0.78968
KUL01-T_AAAGATGGTATAGGGC     0.243564
KUL01-T_AAAGATGGTGGCCCTA     0.8796
KUL01-T_AAAGCAAGTAAACACA     0.6598
Name: malignancy_probability, dtype: float64
```

### Identify malignant regions from spatial transcriptomics

On top of the malignancy probability from [Identify malignant cells from scRNA-seq data](#identify-malignant-cells-from-scRNA-seq-data), malignant regions in spatial transcriptomics data can be further identified by integrating gene signatures and image-based features:
```python
from scMalignantFinder import spatial, utils # assuming spatial module provides relevant methods

# Step 1: Calculate AUCell score using scRNA-seq-derived gene sets
sc_gmt = './model/sc_malignant_deg.gmt'
adata = utils.aucell_cal(adata, sc_gmt)

# Step 2: Extract image-based features from spatial images
adata = spatial.image_cal(adata)

# Step 3: Integrate multi-modal features to segment malignant regions

adata = spatial.region_identification(
    adata,
    
    features=['malignancy_probability', 'Malignant_up', 'image_score'],  
    # A list of feature names from adata.obs used for clustering.
    # These features are used to compute pairwise distances and perform hierarchical clustering.
    
    nclus=3,  
    # Number of clusters to define from hierarchical clustering (default: 3).
    # One of the clusters will be identified as "Malignant", and the others as "Normal".
    
    define_feature='Malignant_up',  
    # The key feature used to determine which cluster corresponds to malignant regions.
    # The cluster with the highest average value of this feature will be labeled as "Malignant".
    
    spatial_nn=True  
    # Whether to refine region labels using spatial neighbor information (default: True).
    # If True, each spot's label may be adjusted by majority vote among its spatial neighbors.
)

# Example output for scMalignantFinder_prediction:
print(adata.obs[['cluster','region_prediction']])
                   cluster region_prediction
AAACAAGTATCTCCCA-1       0            Normal
AAACACCAATAACTGC-1       1         Malignant
AAACAGAGCGACTCCT-1       2            Normal
AAACAGGGTCTATATT-1       0            Normal
AAACAGTGTTCCTGGG-1       1         Malignant

```

### Analyze cancer cell states using curated gene sets

To support downstream functional interpretation of malignant cells, we collected and curated **67 cancer cell state gene sets** from [a pan-cancer study](https://doi.org/10.1038/s41586-023-06130-4). These gene sets represent a wide spectrum of cancer-associated cellular programs (e.g., cell cycle, EMT, immune evasion, hypoxia) across multiple cancer types.

You can quantify the enrichment of these gene sets in individual cells using AUCell scoring:

```python
from scMalignantFinder import utils

# Path to the pan-cancer curated gene sets
pan_cancer_gene_sets = "/path/to/model/Malignant_MPs.Gavish_2023.gmt"

# Compute AUCell scores for each gene set
adata = utils.aucell_cal(adata, pan_cancer_gene_sets, norm_type=False)

# View results
print(adata.obs.loc[:,adata.obs.columns.str.startswith('MP')].iloc[:5,:3]) 
                MP1 Cell Cycle - G2/M	MP2 Cell Cycle - G1/S	MP3 Cell Cylce HMG-rich				
KUL01-T_AAACCT	             0.045819	             0.000000	               0.306887
KUL01-T_AAACGG	             0.155027	             0.078003	               0.227548
KUL01-T_AAAGAT	             0.000000	             0.000000	               0.293480
KUL01-T_AAAGAG	             0.000000	             0.000000	               0.239118
KUL01-T_AAAGCA	             0.068728	             0.000000	               0.272176


```

# Citation

If you use scMalignantFinder in your research, please cite:

Qiaoni, Yu, et al. "scMalignantFinder distinguishes malignant cells in single-cell and spatial transcriptomics by leveraging cancer signatures.", *Communications Biology*, 2025.
DOI: [https://doi.org/10.1038/s42003-025-07942-y](https://doi.org/10.1038/s42003-025-07942-y)

BibTeX:
```bibtex
@article{yu2025scmalignantfinder,
  title={scMalignantFinder distinguishes malignant cells in single-cell and spatial transcriptomics by leveraging cancer signatures},
  author={Yu, Qiaoni and Li, Yuan-Yuan and Chen, Yunqin},
  journal={Communications Biology},
  volume={8},
  number={1},
  pages={504},
  year={2025},
  publisher={Nature Publishing Group UK London}
}


