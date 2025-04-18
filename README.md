# scMalignantFinder: Distinguishing Malignant Cells in Single-Cell and Spatial Transcriptomics by Leveraging Cancer Signatures

**scMalignantFinder** is a Python package designed for analyzing cancer single-cell RNA-seq and spatial transcriptomics datasets to distinguish malignant cells from their normal counterparts. Trained on over 400,000 high-quality single-cell transcriptomes, scMalignantFinder uses curated pan-cancer gene signatures for training set calibration and selects features by taking the union of differentially expressed genes across each dataset.

![workflow](docs/workflow.png)

# Latest updates

## **Version 1.1.5 (2025-04-20)**

**New spatial region identification feature**

- Add the ability to identify malignant regions from spatial transcriptomics data by applying a clustering-based method.

## **Version 1.0.5 (2025-01-11)**

**Enhanced Flexibility for Test Input**

- Test data can now be provided as a path to an .h5ad file or directly as an AnnData object.

**Dynamic Feature Handling**

- During pretrained model usage, any missing features in the test dataset are temporarily filled with zeros during prediction to ensure compatibility.


## Version 1.0.0 2024-12-24

**New features**

- Introduced malignancy probability output

# Installation

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

# Data preparation

A pretrained model and a list of ordered features are provided in the model directory. Users can also download or use the training data for training the model. 

1. **Training data**: Download the training data used in the original study from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/combine_training.h5ad), or use your own dataset to train the model.
2. **Feature file**: The feature list file can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/combined_tumor_up_down_degs.txt).
3. **Example test data**: 
   - Cancer cell line data containing malignant cells can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/test_cancerCellLine.h5ad).
   - Healthy tissue data containing normal epithelial cells can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/test_TabulaSapiens.h5ad).

# User guidance

## Part 1: Identify malignant cells from scRNA-seq data

```python
### Load package
from scMalignantFinder import classifier

# Initialize model
model = classifier.scMalignantFinder(
  # Path to test data or AnnData object
  test_input="path/to/test_data.h5ad",
  # Directory containing the pretrained model and ordered feature list. If None, the model will be trained. Default: None.
  pretrain_dir=None,
  # Path to training data
  train_h5ad_path="/path/to/train_data.h5ad",
  # Path to feature list
  feature_path="/path/to/features.txt",
  model_method="LogisticRegression",              
  norm_type=True,                           
  n_thread=1
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

## Part 2: Identify malignant regions from spatial transcriptomics

On top of the malignancy probability from Part 1, malignant regions in spatial transcriptomics data can be further identified by integrating gene signatures and image-based features:
```python
from scMalignantFinder import spatial # assuming spatial module provides relevant methods

# Step 1: Calculate AUCell score using scRNA-seq-derived gene sets
sc_gmt = '/path/to/sc_malignant_deg.gmt'
adata = spatial.aucell_cal(adata, sc_gmt)

# Step 2: Extract image-based features from spatial images
adata = spatial.image_cal(adata)

# Step 3: Integrate multi-modal features to segment malignant regions

adata = spatial.region_identification(
  adata,
  # features: A list of feature names from `adata.obs` used for clustering.
  # These features are used to compute pairwise distances and perform hierarchical clustering.
  features=['malignancy_probability', 'Malignant_up', 'image_score'],
  # nclus: Number of clusters to define from hierarchical clustering (default: 3).
  # One of the clusters will be identified as "Malignant", the others as "Normal".
  nclus=3,
  # define_feature: The key feature used to determine which cluster corresponds to malignant regions.
  # The cluster with the highest average value of this feature will be labeled as "Malignant".
  define_feature='Malignant_up',
  # spatial_nn: Whether to refine the predicted region labels using spatial neighbor information (default: True).
  # If True, each spot's label may be adjusted by majority vote from its spatial neighbors.
  spatial_nn=True
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
