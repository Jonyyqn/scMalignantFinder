# scMalignantFinder: Malignant cell detection in cancer lineages at single-cell resolution

**scMalignantFinder** is a Python package designed for analyzing cancer single-cell RNA-seq datasets to distinguish malignant cells from their normal counterparts. Trained on over 400,000 high-quality single-cell transcriptomes, scMalignantFinder uses curated pan-cancer gene signatures for calibration and selects features by taking the union of differentially expressed genes across each dataset. More information can be found in the corresponding publication.

![workflow](docs/workflow.png)

# Installation

We recommend using a conda environment to install scMalignantFinder.

1. Create and activate a conda environment

```linux
conda create -n scmalignant python=3.10.10
conda activate scmalignant
```

2. Install `scMalignantFinder` from PyPI:

```shell
pip install scMalignantFinder
```

Optional: scMalignantFinder includes a built-in pan-cancer cell type annotation tool, scATOMIC. If you want to perform basic cell type annotation before identifying malignant cells, follow the [scATOMIC official tutorial](https://github.com/abelson-lab/scATOMIC) to complete its installation in the same conda environment.

# Data preparation

1. **Training data**: Download the training data used in the original study from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/combine_training.h5ad), or use your own dataset to train the model.
2. **Feature file**: The feature list file can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/combined_tumor_up_down_degs.txt).
3. **Example test data**: 
   - Cancer cell line data containing malignant cells can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/test_cancerCellLine.h5ad).
   - Healthy tissue data containing normal epithelial cells can be collected from [here](http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/test_TabulaSapiens.h5ad).

# User guidance

```python
### Load package
from scMalignantFinder import classifier

# Initialize model
model = classifier.scMalignantFinder(
    train_h5ad_path='/path/to/training_data.h5ad',
    feature_path='/path/to/feature_list',
    test_h5ad_path='/path/to/test_data.h5ad', 
    celltype_annotation=False)
# celltype_annotation: If False, the cell type annotation process will not be performed. If True, use scAtomic for cell type annotation

# Model prediction
features = model.fit()
test_adata = model.predict(features)

# View prediction
print(test_adata.obs['scMalignantFinder_prediction'].head())

# Output example:
## Index
## KUL01-T_AAACCTGGTCTTTCAT     Tumor
## KUL01-T_AAACGGGTCGGTTAAC     Tumor
## KUL01-T_AAAGATGGTATAGGGC    Normal
## KUL01-T_AAAGATGGTGGCCCTA     Tumor
## KUL01-T_AAAGCAAGTAAACACA     Tumor
## Name: scMalignantFinder_prediction, dtype: category
## Categories (2, object): ['Tumor', 'Normal']
```



# Citation

