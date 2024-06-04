# scMalignantFinder: Malignant cell detection in cancer lineages at single-cell resolution

scMalignantFinder is a Python package specially designed for analyzing cancer single-cell RNA-seq datasets to distinguish malignant cells from their normal counterparts. scMalignantFinder was trained on >400,000 high-quality single-cell transcriptomes, in which malignant cells were calibrated using curated pan-cancer gene signatures, and feature selection was performed by taking the union of differential expressed genes across each dataset. Please find more information in the corresponding publication.

![workflow](docs/workflow.png)

# Installation

We suggest using a conda environment for installing scMalignantFinder.

1. Create and activate a conda environment

```linux
conda create -n scmalignant python=3.10.10
conda activate scmalignant
```

2. Install `scMalignantFinder` from PyPI:

```shell
pip install scMalignantFinder
```

 Or, install `scMalignantFinder` from github:

``` shell
pip install git+https://github.com/Jonyyqn/scMalignantFinder.git
```

Optional: scMalignantFinder software has a built-in pan-cancer cell type annotation tool scATOMIC. If you want to complete basic annotation of cell types first and then find malignant cells, please follow the official tutorial of scATOMIC (https://github.com/abelson-lab/scATOMIC) to complete its installation in the same conda environment.

# Example

```python
### load package
from scMalignantFinder import classifier

# init model
model = classifier.scMalignantFinder(train_h5ad_path='combine_training.h5ad',
                                     feature_path='combined_tumor_up_down_degs.txt',
                                     test_h5ad_path='/path/to/test_data.h5ad', 
                                     celltype_annotation=False)
# downlaod the training data from http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/combine_training.h5ad
# download the feature file from http://home.ustc.edu.cn/~jonyyqn/scMalignantFinder_data/combined_tumor_up_down_degs.txt
# celltype_annotation: If False, the cell type annotation process will not be performed. If True, use scAtomic for cell type annotation

# model prediction
features = model.fit()
test_adata = model.predict(features)

# view prediction
test_adata.obs['scMalignantFinder_prediction'].head()

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

