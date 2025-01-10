import os
import sys
import joblib
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from sklearn.linear_model import LogisticRegression, Lasso
from sklearn.ensemble import RandomForestClassifier
import warnings
from sklearn.exceptions import ConvergenceWarning
# Ignore only ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)

def init_core_model(model_method, n_thread):
    """Method to set the underlying core model.

    The core model is used to make first predictions 
    based on ranked scores.
    """

    if model_method == "LogisticRegression":
        return LogisticRegression(n_jobs=n_thread)
    elif model_method == "RandomForest":
        return RandomForestClassifier(n_jobs=n_thread)
    elif model_method == "XGBoost":
        from xgboost.sklearn import XGBClassifier
        return XGBClassifier(n_jobs=n_thread)
    else:
        raise NotImplementedError(
            f"core model '{model_method}' is not yet implemented. "
            f"Please try 'LogisticRegression, RandomForest, or XGBoost' instead."
        )

def data_preprocess(h5ad_path, norm=True):
    adata = sc.read_h5ad(h5ad_path)
    if norm==True: sc.pp.normalize_total(adata, target_sum=1e4)
    
    return adata

def load_feature(feature_path):
    features = list(np.loadtxt(feature_path, dtype=str))
    
    return features

def scatomic_celltype_annotation(test_h5ad_path, n_thread):
    tmpdir = os.getcwd()
    # print(tmpdir)
    # Define the path to the R script
    script_dir = os.path.dirname(__file__)
    r_script = os.path.join(script_dir, 'scAtomic.r')
    
    # Call the R script using subprocess
    result = subprocess.run(
        ['Rscript', r_script, test_h5ad_path, str(n_thread), tmpdir],
        capture_output=True,
        text=True
    )
    
    # Check for errors
    if result.returncode != 0:
        raise RuntimeError(f'Error running R script: {result.stderr}')
    
    # print(result.stdout)

    celltype_annotation_df = pd.read_csv(f'{tmpdir}/celltype_annotation.scAtomic.csv', index_col=0)

    return celltype_annotation_df

class scMalignantFinder:
    """scMalignantFinder class
    
    scMalignantFinder holds the main modelling functionalities. It's basic usage
    is based on the typical scikit-learn workflow. That means 
    1. load data, 
    2. initialize a model, 
    3. fit the model and 
    4. make the actual predictions on unknown data. 
    In general, annotated data objects are used as data format (AnnData).
    
    Parameters
    -----------
    test_h5ad_path : str
        Path of the test data.
    pretrain_path: str
        Path of the pretrained model.
        Default: None
    train_h5ad_path : str
        Path of the training data.
        Default: /mnt/home/qnyu/workspace/scOmics/malignantModel/multi_tissue/data/combine_training.h5ad.
    celltype_annotation : bool
        If False, the cell type annotation process will not be performed. If True, use scAtomic for cell type annotation
        Default: False
    output_annotation : str
        Column name of AnnData storing the prediction.
        Default: Malignant_prediction.
    model_method : str
        Core model algorithm. Currently supported: "LogisticRegression, RandomForest, XGBoost".
        Default: LogisticRegression.
    feature_path : str
        Path storing the feature list.
        Default: /mnt/home/qnyu/workspace/scOmics/malignantModel/multi_tissue/file/deg/combined_tumor_up_down_degs.txt.
    norm_type: bool
        If True, Normalize the test data.
        Default: True
    n_thread : int
        Number of multiprocessing threads.
        Default: 1
    probability: bool
        If True, Output the prediction probability for each class.
        Default: False
    out_path: str
        Path storing the prediction probability result.
        Default: None
    """
    
    def __init__(
        self,
        test_h5ad_path,
        pretrain_path=None,
        train_h5ad_path="./combine_training.h5ad",
        celltype_annotation=False,
        output_annotation="scMalignantFinder_prediction",
        model_method="LogisticRegression",
        feature_path="./combined_tumor_up_down_degs.txt",
        norm_type=True,
        n_thread=1,
        probability=False,
        out_path=None
    ):
        self.test_h5ad_path = test_h5ad_path
        self.pretrain_path = pretrain_path
        self.train_h5ad_path = train_h5ad_path
        self.celltype_annotation = celltype_annotation
        self.output_annotation = output_annotation
        self.model_method = model_method
        self.feature_path = feature_path
        self.norm_type = norm_type
        self.n_thread = 1
        self.fitted = False
        self.probability = probability
        self.out_path = out_path
        
    
    def fit(self):
        label_col = 'Raw_annotation'
        label_dict = {'Normal':0,'Tumor':1}
        r_label_dict = dict(zip(label_dict.values(), label_dict.keys()))

        if self.pretrain_path == None:
            train_adata = data_preprocess(self.train_h5ad_path)
            features = load_feature(self.feature_path)
            filter_features = list(set(features) & set(sc.read_h5ad(self.test_h5ad_path).var_names))
            self.core_model = init_core_model(self.model_method, self.n_thread)
            if len(filter_features)>0:
                _ = self.core_model.fit(train_adata[:,filter_features].X, train_adata.obs[label_col].map(label_dict).tolist())
                self.fitted = True
            else:
                print('The number of features is 0. Please check whether the features of the test set are gene symbols.')
                self.fitted = False

            return filter_features
        
        else:
            self.core_model = joblib.load(os.path.join(self.pretrain_path,'model.joblib'))
            features = load_feature(os.path.join(self.pretrain_path,'ordered_feature.tsv'))
            filter_features = list(set(features) & set(sc.read_h5ad(self.test_h5ad_path).var_names))
            if len(filter_features)>0:
                self.fitted = True
            else:
                print('The number of features is 0. Please check whether the features of the test set are gene symbols.')
                self.fitted = False

            return features
            

    def predict(self, features):
        label_dict = {'Normal':0,'Tumor':1}
        r_label_dict = dict(zip(label_dict.values(), label_dict.keys()))
        
        if not self.fitted:
            raise RuntimeError("Model not yet fitted. Please run Model.fit(...) first!")
        
        test_adata = data_preprocess(self.test_h5ad_path, norm=self.norm_type)
        if self.pretrain_path == None:
            y_pred = self.core_model.predict(test_adata[:,features].X)
        else:
            df = test_adata[:,list(set(features) & set(test_adata.var_names))].to_df()
            sub_matrix = pd.DataFrame(index=df.index, columns=features)
            for feature in features:
                if feature in df.columns:
                    sub_matrix[feature] = df[feature]
                else:
                    sub_matrix[feature] = 0
            y_pred = self.core_model.predict(csr_matrix(sub_matrix))
            del df
            del sub_matrix

        pred_df = pd.DataFrame({self.output_annotation: y_pred}, index=test_adata.obs_names)
        pred_df = pred_df.applymap(lambda x: r_label_dict[x])
        
        if self.probability and self.out_path != None:
            pred_prob_df = pd.DataFrame(self.core_model.predict_proba(test_adata[:,features].X), columns=self.core_model.classes_, index=test_adata.obs_names)
            pred_prob_df.columns = pred_prob_df.columns.map(lambda x: r_label_dict[x])
            pred_prob_df.to_csv(self.out_path)
        
        if self.celltype_annotation:
            annotation_df = scatomic_celltype_annotation(self.test_h5ad_path, self.n_thread)
            test_adata.obs[self.output_annotation] = annotation_df.loc[test_adata.obs_names,'scATOMIC_pred']
            
            potential_cancer_prefix = 'Cancer|Normal|Soft Tissue|Brain|Neuroblastoma|Oligodendrocytes|Bile Duct|Bladder|\
                                       Bone|Brain|Breast|Colon|Colorectal|Endometrial|Uterine|Esophageal|Gallbladder|Gastric|\
                                       Glial|Kidney|Liver|Lung|Oligodendrocytes|Ovarian|Pancreatic|Prostate|Skin|Sarcoma|Melanoma|Hepatobiliary'
            potential_tumor_cells = annotation_df[annotation_df['scATOMIC_pred'].str.contains(potential_cancer_prefix)].index.tolist()
            test_adata.obs.loc[potential_tumor_cells,self.output_annotation] = pred_df.loc[potential_tumor_cells,self.output_annotation]
        else:
            test_adata.obs[self.output_annotation] = pd.Categorical(pred_df[self.output_annotation], 
                                                                     categories=['Tumor','Normal'])

        return test_adata

    
