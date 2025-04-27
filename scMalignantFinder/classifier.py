import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
import warnings

# Suppress convergence warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)


def init_core_model(model_method, n_thread):
    """
    Initialize the core classification model.
    """
    if model_method == "LogisticRegression":
        return LogisticRegression(n_jobs=n_thread)
    elif model_method == "RandomForest":
        return RandomForestClassifier(n_jobs=n_thread)
    elif model_method == "XGBoost":
        from xgboost.sklearn import XGBClassifier
        return XGBClassifier(n_jobs=n_thread)
    else:
        raise NotImplementedError(f"Model '{model_method}' is not supported.")


def preprocess_adata(adata, use_raw=False, norm=True):
    """
    Preprocess the AnnData object by applying `.raw` (if specified) and normalization.
    """
    adata.var_names_make_unique()
    if use_raw:
        if adata.raw is not None:
            adata.X = adata.raw.X.copy()
        else:
            raise ValueError("`use_raw` is True, but `raw` attribute is not available in AnnData.")
    if norm:
        sc.pp.normalize_total(adata, target_sum=1e4)
    return adata


def load_feature(feature_path):
    """
    Load the feature list from a file.
    """
    return list(np.loadtxt(feature_path, dtype=str))





class scMalignantFinder:
    """
    Main class for scMalignantFinder.

    Parameters:
    - test_input (str or AnnData): Test data input. It can be:
    - a file path to a `.h5ad`, `.txt`, `.txt.gz`, `.tsv`, `.tsv.gz`, `.csv`, `.csv.gz` file, or
    - an `AnnData` object.
    For `.h5ad` files, the dataset is directly loaded as an AnnData object.  
    For `.txt` (tab-separated), `.tsv` (tab-separated), or `.csv` (comma-separated) files, with rows treated as features and columns as cells.
    - pretrain_dir (str): Directory containing the pretrained model and ordered feature list. If None, the model will be trained. Default: None.
    - train_h5ad_path (str): Path to the training data. Default: "./combine_training.h5ad".
    - model_method (str): Classification algorithm. Supported: 'LogisticRegression', 'RandomForest', 'XGBoost'.
                          Default: "LogisticRegression".
    - feature_path (str): Path to the feature list. Default: "./combined_tumor_up_down_degs.txt".
    - norm_type (bool): Whether to normalize the test data. Default: True.
    - n_thread (int): Number of threads for parallel processing. Default: 1.
    - use_raw (bool): Use `.raw` attribute of AnnData if available. Default: False.
    - celltype_annotation (bool): use scAtomic for cell type annotation. Default: False.
    - r_script_path (str): Path to the scAtomic r script. Default: None.
    - r_env_path (str): Path to the desired R environment. Default is 'Rscript', which uses the system's default R.

    """

    def __init__(
        self,
        test_input,
        pretrain_dir=None,
        train_h5ad_path="./combine_training.h5ad",
        model_method="LogisticRegression",
        feature_path="./combined_tumor_up_down_degs.txt",
        norm_type=True,
        n_thread=1,
        use_raw=False,
        celltype_annotation=False,
        r_script_path=None,
        r_env_path='Rscript',

    ):
        self.test_input = test_input
        self.pretrain_dir = pretrain_dir
        self.train_h5ad_path = train_h5ad_path
        self.model_method = model_method
        self.feature_path = feature_path
        self.norm_type = norm_type
        self.n_thread = n_thread
        self.use_raw = use_raw

        self.test_adata = None
        self.train_adata = None
        self.features = None
        self.core_model = None
        self.fitted = False  # Indicates whether the model is ready for prediction
        self.pretrain_load = False
        self.celltype_annotation = celltype_annotation
        self.r_script_path = r_script_path
        self.r_env_path = r_env_path
        self.missing_feature = []

    def load(self):
        """
        Load and preprocess all necessary data.
        """
        # Load test data
        self._load_test_data()

        # If no pretrain_dir is provided, load and preprocess training data
        if not self.pretrain_dir:
            self._load_train_data()
            self.features = load_feature(self.feature_path)
        else:
            self._load_pretrained_model()

    def _load_test_data(self):
        """
        Load and preprocess test data.
        """
        if isinstance(self.test_input, str):
            if self.test_input.endswith('.h5ad'):
                self.test_adata = sc.read_h5ad(self.test_input)
            elif self.test_input.endswith(('.txt', '.txt.gz', '.tsv', '.tsv.gz')):
                self.test_adata = sc.AnnData(pd.read_csv(self.test_input, sep='\t', index_col=0, 
                                                         low_memory=False, engine='c').T)
                self.test_adata.X = sparse.csr_matrix(self.test_adata.X)
            elif self.test_input.endswith(('.csv', '.csv.gz')):
                self.test_adata = sc.AnnData(pd.read_csv(self.test_input, sep=',', index_col=0, 
                                                         low_memory=False, engine='c').T)
                self.test_adata.X = sparse.csr_matrix(self.test_adata.X)
            else:
                raise ValueError("`test_input` must be a .h5ad, .txt, .tsv, .csv or .gz file path.")
        elif isinstance(self.test_input, sc.AnnData):
            self.test_adata = self.test_input
        else:
            raise ValueError("`test_input` must be a file path (.h5ad, .txt, .tsv, .csv) or an AnnData object.")

        # Apply preprocessing to test data
        self.test_adata = preprocess_adata(self.test_adata, use_raw=self.use_raw, norm=self.norm_type)

    def _load_train_data(self):
        """
        Load and preprocess training data.
        """
        self.train_adata = sc.read_h5ad(self.train_h5ad_path)
        self.train_adata = preprocess_adata(self.train_adata, norm=True)

    def _load_pretrained_model(self):
        """
        Load pretrained model and ensure test data compatibility.
        """
        import joblib

        # Load pretrained model
        model_path = os.path.join(self.pretrain_dir, "model.joblib")
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Pretrained model not found at {model_path}.")
        self.core_model = joblib.load(model_path)

        # Ensure features compatibility with test data
        ordered_features_path = os.path.join(self.pretrain_dir, "ordered_feature.tsv")
        self.features = load_feature(ordered_features_path)

        # Check for missing features (inform user but don't modify test_adata permanently)
        missing_features = set(self.features) - set(self.test_adata.var_names)
        num_missing = len(missing_features)
        total_features = len(self.features)
        missing_percentage = (num_missing / total_features) * 100

        print(f"Model features: {total_features}")
        print(f"Missing features: {num_missing} ({missing_percentage:.2f}%)")

        if missing_percentage > 20:
            print("Warning: Missing feature percentage is high (>20%). Consider retraining the model.")

        self.fitted = True  # Mark model as ready for prediction
        self.pretrain_load = True
        self.missing_feature = list(missing_features)

    def predict(self):
        """
        Predict malignancy for the test data.

        If no pretrained model is provided, train the model before prediction.
        """
        if not self.features:
            raise RuntimeError("Features are not loaded. Please run `load()` first!")

        # Train the model if no pretrained model is provided
        if not self.pretrain_dir:
            self._train_model()

        # Perform predictions only if the model is fitted
        if self.fitted:
            return self._make_predictions()
        else:
            raise RuntimeError("Model is not fitted. Please check the training or loading process.")

    def _train_model(self):
        """
        Train the model using the training data.
        """
        filter_features = list(set(self.features) & set(self.train_adata.var_names) & set(self.test_adata.var_names))
        if not filter_features:
            raise RuntimeError("No overlapping features found between test and training sets.")

        num_cells = self.train_adata.shape[0]
        num_features = len(filter_features)
        print(f"Training model with {num_cells} cells and {num_features} features using {self.model_method}.")

        self.core_model = init_core_model(self.model_method, self.n_thread)
        self.core_model.fit(
            self.train_adata[:, filter_features].X,
            self.train_adata.obs["Raw_annotation"].map({"Normal": 0, "Malignant": 1, 'Tumor': 1}).tolist(),
        )
        self.fitted = True  # Mark model as ready for prediction
        self.features = filter_features

    def _make_predictions(self):
        """
        Perform predictions and add results to the test AnnData object.
        """

        if len(self.missing_feature)==0:
            y_pred = self.core_model.predict(np.array(self.test_adata[:,self.features].X.todense()))
            y_prob = self.core_model.predict_proba(np.array(self.test_adata[:,self.features].X.todense()))[:, 1]
        else:
            missing_df = pd.DataFrame(0, index=self.test_adata.obs_names, columns=self.missing_feature)
            full_matrix = pd.concat([self.test_adata.to_df(), missing_df], axis=1).loc[:,self.features].values
            y_pred = self.core_model.predict(full_matrix)
            y_prob = self.core_model.predict_proba(full_matrix)[:, 1]

            del full_matrix, missing_df

        # Add predictions and probabilities to AnnData (without altering its original structure)
        self.test_adata.obs["scMalignantFinder_prediction"] = pd.Categorical(
            ["Malignant" if pred else "Normal" for pred in y_pred]
        )
        self.test_adata.obs["malignancy_probability"] = y_prob

        # Perform cell type annotation
        if self.celltype_annotation:
            tmpdir = os.getcwd()
            r_script = os.path.join(self.r_script_path, 'scAtomic.r')
            python_path = sys.executable
            
            if isinstance(self.test_input, str):
                adata = sc.read_h5ad(self.test_input)
            elif isinstance(self.test_input, sc.AnnData):
                adata = self.test_input
            if self.use_raw:
                if adata.raw is not None:
                    adata.X = adata.raw.X.copy()
            tmp_adata_path = os.path.join(tmpdir,'input.h5ad')
            adata.write(tmp_adata_path)

            import subprocess
            result = subprocess.run(
                [self.r_env_path, r_script, tmp_adata_path, str(self.n_thread), tmpdir, python_path],
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                raise RuntimeError(f'Error running R script: {result.stderr}')


            annot_df = pd.read_csv(f'{tmpdir}/celltype_annotation.scAtomic.csv', index_col=0)
            self.test_adata.obs['scATOMIC_pred'] = annot_df.loc[self.test_adata.obs_names,'scATOMIC_pred']
            os.system('rm {0} {1}'.format(tmp_adata_path, f'{tmpdir}/celltype_annotation.scAtomic.csv'))

            # potential_prefix = 'Cancer|Normal|Soft Tissue|Brain|Neuroblastoma|Oligodendrocytes|Bile Duct|\
            #                     Bladder|Bone|Brain|Breast|Colon|Colorectal|Endometrial|Uterine|Esophageal|\
            #                     Gallbladder|Gastric|Glial|Kidney|Liver|Lung|Oligodendrocytes|Ovarian|\
            #                     Pancreatic|Prostate|Skin|Sarcoma|Melanoma|Hepatobiliary'

            # potential_cells = annot_df[annot_df['scATOMIC_pred'].str.contains(potential_prefix)].index.tolist()


        return self.test_adata


