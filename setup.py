from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

d = {}
with open("scMalignantFinder/_version.py") as f:
    exec(f.read(), d)
    
setup(
    name='scMalignantFinder',
    version=d["__version__"],
    description='A tool for identitying malignant cells from single-cell RNA-seq data',
    long_description=long_description,
    author='Journey Yu',
    author_email='jonyyqn@mail.ustc.edu.cn',
    long_description_content_type="text/markdown",
    url="https://github.com/jonyyqn/scMalignantFinder",
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=[
        'scanpy==1.9.3',
        'anndata==0.9.1',
        'scikit-learn==1.2.2',
        'numpy==1.23.4',
        'joblib==1.1.1',
        'pandas==2.1.4',
        'matplotlib==3.6.3',
        'scipy==1.13.1',
        'numba==0.56.4',
        'h5py==3.7.0'
        
    ],
)
