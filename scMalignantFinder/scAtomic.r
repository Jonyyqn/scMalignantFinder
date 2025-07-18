args <- commandArgs(trailingOnly = TRUE)

# .libPaths('/path/to/r/library')

library(scATOMIC)
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)

h5ad.path <- args[1]
ncpu <- as.numeric(args[2])
outdir <- args[3]
python_path <- args[4]
use_python(python_path)

input <- anndata::read_h5ad(h5ad.path)
sparse_matrix <- as(as.matrix(t(as.matrix(input$X))), "sparseMatrix")
cell_predictions <- run_scATOMIC(sparse_matrix, mc.cores = ncpu)
results_pancancer <- create_summary_matrix(prediction_list = cell_predictions, 
                                           use_CNVs = F, 
                                           modify_results = T,
                                           mc.cores = ncpu, 
                                           raw_counts = sparse_matrix, 
                                           min_prop = 0.5)
write.csv(results_pancancer, paste0(outdir,'/celltype_annotation.scAtomic.csv'), quote=FALSE)
