# Consensus Clustering Imputation (CCI)

## Inroduction
**CCI (Consensus Clustering Imputation)** is a novel imputation method inspired by consensus clustering techniques. It aims to address the issue of dropouts in single-cell RNA sequencing (scRNA-seq) and other high-dimensional data. CCI borrows the information from similar cells by ensemble clustering to identify robust neighbors, and employs the cell-to-cell similarities to impute gene expression levels. 


## Installation  
CCI is implemented in **R**. You can install it via the following steps:

```R
# Install necessary dependencies (if applicable)
install.packages(c("matrixStats", "ClusterR", "dplyr")) 

# If hosted on GitHub
devtools::install_github("wanlinjuan/CCI")
```

## Quick start

```R
cc_impute = function(input,                         # expression matrix which needs to be imputed
                     num_sampling=50,               # number of repeated samplings
                     prop_sampling=0.8,             # sampling proportion from genes
                     num_clusters=4,                # number of clusters prespecified for clustering on the subsets
                     resolution=NULL,               # resolution prespecified for clustering on the subsets
                     cutoff=0.2,                    # cutoff used to define neighbors
                     select_genes=NULL,             # genes specified to be sampled from
                     clustering_method="K-means",   # clustering method 
                     impute_what=NULL,              # specify which values will be imputed
                     normalize_method="sct"         # normalization method
                     )
```
The function returns an output imputed matrix. 

## Contact
For questions or support, please contact:
Wanlin Juan - wjuan@mcw.edu
