# Consensus Clustering Imputation (CCI)

## Overview  
**CCI (Consensus Clustering Imputation)** is a novel imputation method inspired by consensus clustering techniques. It aims to address the issue of dropouts in single-cell RNA sequencing (scRNA-seq) and other high-dimensional data. CCI enhances imputation accuracy by leveraging ensemble clustering methods to identify robust patterns in the data.

## Key Features  
- **Robust Imputation:** Uses consensus clustering to improve the accuracy of missing value imputation.
- **Scalable:** Can handle large `n × p` matrices where both `n` (samples) and `p` (features) are large.
- **High Performance:** Suitable for single-cell RNA-seq data, reducing the impact of dropout events.
- **Versatile:** Compatible with various types of high-dimensional datasets, not limited to scRNA-seq.

## Installation  
CCI is implemented in **R**. You can install it via the following steps:

```R
# Install necessary dependencies (if applicable)
install.packages(c("matrixStats", "ClusterR", "dplyr")) 

# If hosted on GitHub
devtools::install_github("your_github_username/CCI")

## Method
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

