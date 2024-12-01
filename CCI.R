
cc_impute <- function(input,                         # Expression matrix to be imputed (rows = genes, columns = cells).
                      num_sampling <- 50,            # Number of repeated samplings.
                      prop_sampling <- 0.8,          # Proportion of genes to sample in each iteration.
                      num_clusters <- 4,             # Prespecified number of clusters for clustering subsets.
                      resolution <- NULL,            # Resolution parameter for clustering (alternative to num_clusters).
                      cutoff <- 0.2,                 # Threshold for defining neighbors in the consensus matrix.
                      select_genes <- NULL,          # List of genes to sample from. Defaults to all genes if NULL.
                      clustering_method <- "K-means",# Clustering method, either "K-means" or "SNN".
                      impute_what <- NULL,           # Logical matrix indicating which values to impute (default: all zeros).
                      normalize_method <- "sct"      # Normalization method: "sct" for SCTransform or "log" for log normalization.
){
  # Initialize dimensions and consensus matrix
  n_genes <- dim(input)[1]
  n_cells <- dim(input)[2]
  cm <- matrix(0, n_cells, n_cells)
  
  # Determine genes to sample
  if (is.null(select_genes)) {
    genes2sample <- rownames(input)  # Use all genes if none are specified
  } else {
    genes2sample <- intersect(rownames(input), select_genes)  # Use the intersection of specified genes and input rows
  }
  
  # Perform repeated sampling
  for (i in 1:num_sampling) {
    set.seed(NULL)  # Use random seed for sampling
    genes_sample <- sample(genes2sample, prop_sampling * length(genes2sample), replace = FALSE)
    subset <- input[genes_sample, ]
    subset.seurat <- CreateSeuratObject(counts = subset)
    
    # Apply normalization
    if (normalize_method == "sct") {
      subset.seurat <- SCTransform(subset.seurat, verbose = FALSE)
    } else if (normalize_method == "log") {
      subset.seurat <- NormalizeData(subset.seurat)
      subset.seurat <- FindVariableFeatures(subset.seurat)
      subset.seurat <- ScaleData(subset.seurat)
    }
    
    # Perform PCA and UMAP
    subset.seurat <- RunPCA(subset.seurat, features = VariableFeatures(object = subset.seurat), verbose = FALSE,
                            reduction.name = "pca", reduction.key = "pca_")
    subset.seurat <- RunUMAP(object = subset.seurat, dims = 1:10, verbose = FALSE)
    
    # Clustering using specified method
    if (clustering_method == "K-means") {
      km_data <- kmeans(subset.seurat@reductions$umap@cell.embeddings, centers = num_clusters, iter.max = 10, nstart = 1)
      membership <- as.numeric(km_data$cluster)
    }
    
    # Create membership matrix
    membership_matrix <- dist(membership, diag = TRUE, method = "manhattan", upper = TRUE)
    membership_matrix <- as.matrix(membership_matrix)
    membership_matrix <- (membership_matrix == 0) * 1
    
    # Update consensus matrix
    cm <- cm + membership_matrix
  }
  
  # Normalize consensus matrix and apply cutoff
  cm <- cm / num_sampling
  cm_neighbor <- cm * (cm >= cutoff)
  
  # Function to compute weighted average for imputation
  impute_weighted_ave <- function(i) {
    (cm_neighbor %*% data_input[i, ]) / (cm_neighbor %*% (data_input[i, ] != 0) + 1e-6)
  }
  
  # Determine indices to impute
  if (is.null(impute_what)) {
    impute_index <- (input == 0)
  } else {
    impute_index <- impute_what
  }
  
  # Impute data
  if (is.null(select_genes)) {
    data_input <- input
    data_output <- data_input
    final_data <- data_input
    
    data_output <- t(sapply(1:dim(data_input)[1], impute_weighted_ave))
    replace_indices <- which(impute_index, arr.ind = TRUE)
    final_data[replace_indices] <- data_output[replace_indices]
  } else {
    data_input <- input[genes2sample, ]
    data_output <- data_input
    final_data_select_genes <- data_input
    final_data <- input
    impute_index_genes2sample <- impute_index[intersect(genes2sample, rownames(impute_index)), ]
    
    data_output <- t(sapply(1:dim(data_input)[1], impute_weighted_ave))
    
    replace_indices <- which(impute_index_genes2sample, arr.ind = TRUE)
    final_data_select_genes[replace_indices] <- data_output[replace_indices]
    final_data[genes2sample, ] <- final_data_select_genes
  }
  
  # Finalize and return imputed data
  finaldata <- as.matrix(final_data)
  rownames(finaldata) <- rownames(input)
  colnames(finaldata) <- colnames(input)
  
  return(finaldata)
}
