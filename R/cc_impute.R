#' Consensus Clustering-Based Imputation for scRNA-seq
#'
#' This function performs imputation on scRNA-seq data using a consensus clustering approach,
#' aggregating multiple clustering results over sampled gene subsets.
#'
#' @param input Gene-by-cell expression matrix to impute.
#' @param num_sampling Number of repeated samplings (default = 50).
#' @param prop_sampling Proportion of genes sampled in each iteration (default = 0.8).
#' @param num_clusters Number of clusters for K-means (default = 4).
#' @param resolution Optional Seurat clustering resolution (not used if `num_clusters` is set).
#' @param cutoff Threshold to define neighboring cells (default = 0.2).
#' @param select_genes Vector of gene names to sample from (default = NULL means all genes).
#' @param clustering_method Clustering method to use: "K-means" or "SNN" (default = "K-means").
#' @param impute_what Logical matrix indicating which values to impute (default = zeros).
#' @param normalize_method Normalization method: "sct" or "log" (default = "sct").
#'
#' @return A gene-by-cell matrix with imputed values.
#' @export
#' @importFrom Seurat CreateSeuratObject SCTransform NormalizeData FindVariableFeatures ScaleData RunPCA RunUMAP GetAssayData
#' @importFrom scuttle logNormCounts
#' @importFrom scater runPCA
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("splatter", quietly = TRUE)) {
#'   library(Seurat)
#'   library(splatter)
#'   set.seed(1)
#'   params <- newSplatParams(batchCells = 200, nGenes = 100)
#'   sim <- splatSimulate(params, verbose = FALSE)
#'   sim <- logNormCounts(sim)
#'   seurat_obj <- as.Seurat(sim)
#'   data <- GetAssayData(seurat_obj, slot = "data")
#'   hv_genes <- rownames(data)[1:50]
#'   result <- cc_impute(data, num_sampling = 5, select_genes = hv_genes, num_clusters = 2)
#' }
#' \dontrun{
#' imputed_data <- cc_impute(input = expression_matrix)
#' }
cc_impute = function(input,                         # expression matrix which needs to be imputed. should involve the rownames and colnames.
                     num_sampling=50,               # number of repeated samplings
                     prop_sampling=0.8,             # sampling proportion from genes
                     num_clusters=4,                # number of clusters prespecified for clustering on the subsets. num_clusters or resolution must be specified.
                     resolution=NULL,               # resolution prespecified for clustering on the subsets. num_clusters or resolution must be specified
                     cutoff=0.2,                    # cutoff used to define neighbors
                     select_genes=NULL,             # genes specified to be sampled from. If NULL, sample from all genes.
                     clustering_method="K-means",   # clustering_method including SNN or K-means
                     impute_what=NULL,              # specify which values will be imputed. default is all the zero counts. It should be a logical matrix with the same size
                     normalize_method="log"         # normalization method including SCTranform "sct" or log normalization "log"
){



  # create consensus matrix (cell by cell)
  n_genes = dim(input)[1]
  n_cells = dim(input)[2]
  cm <- matrix(0,n_cells,n_cells)

  # genes to be sampled from
  if(is.null(select_genes)){
    # genes2sample = 1:ngenes
    genes2sample = rownames(input)
  }else{
    # select_genes_i = match(select_genes, rownames(input))
    # genes2sample = na.omit(select_genes_i)
    genes2sample = intersect(rownames(input), select_genes)
  }

  for (i in 1:num_sampling){
    set.seed(NULL)
    genes_sample <- sample(genes2sample, prop_sampling*length(genes2sample), replace=FALSE)
    subset = input[genes_sample,]
    subset.seurat = CreateSeuratObject(counts=subset)
    if(normalize_method=="sct"){
      subset.seurat = SCTransform(subset.seurat,assay="originalexp",verbose=FALSE)
    }else if(normalize_method=="log"){
      subset.seurat = NormalizeData(subset.seurat)
      subset.seurat = FindVariableFeatures(subset.seurat)
      subset.seurat = ScaleData(subset.seurat)
    }


    subset.seurat = RunPCA(subset.seurat, features = VariableFeatures(object = subset.seurat), verbose=FALSE,
                           reduction.name = "pca",reduction.key = "pca_")
    subset.seurat = RunUMAP(object=subset.seurat, dims=1:10, verbose=FALSE)


    newdata = t(GetAssayData(subset.seurat,slot="scale.data"))
    # pc=prcomp(newdata,center = TRUE,scale. = FALSE)
    # km_data=kmeans(pc$x[,1:2],centers=2, iter.max = 10, nstart = 1)
    if(clustering_method=="K-means"){
      km_data=kmeans(subset.seurat@reductions$umap@cell.embeddings,centers=num_clusters, iter.max = 10, nstart = 1)
      membership <- as.numeric(km_data$cluster)
    }


    # membership matrix
    membership_matrix <-dist(membership,diag=TRUE,method="manhattan",upper=TRUE)
    membership_matrix <-as.matrix(membership_matrix)
    membership_matrix <- (membership_matrix==0)*1

    cm <- cm+membership_matrix
  }
  cm <- cm/num_sampling
  cm_neighbor = cm*(cm>=cutoff)

  impute_weighted_ave <- function(i) {
    (cm_neighbor %*% data_input[i,]) / (cm_neighbor %*% (data_input[i,] != 0) + 0.000001)
  }

  if(is.null(impute_what)){
    impute_index = (data==0)
  }else{
    impute_index = (impute_what)
  }

  if(is.null(select_genes)){
    # Case 1: All genes to be considered
    data_input = data
    data_output = data_input
    final_data =data_input

    data_output <- t(sapply(1:dim(data_input)[1], impute_weighted_ave))
    replace_indices <- which(impute_index, arr.ind = TRUE)
    final_data[replace_indices] = data_output[replace_indices]

  }else{
    # Case 2: Subset of genes specified by genes2sample
    data_input = data[genes2sample,]
    data_output = data_input
    final_data_select_genes = data_input
    final_data = data
    impute_index_genes2sample = impute_index[intersect(genes2sample,rownames(impute_index)),]

    data_output <- t(sapply(1:dim(data_input)[1], impute_weighted_ave))

    replace_indices <- which(impute_index_genes2sample, arr.ind = TRUE)
    final_data_select_genes[replace_indices] = data_output[replace_indices]
    final_data[genes2sample,] = final_data_select_genes


  }

  finalcount=as.matrix(final_data)
  newcount = finalcount

  return(finaldata = finalcount)
}


