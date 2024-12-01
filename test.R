library(Seurat)
library(splatter)
library(scater)
library(ggpubr)


set.seed(1)
ngroups=2;batchCells = 2000;nGenes = 500; de.prob = 0.1; de.facLoc= 0.3; de.downProb=0.5; seed=1
dropout.type="experiment"; dropout.mid=2; dropout.shape=-1;de.facScale=0.1


# sim1: simulated data without dropouts
params.groups <- newSplatParams(batchCells = batchCells, nGenes = nGenes, seed=seed)
sim1 <- splatSimulateGroups(params.groups, group.prob = rep(1,ngroups)/ngroups,
                            de.prob = de.prob, de.facLoc= de.facLoc, de.facScale=de.facScale,
                            de.downProb=de.downProb, verbose = FALSE, seed=seed)
sim1 <- logNormCounts(sim1)
sim1 <- runPCA(sim1)

params = newSplatParams(batchCells = batchCells, nGenes = nGenes, seed=seed)
params = setParams(params, update=list(dropout.type="experiment",dropout.mid=-99999, seed=seed))
sim1 = splatter:::splatSimDropout(sim1, params)

# convert singlecellexperiment object to seurat object
sim1.s = as.Seurat(sim1)
sim1.s$seurat_clusters = as.numeric(as.factor(sim1$Group))
Idents(sim1.s) = as.numeric(as.factor(sim1$Group))

sim1.s = NormalizeData(sim1.s)
sim1.s = FindVariableFeatures(sim1.s)
sim1.s = ScaleData(sim1.s)
sim1.s = RunPCA(sim1.s, features=VariableFeatures(sim1.s), verbose = FALSE)
sim1.s = RunUMAP(object=sim1.s, dims=1:10, verbose=FALSE)



# sim2: add dropouts to sim1
params = newSplatParams(seed=seed)
params = setParams(params, nGenes=nGenes, update=list(dropout.type="experiment", dropout.mid=dropout.mid, dropout.shape=dropout.shape,seed=seed))
sim2 = splatter:::splatSimDropout(sim1, params)
sim2 <- logNormCounts(sim2)
sim2 <- runPCA(sim2)

# convert singlecellexperiment object to seurat object
sim2.s = as.Seurat(sim2)
sim2.s$seurat_clusters = as.numeric(as.factor(sim2$Group))
Idents(sim2.s) = as.numeric(as.factor(sim2$Group))

sim2.s = NormalizeData(sim2.s)
sim2.s = FindVariableFeatures(sim2.s)  
sim2.s = ScaleData(sim2.s)
sim2.s = RunPCA(sim2.s, features=VariableFeatures(sim2.s), verbose = FALSE)
sim2.s = RunUMAP(object=sim2.s, dims=1:10, verbose=FALSE)



# Define top100 highly variable genes
hv_genes = VariableFeatures(sim2.s)[1:100]
data = GetAssayData(sim2.s, slot="data")

# Impute top100 genes using CCI
newcount = cc_impute(data, num_sampling=10, prop_sampling=0.8, 
                     num_clusters=4, resolution=NULL,
                     cutoff=0.2, select_genes=hv_genes, clustering_method="K-means",
                     normalize_method="log")


cci_imputed = sim2.s
cci_imputed[["imputed"]] <- CreateAssayObject(data=newcount)
cci_imputed[["imputed"]] <- CreateAssayObject(counts=newcount)
DefaultAssay(object = cci_imputed) <- "imputed"
# cci_imputed = CreateSeuratObject(assay="imputed", counts=newcount,data=newcount)
cci_imputed = FindVariableFeatures(cci_imputed, verbose=FALSE)
cci_imputed = ScaleData(cci_imputed, assay="imputed")
cci_imputed = RunPCA(cci_imputed, features=VariableFeatures(cci_imputed),  verbose = FALSE,
                reduction.name = "pca",reduction.key = "pca_")
cci_imputed = RunUMAP(object=cci_imputed, dims=1:10, verbose=FALSE)


# PCA plots
p1 = DimPlot(sim1.s, reduction="pca",group.by="Group")+ggtitle("Without dropout")
p2 = DimPlot(sim2.s, reduction="pca",group.by="Group")+ggtitle("With dropout")
p3 = DimPlot(cci_imputed, reduction="pca",group.by="Group")+ggtitle("Imputed by CCI")
summary = ggarrange(p1,p2,p3,nrow=3,ncol=1)
print(summary)


# UMAP plots
p1 = DimPlot(sim1.s, reduction="umap",group.by="Group")+ggtitle("Without dropout")
p2 = DimPlot(sim2.s, reduction="umap",group.by="Group")+ggtitle("With dropout")
p3 = DimPlot(cci_imputed, reduction="umap",group.by="Group")+ggtitle("Imputed by CCI")
summary = ggarrange(p1,p2,p3,nrow=3,ncol=1)
print(summary)

