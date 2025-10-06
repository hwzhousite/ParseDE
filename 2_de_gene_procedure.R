source("utility.r")

K = 4
parse.fit = readRDS(paste0(c("islet/parse_islet_K",K,".rds"),collapse = ""))
parse.fit.null = readRDS(paste0(c("islet/parse_islet_null_K",K,".rds"),collapse = ""))

# Contrast Score
origin_score = score_func(mu=parse.fit$mu.hat.best, sigma = (parse.fit$sigma.hat.best), gene = parse.fit$feature.hat.best)
null_score =  score_func(mu=parse.fit.null$mu.hat.best, sigma = diag(parse.fit.null$sigma.hat.best), gene = parse.fit.null$feature.hat.best)
contrast_score =  (origin_score - null_score)


# DE gene
thres = knockoff::knockoff.threshold(contrast_score[parse.fit$feature.hat.best] ,fdr = 0.05, offset = 1)
de_gene = which( contrast_score > thres)



# Seurat Performance
set.seed(123)
# Step 1: Create Seurat object
seurat_obj <- CreateSeuratObject(counts = t(data_norm[, de_gene]))
seurat_obj <- UpdateSeuratObject(seurat_obj) ## Avoid errors
# 1. Normalize and find variable features
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = ncol(data_norm))

# 2. Scale and run PCA
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# 3. Find neighbors & clusters
# Tuning the resolution to get the cluster
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)

# 4. Run t-SNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)

# 5. Plot with DimPlot
#DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters")

saveRDS(list(seurat_obj = seurat_obj, de_gene = de_gene), "islet/parseDE_K4.rds")


