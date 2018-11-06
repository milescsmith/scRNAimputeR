require(reticulate)

dim_from_python <- function(seuratObj, dataframe, reduction.use = "tsne"){
  dim_xfer <- as.matrix(dataframe)
  rownames(dim_xfer) <- rownames(seuratObj@meta.data)
  seuratObj <- SetDimReduction(object = seuratObj,
                               reduction.type = reduction.use,
                               slot = "cell.embeddings",
                               new.data = dim_xfer)
  seuratObj <- SetDimReduction(object = seuratObj,
                               reduction.type = reduction.use,
                               slot = "key",
                               new.data = paste0(reduction.use, "_"))
  return(seuratObj)
}

mctsne <- function(seuratObj,
                          reduction.use = "pca",
                          reduction.save = "tsne",
                          ...){
  cell_embeddings <- GetDimReduction(object = seuratObj,
                                     reduction.type = reduction.use,
                                     slot = "cell.embeddings")
  mtsne_df <- mtsne_for_r(r_data_frame = cell_embeddings, ...)
  seuratObj <- dim_from_python(seuratObj = seuratObj,
                               python_dataframe = mtsne_df,
                               reduction.use = reduction.save)
  return(seuratObj)
}

umap <- function(seuratObj,
                 reduction.use = "pca",
                 reduction.save = "umap",
                 ...){
  cell_embeddings <- GetDimReduction(object = seuratObj,
                                     reduction.type = reduction.use,
                                     slot = "cell.embeddings")
  umap_df <- umap_for_r(r_data_frame = cell_embeddings, ...)
  seuratObj <- dim_from_python(seuratObj = seuratObj,
                               python_dataframe = umap_df,
                               reduction.use = reduction.save)
  return(seuratObj)
}

phate <- function(seuratObj,
                  reduction.use = "pca",
                  reduction.save = "phate",
                  ...){
  cell_embeddings <- GetDimReduction(object = seuratObj,
                                     reduction.type = reduction.use,
                                     slot = "cell.embeddings")
  phate_df <- phate_for_r(r_data_frame = cell_embeddings, ...)
  seuratObj <- dim_from_python(seuratObj = seuratObj,
                               python_dataframe = phate_df,
                               reduction.use = reduction.save)
  return(seuratObj)
}

fiTSNE <- function(seuratObj,
                   reduction.use = "pca",
                   reduction.save = "fitsne",
                   ...){
  cell_embeddings <- GetDimReduction(object = seuratObj,
                                     reduction.type = reduction.use,
                                     slot = "cell.embeddings")
  fitsne_df <- fitsne_for_r(r_data_frame = cell_embeddings, ...)
  seuratObj <- dim_from_python(seuratObj = seuratObj,
                               python_dataframe = fitsne_df,
                               reduction.use = reduction.save)
  return(seuratObj)
}

uncurl <- function(seuratObj,
                   data.slot.use = "data",
                   num_clusters = NULL,
                   ...){
  data_use <- GetAssayData(object = seuratObj,
                           slot = data.slot.use)
  if(!is.double(num_clusters)){
    num_clusters <- length(unique(phenograph(seuratObj)@meta.data$community))
  } else {
    num_clusters <- as.integer(num_clusters)
  }
  data_use <- as.matrix(data_use)
  uncurl_return <- uncurl_for_r(data_use, num_clusters, ...)
  mw_log <- uncurl_return[[1]]
  mw_tsvd <- uncurl_return[[2]]
  seuratObj <- dim_from_python(seuratObj, python_dataframe = mw_log, reduction.use = "uncurl")
  seuratObj <- dim_from_python(seuratObj, python_dataframe = mw_tsvd, reduction.use = "uncurl_tsvd")
  seuratObj <- umap(seuratObj, reduction.use = "uncurl_tsvd", reduction.save = "uncurl_umap")
  seuratObj <- multicoreTSNE(seuratObj, reduction.use = "uncurl_tsvd", reduction.save = "uncurl_tsne")
  return(seuratObj)
}

deep_autoencode_seurat <- function(seuratObj){
  raw_data <- as.matrix(seuratObj@raw.data)
  cell_names <- seuratObj@cell.names
  gene_names <- rownames(raw_data)
  dca_data <- deep_autoencode(raw_data, cell_names = cell_names, gene_names = gene_names)
  dca <- t(dca_data[[1]])
  rownames(dca) <- paste0("DCA_",dca_data[[3]])
  colnames(dca) <- dca_data[[2]]
  seuratObj <- SetAssayData(seuratObj, assay.type = "DCA", slot = "raw.data", new.data = dca)
  seuratObj <- NormalizeData(seuratObj, assay.type = "DCA")
  seuratObj <- ScaleData(seuratObj, assay.type = "DCA", do.par = TRUE, num.cores = parallel::detectCores())
  return(seuratObj)
}

phenograph <- function(seuratObj, reduction.use = "pca", k.param = 30, ...){
  ce <- GetDimReduction(seuratObj, reduction.type = reduction.use, slot = "cell.embeddings")
  communities <- phenograph_for_r(ce, k.param, ...)
  seuratObj@meta.data[,'community'] <- communities
  seuratObj <- SetAllIdent(seuratObj, "community")
  return(seuratObj)
}

python_knn_shunt <- function(seuratObj, reduction.use = "pca",...){
  obj_mat <- GetCellEmbeddings(seuratObj, reduction.type = reduction.use)
  distance.metrics <- get_nn(obj_mat)
  n <- nrow(distance.metrics)
  k.for.nn <- 10
  knn.mat <- matrix(data=0, ncol = k.for.nn, nrow = n)
  knd.mat <- knn.mat
  for (i in 1:n){
    knn.mat[i, ] <- order(distance.metrics[i, ])[1:k.for.nn]
    knd.mat[i, ] <- distance.metrics[i, knn.mat[i, ]]
  }
  nn.ranked <- knn.mat[, 1:10]
  seuratObj@snn <- ComputeSNN(nn_ranked = nn.ranked, prune = 1/15)
  rownames(seuratObj@snn) <- seuratObj@cell.names
  colnames(seuratObj@snn) <- seuratObj@cell.names
  return(seuratObj)
}

# pca <- function(seuratObj,
#                 genes.use = NULL,
#                 num.dims = 50,
#                 reduction.save = "pypca",
#                  ...){
#   if (is.null(genes.use)){
#     expr_data <- as.matrix(FetchData(object = seuratObj, vars.all = rownames(seuratObj@data)))
#   } else {
#     expr_data <- as.matrix(FetchData(object = seuratObj, vars.all = genes.use))
#   }
#
#   pca_results <- fbpca$pca(t(expr_data), k = as.integer(num.dims), n_iter = as.integer(12))
#   rownames(pca_results[[1]]) <- colnames(expr_data)
#   #colnames(pca_results[[1]]) <- paste0(reduction.save, 1:num.dims, sep = "_")
#   colnames(pca_results[[3]]) <- rownames(expr_data)
#   #rownames(pca_results[[3]]) <- paste0(reduction.save, 1:num.dims,sep = "_")
#   seuratObj <- SetDimReduction(object = seuratObj,
#                                reduction.type = reduction.save,
#                                slot = 'cell.embeddings',
#                                new.data = t(pca_results[[3]]))
#   seuratObj <- SetDimReduction(object = seuratObj,
#                                reduction.type = reduction.save,
#                                slot = 'gene.loadings',
#                                new.data = pca_results[[1]])
#   seuratObj <- SetDimReduction(object = seuratObj,
#                                reduction.type = reduction.save,
#                                slot = "key",
#                                new.data = paste0(reduction.save, "_"))
#   return(seuratObj)
# }
