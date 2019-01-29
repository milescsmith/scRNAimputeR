#' Magic
#' 
#' Wrapper around the Markov Affinity-based Graph Imputation of Cells (MAGIC) 
#' imputation algorithm from the KrishnaswamyLab
#'
#' @param genes Default: NULL
#' @param k Default: 10
#' @param alpha Default: 15
#' @param t Default: "auto"
#' @param npca Default: 10
#' @param init Default: NULL
#' @param t.max Default: 20
#' @param knn.dist.method Default: "euclidean"
#' @param verbose Default: 1
#' @param n.jobs Default: 1
#' @param seed Default: NULL
#'
#' @return
#' @export
#'
#' @examples
Magic <- function(object, ...) {
  UseMethod("Magic")
}

#' @rdname Magic
#' @method Magic matrix
#' 
#' @param exprDat
#' 
#' @importFrom reticulate py_module_available import
#' @importFrom glue glue
#' 
#' @return matrix
#' @export
Magic.default <- function(exprDat,
                          genes = NULL, 
                          k = 10, 
                          alpha = 15, 
                          t = "auto",
                          npca = 100, 
                          init = NULL, 
                          t.max = 20,
                          knn.dist.method = "euclidean", 
                          verbose = 1, 
                          n.jobs = 1,
                          seed = NULL){
  
  required_modules <- c("magic")
  for (i in required_modules){
    if(!py_module_available(i)){
      stop(glue("The {i} module is unavailable.
         Please activate the appropriate environment or
         install the module."))
    }
  }
  
  magic.module <- import(module = 'magic', delay_load = TRUE)
  
  exprDat %<>% as.data.frame()
  
  if (is.null(genes)) {
    genes <- colnames(exprDat)
  }
  
  magic_operator <- magic.module$MAGIC(k = as.integer(k), 
                                       a = as.integer(alpha), 
                                       t = t, 
                                       n_pca = as.integer(npca), 
                                       knn_dist = knn.dist.method, 
                                       n_jobs = as.integer(n.jobs), 
                                       random_state = seed, 
                                       verbose = verbose)
  
  magic_result <- magic_operator$fit_transform(X = exprDat, 
                                               genes = genes, 
                                               t_max = t.max)
  return(magic_result)
}

#' @rdname Magic
#' @method Magic Seurat
#' 
#' @param object Seurat object with data to impute by MAGIC.
#' @param assay_use Object assay containing the data to impute. Default: "RNA"
#' @param slot_use Assay slot containing data to impute. Default: "counts"
#' @param assay_store Name to use when storing imputed data. Default: "magic"
#' 
#' @return Seurat object
#' @export
Magic.Seurat <- function(object,
                         assay_use = "RNA",
                         slot_use = "counts",
                         assay_store = "magic",
                         genes = NULL, 
                         k = 10, 
                         alpha = 15, 
                         t = "auto",
                         npca = 100, 
                         init = NULL, 
                         t.max = 20,
                         knn.dist.method = "euclidean", 
                         verbose = 1, 
                         n.jobs = 1,
                         seed = NULL){

  exprDat <- GatherData(object, assay_use, slot_use) %>% t()

  magic_result <- Magic(exprDat, genes, k, alpha, t, npca, init, 
                        t.max,knn.dist.method, verbose, n.jobs,seed)
  
  object <- PlaceData(object, assay_store = assay_store, imputed_data = t(magic_result))
}

#' @rdname Magic
#' @method Magic seurat
#' 
#' @param object Seurat object with data to impute by MAGIC.
#' @param slot_use Data slot containing data to impute. Default: "counts"
#' @param assay_store Name to use when storing imputed data. Default: "magic"
#' 
#' @return seurat object
#' @export
Magic.seurat <- function(object,
                         slot_use = "data",
                         assay_store = "magic",
                         genes = NULL, 
                         k = 10, 
                         alpha = 15, 
                         t = "auto",
                         npca = 100, 
                         init = NULL, 
                         t.max = 20,
                         knn.dist.method = "euclidean", 
                         verbose = 1, 
                         n.jobs = 1,
                         seed = NULL){
  
  exprDat <- GatherData(object, slot_use) %>% t()

  magic_result <- Magic(exprDat, genes, k, alpha, t, npca, init, 
                        t.max,knn.dist.method, verbose, n.jobs,seed)
  
  object <- PlaceData(object, assay_store = assay_store, imputed_data = t(magic_result))
}

#' @rdname Magic
#' @method Magic SingleCellExperiment
#' 
#' @param object SingleCellExperiment object with data to impute by MAGIC.
#' @param assay_use Object assay containing the data to impute. Default: "logcounts"
#' @param assay_store Name to use when storing imputed data. Default: "magic"
#' 
#' @return SingleCellExperiment object
#' @export
Magic.SingleCellExperiment <- function(object,
                                       assay_use = "logcounts",
                                       assay_store = "magic",
                                       genes = NULL, 
                                       k = 10, 
                                       alpha = 15, 
                                       t = "auto",
                                       npca = 100, 
                                       init = NULL, 
                                       t.max = 20,
                                       knn.dist.method = "euclidean", 
                                       verbose = 1, 
                                       n.jobs = 1,
                                       seed = NULL){
  
  exprDat <- GatherData(object, assay_use) %>% t()
  
  magic_result <- Magic(exprDat, genes, k, alpha, t, npca, init, 
                        t.max,knn.dist.method, verbose, n.jobs,seed)
  
  object <- PlaceData(object, assay_store = assay_store, imputed_data = t(magic_result))
}