setGeneric(
  name = "uncurl",
  def = function(x, ...) {
    standardGeneric("uncurl")
  }
)

#' Title
#'
#' @param seuratObj
#' @param data.slot.use
#' @param num_clusters
#' @param ...
#'
#' @importFrom gsubfn list
#'
#' @return
#' @export
#'
#' @examples
setMethod("uncurl",
          signature(x='seurat'),
          function(seuratObj,
                   data.slot.use = "data",
                   num_clusters = NULL,
                   ...){
  uncurl_module <- import(module = 'uncurl', delay_load = TRUE)
  sklearn.decomposition <- import(module = 'sklearn.decomposition', delay_load = TRUE)

  data_use <- GetAssayData(object = seuratObj,
                           slot = data.slot.use)

  if(!is.double(num_clusters)){
    num_clusters <- length(unique(phenograph(seuratObj)@meta.data$community))
  } else {
    num_clusters <- as.integer(num_clusters)
  }
  data_use <- as.matrix(data_use)

  genes <- uncurl_module$max_variance_genes(data_use,
                                            nbins = as.integer(5),
                                            frac = as.numeric(0.2)) %>%
    unlist()
  # translate...
  data_subset = data_use[genes, ]
  list[M, W, ll] <- uncurl_module$run_state_estimation(data_subset,
                                                 num_clusters,
                                                 dist='Poiss',
                                                 disp=FALSE,
                                                 reps=as.integer(6),
                                                 threads = as.integer(4))
  mw <- M %*% W
  mw_log <- log(1 + mw) %>% t()

  mw_tsvd <- irlba(mw_log, nv = 50)
  tsvd = sklearn.decomposition$TruncatedSVD(as.integer(50))
  mw_tsvd = tsvd$fit_transform(mw_log)
  uncurl_return <- uncurl_for_r(data_use, num_clusters, ...)
  mw_log <- uncurl_return[[1]]
  mw_tsvd <- uncurl_return[[2]]
  seuratObj <- dim_from_python(seuratObj, dataframe = mw_log, reduction.use = "uncurl")
  seuratObj <- dim_from_python(seuratObj, dataframe = mw_tsvd, reduction.use = "uncurl_tsvd")
  return(seuratObj)
})
