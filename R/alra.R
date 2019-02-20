#' @title ALRA
#' @description Use Adaptively-thresholded Low Rank Approximation to perform imputation
#'
#' @param object Data object to impute
#' @param assay_use The assay from which to retrieve data. Default: "RNA"
#' @param slot_use The slot within the assay to impute. Default: "data"
#' @param ... additional parameters to pass to ALRA::alra()
#'
#' @importFrom glue glue
#' @importFrom Matrix Matrix
#'
#' @export
#' @return
ALRA <- function(object, 
                 assay_use = "RNA", 
                 slot_use = "data", 
                 ...) {
  datExprs <- GatherData(object = object,
                         assay = assay_use, 
                         slot = slot_use) %>% t()
  alraExprs <- alra(datExprs, 
                    ...) %>% 
    '[['(3) %>%
    t()
  colnames(alraExprs) <- colnames(datExprs)
  rownames(alraExprs) <- rownames(datExprs)
  alraExprs <- Matrix(alraExprs, sparse = T)
  object <- PlaceData(object = object, 
                      assay_store = "ALRA", 
                      imputed_data = alraExprs)
  return(object)
}
