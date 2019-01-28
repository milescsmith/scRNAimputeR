PlaceData <- function(object, ...) {
    UseMethod("PlaceData")
  }

#' PlaceData
#'
#' @param object Seurat object to pull data from
#' @param assay Name of assay slot to place imputated data
#' @param ...
#'
#' @importFrom Seurat CreateAssayObject NormalizeData ScaleData SetAssayData
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom Matrix Matrix
#'
#' @rdname PlaceData
#' @method PlaceData Seurat
#' @return
#' @export
PlaceData.Seurat <- function(object,
                   assay_store = "new_assay",
                   imputed_data,
                   normalize = "FALSE",
                   scale = "FALSE",
                   ...){
            
  imputed_data %<>% as.matrix()
  object[[assay_store]] <- CreateAssayObject(data = imputed_data)
  if (normalize)
  {
    object <- NormalizeData(object, assay = assay_store)
  }
  if (scale){
    object <- ScaleData(object, assay = assay_store, ...)
  }
  return(object)
}

#' @rdname PlaceData
#' @method PlaceData seurat
#' @return
#' @export
PlaceData.seurat <- function(object,
                             assay_store = "new_assay",
                             slot_use = "data",
                             imputed_data,
                             normalize = "FALSE",
                             scale = "FALSE",
                             ...){
            
  data <- Matrix(data = imputed_data, sparse = TRUE)
  object <- SetAssayData(object = seuratObj, 
                         assay = assay_store,
                         slot = slot_use,
                         new.data = imputed_data)
  if (normalize){
    object <- NormalizeData(object, assay.type = assay_store)
  }
  if (scale){
    object <- ScaleData(object, ...)
  }
  return(object)
}

#' @rdname PlaceData
#' @method PlaceData SingleCellExperiment
#' @return
#' @export
PlaceData.SingleCellExperiment <- function(object,
                                           assay_store = "new_assay",
                                           imputed_data,
                                           normalize = "FALSE",
                                           ...){
  
  data <- Matrix(data = imputed_data, sparse = TRUE)
  assay(object = object, i = assay_store) <- data
  if (normalize){
    object <- normalize(object, exprs_values = assay_store)
  }
  return(object)
}

GatherData <- function(object, ...) {
    UseMethod("GatherData")
  }

#' GatherData
#'
#' @param object Seurat object to pull data from
#' @param assay (Default: "RNA")
#' @param slot_use (Default: "data)
#' @param ...
#'
#' @importFrom Seurat GetAssayData
#' @importFrom magrittr %>%
#'
#' @rdname GatherData
#' @method GatherData Seurat
#' @return matrix
#' @export
GatherData.Seurat <- function(object,
                              assay = "RNA",
                              slot_use = "counts",
                              ...){
            
  obj_data <- GetAssayData(object = object, 
                           assay = assay,
                           slot = slot_use) %>%
    as.matrix()
  return(obj_data)
}

#' @rdname GatherData
#' @method GatherData seurat
#' @return matrix
#' @export
GatherData.seurat <- function(object,
                              slot_use = "data",
                              ...){
  obj_data <- GetAssayData(object = object, 
                           assay = assay,
                           slot = slot_use) %>%
    as.matrix()
  return(obj_data)
}

#' @rdname GatherData
#' @method GatherData SingleCellExperiment
#' @return matrix
#' @export
GatherData.SingleCellExperiment <- function(object,
                                            assay_use = "logcounts",
                                            ...){
  if (!assay_use %in% names(assays(object))){
    stop(glue("Sorry, but {assay_use} is not present in the object.  
              Please run the assay or choose a different assay slot."))
  }
  
  obj_data <- assay(object = object, 
                    i = assay_use) %>% 
    as.matrix()
  
  return(obj_data)
}