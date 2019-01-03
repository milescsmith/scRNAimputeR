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
                   assay.store = "new_assay",
                   imputed_data,
                   normalize = "FALSE",
                   scale = "FALSE",
                   ...){
            
  imputed_data %<>% as.matrix()
  object[[assay.store]] <- CreateAssayObject(data = imputed_data)
  if (normalize)
  {
    object <- NormalizeData(object, assay = assay.store)
  }
  if (scale){
    object <- ScaleData(object, assay = assay.store, ...)
  }
  return(object)
}

#' @rdname PlaceData
#' @method PlaceData seurat
#' @return
#' @export
PlaceData.seurat <- function(object,
                   assay = "new_assay",
                   slot.use = "data",
                   imputed_data,
                   normalize = "FALSE",
                   scale = "FALSE",
                   ...){
            
  data <- Matrix(data = data, sparse = TRUE)
  object <- SetAssayData(object = seuratObj, 
                         assay = assay,
                         slot = slot.use,
                         new.data = imputed_data)
  if (normalize){
    object <- NormalizeData(object, assay.type = assay)
  }
  if (scale){
    object <- ScaleData(object, ...)
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
#' @param slot.use (Default: "data)
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
                              slot.use = "data",
                              ...){
            
  obj_data <- GetAssayData(object = object, 
                           assay = assay,
                           slot = slot.use) %>%
    as.matrix()
  return(obj_data)
}

#' @rdname GatherData
#' @method GatherData seurat
#' @return matrix
#' @export
GatherData.seurat <- function(object,
                              assay = "RNA",
                              slot.use = "data",
                              ...){
  obj_data <- GetAssayData(object = object, 
                           assay = assay,
                           slot = slot.use) %>%
    as.matrix()
  return(obj_data)
}