setGeneric(
  name = "PlaceData",
  def = function(object, ...) {
    standardGeneric("PlaceData")
  }
)

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
#' @return matrix
#' @export
#'
#' @examples
setMethod("PlaceData",
          signature(object='Seurat'),
          function(object,
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
          })

setMethod("PlaceData",
          signature(object='seurat'),
          function(object,
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
          })

setGeneric(
  name = "GatherData",
  def = function(object, ...) {
    standardGeneric("GatherData")
  }
)

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
#' @return matrix
#' @export
#'
#' @examples
setMethod("GatherData",
          signature(object='Seurat'),
          function(object,
                   assay = "RNA",
                   slot.use = "data",
                   ...){
            
            obj_data <- GetAssayData(object = object, 
                                     assay = assay,
                                     slot = slot.use) %>%
              as.matrix()
            return(obj_data)
          })

setMethod("GatherData",
          signature(object='seurat'),
          function(object,
                   assay = "RNA",
                   slot.use = "data",
                   ...){
            obj_data <- GetAssayData(object = object, 
                                     assay = assay,
                                     slot = slot.use) %>%
              as.matrix()
            return(obj_data)
          })