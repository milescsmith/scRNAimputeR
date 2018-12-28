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
#' @importFrom Seurat SetAssayData
#' @importFrom magrittr %>%
#' @importFrom Matrix Matrix
#'
#' @return matrix
#' @export
#'
#' @examples
setMethod("PlaceData",
          signature(object='Seurat'),
          function(object,
                   assay = "new_assay",
                   slot.use = "data",
                   imputed_data,
                   normalize = "FALSE",
                   scale = "FALSE",
                   ...){
            
            data <- Matrix(data = data, sparse = TRUE)
            object[[assay]] <- CreateAssayObject(data = new.data = imputed_data)
            if (normalize)
            {
              object <- NormalizeData(object, assay = assay)
            }
            if (scale){
              object <- ScaleData(object, assay = assay, ...)
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
