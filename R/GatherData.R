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