#' @title PlaceData
#'
#' @param object Data object to pull data from
#' @param assay_store Name of assay slot to place imputated data
#' @param slot_use Assay slot to pull data from
#' @param imputed_data Imputed matrix
#' @param normalize Should the imputed data be normalized?
#' @param scale Should the imputed data be scaled?
#' @param ... Additional arguments
#'
#' @return
#' @export
#'
#' @examples
PlaceData <- function(object, ...) {
  UseMethod("PlaceData")
}

#' @rdname PlaceData
#' @method PlaceData Seurat
# @importFrom Seurat CreateAssayObject NormalizeData ScaleData
#' @import Seurat
#' @importFrom Matrix Matrix
#' @return
#' @export
PlaceData.Seurat <- function(object,
                             assay_store = "new_assay",
                             imputed_data,
                             normalize = "FALSE",
                             scale = "FALSE",
                             ...) {
  imputed_data %<>% as.matrix()
  object[[assay_store]] <- CreateAssayObject(data = imputed_data)
  if (normalize) {
    object <- NormalizeData(object, assay = assay_store)
  }
  if (scale) {
    object <- ScaleData(object, assay = assay_store)
  }
  return(object)
}

#' @rdname PlaceData
#' @method PlaceData seurat
# @importFrom Seurat SetAssayData NormalizeData ScaleData
#' @import Seurat
#' @importFrom Matrix Matrix
#' @return
#' @export
PlaceData.seurat <- function(object,
                             assay_store = "new_assay",
                             slot_use = "data",
                             imputed_data,
                             normalize = "FALSE",
                             scale = "FALSE",
                             ...) {
  data <- Matrix(data = imputed_data, sparse = TRUE)
  object <- SetAssayData(
    object = object,
    assay = assay_store,
    slot = slot_use,
    new.data = imputed_data
  )
  if (normalize) {
    object <- NormalizeData(object, assay.type = assay_store)
  }
  if (scale) {
    object <- ScaleData(object)
  }
  return(object)
}

#' @rdname PlaceData
#' @method PlaceData SingleCellExperiment
#' @importFrom Matrix Matrix
#' @importFrom scater normalize
#' @importFrom SummarizedExperiment assay assay<-
#' @return
#' @export
PlaceData.SingleCellExperiment <- function(object,
                                           assay_store = "new_assay",
                                           imputed_data,
                                           normalize = "FALSE",
                                           ...) {
  data <- Matrix(data = imputed_data, sparse = TRUE)
  assay(object = object, i = assay_store) <- data
  if (normalize) {
    object <- normalize(object, exprs_values = assay_store)
  }
  return(object)
}

#' @title GatherData
#'
#' @param object Seurat object to pull data from
#' @param assay (Default: "RNA")
#' @param slot_use (Default: "data)
#' @param ... Additional arguments
GatherData <- function(object, ...) {
  UseMethod("GatherData")
}

#' @rdname GatherData
#' @method GatherData Seurat
#' @importFrom Seurat GetAssayData
#' @return matrix
#' @export
GatherData.Seurat <- function(object,
                              assay = "RNA",
                              slot_use = "counts",
                              ...) {
  obj_data <- GetAssayData(
    object = object,
    assay = assay,
    slot = slot_use
  ) %>%
    as.matrix()
  return(obj_data)
}

#' @rdname GatherData
#' @method GatherData seurat
#' @importFrom Seurat GetAssayData
#' @return matrix
#' @export
GatherData.seurat <- function(object,
                              assay = "RNA",
                              slot_use = "data",
                              ...) {
  obj_data <- GetAssayData(
    object = object,
    assay.type = assay,
    slot = slot_use
  ) %>%
    as.matrix()
  return(obj_data)
}

#' @rdname GatherData
#' @method GatherData SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays
#' @return matrix
#' @export
GatherData.SingleCellExperiment <- function(object,
                                            assay = "logcounts",
                                            ...) {
  if (!assay %in% names(assays(object))) {
    stop(glue("Sorry, but {assay} is not present in the object.  
              Please run the assay or choose a different assay slot."))
  }

  obj_data <- assay(
    object = object,
    i = assay
  ) %>%
    as.matrix()

  return(obj_data)
}


#' @title RetrieveIdents
#'
#' @description Helper function to retrive cell identity information
#'
#' @param object Data object from which to retrieve identiy information
#' @param identity Variable to use for grouping identity
#' @param ... Additional arguments to pass to child functions
#'
#' @return A list corresponding to the identities of each cell.
#' @export
#'
#' @examples
RetrieveIdents <- function(object, ...) {
  UseMethod("RetrieveIdents")
}

#' @rdname RetrieveIdents
#' @method RetrieveIdents Seurat
#' @import Seurat
#' @return matrix
#' @export
RetrieveIdents.Seurat <- function(object,
                                  identity = NULL,
                                  ...) {
  if (is.null(identity)) {
    idents <- Idents(object)
  } else {
    idents <- FetchData(object = object, vars = identity, ...)
  }
  return(idents)
}

#' @rdname RetrieveIdents
#' @method RetrieveIdents seurat
#' @import Seurat
#' @return matrix
#' @export
RetrieveIdents.seurat <- function(object,
                                  identity = NULL,
                                  ...) {
  if (is.null(identity)) {
    idents <- object@ident
  } else {
    idents <- FetchData(object = object, vars.all = identity, ...)
  }
  return(idents)
}

#' @rdname RetrieveIdents
#' @method RetrieveIdents SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @return matrix
#' @export
RetrieveIdents.SingleCellExperiment <- function(object,
                                                identity = NULL,
                                                ...) {
  if (is.null(identity)) {
    if ("ident" %in% colnames(colData(object))) {
      idents <- colData(object)[[identity]]
    } else {
      stop("No identity variable given and no default ident present.")
    }
  } else {
    idents <- colData(object)[[identity]]
  }
  return(idents)
}
