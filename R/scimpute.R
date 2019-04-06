#' @title scImpute
#'
#' @description Impute data using scImpute
#'
#' @param object Data object to impute
#' @param assay Assay from which to retrieve expression data. Default: "RNA"
#' @param slot Slot from which to retrieve expression data. Default: "counts"
#' @param num_clusters Number of clusters to use when imputing data.  If left NULL,
#' will attempt to retrive the current identity information and use a count equal to the
#' number of unique identities.
#' @param ident_use Identity variable to use when imputing the data.
#' @param ensembl_db The Ensembl annotation package to use to compute gene lengths.
#' Default: EnsDb.Hsapiens.v86 (if available)
#' @param ... Additional arguments to pass to \code{\link{scimpute}}
#'
#' @importFrom scImpute scimpute
#'
#' @return Seurat
#' @export
scImpute <- function(object,
                     assay = "RNA",
                     slot = "counts",
                     num_clusters = NULL,
                     ident_use = NULL,
                     ensembl_db = NULL,
                     ...) {
  datExprs <- GatherData(object, assay = assay, slot_use = slot)
  if (is.null(ensembl_db)) {
    if(require(EnsDb.Hsapiens.v86)){
      ensembl_db <- EnsDb.Hsapiens.v86 
    } else {
      stop("Please provide an Ensembl-based annotation package")
    }
  } else {
    gl <- getGeneLengths(genes = rownames(datExprs), edb = ensembl_db)
  }

  if (is.null(num_clusters) & is.null(ident_use)) {
    num_clusters <- RetrieveIdents(object) %>% unique() %>% length()
    idents <- NULL
  } else if (is.null(num_clusters) & !is.null(ident_use)) {
    idents <- RetrieveIdents(object, identity = ident_use)
  }

  sc_data <- scimpute(
    datExpr = datExprs[gl$gene_name, ],
    type = "count",
    Kcluster = num_clusters,
    labels = idents,
    genelen = gl$length,
    ...
  )

  object <- PlaceData(
    object = object,
    imputed_data = sc_data,
    assay_store = "scimpute"
  )
  return(object)
}

#' @importFrom ensembldb cdsBy
#' @importFrom dplyr select mutate group_by top_n distinct
#' @importFrom AnnotationFilter GenenameFilter
#' @importFrom tibble as_tibble
getGeneLengths <- function(genes, edb = NULL) {
  if (is.null(edb)) {
    stop("Must specify the Ensembl annotation package appropriate for your species")
  }

  txs <- cdsBy(edb, filter = GenenameFilter(genes)) %>%
    unlist() %>%
    as_tibble() %>%
    mutate(length = end - start) %>%
    group_by(gene_name) %>%
    top_n(n = 1, wt = length) %>%
    dplyr::select(gene_name, length) %>%
    distinct()
  return(txs)
}
