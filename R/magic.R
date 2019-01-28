Magic <- function(object, ...) {
  UseMethod("magic")
}

Magic.Seurat <- function(object,
                         assay = "RNA",
                         slot.use = "counts",
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
  if (!"pca" %in% names(object)) {
    message("PCA must be computed before running Harmony.  Since it is missing, will not compute PCA...")
    object <- RunPCA(object, assay = assay.use)
  }
  if (ncol(object[['pca']]@cell.embeddings) < npca){
    message("Number of available PCs is less than the npcs argument passed.  Adjusting...")
    npcs <- ncol(object[['pca']]@cell.embeddings)
  }

  required_modules <- c("magic")
  for (i in required_modules){
    if(!py_module_available(i)){
      stop(glue("The {i} module is unavailable.
         Please activate the appropriate environment or
         install the module."))
    }
  }

  exprDat <- GatherData(object, assay, slot.use) %>% t()
  magic_result <- magic(data = exprDat, 
                        genes, 
                        k, 
                        alpha, 
                        t, 
                        npca, 
                        init, 
                        t.max, 
                        knn.dist.method, 
                        verbose, 
                        n.jobs, 
                        seed)
}