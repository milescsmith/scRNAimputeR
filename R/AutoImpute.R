#' AutoImpute
#'
#' AutoImpute is an auto-encoder based gene-expression (sparse) matrix imputation.
#' (https://github.com/divyanshu-talwar/AutoImpute)
#' Some arguments have been removed due to being troublesome in the R-Python
#' interface
#'
#'
#' @param exprData 
#' @param debug 
#' @param debug_display_step 
#' @param hidden_units 
#' @param lambda_val 
#' @param initial_learning_rate 
#' @param iterations 
#' @param threshold 
#' @param masked_matrix_test 
#' @param masking_percentage 
#' @param log_file 
#' @param load_saved 
#'
#' @importFrom reticulate import py_module_available
#' @importFrom magrittr %<>% %>%
#' @importFrom glue glue
#'
#' @return
#' @export
#'
#' @examples
autoimpute <- function(object,
                       assay = 'RNA',
                       slot_use = 'data',
                       assay_store = "RNA",
                       normalize.data = FALSE,
                       scale.data = TRUE,
                       debug = TRUE, 
                       debug_display_step = 1, 
                       hidden_units = 2000, 
                       lambda_val = 1, 
                       initial_learning_rate = 0.0001,
                       iterations = 7000,
                       threshold = 0.0001,
                       masked_matrix_test = FALSE,
                       masking_percentage = 10,
                       log_file = 'log.txt',
                       load_saved = FALSE){
  
  required_modules <- c("AutoImpute", "numpy", "scipy", "tensorflow", "scikit-learn", "typing")
  for (i in required_modules){
    if(!py_module_available(i)){
      stop(glue("The {i} module is unavailable.
         Please activate the appropriate environment or
         install the module."))
    }
  }
  
  exprDat <- GatherData(object, assay, slot_use)
  ai.module <- import(module = 'AutoImpute', delay_load = TRUE)
  
  exprDat %<>% t()
  cell.names <- rownames(exprDat)
  gene.names <- colnames(exprDat)
  
  ai_data <- ai$autoimpute(data = exprData, 
                           debug = debug,
                           debug_display_step = as.integer(debug_display_step),
                           hidden_units = as.integer(hidden_units),
                           lambda_val = as.integer(lambda_val),
                           initial_learning_rate = as.numeric(initial_learning_rate),
                           iterations = as.integer(iterations),
                           threshold = as.numeric(threshold),
                           masked_matrix_test = masked_matrix_test,
                           masking_percentage = as.integer(masking_percentage),
                           log_file = log_file ,
                           load_saved = FALSE) %>% t()
  
  colnames(ai_data) <- cell.names
  rownames(ai_data) <- glue("AI_{gene.names}")
  object <- PlaceData(object = object,
                      assay_store = "autoimpute", 
                      imputed_data = ai_data)
  return(object)
}