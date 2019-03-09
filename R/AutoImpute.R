#' @title AutoImpute
#'
#' @description \href{https://github.com/divyanshu-talwar/AutoImpute}{AutoImpute} is an
#' auto-encoder based gene-expression (sparse) matrix imputation.
#'
#' @param debug To print and save debug statements (loss function value). Default: TRUE
#' @param debug_display_step Number of steps to display loss function value after. Default: TRUE
#' @param hidden_units Size of hidden layer or latent space dimensions. Default: 2000
#' @param lambda_val Regularization coefficient, to control the contribution of
#' regularization term in the cost function. Default: 1
#' @param initial_learning_rate Initial value of learning rate. Default: 0.0001
#' @param iterations Number of iterations to train the model for. Default: 7000
#' @param threshold To stop gradient descent after the change in loss function
#' value in consecutive iterations is less than the threshold,
#' implying convergence. Default: 0.0001
#' @param masked_matrix_test Run the masked matrix recovery test? Default: FALSE
#' @param masking_percentage Percentage of masking required. Like 10, 20, 12.5 etc. Default: 10
#' @param log_file Text file to save training logs. Default: log.txt
#' @param load_saved Flag to indicate if a saved model will be loaded. Default: FALSE
#' @param object Data object to impute
#' @param assay Assay to take data from
#' @param slot_use Slot within assay to take data from
#' @param normalize_data Should the data be normalized after imputation?
#' @param scale_data Should the data be scaled after imputation?
#'
#' @importFrom reticulate import py_module_available
#' @importFrom glue glue
#'
#' @return
#' @export
#'
#' @examples
autoimpute <- function(object,
                       assay = "RNA",
                       slot_use = "data",
                       normalize_data = FALSE,
                       scale_data = TRUE,
                       debug = TRUE,
                       debug_display_step = 1,
                       hidden_units = 2000,
                       lambda_val = 1,
                       initial_learning_rate = 0.0001,
                       iterations = 7000,
                       threshold = 0.0001,
                       masked_matrix_test = FALSE,
                       masking_percentage = 10,
                       log_file = "log.txt",
                       load_saved = FALSE) {
  required_modules <- c("AutoImpute", 
                        "numpy", 
                        "scipy", 
                        "tensorflow", 
                        "sklearn", 
                        "typing")
  for (i in required_modules) {
    if (!py_module_available(i)) {
      stop(glue("The {i} module is unavailable.
         Please activate the appropriate environment or
         install the module."))
    }
  }

  exprDat <- GatherData(object, assay, slot_use)
  ai <- import(module = "AutoImpute", delay_load = TRUE)

  exprDat %<>% t()
  cellnames <- rownames(exprDat)
  genenames <- colnames(exprDat)

  ai_data <- ai$AutoImpute$autoimpute(
    data = exprDat,
    debug = debug,
    debug_display_step = as.integer(debug_display_step),
    hidden_units = as.integer(hidden_units),
    lambda_val = as.integer(lambda_val),
    initial_learning_rate = as.numeric(initial_learning_rate),
    iterations = as.integer(iterations),
    threshold = as.numeric(threshold),
    masked_matrix_test = masked_matrix_test,
    masking_percentage = as.integer(masking_percentage),
    log_file = log_file,
    load_saved = FALSE
  ) %>%
    t() %>%
    `[[`(1) %>%
    `[[`(1) %>%
    `[[`(1)

  colnames(ai_data) <- cellnames
  rownames(ai_data) <- genenames
  object <- PlaceData(
    object = object,
    assay_store = "autoimpute",
    imputed_data = ai_data
  )
  return(object)
}
