.onLoad <- function(libname, pkgname) {
  required_modules <- c("AutoImpute", "numpy", "scipy", "tensorflow", "sklearn", "typing", 'dca', 'scanpy')
  for (i in required_modules){
    if(!py_module_available(i)){
      stop(glue("The {i} module is unavailable.
         Please activate the appropriate environment or
         install the module."))
    }
  }
  
  invisible()
}