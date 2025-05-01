#' Create a BedFileReader Object from the BedFileReader C++ Class
#'
#' later add more documentation
#' 
#' @name BedFileReader

# ^^^^^^^^^^^^^^^^
# Export the "BedFileReader" C++ class by explicitly requesting BedFileReader be
# exported via roxygen2's export tag.
# Also, provide a name for the Rd file.

loadModule(module = "BedFileReader_module", TRUE)
