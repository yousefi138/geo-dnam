library(data.table)

#' gread
#'
#' Similar to fread but handles output from write.table/write.csv
#' where the first column is row names.
#' 
gread <- function(..., fill=TRUE) {
  x = fread(..., fill=TRUE)
  if (tail(colnames(x),1) == paste0("V", ncol(x))) 
    colnames(x) = c("V0", colnames(x)[-ncol(x)])
  x
}

