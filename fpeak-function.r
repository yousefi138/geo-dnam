fpeak = function(file, nrows=10, ...) {
  stopifnot(!is.null(file))
  fread(cmd = paste("gzip -dc", shQuote(file), "| head -n", nrows), ...)
}
