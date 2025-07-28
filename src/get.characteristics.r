get.characteristics <- function(filenames){

	geo <- lapply(filenames, function(i){
	 		  gse <- geograbi.get.samples(i)
	 		  data.frame(cbind(gse, 
	 		  	geograbi.extract.characteristics(gse)), 
				stringsAsFactors = F)
			})
	names(geo) <- basename(filenames)

	chrs <- lapply(geo, geograbi.extract.characteristics)
	chr.tabs <- lapply(chrs, function(i) {
					out <- lapply(na.omit(i), function(i2) {
						ret <- table(i2)
						if (length(ret) <6) return(ret)
				})
				if (length(out) >56 ) return(na.omit(names(i)))
				out
	})

	chr.tabs <- lapply(chr.tabs, function(i) paste(names(i), i))

	## set all chr.names to the max name length
	ncols <- max(sapply(chr.tabs, length))
	chr.tabs <- lapply(chr.tabs, function(i){
		length(i) <- ncols
		i
	})

	chr.tabs <- do.call(rbind, chr.tabs)

	colnames(chr.tabs) <- paste0("chr.fld.",seq(1, ncol(chr.tabs)))
	list(chr.tabs = 
			data.frame(chr.tabs, 
			stringsAsFactors = F),
		chrs = chrs)
}

