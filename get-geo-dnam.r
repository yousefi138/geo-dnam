## ----globals -------------------------------------------------------------
packages <- c("geograbi", # https://github.com/yousefi138/geograbi
            "eval.save",
			"purrr", # for map2)
            "easyPubMed",
			"dplyr")
lapply(packages, require, character.only=T)

# set dirs  
dir <- paths
if(!dir.exists(dir$cache)) dir.create(dir$cache)
eval.save.dir(dir$cache)

## check or make all output dirs
c("series.files", "dnam") |>
	sapply(function(.x){
		if (!dir.exists(file.path(dir$output,.x)))
			dir.create(file.path(dir$output,.x), recursive = T)
		dir.exists(file.path(dir$output,.x))		
	})

## source py written functions
source(file.path(dir$scripts, "src/get.characteristics.r"))
source(file.path(dir$scripts, "src/retrieve.papers.r"))

## ----gses -------------------------------------------------------------
eval.save({
	## get gses and basic info on all methylation array records
	gses <- geograbi.retrieve.datasets(
						gpl = c(
							"GPL13534", # 450k
							"GPL21145", "GPL21145", # EPIC
							"GPL33022")) |> ## 5sec
						subset(organism =="Homo sapiens") |>
						subset(samples >100)

}, "gses", redo=F)
gses <- eval.ret("gses")

## ----series.files -------------------------------------------------------------
eval.save({
	## download the series files, keep filename
	## ## 2 hours
	filenames <- geograbi.download.series.files( 
				path = file.path(dir$output, "series.files"),
				gses = gses$accession, 
				gpls = gses$platform)

}, "filenames", redo=F)
filenames <- eval.ret("filenames")
gses$filenames <- ifelse(filenames == "", NA, filenames)

## ----get.pmids -------------------------------------------------------------
eval.save({
	## get pmids from series files
	series <- lapply(na.omit(gses$filenames), geograbi.get.series) # 2mins
	keep <-sapply(series, function(i) any(grepl("pubmed", names(i))))
	series<- series[keep]
}, "series", redo=F)
series <- eval.ret("series")

series <- lapply(series, function(i){
	keep <- c(grep("accession", names(i)), grep("pubmed", names(i)))[1:2]
	i[keep]
})
series <- as.data.frame(do.call(rbind, series), 
				stringsAsFactors = F)
names(series)[names(series) == "pubmed_id"] <- "pmid"

## add pmids to gses
gses <- merge(gses, series, 
				by.x="accession", 
				by.y="geo_accession",
				all.x =T)

## ----get.chars -------------------------------------------------------------
# eval.save({
# 
# 	chars <- get.characteristics(gses$filenames)
# 
# }, "chars", redo=F)
# chars <- eval.ret("chars")
# 
# gses <- cbind(gses, 
# 			chars$chr.tabs)
# 
## ----get.data -------------------------------------------------------------
eval.save({
	dnam <- 
		purr::map2(gses$filenames, gses$accession, ~ {
						path <- file.path(dir$output, "dnam")		
						ret <-tryCatch(geograbi.get.data(filename = .x,
										gse = .y,
										path = path),
							error=function(e) e)

						if (any(grepl("error", class(ret)))){
							return(c(class(ret), dim(ret), length(ret)))
						}

						if (length(ret)==1 && ret==-1){
							gse <- tolower(.y)
							data.filename <- file.path(path, 
								paste(tolower(gse), "csv.gz", sep="."))
							ret <- as.data.frame(data.table::fread(data.filename))
						}
						c(class(ret)[1], dim(ret), length(ret))
		})
}, "dnam", redo=F)
dnam <- eval.ret("dnam")

# make a nice summry data frame of  the success status of the download
names(dnam) <- gses$accession
dnam <- lapply(dnam, function(i){
				if (any(grepl("Error|numeric", i))) {
					c("error", rep("0",3))
				} else i
			})
dnam <- do.call(rbind, dnam)
identical(gses$accession, rownames(dnam))

colnames(dnam) <- c("class", "nrow", "ncol", "length")
dnam <- data.frame(dnam, stringsAsFactors = F)
dnam[grep("[0-9]", dnam[1,])] <- 
						lapply(dnam[grep("[0-9]", 
							dnam[1,])], as.numeric)

gses <- cbind(gses, dnam)

## ----pubmed -------------------------------------------------------------
eval.save({
	pubmed <- retrieve.papers(gses$pmid) 
}, "pubmed", redo=F)
pubmed <- eval.ret("pubmed")

gses <- pubmed |>
		rename(pub.title = title) |>
		select(all_of(c("pmid", 
					"pub.title", "abstract")))  %>%
		left_join(gses, ., by = c("pmid"))

## ----write.gses -------------------------------------------------------------
data.table::fwrite(gses, 
	file.path(dir$output, "gses.csv"),
	row.names = F,
	na = "NA")
