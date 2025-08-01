## ----globals -------------------------------------------------------------
packages <- c(
  "meffil",
  "geograbi", # https://github.com/yousefi138/geograbi
  "eval.save",
  "purrr", # for map2)
  "easyPubMed",
  "dplyr",
  "parallel")
lapply(packages, require, character.only=T)

options(mc.cores=20)

# set dirs  
dir <- paths
if(!dir.exists(dir$cache)) dir.create(dir$cache)
eval.save.dir(dir$cache)

## check or make all output dirs
c("series.files", "dnam", "samples", "supp.files") |>
  sapply(function(.x){
    if (!dir.exists(file.path(dir$output,.x)))
      dir.create(file.path(dir$output,.x), recursive = T)
    dir.exists(file.path(dir$output,.x))		
  })

## source py written functions
#source(file.path(dir$scripts, "src/get.characteristics.r"))
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

## ----get.data -------------------------------------------------------------
eval.save({
  dnam <- 
    purrr::map2(gses$filenames, gses$accession, ~ {
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
gses$dnam.file <- 
		file.path(dir$output, "dnam", 
          paste(tolower(gses$accession), "csv.gz", sep="."))	

## omit any where the series matrix file contains no DNAm
gses$dnam.file[gses$nrow < 300000] = NA

# ----blood.pred -------------------------------------------------------------
## retrieve sample info for each dataset
eval.save({
  mclapply(1:nrow(gses), function(i) {
    ret = NULL
    	try({
        ret = geograbi.get.samples(gses$filename[i])
      })
      ret
    })
}, "samples", redo=F)
samples <- eval.ret("samples")
names(samples) <- gses$accession

# write all the sample info to file
samples.file <- purrr::imap(samples, ~ {
		if(!is.null(.x)){
			file <- file.path(dir$output, "samples", 
						paste0(tolower(.y), ".csv.gz"))
			if(!file.exists(file)) data.table::fwrite(.x, file)
			file
		} else NA
}) 
gses$samples.file <- unlist(samples.file)

## extract just the info on sample source
source <- mclapply(samples, function(i){
		i$source_name_ch1
})

omit.patterns <- paste0("(",paste(c("tumour","spot","tumor", "cancer", "placenta","tissue"),collapse="|"),")")
keep.patterns <- paste0("(",paste(c("blood", "WB", "PBMC"),collapse="|"),")")
 
## uses patterns to guess if the source is likely to be  peripheral blood
blood.pred = sapply(source, function(samp) {
	keep = NA
  if (length(samp) > 0) {
    keep = !grepl(omit.patterns,samp,ignore.case=T) & grepl(keep.patterns,samp,ignore.case=T)
  }
  keep
})
gses$blood.pred = sapply(blood.pred, all)

table(gses$blood.pred, useNA = "ifany")
#FALSE  TRUE  <NA> 
#  355   120   201 

## adhoc additions
gses$blood.pred <- (
  gses$blood.pred
  | gses$accession %in% c("GSE40279","GSE42861","GSE51032","GSE51057","GSE74548"))

## corrections
gses$ncol[gses$accession == "GSE42861"] = 689
gses$ncol[gses$accession == "GSE51032"] = 845

## adhoc removals
gses$blood.pred[gses$accession %in% readLines("omit.txt")] = FALSE

table(gses$blood.pred, useNA = "ifany")
## FALSE  TRUE  <NA> 
##   362   115   199 

#str(source[which(gses$blood.pred==F)])
## check on count
sum(gses[gses$blood.pred == TRUE & !is.na(gses$blood.pred), "ncol"])
# [1] 40953

## ----pubmed -------------------------------------------------------------
eval.save({
  pubmed <- retrieve.papers(gses$pmid) 
}, "pubmed", redo=F)
pubmed <- eval.ret("pubmed")

gses <- pubmed |>
  rename(c("pub.title" = "title")) |>
  select(all_of(c("pmid", "pub.title", "abstract")))  %>%
  left_join(gses, ., by = c("pmid"))


## ----blood.gses -------------------------------------------------------

blood.gses <- gses |>
  subset(blood.pred == TRUE & !is.na(blood.pred))

## ----supp.files------------------------------------------------------

## retrieve supplementary lists for each dataset
eval.save({
  mclapply(1:nrow(blood.gses), function(i) {
    ret = NULL
    if (blood.gses$nrow[i] < 300000) {
      try({
        ret = geograbi.list.supplementary.files(blood.gses$accession[i])
      })
      ret
    }})
}, "supplementary-files", redo=F)
supp.files <- eval.ret("supplementary-files")

## filename clues that a supplementary file contains normalized DNA methylation data
## or, if none exists, the raw IDAT files
omit.patterns <- paste0("(",paste(c("intensities","unmeth","unprocessed","non[-_]+normalized","unnormalized","unnormalised","raw","xlsx","tumor"),collapse="|"),")")
keep.patterns <- paste0("(",paste(c("normalized","normalised","methylation","processed","beta"),collapse="|"),")")

## uses patterns to guess appropriate supplementary file
## for datasets without series matrix files
blood.gses$supp.file = sapply(supp.files, function(files) {
  keep = files
  if (length(files) > 0) {
    keep = files[!grepl(omit.patterns,files,ignore.case=T) & grepl(keep.patterns,files,ignore.case=T)]
  }
  if (length(keep) == 0)
    keep = files[grepl("raw.tar$",files,ignore.case=T)]
  keep
})

## adhoc corrections/additions
additions = c(
  "GSE191082_Normalized_data_Blood.txt.gz",
  "GSE74414_all_beta_values.txt.gz",
  "GSE54882_processed_methylation_matrix.txt.gz",
  "GSE89218_Meth_163_matrix_processed.txt.gz",
  "GSE42861_processed_methylation_matrix.txt.gz",
  "GSE201322_RAW.tar",
  "GSE163970_RAW.tar",
  "GSE198904_DHRC_Matrix_processed.txt.gz;GSE198904_OBSERVEMDD0001_Matrix_processed.txt.gz")
names(additions) = sub("_.*", "", additions)

blood.gses = rbind(
  blood.gses[!blood.gses$accession %in% names(additions),],
  cbind(gses[gses$accession %in% names(additions),], supp.file=""))

blood.gses$supp.file[match(names(additions), blood.gses$accession)] = additions

stopifnot(all(sapply(blood.gses$supp.file,length) <= 1))

blood.gses$supp.file <- sapply(blood.gses$supp.file, function(filename) ifelse(length(filename) > 0, unlist(filename), ""))

stopifnot(all(blood.gses$supp.file != "" | blood.gses$nrow > 300000))

is.supp = blood.gses$supp.file != ""
supp.dir <- file.path(dir$output, "supp.files")
blood.gses$dnam.file[is.supp] = file.path(supp.dir, "clean", paste0(tolower(blood.gses$accession[is.supp]), ".csv.gz"))

## ----download.supp.files ----------------------------------------------------------

supp.indices <- which(blood.gses$supp.file != "")
dir.create(download.dir <- file.path(supp.dir, "downloads"))
for (i in supp.indices) {
  filenames = unlist(strsplit(blood.gses$supp.file[i], ";"))
  for (filename in filenames) {
    geograbi.download.supplementary.file(
      path=download.dir,
      gse=blood.gses$accession[i],
      filename=filename)
  }
}

## ----clean.supp.files -------------------------------------------------

## convert supplementary files to csv files where rows=cpgs and columns=samples
supp.filenames <- with(blood.gses, supp.file[grepl("(csv|tsv|txt).gz$", supp.file)])
supp.filenames <- unlist(strsplit(supp.filenames,";"))
supp.gses = unique(sub("_.*", "", supp.filenames))
dir.create(save.dir <- file.path(supp.dir, "clean"))

for (supp.gse in supp.gses) {
  out.filename <- file.path(save.dir, paste0(tolower(supp.gse),".csv.gz"))
  if (!file.exists(out.filename)) {
    cat(date(), " cleaning supplementary file and saving to ", out.filename, "\n")
    is.filename = grepl(paste0(supp.gse,"_"), supp.filenames)
    dats <- lapply(supp.filenames[is.filename], function(filename) {  
      dat <- fread(file.path(download.dir, filename), data.table=F)
      cg.idx = which(colnames(dat) == "ID_REF")
      cg.idx = ifelse(length(cg.idx) == 0, 1, cg.idx[1])
      rownames(dat) <- dat[[cg.idx]]
      dat[[cg.idx]] = NULL
      dat <- as.matrix(dat[,!grepl("(cpg|pval|id_ref)",colnames(dat),ignore.case=T)])
      rows <- rownames(dat)
      dat <- apply(dat,2,as.numeric)
      rownames(dat) = rows
      vs <- apply(dat,2,var,na.rm=T)
      stopifnot(!all(vs < 0.0001)) ## stop if all sample variance is tiny
      dat[,vs > 0.0001]            ## omit samples with little variance
    })
    if (length(dats) > 1) {
      rows <- unlist(lapply(dats,rownames))
      rows <- names(which(table(rows)==length(dats)))
      dats <- lapply(dats, function(dat) dat[match(rows,rownames(dat)),])
      dat <- do.call(cbind, dats)
    }
    else dat <- dats[[1]]
    fwrite(data.frame(cg=rownames(dat), dat, check.names=F), file=out.filename)
  }
}

## ----idats.extraction --------------------------------------------------------------------

arc.filenames <- with(blood.gses, supp.file[grepl("_RAW.tar$", supp.file)])
dir.create(idats.dir <- file.path(supp.dir, "idats"))
for (filename in arc.filenames) {
  contents.out <- file.path(idats.dir, sub(".tar",".out",filename))
  if (!file.exists(contents.out)) {
    cat(date(), " extracting idats from ", filename, "\n")
    system(paste("tar -tf", file.path(download.dir, filename), ">", contents.out))
    system(paste("tar -xf", file.path(download.dir, filename), "-C", idats.dir))
  }
}

## ----normalize.idats ----------------------------------------------------------

## in:  paths$output, paths$cache,
##      idats.dir ([paths$output]/supp.files/idats),
##      save.dir ([paths$output]/supp.files/clean)
## out: normalized dna methylation files: ([idats.dir]/gse*.csv.gz),
##      normalization reports: [paths$output]/normalization-reports/*
source("normalize-idat-files.r",echo=T)
## ~2 hours to run

## ----check.sample.counts -------------------------------------------------------------

## compare samples counts from blood.gses$ncol to sample counts from data files

## datasets with a series matrix data
sm.counts = sapply(blood.gses$accession, function(gse) {
  filename = file.path(dir$output, "dnam", paste0(tolower(gse), ".csv.gz"))
  if (file.exists(filename)) {
    ret = system(paste("gunzip -c", filename, "| head -n1 | tr -cd , | wc -c"), intern=TRUE)
    as.integer(ret) + 1
  } else 0
})
sum(sm.counts)
## [1] 39622

## datasets with supplementary data
supp.counts = system(
  paste(
    "for i in", file.path(supp.dir, "clean", "*.gz"), "; do",
    "echo `basename $i .csv.gz` `gunzip -c $i | head -n1 | tr -cd , | wc -c`; done"),
  intern=TRUE)
names(supp.counts) = toupper(sub(" .*", "", supp.counts))
supp.counts = sub(".* ", "", supp.counts)
supp.counts = sapply(supp.counts, as.integer)
sum(supp.counts)
## [1] 36089

## datasets with idat files
idat.counts = sapply(list.files(file.path(supp.dir, "idats"), "out$", full.names=T), function(filename) {
  sum(grepl("_red.idat.gz", readLines(filename), ignore.case=T))
})
names(idat.counts) = sub("_RAW.out", "", basename(names(idat.counts)))
sum(idat.counts)
## 4465

blood.gses$sm.samples = 0
blood.gses$sm.samples[match(names(sm.counts), blood.gses$accession)] = sm.counts
blood.gses$supp.counts = 0
blood.gses$supp.counts[match(names(supp.counts), blood.gses$accession)] = supp.counts
blood.gses$idat.counts = 0
blood.gses$idat.counts[match(names(idat.counts), blood.gses$accession)] = idat.counts

stopifnot(
  with(blood.gses, {
    all(ncol == 0 | sm.samples == 0 | abs(ncol-sm.samples)<2)
  }))

stopifnot(
  with(blood.gses, {
    all(accession == "GSE191082" | ncol == 0 | supp.counts == 0 | abs(ncol-supp.counts) < 10)
  }))

stopifnot(
  with(blood.gses, {
    all(ncol == 0 | idat.counts == 0 | abs(ncol-idat.counts) < 2)
  }))

stopifnot(all(file.exists(blood.gses$dnam.file)))

## ----common.cpgs ------------------------------------------------------------

cpgs = eval.save({
  extract.cpgs = function(filename) {
    cat(date(), " extracting cpg sites from ", filename, "\n") 
    system(paste("gunzip -c", filename, "| cut -d ',' -f1"), intern=T)[-1]
  }    
  ret = sapply(blood.gses$dnam.file, extract.cpgs, simplify=F)
  names(ret) = sub(".csv.gz$", "", basename(names(ret)))
  ret
}, "cpgs", redo=F)
## ~1 hour

quantile(sapply(cpgs, length))
##       0%      25%      50%      75%     100% 
## 369696.0 458636.2 485577.0 785134.5 937690.0 

cpg.freqs = table(unlist(sapply(cpgs, function(cpgs) unique(sub("_.*","", cpgs)), simplify=F)))
## beware: cg00000029_TC21

underscores = sapply(cpgs, function(cpgs) sum(grepl("_", cpgs)))

common.cpgs = names(cpg.freqs)[which(cpg.freqs > nrow(blood.gses)*0.75)]
                 
length(common.cpgs)
## [1] 401709

writeLines(common.cpgs, con=file.path(paths$output, "common-cpgs.txt"))

## -----harmonize.cpgs --------------------------------------------------------------

## harmonize cpg identifiers in one dataset
acc = "GSE246337"
current.filename = blood.gses$dnam.file[which(blood.gses$accession == acc)]
new.filename = file.path(dirname(current.filename), paste0(tolower(acc), "-fixed.csv.gz"))
if (!file.exists(new.filename)) {
  meth = fread(current.filename, data.table=F)
  cg = sub("_.*$", "", meth$cg)
  meth = meth[match(unique(cg), cg),]
  meth$cg = sub("_.*$", "", meth$cg)
  fwrite(meth, file=new.filename)
}
blood.gses$dnam.file[which(blood.gses$accession == acc)] = new.filename

## -----remove.bad.samples --------------------------------------------------------------

## two datasets have samples with entirely missing data
accs = c("GSE59592", "GSE69138")
for (acc in accs) {
  current.filename = blood.gses$dnam.file[which(blood.gses$accession == acc)]
  new.filename = file.path(dirname(current.filename), paste0(tolower(acc), "-fixed.csv.gz"))
  if (!file.exists(new.filename)) {
    meth = fread(current.filename, data.table=F)
    not.missing = sapply(meth, function(v) sum(!is.na(v)))
    fwrite(meth[,which(not.missing > 0)], file=new.filename)
  }
  blood.gses$dnam.file[which(blood.gses$accession == acc)] = new.filename
}  

## -----cell.counts -----------------------------------------------------------

dir.create(counts.dir <- file.path(paths$output, "cell-counts"))
ret = lapply(1:nrow(blood.gses), function(i) {
  meth.filename = blood.gses$dnam.file[i]
  counts.filename = file.path(counts.dir, basename(meth.filename))
  if (!file.exists(counts.filename)) {
    cat(date(), " estimating cell counts for ", basename(counts.filename), "\n")
    try({
      meth = fread(meth.filename, data.table=F)    
      rownames(meth) = meth[[1]]
      meth[[1]] = NULL
      counts = meffil.estimate.cell.counts.from.betas(as.matrix(meth), "blood gse35069 complete")
      fwrite(data.frame(sample=rownames(counts), counts), file=counts.filename)
      return(TRUE)
    })
    return(FALSE)
  }
  return(TRUE)
}) ## ~3 hours

blood.gses$counts.file = file.path(counts.dir, basename(blood.gses$dnam.file))
stopifnot(all(sapply(blood.gses$counts.file, file.exists)))

## ----write.gses -------------------------------------------------------------

# full query results
data.table::fwrite(
  gses, 
  file.path(dir$output, "gses.csv"),
  row.names = F,
  na = "NA")

## datasets derived from blood
data.table::fwrite(
  blood.gses,
  file.path(dir$output, "blood.gses.csv"),
  row.names = F,
  na = "NA")


             



