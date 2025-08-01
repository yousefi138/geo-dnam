#require(meffil)
#require(parallel)
#options(mc.cores=20)
#require(eval.save)
#eval.save.dir(paths$cache)
#idats.dir = file.path(paths$output, "supp.files", "idats")
#save.dir = file.path(paths$output, "supp.files", "clean")
dir.create(reports.dir <- file.path(paths$output, "normalization-reports"))

samplesheet = meffil.create.samplesheet(idats.dir, recursive=T)

dataset = sapply(list.files(idats.dir, "\\.out$", full.names=T), readLines, simplify=F)
names(dataset) = sub("_RAW.out", "", basename(names(dataset)))
dataset = lapply(names(dataset), function(gse) data.frame(accession=gse, file=dataset[[gse]]))
dataset = do.call(rbind, dataset)

samplesheet$dataset = dataset$accession[match(samplesheet$Sample_Name, sub("_.*", "", dataset$file))]


##################################################################
## QC 
qc.objects = eval.save({
  mclapply(1:nrow(samplesheet), function(i) {
    try({
      meffil.create.qc.object(
        samplesheet[i,],
        cell.type.reference="blood gse35069 complete",
        featureset="450k:epic:epic2",
        verbose=TRUE)
    })
  })
}, "qc-objects", redo=F)
## 1 hour

table(sapply(qc.objects, class))
## list try-error 
## 5181         5

qc.objects = qc.objects[sapply(qc.objects, class)=="list"]

names(qc.objects) = sapply(qc.objects, function(x) x$sample.name)

report.filename = file.path(reports.dir, "qc-report.html")
if (!file.exists(report.filename)) {
  qc.summary = meffil.qc.summary(qc.objects)        
  meffil.qc.report(qc.summary, output.file=report.filename)
}

with(qc.summary$meth.unmeth.summary, quantile(tab$methylated))
##       0%       25%       50%       75%      100% 
## 827.9715 2976.1799 3205.6788 3382.1630 4869.9397 

with(qc.summary$meth.unmeth.summary, tab$sample.name[tab$methylated < 2000])
## [1] "GSM1236218" "GSM1236317" "GSM2334615" "GSM2656366"

with(qc.summary$controlmeans.summary,{
    sort(table(tab$sample.name[tab$outliers]),decreasing=T)
})
## GSM7020894 GSM7021250 GSM2656366 GSM1235872 GSM1236037 GSM1236194 GSM1236200 
##          4          4          3          2          2          2          2 
## GSM1235689 GSM1235912 GSM1236008 GSM1236056 GSM1236077 GSM1236139 GSM2334289 
##          1          1          1          1          1          1          1 
## GSM3853565 GSM6047256 
##          1          1 

with(qc.summary$sample.detectionp.summary, {
    tab[tab$prop.badprobes > 0.05,]
})
##      sample.name prop.badprobes colour.code   id outliers
## 1772  GSM2656366      0.6521804           1 1772     TRUE
## 4730  GSM7020894      0.9339097           1 4730     TRUE
## 5085  GSM7021250      0.9339097           1 5085     TRUE

qc.objects[["GSM2656366"]] = NULL
qc.objects[["GSM7020894"]] = NULL
qc.objects[["GSM7021250"]] = NULL

samplesheet = samplesheet[match(names(qc.objects), samplesheet$Sample_Name),]

##################################################################
## normalize quantiles
pcfit = eval.save(meffil.plot.pc.fit(qc.objects), "pcfit")
ggsave(file.path(reports.dir, "pcfit.pdf"), pcfit$plot)

number.pcs = 22

norm.objects = eval.save({
  meffil.normalize.quantiles(
    qc.objects,
    fixed.effects="dataset",
    number.pcs=number.pcs)
}, "norm-objects", redo=F)

## 2 minutes

##################################################################
## normalize dna methylation levels
for (dname in unique(samplesheet$dataset)) {
  is.dat = samplesheet$dataset == dname
  meth.filename = file.path(save.dir, paste0(tolower(dname), ".csv.gz"))
  report.filename = file.path(reports.dir, paste0(dname, "-norm.html"))
  if (!file.exists(meth.filename)) {
    meth = meffil.normalize.samples(
      norm.objects[is.dat],
      cpglist.remove=qc.summary$bad.cpgs$name,
      verbose=T)
    
    pcs = meffil.methylation.pcs(meth)
    norm.summary = meffil.normalization.summary(
      norm.objects[is.dat],
      pcs=pcs,
      parameters=meffil.normalization.parameters(
        norm.objects[is.dat],
        variables=c("Slide","sentrix_row")))
    meffil.normalization.report(
      norm.summary,
      output.file=report.filename)
    
    fwrite(
      data.frame(cg=rownames(meth), meth, check.names=F),
      file=meth.filename)
  }
}

## 1 hour
