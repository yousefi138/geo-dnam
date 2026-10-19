normalize.idat.archive = function(archive,norm.dir,output) {
  require(meffil)
  require(parallel)
  options(mc.cores=20)

  idats.dir <- tempfile(pattern="idats", tmpdir=norm.dir)
  dir.create(idats.dir, recursive=TRUE)
  on.exit(unlink(idats.dir, recursive = TRUE), add = TRUE)

  system(paste("tar -xf", archive, "-C", idats.dir))

  samplesheet = meffil.create.samplesheet(idats.dir, recursive=T)

  qc.objects.filename = file.path(norm.dir, "qc-objects.rds")

  if (!file.exists(qc.objects.filename)) { 
    qc.objects=mclapply(1:nrow(samplesheet), function(i) {
      try({
        meffil.create.qc.object(
          samplesheet[i,],
          cell.type.reference="blood gse35069 complete",
          featureset="450k:epic:epic2",
          verbose=TRUE)
      })
    })
    
    qc.objects = qc.objects[sapply(qc.objects, class)=="list"]
    names(qc.objects) = sapply(qc.objects, function(x) x$sample.name)
    
    report.filename = file.path(norm.dir, "qc-report.html")
    qc.summary = meffil.qc.summary(qc.objects)        
    meffil.qc.report(qc.summary, output.file=report.filename)
    
    bad.samples = with(qc.summary$sample.detectionp.summary, {
      tab[which(tab$prop.badprobes > 0.2),]
    })  
    if (nrow(bad.samples) > 0) 
      qc.objects = qc.objects[!names(qc.objects) %in% bad.samples$sample.name]

    saveRDS(qc.objects, file=qc.objects.filename)
  } else {
    qc.objects = readRDS(qc.objects.filename)
    qc.summary = meffil.qc.summary(qc.objects)        
  }
  
  samplesheet$failed = !(samplesheet$Sample_Name %in% names(qc.objects))
  fwrite(samplesheet, file=file.path(norm.dir, "samplesheet.csv"), row.names=F)

  pcfit.filename = file.path(norm.dir, "pc-fit.csv")
  if (!file.exists(pcfit.filename)) {
    pcfit = meffil.plot.pc.fit(qc.objects)
    fwrite(pcfit$data, file=pcfit.filename, row.names=F)
  } else
    pcfit = list(data=fread(pcfit.filename, data.table=F))
  number.pcs = floor(nrow(pcfit$data)/2)

  norm.objects.filename = file.path(norm.dir, "norm-objects.rds")
  if (!file.exists(norm.objects.filename)) {
    norm.objects = meffil.normalize.quantiles(
      qc.objects,
      number.pcs=number.pcs)
    saveRDS(norm.objects, file=norm.objects.filename)
  } else {
    norm.objects = readRDS(norm.objects.filename)
    basenames = samplesheet$Basename[match(names(norm.objects), samplesheet$Sample_Name)]
    for (i in 1:length(norm.objects))
      norm.objects[[i]]$basename = basenames[i]
  }

  if (!file.exists(output)) {
    meth = meffil.normalize.samples(
      norm.objects,
      cpglist.remove=qc.summary$bad.cpgs$name,
      verbose=T)
    fwrite(
      data.frame(cg=rownames(meth), meth, check.names=F),
      file=output)
  } else
    meth = fread(output, data.table=F, check.names=F, tmpdir=".")
}
