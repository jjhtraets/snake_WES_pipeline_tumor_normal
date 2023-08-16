## adapted from cnv_facets (https://github.com/dariober/cnv_facets/blob/master/bin/cnv_facets.R#L333) ##
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library in container
.libPaths(Sys.getenv("R_LIBS_USER"))

#install/load the libraries, TODO create docker file!
library(facets)
install.packages(c("data.table","optparse","BiocManager"),repos = "http://cran.us.r-project.org",quiet = TRUE,verbose=FALSE)
library(data.table)
library(optparse)
BiocManager::install("GenomicRanges",quietly = TRUE)
library(GenomicRanges)

set.seed(1234)

VERSION= sprintf('facets=%s', packageVersion('facets'))

filter_rcmat<- function(rcmat, min_ndepth, max_ndepth, target_bed){
  rcmat_flt<- rcmat[NOR.DP >= min_ndepth & NOR.DP < max_ndepth]
  if(is.null(target_bed) || nrow(rcmat_flt) == 0){
    return(rcmat_flt)
  }
  targets<- makeGRangesFromDataFrame(target_bed, seqnames.field= 'V1', 
                                     start.field= 'V2', end.field= 'V3', starts.in.df.are.0based= TRUE)
  
  # Convert rcmat to GRanges
  rcmat_flt<- makeGRangesFromDataFrame(rcmat_flt, seqnames.field= 'Chromosome', 
                                       start.field= 'Position', end.field= 'Position', keep.extra.columns= TRUE)
  
  hits<- findOverlaps(query= rcmat_flt, subject= targets, ignore.strand= TRUE)
  
  rcmat_flt<- as.data.table(rcmat_flt[unique(hits@from)])
  rcmat_flt[, start := NULL]
  rcmat_flt[, width := NULL]
  rcmat_flt[, strand := NULL]
  setnames(rcmat_flt, c('seqnames', 'end'), c('Chromosome', 'Position'))
  rcmat_flt[]
  return(rcmat_flt)
}

readSnpMatrix2<- function(pileup, gbuild){
  xf<- file(pileup, open= 'r')
  if(summary(xf)$class == 'gzfile'){
    conn<- sprintf('gzip -d -c %s', pileup)
  } else {
    conn<- pileup
  }
  close(xf)
  
  options(datatable.fread.input.cmd.message= FALSE)
  rcmat<- fread(conn, select= c('Chromosome', 'Position', 'File1R', 'File1A', 'File2R', 'File2A'))
  options(datatable.fread.input.cmd.message= TRUE)
  
  setnames(rcmat, c('File1R', 'File1A', 'File2R', 'File2A'),
           c('NOR.RD', 'NOR.DP', 'TUM.RD', 'TUM.DP'))
  rcmat[, NOR.DP := NOR.DP + NOR.RD]
  rcmat[, TUM.DP := TUM.DP + TUM.RD]
  
  chr_prefix<- any(rcmat$Chromosome %in% c(paste0('chr', 1:22), 'chrX'))
  
  # Unfortunately, facets needs numeric chromsomes. X will be converted later by facets
  rcmat[, Chromosome := sub("^chr", "", Chromosome)]
  setcolorder(rcmat, c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD', 'TUM.DP', 'TUM.RD'))
  
  # We only keep the major chromosomes 
  if(gbuild %in% c("hg19", "hg38")){
    rcmat<- rcmat[Chromosome %in% c(1:22, 'X')]
  } else if(gbuild %in% c("mm9", "mm10")){
    rcmat<- rcmat[Chromosome %in% c(1:19, 'X')]
  } else {
    write(sprintf('Invalid genome build: %s', gbuild), stderr())
    quit(status= 1)
  }
  return(list(pileup= rcmat, chr_prefix= chr_prefix))
}

run_facets<- function(
    pre_rcmat, pre_gbuild, pre_snp.nbhd, pre_het.thresh, pre_cval, pre_ndepth, pre_ndepthmax,
    proc_cval, proc_min.nhet,
    emcncf_unif, emcncf_min.nhet
){
  # Run the core functions of facets for segmentation, purity etc.
  # Here is where the actual CNV discovery happen.
  # Param prefix matches the facets function they go to.
  write(sprintf('[%s] Preprocessing sample...', Sys.time()), stderr())
  
  xx<- preProcSample(
    rcmat=      pre_rcmat,
    gbuild=     pre_gbuild, 
    snp.nbhd=   pre_snp.nbhd, 
    het.thresh= pre_het.thresh, 
    cval=       pre_cval, 
    ndepth=     pre_ndepth,    
    ndepthmax=  pre_ndepthmax
  )
  #rm(pre_rcmat)
  x_ <- gc(verbose= FALSE)
  
  write(sprintf('[%s] Processing sample...', Sys.time()), stderr())
  proc_out<- procSample(xx, 
                        cval=     proc_cval, 
                        min.nhet= proc_min.nhet)
  proc_out$jointseg<- data.table(proc_out$jointseg)
  proc_out$out<- data.table(proc_out$out)
  
  write(sprintf('[%s] Fitting model...', Sys.time()), stderr())
  emcncf_fit<- emcncf(x=        proc_out, 
                      unif=     emcncf_unif, 
                      min.nhet= emcncf_min.nhet)
  names(emcncf_fit$purity)<- NULL
  names(emcncf_fit$ploidy)<- NULL
  emcncf_fit[['cncf']]<- data.table(emcncf_fit[['cncf']])[order(chrom, start)]
  return(list(proc_out= proc_out, emcncf_fit= emcncf_fit))
}


classify_cnv<- function(dat){
  # Classify CNV. See also https://github.com/mskcc/facets/issues/62
  # dat is a data.table modifed in-place
  dat[, type := NA]
  dat[, type := ifelse((tcn.em == 2 & (lcn.em == 1 | is.na(lcn.em))), 'NEUTR', type)]
  dat[, type := ifelse(is.na(type) & tcn.em == 2 & lcn.em == 2, 'DUP', type)]
  dat[, type := ifelse(is.na(type) & tcn.em == 0, 'DEL', type)]
  dat[, type := ifelse(is.na(type) & tcn.em > 2 & (lcn.em > 0 | is.na(lcn.em)), 'DUP', type)]
  dat[, type := ifelse(is.na(type) & tcn.em == 1, 'HEMIZYG', type)]
  dat[, type := ifelse(is.na(type) & tcn.em == 2 & lcn.em == 0, 'LOH', type)]
  dat[, type := ifelse(is.na(type) & tcn.em > 2 & lcn.em == 0, 'DUP-LOH', type)]
  stopifnot(all(!is.na(dat$type))) # Everything has been classified
  
  return(dat)
}


## arguments ##
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="SNP file name, full path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output folder/name", metavar="character"),
  make_option(c("-t", "--targets"), type="character", default=NULL, 
              help="targets file, bed format, full path", metavar="character"),
  make_option(c("-m", "--min_reads"), type="numeric", default=20, 
              help="min reads", metavar="character"),
  make_option(c("-b", "--gbuild"), type="character", default="hg38", 
              help="genome build, chose between hg19 or hg38", metavar="character"),
  make_option(c("-c", "--cvalue_1"), type="numeric", default= 25,  # WES
              help="pre_cval", metavar="character"),
  make_option(c("-C", "--cvalue_2"), type="numeric", default= 150,  # WES
              help="proc_cval", metavar="character")
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$input)
print(paste0("Reading file: ", opt$input))
rcmat = readSnpMatrix2(opt$input,"hg38")

# reorder, numerically
rcmat$pileup[rcmat$pileup$Chromosome == "X",]$Chromosome = 23
rcmat$pileup[rcmat$pileup$Chromosome == "Y",]$Chromosome = 24
rcmat$pileup = rcmat$pileup[order(as.numeric(rcmat$pileup$Chromosome)),]
rcmat$pileup[rcmat$pileup$Chromosome == 23,]$Chromosome = "X"
rcmat$pileup[rcmat$pileup$Chromosome == 24,]$Chromosome = "Y"

# filter with target file
targets<- data.table(read.table(opt$targets, comment.char= '#', header= FALSE, sep= '\t'))

rcmat_flt<- filter_rcmat(rcmat= rcmat[['pileup']], min_ndepth= opt$min_reads, 
                         max_ndepth= 10000000, target_bed= targets)
nbhd_snp <- 250 # default
# run facets
facets<- run_facets(
  pre_rcmat=       rcmat_flt,
  pre_gbuild=      opt$gbuild, 
  pre_snp.nbhd=    nbhd_snp, 
  pre_het.thresh=  0.25, 
  pre_cval=        opt$cvalue_1, 
  pre_ndepth=      25,
  pre_ndepthmax=   1e8, # options
  proc_cval=       opt$cvalue_2, 
  proc_min.nhet=   15, 
  emcncf_unif=     FALSE, 
  emcncf_min.nhet= 15
) 



write(sprintf('[%s] Writing output', Sys.time()), stderr())
cncf<- copy(facets$emcncf_fit$cncf)

cncf = as.data.frame(cncf)
cncf= cncf[order(as.numeric(cncf$chrom)),]
cncf[cncf$chrom == 23,]$chrom = "X"
#cncf[cncf$chrom == 24,]$chrom = "Y"

cncf = as.data.table(cncf)
cncf[,"Purity" := facets$emcncf_fit$purity]
cncf[,"Ploidy" := facets$emcncf_fit$ploidy]
cncf[,"dipLogR" := facets$emcncf_fit$dipLogR]
cncf[,"emflags" := facets$emcncf_fit$emflags]

classify_cnv(cncf)

fwrite(cncf, opt$output)
print(paste(opt$output, "written"))

print("Drawing diagnostic plots")

sname<- sprintf('%s; ploidy= %.2f; purity= %.2f', basename(opt$output), facets$emcncf_fit$ploidy, facets$emcncf_fit$purity)
png(paste0(opt$output, opt$cvalue_1, 15, "_fitted.png"), width = 3.25, height = 3.25, units = "in", res=1200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i",cex.axis = 2,cex.lab  = 2)
print("Plotting Sample")
plotSample(x= facets$proc_out, emfit= facets$emcncf_fit, sname= sname)
while(!is.null(dev.list())){dev.off()}

print("Ploting Spider")
png(paste0(opt$output, opt$cvalue_1, 15, "_diagnostic.png"), width = 3.25, height = 3.25, units = "in", res=1200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i",cex.axis = 2,cex.lab  = 2)
logRlogORspider(facets$proc_out$out, facets$proc_out$dipLogR)
while(!is.null(dev.list())){dev.off()}

si<- capture.output(sessionInfo())
write(si, stderr())




