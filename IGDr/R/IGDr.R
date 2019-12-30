# S4 class of iGD
setClass("IGDr",
         representation(ref="externalptr")
)

setMethod("close", "IGDr",
          function(con, ...) {
           .Call("iGD_free", con@ref, PACKAGE = "IGDr")
           con
         })

IGDr <- function(igd_file)
{
  #check if igd file exist
  if(!file.exists(igd_file))
     stop("File '", igd_file, "' is not found. ")
  tsv_file = paste(substr(igd_file, 1, nchar(igd_file)-4),"_index.tsv", sep="")
  if(!file.exists(tsv_file))
    stop("IGD tsv file '", tsv_file, "' not found. ")

  ans <- .Call("iGD_new", igd_file, PACKAGE = "IGDr")
}

search_1r <- function(igdr, chrm, qs, qe)
{
  hits <- .Call("search_1r", igdr@ref, as.character(chrm), as.integer(qs), as.integer(qe), PACKAGE = "IGDr")
}

search_nr <- function(igdr, n, chrm, qs, qe)
{
  hits <- .Call("search_nr", igdr@ref, as.integer(n), as.character(chrm), as.integer(qs), as.integer(qe), PACKAGE = "IGDr")
}

search_qfile <- function(igdr, qfile) { #int32 for counts
  if(!file.exists(qfile))
    stop("File '", qfile, "' is not found. ")
  qinfo <- read.csv(file=qfile, head=FALSE, sep="\t") #Index(0-based), File, Number of regions, Avg size
  nfiles <- length(qinfo[,1])
  hits <- .Call("search_nr", igdr@ref, as.integer(nfiles), as.character(qinfo[,1]), as.integer(qinfo[,2]), as.integer(qinfo[,3]), PACKAGE = "IGDr")
}

get_ctgId <- function(igdr, chrm)
{
  ctgId <- .Call("get_cid", igdr@ref, as.character(chrm), PACKAGE = "IGDr")
}

get_nCtgs <- function(igdr)
{
  nCtgs <- .Call("get_nCtgs", igdr@ref, PACKAGE = "IGDr")
}

get_nFiles <- function(igdr)
{
  nFiles <- .Call("get_nFiles", igdr@ref, PACKAGE = "IGDr")
}

get_binSize <- function(igdr)
{
  binSize <- .Call("get_nbp", igdr@ref, PACKAGE = "IGDr")
}

get_binLen <- function(igdr, ichr, binID)
{
  binLen <- .Call("get_binLen", igdr@ref, as.integer(ichr), as.integer(binID), PACKAGE = "IGDr")
}

get_binData <- function(igdr, ichr, binID)
{
  binData <- .Call("get_binData", igdr@ref, as.integer(ichr), as.integer(binID), PACKAGE = "IGDr")
}
