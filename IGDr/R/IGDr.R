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

get_ctgId <- function(igdr, chrm)
{
  ctgId <- .Call("get_cid", igdr@ref, as.character(chrm), PACKAGE = "IGDr")
}

get_nCtgs <- function(igdr)
{
  nCtgs <- .Call("get_nCtgs", igrd@ref, PACKAGE = "IGDr")
}

get_nFiles <- function(igdr)
{
  nFiles <- .Call("get_nFiles", igrd@ref, PACKAGE = "IGDr")
}

get_binSize <- function(igdr)
{
  binSize <- .Call("get_nbp", igrd@ref, PACKAGE = "IGDr")
}

get_igdInfo <- function(igdr)
{
  igdInfo <- .Call("get_fInfo", igdr@ref, PACKAGE = "IGDr")
}

get_bin <- function(igdr, binNumber)
{
  binCount <- .Call("get_nCnt", igdr@ref, as.integer(binNumber), PACKAGE = "IGDr")
}

get_binData <- function(igdr, binNumber)
{
  binData <- .Call("get_gData", igdr@ref, as.integer(binNumber), PACKAGE = "IGDr")
}
