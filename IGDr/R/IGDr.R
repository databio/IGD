# S4 class of iGD
setClass("IGDr",
         representation(ref="externalptr")
)

setMethod("close", "IGDr",
          function(con, ...) {
           .Call("iGD_free", con@ref, PACKAGE = "IGDr")
           con
         })

#' Function to open/load an igd database for search 
#'
#' @param igd_file	the path to the igd database file
#' @return an IGDr object
#' @export
#' @examples
#' library(IGDr)
#' igdr <- IGDr("roadmap_b14.igd")
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

#' Function to search the igd database for a single query
#'
#' @param igdr	an igd database object (loaded) 
#' @param chrm the chromosome name of the query in format of chr1, chrX, chrY, chrM, ...
#' @param qs the start location of the query
#' @param qe the end location of the query
#' @return hits: number of intersections to each database source file
#' @export
#' @examples
#' library(IGDr)
#' igdr <- IGDr("roadmap_b14.igd")
#' search_1r(igdr, chrX, 1000000, 1100000)
search_1r <- function(igdr, chrm, qs, qe)
{
  hits <- .Call("search_1r", igdr@ref, as.character(chrm), as.integer(qs), as.integer(qe), PACKAGE = "IGDr")
}

#' Function to search the igd database for multiple queries 
#'
#' @param igdr	an igd database object (loaded) 
#' @param n number of queries to be searched
#' @param chrm vector of chromosome names
#' @param qs vector of the start locations of the queries
#' @param qe vector of the end locations of the queries
#' @return hits: number of intersections to each database source file
#' @export
#' @examples
#' library(IGDr)
#' igdr <- IGDr("roadmap_b14.igd")
#' search_nr(igdr, chrm, qs, qe)
search_nr <- function(igdr, n, chrm, qs, qe)
{
  hits <- .Call("search_nr", igdr@ref, as.integer(n), as.character(chrm), as.integer(qs), as.integer(qe), PACKAGE = "IGDr")
}

#' Function to search the igd database for a query set from a file 
#'
#' @param igdr	an igd database object (loaded) 
#' @param qfile path to the query file (file type of .bed or .bed.gz)
#' @return hits: number of intersections to each database source file
#' @export
#' @examples
#' search_qfile(igdr, "q1000.bed")
search_qfile <- function(igdr, qfile) { #int32 for counts
  if(!file.exists(qfile))
    stop("File '", qfile, "' is not found. ")
  qinfo <- read.csv(file=qfile, head=FALSE, sep="\t") #Index(0-based), File, Number of regions, Avg size
  nfiles <- length(qinfo[,1])
  hits <- .Call("search_nr", igdr@ref, as.integer(nfiles), as.character(qinfo[,1]), as.integer(qinfo[,2]), as.integer(qinfo[,3]), PACKAGE = "IGDr")
}

#' Function to get the contig id of a chromosome name 
#'
#' @param igdr	an igd database object (loaded) 
#' @param chrm chromosome name ("chr1", "chrX", ...)
#' @return ichr (0 if not exist)
#' @export
#' @examples
#' get_ctgId(igdr, "chrX")
get_ctgId <- function(igdr, chrm)
{
  ichr <- .Call("get_cid", igdr@ref, as.character(chrm), PACKAGE = "IGDr")
}

#' Function to get the number of contigs in an igd database 
#'
#' @param igdr	an igd database object (loaded) 
#' @return nCtgs: number of contigs
#' @export
#' @examples
#' get_ctgId(igdr)
get_nCtgs <- function(igdr)
{
  nCtgs <- .Call("get_nCtgs", igdr@ref, PACKAGE = "IGDr")
}

#' Function to get the number of source files in an igd database 
#'
#' @param igdr	an igd database object (loaded) 
#' @return nCtgs: number of source files
#' @export
#' @examples
#' get_nFiles(igdr)
get_nFiles <- function(igdr)
{
  nFiles <- .Call("get_nFiles", igdr@ref, PACKAGE = "IGDr")
}

#' Function to get the bin size of an igd database 
#'
#' @param igdr	an igd database object (loaded) 
#' @return binSize
#' @export
#' @examples
#' get_binSize(igdr)
get_binSize <- function(igdr)
{
  binSize <- .Call("get_nbp", igdr@ref, PACKAGE = "IGDr")
}

#' Function to get the number of regions in a given bin and a given contig
#'
#' @param igdr	an igd database object (loaded) 
#' @param ichr contig id
#' @param binID the bin number
#' @return binLen
#' @export
#' @examples
#' get_binLen(igdr, 10, 123)
get_binLen <- function(igdr, ichr, binID)
{
  binLen <- .Call("get_binLen", igdr@ref, as.integer(ichr), as.integer(binID), PACKAGE = "IGDr")
}

#' Function to get region data in the given bin of a given contig
#'
#' @param igdr	an igd database object (loaded) 
#' @param ichr contig id
#' @param binID the bin number
#' @return binData: vector of regions start, end, source_id
#' @export
#' @examples
#' get_binLen(igdr, 10, 123)
get_binData <- function(igdr, ichr, binID)
{
  binData <- .Call("get_binData", igdr@ref, as.integer(ichr), as.integer(binID), PACKAGE = "IGDr")
}
