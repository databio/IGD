#==========================iGD in R========================================================
# Get igd data source info
# Dec, 2019: Jianglin Feng
#-----------------------------------------------------------------------------------------
FInfo <- function(...)
  UseMethod("FInfo")

FInfo.default <- function (nfiles, finfo)
{
  FInfo <- list(nFiles = nfiles, fInfo = finfo)
  class(FInfo) <- c("FInfo","list")
  return(FInfo)
}

is.FInfo <- function(x)
{
  rtn <- any(attr(x,which="class") == "FInfo")
  return(rtn)
}

getFInfo <- function(igd_file) {
  #assume igdFile: .../xxx.igd
  tsv_file = paste(substr(igd_file, 1, nchar(igd_file)-4),"_index.tsv", sep="")
  if(!file.exists(tsv_file))
    stop("File '", igd_file, "' is not found. ")
  finfo <- read.csv(file=tsv_file, head=TRUE, sep="\t") #Index(0-based), File, Number of regions, Avg size
  nfiles <- length(finfo[,1])
  return(FInfo(nfiles, finfo))
}
