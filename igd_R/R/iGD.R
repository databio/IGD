#define class IGD and create a IGD object
IGD <- function(...)
  UseMethod("IGD")

IGD.default <- function (igd_file, nfiles, nbp, gtype, nctg, cName, pre_idx, pre_chr, ntile, ncnt, tidx, gdata, fp, hc, finfo)
{
  IGD <- list(igdFile = igd_file,
              nFiles = nfiles,
              nBp = nbp, gType = gtype,
              nCtg = nctg, cName = cname,
              preIdx = pre_idx,
              preChr = pre_chr,
              nTile = ntile,
              nCnt = ncnt,
              tIdx = tIdx,
              gData = gdata,
              fP = fp,
              hC = hc,
              fInfo = finfo)
  class(IGD) <- c("IGD","list")
  return(IGD)
}

is.IGD <- function(x)
{
  rtn <- any(attr(x,which="class") == "IGD")
  return(rtn)
}
