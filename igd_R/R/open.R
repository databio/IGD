#====================================iGD in R=========================================
# Open iGD database
# Dec, 2019: Jianglin Feng
#-------------------------------------------------------------------------------------
#library(RCurl) #for CFILE fp@ref-- C file pointer
openIGD <- function(igd_file) {
  #fp = CFILE(igd_file, "rb")
  finfo = getFInfo(igd_file)
  fconn = file(igd_file, "rb")
  nbp = readBin(fconn, integer(), 1)
  gtype = readBin(fconn, integer(), 1)
  nctg = readBin(fconn, integer(), 1)
  ntile = readBin(fconn, integer(), nctg)  #vector, not list
  preidx = -1
  prechr = -2
  #---------num of gdata (ncnt) for each tile
  n_t = 0
  for(i in 1:nctg){
    n_t = n_t + ntile[i];
  }
  chr_loc = integer64(1)
  chr_loc = n_t*4 + 12 + 44*nctg
  #---------read in ncnt (list of vectors) for each tile in each ctg
  #---------calc 64-bit tile position: tidx
  ncnt = list()
  tidx = list()
  for(i in 1:nctg){
    k = ntile[i]
    ncnt[i] = readBin(fconn, integer(), k)
    tidx[i] = integer64(k)
    tidx[[i]][1] = chr_loc
    for(j in 2:k){
      tidx[[i]][j] = tidx[[i]][j-1] + ncnt[[i]][j-1]*16
    }
    chr_loc = tidx[[i]][k] + ncnt[[i]][k]*16
  }
  cname = readBin(fconn, character(), nctg, 40)
  #setup hc?
  gdata = vector(mode='integer', length=4) #a single gdata_t: idx, start, end, value
  #---------construct IGD obj
  igd <- IGD(igdFile = igd_file,
              nBp = nbp, gType = gtype,
              nCtg = nctg, cName = cname,
              preIdx = pre_idx,
              preChr = pre_chr,
              nTile = ntile,
              nCnt = ncnt,
              tIdx = tidx,
              gData = gdata,
              fP = fconn,
              fInfo = finfo)
  return(igd)
}
