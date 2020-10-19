//=====================================================================================
//Read igd region data and query data, and then find all overlaps
//by Jianglin Feng  05/12/2018
//
//database intervals sorted by _start: 8/12/2019
//Reorg--simplify: 11/06/19
//-------------------------------------------------------------------------------------
#include "igd_search.h"
//-------------------------------------------------------------------------------------
int search_help(int exit_code)
{
    printf(
"%s, v%s\n"
"usage:   %s search <igd database file> [options]\n"
"         options:\n"
"             -q <query file>\n"
"             -r <a region: chrN start end>\n"
"             -v <signal value 0-1000>\n"
"             -o <output file Name>\n"
"             -c display all intersects\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

void get_overlaps(iGD_t *iGD, char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{
	int ichr = get_id(iGD, chrm);
	if(ichr<0)
		return;
	int i, j, n1 = qs/iGD->nbp, n2 = (qe-1)/iGD->nbp;	//define boundary!
	int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = iGD->nTile[ichr]-1;
	if(n1>mTile)
		return;
	n2 = MIN(n2, mTile);
	tmpi = iGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	long rtn;
	if(tmpi>0){
		if(n1!=iGD->preIdx || ichr!=iGD->preChr){
			fseek(iGD->fP, iGD->tIdx[ichr][n1], SEEK_SET);
			free(iGD->gData);
			iGD->gData = malloc(tmpi*sizeof(gdata_t));
			rtn = fread(iGD->gData, sizeof(gdata_t)*tmpi, 1, iGD->fP);
			iGD->preIdx = n1;
			iGD->preChr = ichr;
		}
		if(qe>iGD->gData[0].start){						//sorted by start
			//find the 1st rs < qe
			tL = 0, tR=tmpi1;
			while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
				tM = (tL+tR)/2;
				if(iGD->gData[tM].start < qe)	//right side:
				    tL = tM;
				else
				    tR = tM;				//left side
			}
			if(iGD->gData[tR].start<qe)tL = tR;
			//-----------------------------------------
			for(i=tL; i>=0; i--){
				if(iGD->gData[i].end>qs){
					hits[iGD->gData[i].idx]++;
				}
			}
		}
		if(n2>n1){									//n2>n1
			int32_t bd = iGD->nbp*(n1+1);			//only keep the first
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
				tmpi = iGD->nCnt[ichr][j];
				tmpi1 = tmpi-1;
				if(tmpi>0){
					if(j!=iGD->preIdx || ichr!=iGD->preChr){
						fseek(iGD->fP, iGD->tIdx[ichr][j], SEEK_SET);
						free(iGD->gData);
						iGD->gData = malloc(tmpi*sizeof(gdata_t));
						rtn = fread(iGD->gData, sizeof(gdata_t)*tmpi, 1, iGD->fP);
						iGD->preIdx = j;
						iGD->preChr = ichr;
					}
					if(qe>iGD->gData[0].start){
						tS = 0;
						while(tS<tmpi && iGD->gData[tS].start<bd)tS++;	//qs<bd
						tL = 0, tR=tmpi1;
						while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
							tM = (tL+tR)/2;
							if(iGD->gData[tM].start < qe)	//right side:
								tL = tM;
							else
								tR = tM;				//left side
						}
						if(iGD->gData[tR].start<qe)tL = tR;
						//-----------------------------------------
						for(i=tL; i>=tS; i--){
							if(iGD->gData[i].end>qs){
								hits[iGD->gData[i].idx]++;
							}
						}
					}
				}
				bd+=iGD->nbp;
			}
		}
	}
}

void get_overlaps32(iGD_t *iGD, char *chrm, int32_t qs, int32_t qe, int32_t *hits)
{
  //printf("%i\t%i\t%s\n", qs, qe, chrm);
  //printf("%i\t%i\t%i\t%i\n", hits[0], hits[1], hits[2], hits[3]);
  int ichr = get_id(iGD, chrm);
  if(ichr<0)
    return;
  int i, j, n1 = qs/iGD->nbp, n2 = (qe-1)/iGD->nbp;	//define boundary!
  int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = iGD->nTile[ichr]-1;
  if(n1>mTile)
    return;
  n2 = MIN(n2, mTile);
  tmpi = iGD->nCnt[ichr][n1];
  tmpi1 = tmpi-1;
  long rtn;
  if(tmpi>0){
    if(n1!=iGD->preIdx || ichr!=iGD->preChr){
      fseek(iGD->fP, iGD->tIdx[ichr][n1], SEEK_SET);
      free(iGD->gData);
      iGD->gData = malloc(tmpi*sizeof(gdata_t));
      rtn = fread(iGD->gData, sizeof(gdata_t)*tmpi, 1, iGD->fP);
      iGD->preIdx = n1;
      iGD->preChr = ichr;
    }
    if(qe>iGD->gData[0].start){						//sorted by start
      //find the 1st rs < qe
      tL = 0, tR=tmpi1;
      while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
        tM = (tL+tR)/2;
        if(iGD->gData[tM].start < qe)	//right side:
          tL = tM;
        else
          tR = tM;				//left side
      }
      if(iGD->gData[tR].start<qe)tL = tR;
      //-----------------------------------------
      for(i=tL; i>=0; i--){
        if(iGD->gData[i].end>qs){
          hits[iGD->gData[i].idx]++;
        }
      }
    }
    if(n2>n1){									//n2>n1
      int32_t bd = iGD->nbp*(n1+1);			//only keep the first
      for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
        tmpi = iGD->nCnt[ichr][j];
        tmpi1 = tmpi-1;
        if(tmpi>0){
          if(j!=iGD->preIdx || ichr!=iGD->preChr){
            fseek(iGD->fP, iGD->tIdx[ichr][j], SEEK_SET);
            free(iGD->gData);
            iGD->gData = malloc(tmpi*sizeof(gdata_t));
            rtn = fread(iGD->gData, sizeof(gdata_t)*tmpi, 1, iGD->fP);
            iGD->preIdx = j;
            iGD->preChr = ichr;
          }
          if(qe>iGD->gData[0].start){
            tS = 0;
            while(tS<tmpi && iGD->gData[tS].start<bd)tS++;	//qs<bd
            tL = 0, tR=tmpi1;
            while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
              tM = (tL+tR)/2;
              if(iGD->gData[tM].start < qe)	//right side:
                tL = tM;
              else
                tR = tM;				//left side
            }
            if(iGD->gData[tR].start<qe)tL = tR;
            //-----------------------------------------
            for(i=tL; i>=tS; i--){
              if(iGD->gData[i].end>qs){
                hits[iGD->gData[i].idx]++;
              }
            }
          }
        }
        bd+=iGD->nbp;
      }
    }
  }
  //printf("%i\t%i\t%i\t%i\n", hits[0], hits[1], hits[2], hits[3]);
}


void search_1(char **igdFile, char **qchr, int32_t *qs, int32_t *qe, int64_t *hits)
{
  iGD_t *iGD = open_iGD(*igdFile);
  get_overlaps(iGD, *qchr, *qs, *qe, hits);
  close_iGD(iGD);
}

void getOverlaps(char **igdFile, char **qFile, int64_t *hits)
{
	iGD_t *iGD = open_iGD(*igdFile);
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	if ((fp = gzopen(*qFile, "r")) == 0)
		return;
	ks = ks_init(fp);
    char *chrm;
	int32_t st, en, nl=0;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		chrm = parse_bed(str.s, &st, &en);
		if (chrm) {
			get_overlaps(iGD, chrm, st, en, hits);
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	close_iGD(iGD);
}

//---Hash table for R
int32_t hash(const char *key, int32_t htsize)
{
  int32_t i=0, hash=0;
  while(key && key[i]){
    hash = (hash+key[i])%htsize;
    i++;
  }
  return hash;
}

hTable *ht_init(int32_t htsize)
{
  hTable *ht;
  if(htsize<1)
    return NULL;
  ht = malloc(sizeof(hTable));
  ht->nodes = malloc(htsize*sizeof(htNode));
  memset(ht->nodes, 0, htsize*sizeof(htNode));
  ht->size = htsize;
  return ht;
}

int ht_put(hTable *ht, const char *key, int32_t value)
{
  htNode *node = malloc(sizeof(htNode));
  node->key = strdup(key);
  node->value = value;
  int32_t i = hash(key, value);
  htNode *tmp = ht->nodes[i];//linkList[i]
  if(tmp!=NULL){
    while(tmp!=NULL){
      if(strcmp(tmp->key, node->key)==0)
        break;
      tmp = tmp->next;
    }
    if(tmp==NULL){  //already filled
      node->next = ht->nodes[i];
      ht->nodes[i] = node;
    }
    else{           //alrady exist
      tmp->value = node->value;
      free(node->key);
      free(node);
    }
  }
  else{
    node->next = NULL;
    ht->nodes[i] = node;
  }
}

int32_t ht_get(hTable *ht, const char *key)
{
  char *key1 = strdup(key);
  int32_t i = hash(key, ht->size);
  htNode *tmp = ht->nodes[i];
  while(tmp!=NULL){
    if(strcmp(tmp->key, key1)==0)
      break;
    tmp = tmp->next;
  }
  free(key1);
  if(tmp==NULL)
    return -1;
  return tmp->value;
}

void ht_free(hTable* ht)
{
  htNode *tmp=NULL;
  if(ht==NULL)return;
  for(int i=0;i<ht->size;i++){
    if(ht->nodes[i]!=NULL){
      while(ht->nodes[i]!=NULL){
        tmp = ht->nodes[i]->next;
        free(ht->nodes[i]->key);
        free(ht->nodes[i]);
        ht->nodes[i] = tmp;
      }
      free(ht->nodes[i]);
    }
  }
  free(ht->nodes);
  free(ht);
}

//---------------------------------------------------------------------------
SEXP search_1r(SEXP igdr, SEXP qchrm, SEXP qs, SEXP qe)
{ //NO need to supply output vector!!!
  iGD_t *iGD = (iGD_t *) R_ExternalPtrAddr(igdr);
  if(iGD==NULL)
    error("iGD_free: iGDr external pointer is NULL");
  const char *chrm = CHAR(STRING_ELT(qchrm, 0));
  SEXP hits;
  PROTECT(hits = allocVector(INTSXP, iGD->nFiles));//not initialized
  memset(INTEGER(hits), 0, iGD->nFiles * sizeof(int));
  get_overlaps32(iGD, chrm, INTEGER(qs)[0], INTEGER(qe)[0], INTEGER(hits));
  UNPROTECT(1);
  return(hits);
}

/*SEXP search_nr(SEXP igdr, SEXP n, SEXP qchrm, SEXP qs, SEXP qe)
{ //NO need to supply output vector!!!
  iGD_t *iGD = (iGD_t *) R_ExternalPtrAddr(igdr);
  if(iGD==NULL)
    error("iGD_free: iGDr external pointer is NULL");
  SEXP hits;
  const char *chrm;
  int32_t *tmp = calloc(iGD->nFiles, sizeof(int32_t));
  for(int i=0; i<INTEGER(n)[0]; i++){
    chrm = CHAR(STRING_ELT(qchrm, i));
    get_overlaps32(iGD, chrm, INTEGER(qs)[i], INTEGER(qe)[i], tmp);
  }
  PROTECT(hits = allocVector(INTSXP, iGD->nFiles));//not initialized
  memcpy(INTEGER(hits), tmp, iGD->nFiles * sizeof(int32_t));
  UNPROTECT(1);
  free(tmp);
  return(hits);
}*/

SEXP search_nr(SEXP igdr, SEXP n, SEXP qchrm, SEXP qs, SEXP qe)
{ //NO need to supply output vector!!!
  iGD_t *iGD = (iGD_t *) R_ExternalPtrAddr(igdr);
  if(iGD==NULL)
    error("iGD_free: iGDr external pointer is NULL");
  SEXP hits;
  const char *chrm;
  PROTECT(hits = allocVector(INTSXP, iGD->nFiles));//not initialized
  memset(INTEGER(hits), 0, iGD->nFiles * sizeof(int));
  for(int i=0; i<INTEGER(n)[0]; i++){
    chrm = CHAR(STRING_ELT(qchrm, i));
    get_overlaps32(iGD, chrm, INTEGER(qs)[i], INTEGER(qe)[i], INTEGER(hits));
  }
  UNPROTECT(1);
  return(hits);
}

SEXP get_binData(SEXP igdr, SEXP ichr, SEXP bin)
{
  iGD_t *iGD = (iGD_t *) R_ExternalPtrAddr(igdr);
  if(iGD==NULL)
    error("iGD_free: iGDr external pointer is NULL");
  int ichr0 = INTEGER(ichr)[0]-1;
  int j = INTEGER(bin)[0]-1;
  if(ichr0<0 || j>=iGD->nTile[ichr0] || j<0){
    printf("Max bin number is %i\n", iGD->nTile[ichr0]);
    return(R_NilValue);
  }
  int ncnt = iGD->nCnt[ichr0][j];
  if(ncnt<1){
    printf("No records in bin %i \n", j);
    return(R_NilValue);
  }
  SEXP starts = PROTECT(allocVector(INTSXP, ncnt));
  SEXP ends = PROTECT(allocVector(INTSXP, ncnt));
  SEXP idx = PROTECT(allocVector(INTSXP, ncnt));
  //--------------------------------------------
  gdata_t *gd = malloc(ncnt*sizeof(gdata_t));
  fseek(iGD->fP, iGD->tIdx[ichr0][j], SEEK_SET);
  long rtn = fread(gd, sizeof(gdata_t)*ncnt, 1, iGD->fP);
  //--------------------------------------------
  for(int i=0;i<ncnt;i++){
    INTEGER(idx)[i] = gd[i].idx;
    INTEGER(starts)[i] = gd[i].start;
    INTEGER(ends)[i] = gd[i].end;
  }
  free(gd);
  SEXP gdata = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(gdata, 0, idx);
  SET_VECTOR_ELT(gdata, 1, starts);
  SET_VECTOR_ELT(gdata, 2, ends);
  UNPROTECT(4);
  return(gdata);
}
