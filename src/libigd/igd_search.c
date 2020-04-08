//=====================================================================================
//Read igd region data and query data, and then find all overlaps
//by Jianglin Feng  05/12/2018
//updated bvy Nathan Sheffield 04/2020
//
//database intervals sorted by _start: 8/12/2019
//Reorg--simplify: 11/06/19
//-------------------------------------------------------------------------------------
#include "igd_search.h"
//-------------------------------------------------------------------------------------

void get_overlaps(IGD_t *IGD, char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{
	int ichr = get_id(IGD, chrm);
	if(ichr<0)
		return;
	int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
	if(n1>mTile)
		return;
	n2 = MIN(n2, mTile);
	tmpi = IGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	long rtn;
	if(tmpi>0){
		if(n1!=IGD->preIdx || ichr!=IGD->preChr){
			fseek(IGD->fP, IGD->tIdx[ichr][n1], SEEK_SET);
			free(IGD->gData);
			IGD->gData = malloc(tmpi*sizeof(gdata_t));
			rtn = fread(IGD->gData, sizeof(gdata_t)*tmpi, 1, IGD->fP);
			IGD->preIdx = n1;
			IGD->preChr = ichr;
		}
		if(qe>IGD->gData[0].start){						//sorted by start
			//find the 1st rs < qe
			tL = 0, tR=tmpi1;
			while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
				tM = (tL+tR)/2;
				if(IGD->gData[tM].start < qe)	//right side:
				    tL = tM;
				else
				    tR = tM;				//left side
			}
			if(IGD->gData[tR].start<qe)tL = tR;
			//-----------------------------------------
			for(i=tL; i>=0; i--){
				if(IGD->gData[i].end>qs){
					hits[IGD->gData[i].idx]++;
				}
			}
		}
		if(n2>n1){									//n2>n1
			int32_t bd = IGD->nbp*(n1+1);			//only keep the first
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
				tmpi = IGD->nCnt[ichr][j];
				tmpi1 = tmpi-1;
				if(tmpi>0){
					if(j!=IGD->preIdx || ichr!=IGD->preChr){
						fseek(IGD->fP, IGD->tIdx[ichr][j], SEEK_SET);
						free(IGD->gData);
						IGD->gData = malloc(tmpi*sizeof(gdata_t));
						rtn = fread(IGD->gData, sizeof(gdata_t)*tmpi, 1, IGD->fP);
						IGD->preIdx = j;
						IGD->preChr = ichr;
					}
					if(qe>IGD->gData[0].start){
						tS = 0;
						while(tS<tmpi && IGD->gData[tS].start<bd)tS++;	//qs<bd
						tL = 0, tR=tmpi1;
						while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
							tM = (tL+tR)/2;
							if(IGD->gData[tM].start < qe)	//right side:
								tL = tM;
							else
								tR = tM;				//left side
						}
						if(IGD->gData[tR].start<qe)tL = tR;
						//-----------------------------------------
						for(i=tL; i>=tS; i--){
							if(IGD->gData[i].end>qs){
								hits[IGD->gData[i].idx]++;
							}
						}
					}
				}
				bd+=IGD->nbp;
			}
		}
	}
}

void get_overlaps32(IGD_t *IGD, char *chrm, int32_t qs, int32_t qe, int32_t *hits)
{
  //printf("%i\t%i\t%s\n", qs, qe, chrm);
  //printf("%i\t%i\t%i\t%i\n", hits[0], hits[1], hits[2], hits[3]);
  int ichr = get_id(IGD, chrm);
  if(ichr<0)
    return;
  int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
  int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
  if(n1>mTile)
    return;
  n2 = MIN(n2, mTile);
  tmpi = IGD->nCnt[ichr][n1];
  tmpi1 = tmpi-1;
  long rtn;
  if(tmpi>0){
    if(n1!=IGD->preIdx || ichr!=IGD->preChr){
      fseek(IGD->fP, IGD->tIdx[ichr][n1], SEEK_SET);
      free(IGD->gData);
      IGD->gData = malloc(tmpi*sizeof(gdata_t));
      rtn = fread(IGD->gData, sizeof(gdata_t)*tmpi, 1, IGD->fP);
      IGD->preIdx = n1;
      IGD->preChr = ichr;
    }
    if(qe>IGD->gData[0].start){						//sorted by start
      //find the 1st rs < qe
      tL = 0, tR=tmpi1;
      while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
        tM = (tL+tR)/2;
        if(IGD->gData[tM].start < qe)	//right side:
          tL = tM;
        else
          tR = tM;				//left side
      }
      if(IGD->gData[tR].start<qe)tL = tR;
      //-----------------------------------------
      for(i=tL; i>=0; i--){
        if(IGD->gData[i].end>qs){
          hits[IGD->gData[i].idx]++;
        }
      }
    }
    if(n2>n1){									//n2>n1
      int32_t bd = IGD->nbp*(n1+1);			//only keep the first
      for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
        tmpi = IGD->nCnt[ichr][j];
        tmpi1 = tmpi-1;
        if(tmpi>0){
          if(j!=IGD->preIdx || ichr!=IGD->preChr){
            fseek(IGD->fP, IGD->tIdx[ichr][j], SEEK_SET);
            free(IGD->gData);
            IGD->gData = malloc(tmpi*sizeof(gdata_t));
            rtn = fread(IGD->gData, sizeof(gdata_t)*tmpi, 1, IGD->fP);
            IGD->preIdx = j;
            IGD->preChr = ichr;
          }
          if(qe>IGD->gData[0].start){
            tS = 0;
            while(tS<tmpi && IGD->gData[tS].start<bd)tS++;	//qs<bd
            tL = 0, tR=tmpi1;
            while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
              tM = (tL+tR)/2;
              if(IGD->gData[tM].start < qe)	//right side:
                tL = tM;
              else
                tR = tM;				//left side
            }
            if(IGD->gData[tR].start<qe)tL = tR;
            //-----------------------------------------
            for(i=tL; i>=tS; i--){
              if(IGD->gData[i].end>qs){
                hits[IGD->gData[i].idx]++;
              }
            }
          }
        }
        bd+=IGD->nbp;
      }
    }
  }
  //printf("%i\t%i\t%i\t%i\n", hits[0], hits[1], hits[2], hits[3]);
}

int32_t* getOverlapsFile(IGD_t* IGD, char *qFile, int32_t *hits) {
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	if ((fp = gzopen(qFile, "r")) == 0)
		return 0;
	ks = ks_init(fp);
    char *chrm;
	int32_t st, en, nl=0;
	uint64_t ols = 0;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		chrm = parse_bed(str.s, &st, &en);
		if (chrm) {
			get_overlaps32(IGD, chrm, st, en, hits);
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return hits;
}

int32_t* getOverlapsVectors(IGD_t* IGD, int qLength, char **qChrs,
    					    int32_t *qStarts, int32_t *qEnds, int32_t *hits) {
	for (int i=0; i<qLength; i++){
		get_overlaps32(IGD, qChrs[i], qStarts[i], qEnds[i], hits);
	}
	return hits;
}



// TODO: delete this function
void search_1(char **igdFile, char **qchr, int32_t *qs, int32_t *qe, int64_t *hits)
{
  IGD_t *IGD = open_IGD(*igdFile);
  get_overlaps(IGD, *qchr, *qs, *qe, hits);
  close_IGD(IGD);
}


// TODO: delete this function as it has been decoupled
void getOverlaps_old(char **igdFile, char **qFile, int64_t *hits)
{
	IGD_t *IGD = open_IGD(*igdFile);
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
			get_overlaps(IGD, chrm, st, en, hits);
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	close_IGD(IGD);
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
