//===================================================================================
//Common igd struct, parameters, functions
//by Jianglin Feng  05/12/2018
//database intervals sorted by _start: 8/12/2019
//-----------------------------------------------------------------------------------
#include "igd_base.h"

#define gdata_t_key(r) ((r).start)
KRADIX_SORT_INIT(intv, gdata_t, gdata_t_key, 4)
KHASH_MAP_INIT_STR(str, int32_t)
typedef khash_t(str) strhash_t;

void str_splits( char* str, int *nmax, char **splits)
{   //tsv
    splits[*nmax] = NULL;
    splits[0] = str;
    char *ch = str;
    int ns = 1;
    do {
        if (*ch == '\t'){
            splits[ns++] = &ch[1];
            *ch = '\0';
        }
        ch++;
    } while (*ch != '\0' && ns < *nmax+1);
    *nmax = ns;
}

char *parse_bed(char *s, int32_t *st_, int32_t *en_)
{
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1;
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 0) ctg = p;
			else if (i == 1) st = atol(p);
			else if (i == 2) en = atol(p);
			++i, p = q + 1;
			if (c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 3? ctg : 0;
}

int32_t bSearch(gdata_t *gdata, int32_t t0, int32_t tc, int32_t qe)
{   //find tE: index of the last item satisfying .start < qe from right
	//assuming gdata sorted by start
    int32_t tL=t0, tR=tc, tM, tE = -1;
    if(gdata[tR].start < qe)
    	return tR;
    else if(gdata[tL].start >= qe)
    	return -1;
    while(tL<tR-1){
        tM = (tL+tR)/2;
        if(gdata[tM].start >= qe)
            tR = tM-1;
        else
            tL = tM;
    }
    if(gdata[tR].start < qe)
        tE = tR;
    else if(gdata[tL].start < qe)
        tE = tL;
  	return tE;
}

/**
 * @brief Primary function for adding regions to an IGD object
 *
 * Adds a region to an IGD database?
 *
 * @param *igd a pointer to the IGD object (of type igd_t)
 * @param *chrm Chromosome ???
 * @param s Start coordinate
 * @param e End coordinate
 * @param v Value ???
 * @param idx ???
 * @return Null
 */
void igd_add(igd_t *igd, const char *chrm, int32_t s, int32_t e, int32_t v, int32_t idx)
{	//layers: igd->ctg->gTile->gdata(list)
	if(s >= e)return;
	int absent;  //return value from kh_put function
	khint_t k;
	strhash_t *h = (strhash_t*)igd->hc;
	k = kh_put(str, h, chrm, &absent);
	int32_t n1 = s/igd->nbp;
	int32_t n2 = (e-1)/igd->nbp;
	if (absent) {
		//printf("%s %i %i %i\n", chrm, n1, n2, k);
		//igd
		if (igd->nctg == igd->mctg)
			EXPAND(igd->ctg, igd->mctg);
		kh_val(h, k) = igd->nctg;
		//ctg: initialize
		ctg_t *p = &igd->ctg[igd->nctg++];
		p->name = strdup(chrm);
		p->mTiles= 1 + n2;
		p->gTile = malloc(p->mTiles*sizeof(tile_t));
		kh_key(h, k) = p->name;
		//tile: initialize
		for(int i=0;i<p->mTiles;i++){
			tile_t *tile = &p->gTile[i];
			tile->ncnts = 0;	//each batch
			tile->nCnts = 0;	//total
			tile->mcnts = 4;
			tile->gList = malloc(tile->mcnts*sizeof(gdata_t));
		}
	}
	int32_t kk = kh_val(h, k);
	ctg_t *p = &igd->ctg[kk];
	if (n2+1>=p->mTiles){
		int32_t tt = p->mTiles;
		p->mTiles = n2+1;
	    p->gTile = realloc(p->gTile, p->mTiles*sizeof(tile_t));
	    //initialize new tiles
		for(int i=tt;i<p->mTiles;i++){
			tile_t *tile = &p->gTile[i];
			tile->ncnts = 0;	//each batch
			tile->nCnts = 0;	//total
			tile->mcnts = 16;
			tile->gList = malloc(tile->mcnts*sizeof(gdata_t));
		}
	}
	//add data elements
	for(int i=n1;i<=n2;i++){
		tile_t *tile = &p->gTile[i];
		if(tile->ncnts == tile->mcnts)
			EXPAND(tile->gList, tile->mcnts);
		gdata_t *gdata = &tile->gList[tile->ncnts++];
		gdata->start = s;
		gdata->end   = e;
		gdata->value = v;
		gdata->idx   = idx;
		igd->total++;
	}
	return;
}

info_t* get_fileinfo(char *ifName, int32_t *nFiles)
{   //read head file __index.tsv to get info
    FILE *fp = fopen(ifName, "r");
    if(fp==NULL){
        printf("file not found:%s\n", ifName);
        return NULL;
    }
    char buf[1024], *s0, *s1, *s2, *s3;
    int nfiles=0;
    char* rtn = fgets(buf, 1024, fp);
    while(fgets(buf, 1024, fp)!=NULL)
		nfiles++;

    info_t *fi = (info_t*)malloc(nfiles*sizeof(info_t));
    fseek(fp, 0, SEEK_SET);
    int i=0;
    rtn = fgets(buf, 1024, fp);   //header
    while(fgets(buf, 1024, fp)!=NULL){
        s0 = strtok(buf, "\t");
        s1 = strtok(NULL, "\t");
        fi[i].fileName = strdup(s1);
        s2 = strtok(NULL, "\t");
        fi[i].nr = atol(s2);
        //s3 = strtok(NULL, "\t");
        //fi[i].md = (double)atol(s3);
        i++;
    }
    *nFiles = (int32_t)nfiles;
    fclose(fp);
    return fi;
}

IGD_t* open_IGD(char *igdFile)
{
	IGD_t* IGD = IGD_init();
    char tmp[256];
    strcpy(tmp, igdFile);
    tmp[strrchr(tmp, '.')-tmp] = '\0';
    strcpy(IGD->fname, tmp);
    char *idFile = tmp;					//str_split(tmp, '.', &nCols)[0];
    strcat(idFile, "_index.tsv");
    IGD->finfo = get_fileinfo(idFile, &IGD->nFiles);
    FILE *fp = fopen(igdFile, "rb");
    if(fp == NULL)
        printf("Can't open file %s", igdFile);
    long rtn = fread(&IGD->nbp, sizeof(int32_t), 1, fp);
    rtn = fread(&IGD->gType, sizeof(int32_t), 1, fp);
    rtn = fread(&IGD->nCtg, sizeof(int32_t), 1, fp);
   	int i, k;
   	int32_t gdsize;
   	gdsize = sizeof(gdata_t);
    int32_t tileS, m = IGD->nCtg;		//the idx of a tile in the chrom
    //------------------------------------------
    IGD->nTile = malloc(m*sizeof(int32_t));
    rtn = fread(IGD->nTile, sizeof(int32_t)*m, 1, fp);
    int64_t chr_loc = 12 + 44*m;		//header size in bytes
    for(i=0;i<m;i++) chr_loc += IGD->nTile[i]*4;
    //------------------------------------------
    IGD->nCnt = malloc(m*sizeof(int32_t*));
    IGD->tIdx = malloc(m*sizeof(int64_t*));
    for(i=0;i<m;i++){
    	k = IGD->nTile[i];
    	IGD->nCnt[i] = calloc(k, sizeof(int32_t));
    	rtn = fread(IGD->nCnt[i], sizeof(int32_t)*k, 1, fp);
    	//--------------------------------------
    	IGD->tIdx[i] = calloc(k, sizeof(int64_t));
    	IGD->tIdx[i][0] = chr_loc;
    	for(int j=1; j<k; j++)
    		IGD->tIdx[i][j] = IGD->tIdx[i][j-1]+IGD->nCnt[i][j-1]*gdsize;
    	chr_loc = IGD->tIdx[i][k-1]+IGD->nCnt[i][k-1]*gdsize;
    }

	IGD->cName = malloc(m*sizeof(char*));
    for(i=0;i<m;i++){
		IGD->cName[i] = malloc(40*sizeof(char));
		  rtn = fread(IGD->cName[i], 40, 1, fp);
    }
    IGD->fP = fp;

    //setup hc
	IGD->hc = kh_init(str);
	int absent;
	for(i=0;i<IGD->nCtg;i++){
		khint_t k;
		strhash_t *h = (strhash_t*)IGD->hc;
		k = kh_put(str, h, IGD->cName[i], &absent);
		kh_val(h, k) = i;
		kh_key(h, k) = IGD->cName[i];
	}
	IGD->gData = malloc(1*sizeof(gdata_t));
	IGD->preIdx = -1;
	IGD->preChr = -1;
    return IGD;
}

int32_t get_id(IGD_t *IGD, const char *chrm)
{	//for search
	khint_t k;
	strhash_t *h = (strhash_t*)IGD->hc;
	k = kh_get(str, h, chrm);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

int32_t get_nFiles(IGD_t *IGD)
{
	return IGD->nFiles;
}

void igd_saveT(igd_t *igd, char *oPath)
{	//Save/append tiles to disc, add cnts tp Cnts
	char idFile[256];
	for (int i = 0; i < igd->nctg; i++){
		ctg_t *ctg = &igd->ctg[i];
		for(int j=0; j< ctg->mTiles; j++){
			tile_t *tile = &ctg->gTile[j];
			//---------------------------------------
			if(tile->ncnts>0){
		        sprintf(idFile, "%s%s%s_%i", oPath, "data0/", ctg->name, j);
		        FILE *fp = fopen(idFile, "ab");
		        if(fp==NULL)
		            printf("Can't open file %s", idFile);
		        fwrite(tile->gList, sizeof(gdata_t), tile->ncnts, fp);
		        fclose(fp);
		    }
		    tile->nCnts += tile->ncnts;
			tile->ncnts = 0;
			free(tile->gList);
		    tile->mcnts = 16;//MAX(16, tile->mcnts/16);
		    tile->gList = malloc(tile->mcnts*sizeof(gdata_t));
		    //tile->gList = realloc(tile->gList, tile->mcnts*sizeof(gdata_t));?
		}
	}
	igd->total = 0;	//batch total
}

void igd_save(igd_t *igd, char *oPath, char *igdName)
{
	char idFile[256], iname[256];
	//1. Save IGD data info: ctg string length 40
    int32_t i, j, n, m  = igd->nctg;
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");
    FILE *fp = fopen(idFile, "wb");
    if(fp==NULL)
        printf("Can't open file %s", idFile);
	fwrite(&igd->nbp, sizeof(int32_t), 1, fp); 		//4 bytes
	fwrite(&igd->gType, sizeof(int32_t), 1, fp); 	//4
	fwrite(&m, sizeof(int32_t), 1, fp); 			//4
	//-----------------
	for(i=0;i<m;i++)								//m*4
		fwrite(&igd->ctg[i].mTiles, sizeof(int32_t), 1, fp);
	for(i=0;i<m;i++){								//sum(mTiles)
		ctg_t *p = &igd->ctg[i];
		n = p->mTiles;
		for(j=0;j<n;j++)
			fwrite(&p->gTile[j].nCnts, sizeof(int32_t), 1, fp);
	}
	//write string array
	for(i=0;i<m;i++)								//m*40
		fwrite(igd->ctg[i].name, 40, 1, fp);

	//2. Sort and save tiles data
	for(i=0;i<m;i++){
		ctg_t *p = &igd->ctg[i];
		n = p->mTiles;
		for(j=0;j<n;j++){
			tile_t *q = &p->gTile[j];
			int32_t nrec = q->nCnts, gdsize;
		    if(nrec>0){
		    	sprintf(iname, "%s%s%s_%i", oPath, "data0/", p->name, j);
				FILE *fp0 = fopen(iname, "rb");
				if(fp0 == NULL)
					printf("Can't open file %s", iname);
	    		gdsize = nrec*sizeof(gdata_t);
			    gdata_t *gdata = malloc(gdsize);
			    long rtn = fread(gdata, gdsize, 1, fp0);
			    fclose(fp0);
			    radix_sort_intv(gdata, gdata+nrec);
			    fwrite(gdata, gdsize, 1, fp);
			    free(gdata);
		        remove(iname);
		    }
		}
    }
    fclose(fp);
}

igd_t *igd_init(int tile_size)
{
	igd_t *igd = malloc(1*sizeof(igd_t));
	igd->gType = 1;
	igd->nbp = tile_size;
	igd->hc = kh_init(str);
	igd->nctg = 0;
	igd->mctg = 32;
	igd->ctg = malloc(igd->mctg*sizeof(ctg_t));
	igd->total = 0;
	return igd;
}

void igd_destroy(igd_t *igd)
{
	if (igd == 0) return;
	for (int i = 0; i < igd->nctg; ++i){
		free(igd->ctg[i].name);
		for(int j=0; j< igd->ctg[i].mTiles; j++)
			free(igd->ctg[i].gTile[j].gList);
	}
	free(igd->ctg);
	kh_destroy(str, (strhash_t*)igd->hc);
	free(igd);
}

SearchTask_t *SearchTask_init()
{
    printf("Creating SearchTask\n"); 
    SearchTask_t *sTask = malloc(sizeof(SearchTask_t));
    sTask->datamode = 0;
    sTask->checking = 0;
    sTask->stat2 = 1;
    sTask->status = INIT;
    printf("Finished initializing SearchTask\n"); 
    return sTask;
}

SearchTask_t SearchTask_init2()
{
    printf("Creating SearchTask\n"); 
    SearchTask_t sTask;// = (SearchTask_t *) malloc(1*sizeof(SearchTask_t));
    sTask.datamode = 0;
    sTask.checking = 0;
    sTask.status = FAILED;
    sTask.stat2 = 0;
    printf("Creating SearchTask\n"); 
    return sTask;
}

SearchTask_t *SearchTask_init_old()
{
    printf("Creating SearchTask\n"); 
    SearchTask_t *sTask = (SearchTask_t *) malloc(1*sizeof(SearchTask_t));
    sTask->datamode = 0;
    sTask->checking = 0;
    sTask->status = FAILED;
    sTask->stat2 = 0;
    printf("Creating SearchTask\n"); 
    return sTask;
}


CreateTask_t *CreateTask_init()
{
    CreateTask_t *cTask = (CreateTask_t *) malloc(1*sizeof(CreateTask_t));
    return cTask;
}


IGD_t *IGD_init()
{
    IGD_t *IGD = (IGD_t *) malloc(1*sizeof(IGD_t));
    IGD->nbp = 16384;
    IGD->gType = 1;
    IGD->nCtg = 24;
    return IGD;
}

void close_IGD(IGD_t *IGD)
{
	if(IGD==0) return;
	fclose(IGD->fP);
	free(IGD->gData);
    free(IGD->nTile);
    kh_destroy(str, (strhash_t*)IGD->hc);
    for(int i=0;i<IGD->nCtg;i++){
     	free(IGD->nCnt[i]);
     	free(IGD->tIdx[i]);
    }
    free(IGD->nCnt);
    free(IGD->tIdx);
    free(IGD->cName);
    free(IGD->finfo);
    free(IGD);
}
