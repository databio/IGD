//===================================================================================
//Common igd struct, parameters, functions
//by Jianglin Feng  05/12/2018
//database intervals sorted by _start: 8/12/2019
//-----------------------------------------------------------------------------------
#include "igd_base.h"
#define gdata_t_key(r) ((r).start)
#define gdata0_t_key(r) ((r).start)
KRADIX_SORT_INIT(intv, gdata_t, gdata_t_key, 4)
KRADIX_SORT_INIT(intv0, gdata0_t, gdata0_t_key, 4)
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

int32_t bSearch0(gdata0_t *gdata, int32_t t0, int32_t tc, int32_t qe)
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

void igd_add(igd_t *igd, const char *chrm, int32_t s, int32_t e, int32_t v, int32_t idx)
{	//layers: igd->ctg->gTile->gdata(list)
	if(s >= e)return;
	int absent;
	khint_t k;
	strhash_t *h = (strhash_t*)hc;
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

void igd0_add(igd0_t *igd, const char *chrm, int32_t s, int32_t e, int32_t idx)
{	//layers: igd->ctg->gTile->gdata(list)
	if(s >= e)return;
	int absent;
	khint_t k;
	strhash_t *h = (strhash_t*)hc;
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
		ctg0_t *p = &igd->ctg[igd->nctg++];		
		p->name = strdup(chrm);
		p->mTiles= 1 + n2;
		p->gTile = malloc(p->mTiles*sizeof(tile_t));		
		kh_key(h, k) = p->name;
		//tile: initialize
		for(int i=0;i<p->mTiles;i++){
			tile0_t *tile = &p->gTile[i];
			tile->ncnts = 0;	//each batch 
			tile->nCnts = 0;	//total
			tile->mcnts = 4;	
			tile->gList = malloc(tile->mcnts*sizeof(gdata0_t));
		}	
	}
	int32_t kk = kh_val(h, k);
	ctg0_t *p = &igd->ctg[kk];
	if (n2+1>=p->mTiles){
		int32_t tt = p->mTiles;
		p->mTiles = n2+1;
	    p->gTile = realloc(p->gTile, p->mTiles*sizeof(tile0_t));
	    //initialize new tiles
		for(int i=tt;i<p->mTiles;i++){
			tile0_t *tile = &p->gTile[i];
			tile->ncnts = 0;	//each batch 
			tile->nCnts = 0;	//total
			tile->mcnts = 16;	
			tile->gList = malloc(tile->mcnts*sizeof(gdata0_t));
		}
	}
	//add data elements
	for(int i=n1;i<=n2;i++){
		tile0_t *tile = &p->gTile[i];
		if(tile->ncnts == tile->mcnts)
			EXPAND(tile->gList, tile->mcnts);		
		gdata0_t *gdata = &tile->gList[tile->ncnts++];	
		gdata->start = s;
		gdata->end   = e;
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
    fgets(buf, 1024, fp);//head line
    while(fgets(buf, 1024, fp)!=NULL)
		nfiles++;

    info_t *fi = (info_t*)malloc(nfiles*sizeof(info_t));
    fseek(fp, 0, SEEK_SET);
    int i=0;    
    fgets(buf, 1024, fp);   //header
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

iGD_t *get_igdinfo(char *igdFile)
{
    FILE *fp = fopen(igdFile, "rb");
    if(fp == NULL)
        printf("Can't open file %s", igdFile);
    iGD_t *iGD = malloc(1*sizeof(iGD_t));
    fread(&iGD->nbp, sizeof(int32_t), 1, fp);
    fread(&iGD->gType, sizeof(int32_t), 1, fp);  
    fread(&iGD->nCtg, sizeof(int32_t), 1, fp);    
   	int i, k;
   	int32_t gdsize;
   	if(iGD->gType==0)
   		gdsize = sizeof(gdata0_t);
   	else
   		gdsize = sizeof(gdata_t);       
    int32_t tileS, m = iGD->nCtg;	//the idx of a tile in the chrom 
    //------------------------------------------
    iGD->nTile = malloc(m*sizeof(int32_t));        
    fread(iGD->nTile, sizeof(int32_t)*m, 1, fp);    
    int64_t chr_loc = 12 + 44*m;	//header size in bytes
    for(i=0;i<m;i++) chr_loc += iGD->nTile[i]*4;
    //------------------------------------------
    iGD->nCnt = malloc(m*sizeof(int32_t*));
    iGD->tIdx = malloc(m*sizeof(int64_t*)); 
    for(i=0;i<m;i++){
    	k = iGD->nTile[i];    	
    	iGD->nCnt[i] = calloc(k, sizeof(int32_t));
    	fread(iGD->nCnt[i], sizeof(int32_t)*k, 1, fp);
    	//--------------------------------------     	
    	iGD->tIdx[i] = calloc(k, sizeof(int64_t)); 
    	iGD->tIdx[i][0] = chr_loc;
    	for(int j=1; j<k; j++)
    		iGD->tIdx[i][j] = iGD->tIdx[i][j-1]+iGD->nCnt[i][j-1]*gdsize;
    	chr_loc = iGD->tIdx[i][k-1]+iGD->nCnt[i][k-1]*gdsize;
    }

	iGD->cName = malloc(m*sizeof(char*));     
    for(i=0;i<m;i++){
		iGD->cName[i] = malloc(40*sizeof(char));   	
		fread(iGD->cName[i], 40, 1, fp);     
    }   
    fclose(fp);
   
    //setup hc
	hc = kh_init(str); 
	int absent;
	for(i=0;i<iGD->nCtg;i++){ 
		khint_t k;
		strhash_t *h = (strhash_t*)hc;
		k = kh_put(str, h, iGD->cName[i], &absent);							
		kh_val(h, k) = i;			
		kh_key(h, k) = iGD->cName[i];
	}    
    return iGD;	
}

int32_t get_id(const char *chrm)
{
	khint_t k;
	strhash_t *h = (strhash_t*)hc;
	k = kh_get(str, h, chrm);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

void igd_saveT(igd_t *igd, char *oPath)
{	//Temporarily Save/append tiles to disc, add cnts tp Cnts; reset tile->gList 
	char idFile[128];
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
		        free(tile->gList);
		    }			
		    tile->nCnts += tile->ncnts;
			tile->ncnts = 0;
		    tile->mcnts = 16;//initialize: MAX(16, tile->mcnts/16);
		    tile->gList = malloc(tile->mcnts*sizeof(gdata_t));
		    //tile->gList = realloc(tile->gList, tile->mcnts*sizeof(gdata_t));?		    
		}
	}
	igd->total = 0;	//batch total
}

void igd0_saveT(igd0_t *igd, char *oPath)
{	//Save/append tiles to disc, add cnts tp Cnts 
	char idFile[128];
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
		        fwrite(tile->gList, sizeof(gdata0_t), tile->ncnts, fp);
		        fclose(fp); 			
		        free(tile->gList);
		    }			
		    tile->nCnts += tile->ncnts;
			tile->ncnts = 0;
		    tile->mcnts = 16;//MAX(16, tile->mcnts/16);
		    tile->gList = malloc(tile->mcnts*sizeof(gdata0_t));
		    //tile->gList = realloc(tile->gList, tile->mcnts*sizeof(gdata_t));?		    
		}
	}	
	igd->total = 0;	//batch total
}

void igd_save(igd_t *igd, char *oPath, char *igdName)
{
	char idFile[128], iname[128];
	//1. Save iGD data info: ctg string length 40	
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
			    fread(gdata, gdsize, 1, fp0);
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

void igd0_save(igd0_t *igd, char *oPath, char *igdName)
{
	char idFile[128], iname[128];
	//1. Save iGD data info: ctg string length 40	
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
		ctg0_t *p = &igd->ctg[i];    
		n = p->mTiles;         
		for(j=0;j<n;j++){
			tile0_t *q = &p->gTile[j];	
			int32_t nrec = q->nCnts, gdsize;		
		    if(nrec>0){				    
		    	sprintf(iname, "%s%s%s_%i", oPath, "data0/", p->name, j);
				FILE *fp0 = fopen(iname, "rb");
				if(fp0 == NULL)
					printf("Can't open file %s", iname);
	    		gdsize = nrec*sizeof(gdata0_t);
			    gdata0_t *gdata = malloc(gdsize);
			    fread(gdata, gdsize, 1, fp0);
			    fclose(fp0);
			    radix_sort_intv0(gdata, gdata+nrec); 
			    fwrite(gdata, gdsize, 1, fp);
			    free(gdata);
		        remove(iname);              
		    }
		}
    }
    fclose(fp); 	
}

igd_t *igd_init(void)
{
	igd_t *igd = malloc(1*sizeof(igd_t));
	igd->gType = 1;
	igd->nbp = tile_size;	
	hc = kh_init(str);
	igd->nctg = 0;
	igd->mctg = 32;
	igd->ctg = malloc(igd->mctg*sizeof(ctg_t));
	igd->total = 0;
	return igd;
}

igd0_t *igd0_init(void)
{
	igd0_t *igd = malloc(1*sizeof(igd0_t));
	igd->gType = 0;
	igd->nbp = tile_size;	
	hc = kh_init(str);
	igd->nctg = 0;
	igd->mctg = 32;
	igd->ctg = malloc(igd->mctg*sizeof(ctg0_t));
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
	kh_destroy(str, (strhash_t*)hc);
	free(igd);
}

void igd0_destroy(igd0_t *igd)
{
	if (igd == 0) return;
	for (int i = 0; i < igd->nctg; ++i){
		free(igd->ctg[i].name);
		for(int j=0; j< igd->ctg[i].mTiles; j++)
			free(igd->ctg[i].gTile[j].gList);			
	}	
	free(igd->ctg);
	kh_destroy(str, (strhash_t*)hc);
	free(igd);
}

