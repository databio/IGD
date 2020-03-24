//=====================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//
//time ./igd_search Test110000.bed /media/john/Extra/ucsc_igd/ucsc.igd
//database intervals sorted by _start: 8/12/2019
//flanking mapping: 9/18/2019--x
//multiple queries: 10/17/2019--ma
//Reorg--simplify: 11/06/19
//-------------------------------------------------------------------------------------
#include "igd_search.h"
//-------------------------------------------------------------------------------------
int search_help(int exit_code)
{
    fprintf(stderr,
"%s, v%s\n"
"usage:   %s search <igd database file> [options]\n"
"         options:\n"
"             -q <query file>\n"
"             -r <a region: chrN start end>\n"
"             -v <signal value 0-1000>\n"
"             -o <output file Name>\n"
"             -m heatmap of igd self\n"
"             -c display all intersects\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

int32_t get_overlaps0(char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{   
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;
	int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
	int32_t nols = 0;
	if(n1>mTile) 
		return 0;
	n2 = MIN(n2, mTile);	
	tmpi = IGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	if(tmpi>0){
		if(n1!=preIdx || ichr!=preChr){
			fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
			free(gData0);					
			gData0 = malloc(tmpi*sizeof(gdata0_t));
			fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}		
		if(qe>gData0[0].start){	//sorted by start
			//find the 1st rs < qe from right
			tL = 0, tR=tmpi1;
			while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
				tM = (tL+tR)/2; 			//possible case: tmpi=2 or 1 or [tmpi1].s<qe
				if(gData0[tM].start < qe)	 
				    tL = tM;
				else
				    tR = tM;				
			}
			if(gData0[tR].start<qe)tL = tR;
			//-----------------------------------------
			for(i=tL; i>=0; i--){			
				if(gData0[i].end>qs){
					hits[gData0[i].idx]++;
				} 
			}			
		}
		if(n2>n1){									//n2>n1
			int32_t bd = IGD->nbp*(n1+1);			//only keep the first		
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
				tmpi = IGD->nCnt[ichr][j];
				tmpi1 = tmpi-1;
				if(tmpi>0){
					if(j!=preIdx || ichr!=preChr){
						fseek(fP, IGD->tIdx[ichr][j], SEEK_SET);			
						free(gData0);					
						gData0 = malloc(tmpi*sizeof(gdata0_t));
						fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
						preIdx = j;
						preChr = ichr;
					}	
					if(qe>gData0[0].start){	
						//find the 1st rs < qe						
						tS = 0;
						while(tS<tmpi && gData0[tS].start<bd)tS++;	//qs<bd	
						tL = 0, tR=tmpi1;
						while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
							tM = (tL+tR)/2; 
							if(gData0[tM].start < qe)	 
								tL = tM;
							else
								tR = tM;				
						}
						if(gData0[tR].start<qe)tL = tR;
						//-----------------------------------------
						for(i=tL; i>=tS; i--){ 		
							if(gData0[i].end>qs){
								hits[gData0[i].idx]++;
							} 
						}
					}
				}
				bd+=IGD->nbp;		
			}	
		}
	}
	//-----------------------------------------------------	
    return nols;
}

int64_t getOverlaps0(char *qFile, int64_t *hits)
{	//for gdata0_t
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	if ((fp = gzopen(qFile, "r")) == 0)
		return 0;
	ks = ks_init(fp);                              
    uint64_t ols = 0; 		
    char *chrm;
	int32_t st, en, nl;
	preChr = -6, preIdx=-8;     			
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		chrm = parse_bed(str.s, &st, &en);
		if (chrm) {
			nl = get_overlaps0(chrm, st, en, hits);
			ols += nl;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);	
	return ols;
}

//--------------------------------------------------------------------------------------
void seq_overlaps(char *chrm, int32_t qs, int32_t qe, overlap_t *hits, uint32_t *n, uint32_t *m)
{   
	uint32_t nt=0, mt=*m;
	float qlen=qe-qs, st, rlen;
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;
	int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t re, rs, tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
	if(n1>mTile) 
		return 0;
	n2 = MIN(n2, mTile);	
	tmpi = IGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	if(tmpi>0){
		if(n1!=preIdx || ichr!=preChr){
			fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
			free(gData);					
			gData = malloc(tmpi*sizeof(gdata_t));
			fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}		
		if(qe>gData[0].start){				//sorted by start
			//find the 1st rs < qe from right
			tL = 0, tR=tmpi1;
			while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
				tM = (tL+tR)/2; 			//possible case: tmpi=2 or 1 or [tmpi1].s<qe
				if(gData[tM].start < qe)	 
				    tL = tM;
				else
				    tR = tM;				
			}
			if(gData[tR].start<qe)tL = tR;
			//-----------------------------------------
			for(i=tL; i>=0; i--){			
				if(gData[i].end>qs){
					if(nt==mt)EXPAND(hits, mt);
					re = gData[i].end;
					rs = gData[i].start;					
					st = MIN(qe, re)-MAX(qs, rs);
					rlen = re-rs;//this is necessary: float range small		
					hits[nt].idx_t = n1;
					hits[nt].idx_g = i;
					hits[nt].idx_f = gData[i].idx;
					hits[nt].sm = st/(qlen+rlen-st);
					nt++;
				} 
			}			
		}
		if(n2>n1){									//n2>n1
			int32_t bd = IGD->nbp*(n1+1);			//only keep the first		
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
				tmpi = IGD->nCnt[ichr][j];
				tmpi1 = tmpi-1;
				if(tmpi>0){
					if(j!=preIdx || ichr!=preChr){
						fseek(fP, IGD->tIdx[ichr][j], SEEK_SET);			
						free(gData);					
						gData = malloc(tmpi*sizeof(gdata_t));
						fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
						preIdx = j;
						preChr = ichr;
					}	
					if(qe>gData0[0].start){	
						//find the 1st rs < qe						
						tS = 0;
						while(tS<tmpi && gData[tS].start<bd)tS++;	//qs<bd	
						tL = 0, tR=tmpi1;
						while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
							tM = (tL+tR)/2; 
							if(gData[tM].start < qe)	 
								tL = tM;
							else
								tR = tM;				
						}
						if(gData[tR].start<qe)tL = tR;
						//-----------------------------------------
						for(i=tL; i>=tS; i--){ 		
							if(gData[i].end>qs){
								if(nt==mt)EXPAND(hits, mt);
								re = gData[i].end;
								rs = gData[i].start;					
								st = MIN(qe, re)-MAX(qs, rs);
								rlen = re-rs;//this is necessary: float range small		
								hits[nt].idx_t = j;
								hits[nt].idx_g = i;
								hits[nt].idx_f = gData[i].idx;
								hits[nt].sm = st/(qlen+rlen-st);
								nt++;				
							} 
						}
					}
				}
				bd+=IGD->nbp;		
			}	
		}
	}
	//-----------------------------------------------------	
	*n = nt;
	*m = mt;
    return;
}

void seqOverlaps(char *qFile, double *sm)
{	//return similarities
	//-----------------------------------------------------
	//1. get query set and sort
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	if ((fp = gzopen(qFile, "r")) == 0)
		return 0;
	ks = ks_init(fp); 		
    char *chrm;
	int32_t st, en, nl;
	preChr = -6, preIdx=-8;     			
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		chrm = parse_bed(str.s, &st, &en);
		if (chrm)
			nl = seq_overlaps(chrm, st, en, hits);
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	//-----------------------------------------------------
	//2. calculate overlap for each chromosome 
	
	
	//-----------------------------------------------------
	//3. calculate final sm for each dataset
	
	return;
}

//--------------------------------------------------------------------------------
int32_t get_overlaps(char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{   
	int ichr = get_id(chrm);	
	if(ichr<0)
		return 0;
	int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
	int32_t nols = 0;
	if(n1>mTile) 
		return 0;
	n2 = MIN(n2, mTile);

	tmpi = IGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	if(tmpi>0){
		if(n1!=preIdx || ichr!=preChr){
			fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
			free(gData);					
			gData = malloc(tmpi*sizeof(gdata_t));
			fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}		
		if(qe>gData[0].start){						//sorted by start
			//find the 1st rs < qe
			tL = 0, tR=tmpi1;
			while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
				tM = (tL+tR)/2; 
				if(gData[tM].start < qe)	//right side: 
				    tL = tM;
				else
				    tR = tM;				//left side
			}
			if(gData[tR].start<qe)tL = tR;
			//-----------------------------------------
			for(i=tL; i>=0; i--){			
				if(gData[i].end>qs){
					hits[gData[i].idx]++;
				} 
			}			
		}
		if(n2>n1){									//n2>n1
			int32_t bd = IGD->nbp*(n1+1);			//only keep the first		
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
				tmpi = IGD->nCnt[ichr][j];
				tmpi1 = tmpi-1;
				if(tmpi>0){
					if(j!=preIdx || ichr!=preChr){
						fseek(fP, IGD->tIdx[ichr][j], SEEK_SET);			
						free(gData);					
						gData = malloc(tmpi*sizeof(gdata_t));
						fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
						preIdx = j;
						preChr = ichr;
					}	
					if(qe>gData[0].start){						
						tS = 0;
						while(tS<tmpi && gData[tS].start<bd)tS++;	//qs<bd								
						tL = 0, tR=tmpi1;
						while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
							tM = (tL+tR)/2; 
							if(gData[tM].start < qe)	//right side: 
								tL = tM;
							else
								tR = tM;				//left side
						}
						if(gData[tR].start<qe)tL = tR;
						//-----------------------------------------
						for(i=tL; i>=tS; i--){ 		
							if(gData[i].end>qs){
								hits[gData[i].idx]++;
							} 
						}
					}
				}
				bd+=IGD->nbp;		
			}	
		}
	}
	//-----------------------------------------------------	
    return nols;
}

//for gData with v
int32_t get_overlaps_v(char *chrm, int32_t qs, int32_t qe, int32_t v, int64_t *hits)
{   //no need to store every overlaps, only get the number of hits
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;
	int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tE, tS, tmpi, mTile = IGD->nTile[ichr]-1;
	int32_t nols = 0;
	if(n1>mTile) 
		return 0;
	n2 = MIN(n2, mTile);	
	tmpi = IGD->nCnt[ichr][n1];
	if(tmpi>0){
		if(n1!=preIdx || ichr!=preChr){
			fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
			free(gData);					
			gData = malloc(tmpi*sizeof(gdata_t));
			fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}
		if(qe>gData[0].start){					//sorted by start
			if(tmpi<16){
				tE = tmpi-1;
				while(gData[tE].start>=qe)tE--;
			}
			else
				tE = bSearch(gData, 0, tmpi-1, qe);	//idx
			for(i=tE; i>=0; i--){
				if(gData[i].end>qs && gData[i].value>=v){
					nols++;
					hits[gData[i].idx]++;
				} 
			}
		}
	}
	if(n2>n1){									//n2>n1
		int32_t bd = IGD->nbp*(n1+1);				//only keep the first
		for(j=n1+1; j<=n2; j++){
			tmpi = IGD->nCnt[ichr][j];
			if(tmpi>0){
				if(j!=preIdx || ichr!=preChr){
					fseek(fP, IGD->tIdx[ichr][j], SEEK_SET);			
					free(gData);					
					gData = malloc(tmpi*sizeof(gdata_t));
					fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
					preIdx = j;
					preChr = ichr;
				}
				if(qe>gData[0].start){
					if(tmpi<16){
						tE = tmpi-1;
						while(gData[tE].start>=qe)tE--;
					}
					else
						tE = bSearch(gData, 0, tmpi-1, qe);	//idx
					tS = 0;
					while(tS<tmpi && gData[tS].start<bd)tS++;		//exclude <left boundary
					for(i=tE; i>=tS;i--){ 
						if(gData[i].end>qs && gData[i].value>=v){
							nols++;
							hits[gData[i].idx]++;
						}
					} 
				}
			}
			bd+=IGD->nbp;		
		}	
	}
	//-----------------------------------------------------	
    return nols;
}

int64_t getOverlaps(char *qFile, int64_t *hits)
{	
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	if ((fp = gzopen(qFile, "r")) == 0)
		return 0;
	ks = ks_init(fp);                              
    uint64_t ols = 0; 		
    char *chrm;
	int32_t st, en, nl;
	preChr = -6, preIdx=-8;     			
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		chrm = parse_bed(str.s, &st, &en);
		if (chrm) {
			nl = get_overlaps(chrm, st, en, hits);
			ols += nl;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);	
	return ols;
}

int64_t getOverlaps_v(char *qFile, int64_t *hits, int32_t v)
{	//for gadta1_t
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	if ((fp = gzopen(qFile, "r")) == 0)
		return 0;
	ks = ks_init(fp);                              
    uint64_t ols = 0; 		
    char *chrm;
	int32_t st, en, nl;
	preChr = -6, preIdx=-8;     			
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		chrm = parse_bed(str.s, &st, &en);
		if (chrm) {
			nl = get_overlaps_v(chrm, st, en, v, hits);
			ols += nl;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);	
	return ols;
}

//using AIList: no decomp
int64_t getMap(uint32_t **hitmap)
{	//load igd tile one by one
	int i, j, jj, ichr, n1, m=0;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			m++;
			if(m%1000==0)
				printf("%i\n", m);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//construct ailist--------------------------------
				maxE = malloc(tmpi*sizeof(int32_t));
				tmax = gData[0].end;
				for(i=0;i<tmpi;i++){
					if(gData[i].end>tmax)tmax = gData[i].end;
					maxE[i]=tmax;
				}
				for(j=0;j<tmpi;j++){
					qe = gData[j].end;
					qs = gData[j].start;
					if(qe>gData[0].start){					
						jj = gData[j].idx;						
						tS = 0;
						if(qs<bd)
							while(tS<tmpi && gData[tS].start<bd)tS++;		//exclude 
						//---------------------------------------------------
						if(tmpi<16){
							i = tmpi-1;
							while(gData[i].start>=qe)i--;
						}
						else
							i = bSearch(gData, tS, tmpi-1, qe);	//idx
						while(i>=tS && maxE[i]>qs){
							if(gData[i].end>qs){
								nols++;
								hitmap[jj][gData[i].idx]++;
							}
							i--;
						} 
					}
				}
				free(maxE);	
			}
		}
	}
    return nols;
}

//using AIList: no decomp
int64_t getMap_v(uint32_t **hitmap, int32_t v)
{	//load igd tile one by one
	int i, j, jj, ichr, n1, m=0;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			m++;
			if(m%1000==0)
				printf("%i\n", m);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//construct ailist---------------------------
				maxE = malloc(tmpi*sizeof(int32_t));
				tmax = gData[0].end;
				for(i=0;i<tmpi;i++){
					if(gData[i].end>tmax)tmax = gData[i].end;
					maxE[i]=tmax;
				}
				//-------------------------------------------
				for(j=0;j<tmpi;j++){
					if(gData[j].value>v){
						qe = gData[j].end;
						qs = gData[j].start;
						if(qe>gData[0].start){					
							jj = gData[j].idx;						
							tS = 0;
							if(qs<bd)
								while(tS<tmpi && gData[tS].start<bd)tS++;		//exclude 
							//---------------------------------------------------
							if(tmpi<16){
								i = tmpi-1;
								while(gData[i].start>=qe)i--;
							}
							else
								i = bSearch(gData, tS, tmpi-1, qe);	//idx
							while(i>=tS && maxE[i]>qs){
								if(gData[i].end>qs && gData[i].value>v){
									nols++;
									hitmap[jj][gData[i].idx]++;
								}
								i--;
							} 
						}
					}
				}
				free(maxE);	
			}
		}
	}
    return nols;
}

//the following functions are based on Seqpare (S Feng 2020)
void seq_compare_1n(char *iPath, char *qfile, int32_t *nf, seqpare_t **sps)
{   
    //1. Get the files  
    glob_t gResult;
    int rtn = glob(iPath, 0, NULL, &gResult);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath);
        return;
    }
    char** file_ids = gResult.gl_pathv;
    *nf = gResult.gl_pathc; 
    if(*nf<1)   
        printf("Too few files (add to path /*): %i\n", *nf);  
    seqpare_t *sm = malloc(*nf*sizeof(seqpare_t));
    
	//2. Calculate the s index
	ailist0_t *ail2 = readBED0(qfile);			//as query
	ailist0_construct(ail2, 20);
	//-------------------------
	for(int i=0;i<*nf;i++){
		ailist_t *ail1 = readBED(file_ids[i]);	//as database: fresh q_idx, s_max
		ailist_construct(ail1, 20);
		seq_compare(ail1, ail2, &sm[i]);
		ailist_destroy(ail1);
	}
	//-------------------------
	ailist0_destroy(ail2);   
	globfree(&gResult);
	*sps = sm;
}

void seq_compare_11(char *file1, char *file2, seqpare_t *sp)
{ 
	ailist_t *ail1 = readBED(file1);	//as Database
	ailist_construct(ail1, 20);
	//---------------------------------
	ailist0_t *ail2 = readBED0(file2);	//file2 as query
	ailist0_construct(ail2, 20);
	//---------------------------------
	seq_compare(ail1, ail2, sp);
	//printf("s index = %12.8f\n", *sm);
	ailist_destroy(ail1);
	ailist0_destroy(ail2);
}

void seq_compare(ailist_t *ail1, ailist0_t *ail2, seqpare_t *sp)
{   //need to refresh ail1: database--q_idx, s_max
	float ERR = 0.0000001;
	int32_t i, j, k;
	int32_t N1=0, N2=0;
	for(i=0;i<ail1->nctg;i++)
		N1 += (&ail1->ctg[i])->nr;
	for(i=0;i<ail2->nctg;i++)
		N2 += (&ail2->ctg[i])->nr;	
	sp->N1 = N1;
	sp->N2 = N2;		
	//-------------------------------------------------------------------------
	int32_t qs, qe, rs, re, gid, cs, ce, t, iq, tmax;
	float qlen, s, rlen, smax, st;
	ctg_t *p1;
	ctg0_t *p2;
	for(j=0;j<ail2->nctg;j++){						//ail2: query	
		p2 = &ail2->ctg[j];
		gid = get_ctg(ail1, p2->name);
		if(gid>=ail1->nctg || gid<0)continue;			//get gid in ail1: db					
		p1 = &ail1->ctg[gid];
		for(i=0;i<p2->nr;i++){	
			qs = p2->glist[i].start;
			qe = p2->glist[i].end;				    	
			qlen = qe-qs; 
			tmax=-1;
			smax=0.0;		
			//----------------------------------------
			for(k=0; k<p1->nc; k++){						//search each component
				cs = p1->idxC[k];
				ce = cs + p1->lenC[k];			
			    t = bSearch(p1->glist, cs, ce, qe); 		//rs<qe: inline not better 
				if(t>=cs){
					while(t>=cs && p1->maxE[t]>qs){
					    if((re=p1->glist[t].end)>qs){ 
			    			p1->glist[t].cnts++;					     
					    	rs = p1->glist[t].start;             	
					        s = MIN(qe, re)-MAX(qs, rs);
					        rlen = re-rs;					//this is necessary
					        s = s/(qlen+rlen-s);
					        //------------------------------skip: matched but smaller
					        if(s>smax){						//not matched or matched & larger
					        	st = p1->glist[t].s_max;
					        	if(st<ERR || (st>ERR && s>st)){
							    	smax = s;
							    	tmax = t;
					        	}
					        }
						}
					    t--;
					}
			    }
			}  
			//if(smax<ERR)continue;                 		 
			if(smax>ERR){
				st = p1->glist[tmax].s_max;
				if(st<ERR){								//glist[tmax] not taken, mark it
					p1->glist[tmax].s_max = smax;
					p1->glist[tmax].q_idx = i;
				}
				else if(st<smax){						//mutual max[i, tmax]/rerun .q_idx			
					iq = p1->glist[tmax].q_idx;
					p1->glist[tmax].s_max = smax;
					p1->glist[tmax].q_idx = i;				
					//----rerun on iq to get max other than [tmax]: will be skipped					
					seqpare0(ail1, ail2, gid, j, iq);	
				}
			}
			//else{									//search i again without [tmax]
			//}
		}// glist[i]
	}// ctg[j]	   
	//-------------------------------------------------------------------------
	//calculate si: similarity index
	float sm = 0.0;
	int64_t cnt = 0;
	for(j=0;j<ail1->nctg;j++){
		p1 = &ail1->ctg[j];
		for(i=0;i<p1->nr;i++){
			cnt += p1->glist[i].cnts;
			if((st=p1->glist[i].s_max)>ERR)
				sm += st;
		}
	}	
	//printf("%i\t %i\t %lld\t %12.8f\n", N1, N2, (long long)cnt, sm);
	sp->teo = sm;    
	sp->tc = cnt;
	return;                   
}

void seqpare0(ailist_t *ail1, ailist0_t *ail2, int gid1, int gid2, int iq)
{	//search iq for max other than ir: may recursive
	float ERR = 0.0000001;
	ctg_t *p1 = &ail1->ctg[gid1];
	ctg0_t *p2 = &ail2->ctg[gid2];			
	//------------------------------------
	int32_t qs = p2->glist[iq].start;
	int32_t qe = p2->glist[iq].end;				    	
	int32_t t, cs, ce, rs, re, iq0=iq, tmax=-1; 
	float smax=0.0, s, st, qlen = qe-qs, rlen;
	for(int k=0; k<p1->nc; k++){						//search each component
		cs = p1->idxC[k];
		ce = cs + p1->lenC[k];			
	    t = bSearch(p1->glist, cs, ce, qe); 			//rs<qe: inline not better 
	    if(t>=cs){
			while(t>=cs && p1->maxE[t]>qs){
			    if((re=p1->glist[t].end)>qs){  
			    	rs = p1->glist[t].start;             	
			        s = MIN(qe, re)-MAX(qs, rs);
			        rlen = re-rs;
			        s = s/(qlen+rlen-s);
		        	st = p1->glist[t].s_max;
		        	if(s>smax && (st<ERR || (st>ERR && s>st))){
						smax = s;
						tmax = t;
		        	}
				}
			    t--;
			}
	    }
	} 
	if(smax<ERR)return;                  		 
	//------------------------------------
	st = p1->glist[tmax].s_max;
	if(st<ERR){								//glist[tmax] not taken, mark it
		p1->glist[tmax].s_max = smax;
		p1->glist[tmax].q_idx = iq0;
	}
	else if(st<smax){						//mutual max[i, tmax]/rerun .q_idx			
		iq = p1->glist[tmax].q_idx;
		p1->glist[tmax].s_max = smax;
		p1->glist[tmax].q_idx = iq0;				
		//----rerun on iq to get max other than [tmax]	
		//printf("recursive: %i\t %i\t %i\n", gid1, gid2, iq);				
		seqpare0(ail1, ail2, gid1, gid2, iq);	
	} 
} 


//-------------------------------------------------------------------------------------
int igd_search(int argc, char **argv)
{   //igd[0] search[1] home/john/iGD/rme_igd/roadmap.igd[2] -q[3] query100.bed[4]
    if(argc<4)
        return search_help(EX_OK);     
    
    char *igdName = argv[2]; 
    char *ftype = igdName + strlen(igdName) - 4;
    if(strcmp(".igd", ftype)!=0){
        printf("%s is not an igd database", igdName);
        return EX_OK;
    }    
    FILE* fi = fopen(igdName, "rb");
    if(!fi){
        printf("%s does not exist", igdName);
        return EX_OK;
    }
    fclose(fi); 
        
    int32_t v = 0, qs=1, qe=2;
    int i, i1, j, checking=0, mode=-1, mt=0, ext=0, xlen=0,  mv=0, mx=0, ichr, k;
    char out[64]="";  
    char *chrm;    
    char *qfName = "";  
    uint64_t ols;
    
    //----------------------------------------------------- 
    char tmp[128];  
    IGD = get_igdinfo(igdName);     
    strcpy(tmp, igdName);
    tmp[strrchr(tmp, '.')-tmp] = '\0';
    strcpy(IGD->fname, tmp);
    char *idFile = tmp;	//str_split(tmp, '.', &nCols)[0];
    strcat(idFile, "_index.tsv");            
    IGD->finfo = get_fileinfo(idFile, &IGD->nFiles);      
    int32_t nfiles = IGD->nFiles; 
    int64_t *hits = calloc(nfiles, sizeof(int64_t));  
    //-----------------------------------------------------
    //for(i=0;i<nfiles;i++){
    //	printf("%i\t%i\t%i\n", i, IGD->finfo[i].nr, IGD->finfo[i].md);
    //}
    for(i=3; i<argc; i++){
        if(strcmp(argv[i], "-q")==0){
            if(i+1<argc){
                qfName = argv[i+1];
                mode = 1;
            }
            else{
            	printf("No query file.\n");
				return EX_OK;
			}	
        }
        else if(strcmp(argv[i], "-r")==0){
            if(i+3<argc){
                mode = 2;
                chrm = argv[i+1];
                qs = atoi(argv[i+2]);
                qe = atoi(argv[i+3]);
            }        
        }
        else if(strcmp(argv[i], "-v")==0){//>=v
        	mv = 1;
            if(i+1<argc)
                v = atoi(argv[i+1]);
        }
        else if(strcmp(argv[i], "-m")==0){
            mode = 0;
        }         
        else if(strcmp(argv[i], "-o")==0){
            if(i+1<argc)
                strcpy(out, argv[i+1]);
        }   
        else if(strcmp(argv[i], "-c")==0){
            checking = 1;
        }                              
    }  
     
    //----------------------------------------------------------
	fP = fopen(igdName, "rb");				//share
	if(mode==0){
    	uint32_t **hitmap = malloc(nfiles*sizeof(uint32_t*));
    	for(i=0;i<nfiles;i++)
    		hitmap[i] = calloc(nfiles, sizeof(uint32_t));
    	if(v>0)
			getMap_v(hitmap, v);
		else
			getMap(hitmap);    		   		
    	FILE *fp;
		if(strlen(out)<2)strcpy(out,"Hitmap");
		fp = fopen(out, "w");
	    if(fp==NULL)
	        printf("Can't open file %s\n", out);
	    else{
	        fprintf(fp, "%u\t%u\t%u\n", nfiles, nfiles, v);
	        for(i=0;i<nfiles;i++){
	            for(j=0;j<nfiles;j++)
	                fprintf(fp, "%u\t", hitmap[i][j]); 
	            fprintf(fp, "\n");
	        } 
	        fclose(fp);
	    }     
	
    	for(i=0;i<nfiles;i++)
    		free(hitmap[i]);
    	free(hitmap);  	
	}
	else if(mode==1){//for a query dataset (file)  
		if(IGD->gType==0)
			getOverlaps0(qfName, hits);
		else{
			if(v>0)
				getOverlaps_v(qfName, hits, v);
			else
				getOverlaps(qfName, hits);
		}
		printf("index\t File_name\t number of regions\t number of hits\n"); 
		int64_t total = 0;       
	    for(i=0;i<nfiles;i++){
	    	if(hits[i]>0)
	        	printf("%i\t%i\t%lld\t%s\n", i, IGD->finfo[i].nr, (long long)hits[i], IGD->finfo[i].fileName); 
	    	total += hits[i];
	    }
	    printf("Total: %lld\n", (long long)total);
    }
    else if(mode==2){//mode 2 for a single region
    	if(IGD->gType==0)
    		ols = get_overlaps0(chrm, qs, qe, hits); 
    	else{
    		if(v>0)
    			ols = get_overlaps_v(chrm, qs, qe, v, hits);
    		else
     			ols = get_overlaps(chrm, qs, qe, hits);    			
    	}
    	printf("index\t File_name\t number of regions\t number of hits\n");        
        for(i=0;i<nfiles;i++)
            printf("%i\t%i\t%lld\t%s\n", i, IGD->finfo[i].nr, (long long)hits[i], IGD->finfo[i].fileName);
    }
    else if(mode==3){//output seqpare index
    	double *sm = calloc(nfiles, sizeof(double));
		seqOverlaps(qfName, sm);
		printf("index\t File_name\t number of regions\t similarity\t dataset name\n"); 
		int64_t total = 0;       
	    for(i=0;i<nfiles;i++){
	    	if(hits[i]>0)
	        	printf("%i\t%i\t%lld\t%s\n", i, IGD->finfo[i].nr, sm[i], IGD->finfo[i].fileName); 
	    	total += hits[i];
	    }  
	    free(sm);
    }
    else
        return search_help(EX_OK);      

	fclose(fP);
    free(IGD->nTile);
    for(i=0;i<IGD->nCtg;i++){
     	free(IGD->nCnt[i]);
     	free(IGD->tIdx[i]);
    }
    free(IGD->nCnt);
    free(IGD->tIdx);
    free(IGD->cName);
    free(IGD->finfo);
    free(IGD);
    free(hits);
    return EX_OK;
}

