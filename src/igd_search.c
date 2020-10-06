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
"             -s output Seqpare similarity\n"
"             -f output full overlaps (for -q and -r only)\n"
"             -m hitsmap of igd datasets\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

int32_t get_overlaps0(char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{   
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;		
		
	int i, j, ni, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
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
			ni=fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
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
						ni=fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
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

int32_t get_overlaps_f0(char *chrm, int32_t qs, int32_t qe)
{   //Directly print the overlaps on the terminal < output file by user
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;		
		
	int i, j, ni, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
	int32_t nols = 0;
	if(n1>mTile) 
		return 0;
	//-------------------------------------------------
	printf("Query %s, %i, %i: \n", chrm, qs, qe);
	n2 = MIN(n2, mTile);	
	tmpi = IGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	if(tmpi>0){
		if(n1!=preIdx || ichr!=preChr){
			fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
			free(gData0);					
			gData0 = malloc(tmpi*sizeof(gdata0_t));
			ni=fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
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
					printf("%i\t %i\t %i\t %s\n", nols++, gData0[i].start, gData0[i].end, IGD->finfo[gData0[i].idx].fileName);					
					//hits[gData0[i].idx]++;
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
						ni=fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
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
								printf("%i\t %i\t %i\t %s\n", nols++, gData0[i].start, gData0[i].end, IGD->finfo[gData0[i].idx].fileName);					
								//hits[gData0[i].idx]++;
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

int64_t getOverlaps_f0(char *qFile)
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
			nl = get_overlaps_f0(chrm, st, en);
			ols += nl;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);	
	return ols;
}

//--------------------------------------------------------------------------------------
void seq_overlaps(char *chrm, int32_t qs, int32_t qe, overlaps_t *olp)
{   
	float qlen=qe-qs, st, rlen;
	int ichr = get_id(chrm);
	if(ichr<0)
		return;
	int i, j, ni, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t re, rs, tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
	if(n1>mTile) 
		return;
	n2 = MIN(n2, mTile);	
	tmpi = IGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	if(tmpi>0){
		if(n1!=preIdx || ichr!=preChr){
			fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);	
			free(gData);								
			gData = malloc(tmpi*sizeof(gdata_t));			
			ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}		
		if(qe>gData[0].start){			//sorted by start
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
					if(olp->nn==olp->mm)
						EXPAND(olp->olist, olp->mm);					
					re = gData[i].end;
					rs = gData[i].start;					
					st = MIN(qe, re)-MAX(qs, rs);
					rlen = re-rs;			//this is necessary	
					overlap_t *p = &olp->olist[olp->nn++];						 
					p->idx_g = i;
					p->idx_f = gData[i].idx;				
					p->idx_t = n1;						
					p->sm = st/(qlen+rlen-st);
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
						ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
						preIdx = j;
						preChr = ichr;
					}	
					if(qe>gData[0].start){	
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
								if(olp->nn==olp->mm)
									EXPAND(olp->olist, olp->mm);					
								re = gData[i].end;
								rs = gData[i].start;					
								st = MIN(qe, re)-MAX(qs, rs);
								rlen = re-rs;			//this is necessary	
								overlap_t *p = &olp->olist[olp->nn++];						 
								p->idx_g = i;
								p->idx_f = gData[i].idx;				
								p->idx_t = n1;						
								p->sm = st/(qlen+rlen-st);				
							} 
						}
					}
				}
				bd+=IGD->nbp;		
			}	
		}
	}
    return;
}

void seqOverlaps(char *qFile, double *sm)
{	//calculate similarities
	//-----------------------------------------------------
	ailist_t *ail = readBED(qFile);
	//-----------------------------------------------------
	//calculate overlap for each chromosome 
	int i, j, k, m, nj, nk, ig, it, idx, maxk, maxj, nq, Nq=0;
	float maxf;
	int nfiles = IGD->nFiles;
	preChr = -6, preIdx=-8;  
	for(m=0;m<nfiles;m++)
		sm[m] = 0.0;			
	for(i=0; i<ail->nctg; i++){
		chrom_t *p  = &ail->ctg[i];
		gdata0_t *L1 = p->glist;						
		nq 			= p->nr;
		Nq += nq;		
		//radix_sort_intv(L1, L1+nq);
		qsort(L1, nq, sizeof(gdata0_t), compare_qstart);
		overlap_t **olps = malloc(nq*sizeof(overlap_t*));
		int *nh = malloc(nq*sizeof(int)); 
		//printf("%i\t %i\t %s\n",i, nq, p->name);			
		overlaps_t *olp = (overlaps_t *)malloc(1*sizeof(overlaps_t));
		olp->mm = 1000000;
		olp->olist = malloc(olp->mm*sizeof(overlap_t));		
		for(j=0; j<nq; j++){	
			olp->nn = 0;
			seq_overlaps(p->name, L1[j].start, L1[j].end, olp);
			if(olp->nn > 0){
				olps[j] = (overlap_t *)malloc(olp->nn*sizeof(overlap_t));				
				qsort(olp->olist, olp->nn, sizeof(overlap_t), compare_fidx);				
				memcpy(olps[j], olp->olist, olp->nn*sizeof(overlap_t));
			}
			nh[j]=olp->nn;
		}	
		free(olp->olist);
		free(olp);

		//---deal with one dataset a time
		int *kst0 = calloc(nq, sizeof(int));	//k-start for current idx_f
		int *kst = calloc(nq, sizeof(int));		//k-start for next idx_f
		int *nst0 = calloc(nq, sizeof(int));	//length
		for(m=0; m<nfiles;  m++){
			//1. Find the max
			maxf = 0.0;
			for(j=0;j<nq;j++){
				k=kst[j];
				while(k<nh[j] && olps[j][k].idx_f<m)k++;
				kst0[j] = k;//start of the current m: for section 2
				while(k<nh[j] && olps[j][k].idx_f==m){
					if(olps[j][k].sm>maxf){
						maxf = olps[j][k].sm;
						maxk = k;
						maxj = j;
					}
					k++;	
				}
				kst[j]=k;	//start of the next m
				nst0[j] = k-kst0[j];
				//nols[m] += nst0[j];
			}
		
			//2. Record and Remove the j-th row (set nst0=0) and column(set sm=0.0)
			while(maxf>0.0){
				sm[m] += maxf;
				nst0[maxj] = 0;
				it = olps[maxj][maxk].idx_t;
				ig = olps[maxj][maxk].idx_g;
				maxf = 0.0;
				for(j=0;j<nq;j++){
					if(nst0[j]>0){
						for(k=kst0[j]; k<kst0[j]+nst0[j]; k++){
							if(olps[j][k].idx_g==ig && olps[j][k].idx_t==it)
								olps[j][k].sm = 0.0;
							else if(olps[j][k].sm>maxf){
								maxf = olps[j][k].sm;
								maxk = k;
								maxj = j;
							}
						}
					}
				}			
			}
		}
		free(nst0);
		free(kst);
		free(kst0);	
		free(nh);
		free(olps);
	}
	for(m=0;m<nfiles;m++){
		//printf("%i\t %i\t %i\n", m, nols0[m], nols[m]);
		sm[m] = sm[m]/(Nq+IGD->finfo[m].nr-sm[m]);
	}	
	ailist_destroy(ail);
	return;
}

//--------------------------------------------------------------------------------
int32_t get_overlaps(char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{   
	int ichr = get_id(chrm);	
	if(ichr<0)
		return 0;
	int i, j, ni, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
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
			ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}		
		if(qe>gData[0].start){				//sorted by start
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
						ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
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
    return nols;
}

//--------------------------------------------------------------------------------
int32_t get_overlaps_f1(char *chrm, int32_t qs, int32_t qe)
{   
	int ichr = get_id(chrm);	
	if(ichr<0)
		return 0;
	int i, j, ni, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
	int32_t nols = 0;
	if(n1>mTile) 
		return 0;
	//-------------------------------------------------
	printf("Query %s, %i, %i: \n", chrm, qs, qe);		
		
	n2 = MIN(n2, mTile);
	
	tmpi = IGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	if(tmpi>0){
		if(n1!=preIdx || ichr!=preChr){
			fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
			free(gData);					
			gData = malloc(tmpi*sizeof(gdata_t));
			ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}		
		if(qe>gData[0].start){				//sorted by start
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
					printf("%i\t %i\t %i\t %s\n", nols++, gData[i].start, gData[i].end, IGD->finfo[gData[i].idx].fileName);					
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
						ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
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
								printf("%i\t %i\t %i\t %s\n", nols++, gData[i].start, gData[i].end, IGD->finfo[gData[i].idx].fileName);		
							} 
						}
					}
				}
				bd+=IGD->nbp;		
			}	
		}
	}
    return nols;
}

//for gData with v
int32_t get_overlaps_v(char *chrm, int32_t qs, int32_t qe, int32_t v, int64_t *hits)
{   //no need to store every overlaps, only get the number of hits
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;
	int i, ni, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
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
			ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}
		if(qe>gData[0].start){							//sorted by start
			if(tmpi<16){
				tE = tmpi-1;
				while(gData[tE].start>=qe)tE--;
			}
			else
				tE = bSearch(gData, 0, tmpi-1, qe);		//idx
			for(i=tE; i>=0; i--){
				if(gData[i].end>qs && gData[i].value>=v){
					nols++;
					hits[gData[i].idx]++;
				} 
			}
		}
	}
	if(n2>n1){											//n2>n1
		int32_t bd = IGD->nbp*(n1+1);					//only keep the first
		for(j=n1+1; j<=n2; j++){
			tmpi = IGD->nCnt[ichr][j];
			if(tmpi>0){
				if(j!=preIdx || ichr!=preChr){
					fseek(fP, IGD->tIdx[ichr][j], SEEK_SET);			
					free(gData);					
					gData = malloc(tmpi*sizeof(gdata_t));
					ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
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

int64_t getOverlaps_f1(char *qFile)
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
			nl = get_overlaps_f1(chrm, st, en);
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
	int i, j, jj, ichr, ni, n1, m=0;	
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
				ni = fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
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
	int i, j, jj, ichr, ni, n1, m=0;	
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
				ni=fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
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
    int i, i1, j, checking=0, mode=-1, p_mode=0;
    int mt=0, ext=0, xlen=0,  mv=0, mx=0, ichr, k;
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
    //finfo[i].md is necessary for search -r option
    //for(i=0;i<nfiles;i++){
    //	printf("%i\t%i\t%12.3f\n", i, IGD->finfo[i].nr, IGD->finfo[i].md);
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
        else if(strcmp(argv[i], "-s")==0 && mode!=2){
            mode = 3;	//seqpare
        }                  
        else if(strcmp(argv[i], "-o")==0){
            if(i+1<argc)
                strcpy(out, argv[i+1]);
        } 
        else if(strcmp(argv[i], "-f")==0){
            p_mode = 1;	//output full overlaps
        }          
        else if(strcmp(argv[i], "-c")==0){
            checking = 1;
        }                              
    }  
     
    //----------------------------------------------------------
	fP = fopen(igdName, "rb");				//share
	if(p_mode==1){
		if(mode==1){
			int64_t total = 0;
			if(IGD->gType==0)
				total = getOverlaps_f0(qfName);
			else
				total = getOverlaps_f1(qfName);		       
			printf("Total overlaps: %lld\n", (long long)total);
		}
		else if(mode==2){
			if(IGD->gType==0)
				ols = get_overlaps_f0(chrm, qs, qe); 
			else
				ols = get_overlaps_f1(chrm, qs, qe);    			      
			printf("Total overlaps: %lld\n", (long long)ols);						
		}
		else{
			printf("Not supported -f option\n");
			return EX_OK;
		}
	}
	else if(mode==0){
    	uint32_t **hitmap = malloc(nfiles*sizeof(uint32_t*));
    	for(i=0;i<nfiles;i++)
    		hitmap[i] = calloc(nfiles, sizeof(uint32_t));
    	if(v>0)
			getMap_v(hitmap, v);
		else
			getMap(hitmap);    		   		
    	FILE *fp;
		if(strlen(out)<2)strcpy(out,"Hitsmap");
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
		printf("index\t number of regions\t number of hits\t File_name\n"); 
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
    	printf("index\t number of regions\t number of hits\t File_name\n");        
        for(i=0;i<nfiles;i++)
            printf("%i\t%i\t%lld\t%s\n", i, IGD->finfo[i].nr, (long long)hits[i], IGD->finfo[i].fileName);
    }
    else if(mode==3){//output seqpare index
    	double *sm = malloc(nfiles*sizeof(double));
		seqOverlaps(qfName, sm);
		printf("index\t number of regions\t similarity\t dataset name\n");       
	    for(i=0;i<nfiles;i++){
	        printf("%i\t%i\t%10.6f\t%s\n", i, IGD->finfo[i].nr, sm[i], IGD->finfo[i].fileName); 
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

