//=====================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//
//time ./igd_search Test110000.bed /media/john/Extra/ucsc_igd/ucsc.igd
//database intervals sorted by _start: 8/12/2019
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
"             -x <flank length on query ends>\n"
"             -m heatmap with igd database itself\n"
"             -c display all intersects\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

void constructNew(gdata_t *glist, int32_t nr, int32_t *nc, int32_t *idxC, int32_t *lenC, int32_t *maxE, int cLen)
{   //input list: sorted by starts 
	//sort it be end-->get the idx list: 
	int i, j, t, k, k0, minL=MAX(64, cLen);
	cLen += cLen/2;
    gdata_t *L1 = glist;
	gdata_t *L0 = malloc(nr*sizeof(gdata_t)); 	
    gdata_t *L2 = malloc(nr*sizeof(gdata_t));
    int32_t *di = malloc(nr*sizeof(int32_t));//int64_t?
    gdata_t *D0 = malloc(nr*sizeof(gdata_t)); 
	if(nr <= minL){        
        *nc = 1, lenC[0] = nr, idxC[0] = 0;              
    }
    else{     
		memcpy(L0, L1, nr*sizeof(gdata_t));	
		//assign components
		int iter=0, len, lenT=nr;
		k=0, k0=0;
		while(iter<MAXC && lenT>minL){	
			//---setup di---------------------------------------------------------	
			for(i=0;i<lenT;i++){			//D0:{.start = end, .end=idx}
				D0[i].start =L0[i].end;
				D0[i].end = i;
			}
			radix_sort_intv(D0, D0+lenT);
			for(i=0;i<lenT;i++){			//assign i=29 to L0[i].end=2
				t = D0[i].end;
				di[t] = i-t;				//>0: containment
			}
			//-------------------------------------------------------------------
			len = 0;
			for(t=0;t<lenT-cLen;t++){
				if(di[t]>cLen)
		            memcpy(&L2[len++], &L0[t], sizeof(gdata_t)); 		
				else
					memcpy(&L1[k++], &L0[t], sizeof(gdata_t)); 
			}
			memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata_t)); 
		    k += cLen, lenT = len;                      
		    idxC[iter] = k0;
		    lenC[iter] = k-k0;
		    k0 = k, iter++;
			if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
		        if(lenT>0){
		            memcpy(&L1[k], L2, lenT*sizeof(gdata_t));
		            idxC[iter] = k;
		            lenC[iter] = lenT;
		            lenT = 0;					//exit!!!
		            iter++;
		        }
		       	*nc = iter;                   
		    }
		    else memcpy(L0, L2, lenT*sizeof(gdata_t));        	
		}
	}    
	free(L2),free(L0), free(di), free(D0); 
    //2. Augmentation
    //maxE = malloc(nr*sizeof(int32_t)); 
    for(j=0; j<*nc; j++){ 
        k0 = idxC[j];
        k = k0 + lenC[j];
        int32_t tt = L1[k0].end;
        maxE[k0]=tt;
        for(t=k0+1; t<k; t++){
            if(L1[t].end>tt)tt = L1[t].end;
            maxE[t] = tt;  
        }             
    } 	
}

void construct(gdata_t *glist, int32_t nr, int32_t *nc, int32_t *idxC, int32_t *lenC, int32_t *maxE, int cLen)
{   //input list: sorted by starts  
    int cLen1=cLen/2, j1, minL = MAX(64, cLen);     
    cLen += cLen1;      
    int lenT, len, iter, j, k, k0, t;            	
	//1. Decomposition     	
	gdata_t *L1 = glist;             		               
    if(nr<=minL){        
        *nc = 1, lenC[0] = nr, idxC[0] = 0;                
    }
    else{ 
    	gdata_t *L0 = malloc(nr*sizeof(gdata_t)); 	//L0: serve as input list
        gdata_t *L2 = malloc(nr*sizeof(gdata_t));   //L2: extracted list 
        memcpy(L0, L1, nr*sizeof(gdata_t));			
        iter = 0;	k = 0;	k0 = 0;
        lenT = nr;
        while(iter<MAXC && lenT>minL){   
            len = 0;            
            for(t=0; t<lenT-cLen; t++){
                uint32_t tt = L0[t].end;
                j=1;    j1=1;
                while(j<cLen && j1<cLen1){
                    if(L0[j+t].end>=tt) j1++;
                    j++;
                }
                if(j1<cLen1) 
                	memcpy(&L2[len++], &L0[t], sizeof(gdata_t));
                else 
                	memcpy(&L1[k++], &L0[t], sizeof(gdata_t));                 
            } 
            memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata_t));   
            k += cLen, lenT = len;                
            idxC[iter] = k0;
            lenC[iter] = k-k0;
            k0 = k, iter++;
            if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                if(lenT>0){
                    memcpy(&L1[k], L2, lenT*sizeof(gdata_t));
                    idxC[iter] = k;
                    lenC[iter] = lenT;
                    iter++;
                }
               	*nc = iter;                   
            }
            else memcpy(L0, L2, lenT*sizeof(gdata_t));
        }
        free(L2),free(L0);     
    }
    //2. Augmentation
    //maxE = malloc(nr*sizeof(int32_t)); 
    for(j=0; j<*nc; j++){ 
        k0 = idxC[j];
        k = k0 + lenC[j];
        int32_t tt = L1[k0].end;
        maxE[k0]=tt;
        for(t=k0+1; t<k; t++){
            if(L1[t].end>tt)tt = L1[t].end;
            maxE[t] = tt;  
        }             
    } 
}

//=============================================================================
//For gdata0_t ----------------------------------------------------------------
void construct0(gdata0_t *glist, int32_t nr, int32_t *nc, int32_t *idxC, int32_t *lenC, int32_t *maxE, int cLen)
{   //input list: sorted by starts  
    int cLen1=cLen/2, j1, minL = MAX(64, cLen);     
    cLen += cLen1;      
    int lenT, len, iter, j, k, k0, t;            	
	//1. Decomposition     	
	gdata0_t *L1 = glist;             		               
    if(nr<=minL){        
        *nc = 1, lenC[0] = nr, idxC[0] = 0;                
    }
    else{ 
    	gdata0_t *L0 = malloc(nr*sizeof(gdata0_t)); 	//L0: serve as input list
        gdata0_t *L2 = malloc(nr*sizeof(gdata0_t));   //L2: extracted list 
        memcpy(L0, L1, nr*sizeof(gdata0_t));			
        iter = 0;	k = 0;	k0 = 0;
        lenT = nr;
        while(iter<MAXC && lenT>minL){   
            len = 0;            
            for(t=0; t<lenT-cLen; t++){
                uint32_t tt = L0[t].end;
                j=1;    j1=1;
                while(j<cLen && j1<cLen1){
                    if(L0[j+t].end>=tt) j1++;
                    j++;
                }
                if(j1<cLen1) 
                	memcpy(&L2[len++], &L0[t], sizeof(gdata0_t));
                else 
                	memcpy(&L1[k++], &L0[t], sizeof(gdata0_t));                 
            } 
            memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata0_t));   
            k += cLen, lenT = len;                
            idxC[iter] = k0;
            lenC[iter] = k-k0;
            k0 = k, iter++;
            if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                if(lenT>0){
                    memcpy(&L1[k], L2, lenT*sizeof(gdata0_t));
                    idxC[iter] = k;
                    lenC[iter] = lenT;
                    iter++;
                }
               	*nc = iter;                   
            }
            else memcpy(L0, L2, lenT*sizeof(gdata0_t));
        }
        free(L2),free(L0);     
    }
    //2. Augmentation
    //maxE = malloc(nr*sizeof(int32_t)); 
    for(j=0; j<*nc; j++){ 
        k0 = idxC[j];
        k = k0 + lenC[j];
        int32_t tt = L1[k0].end;
        maxE[k0]=tt;
        for(t=k0+1; t<k; t++){
            if(L1[t].end>tt)tt = L1[t].end;
            maxE[t] = tt;  
        }             
    } 
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
/*
int32_t get_overlaps0(char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{   //for gdat0_t
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
			free(gData0);					
			gData0 = malloc(tmpi*sizeof(gdata0_t));
			fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}		
		if(qe>gData0[0].start){					//sorted by start
			if(tmpi<16){
				tE = tmpi-1;
				while(gData0[tE].start>=qe)tE--;
			}
			else
				tE = bSearch0(gData0, 0, tmpi-1, qe);	//idx
			for(i=tE; i>=0; i--){
				if(gData0[i].end>qs){
					//nols++;
					hits[gData0[i].idx]++;
				} 
			}
		}
		if(n2>n1){									//n2>n1
			int32_t bd = IGD->nbp*(n1+1);			//only keep the first		
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
				tmpi = IGD->nCnt[ichr][j];
				if(tmpi>0){
					if(j!=preIdx || ichr!=preChr){
						fseek(fP, IGD->tIdx[ichr][j], SEEK_SET);			
						free(gData0);					
						gData0 = malloc(tmpi*sizeof(gdata0_t));
						fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
						preIdx = j;
						preChr = ichr;
					}				
					if(qe>gData0[0].start){//bSearch!=-1		
						if(tmpi<16){
							tE = tmpi-1;
							while(gData0[tE].start>=qe)tE--;
						}
						else
							tE = bSearch0(gData0, 0, tmpi-1, qe);	//idx					
						tS = 0;
						while(tS<tmpi && gData0[tS].start<bd)tS++;		//qs<bd
						for(i=tE; i>=tS; i--){ 
							if(gData0[i].end>qs){		
								//nols++;
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
}*/

int32_t get_overlaps0_f(char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{   //replace bSearch with forward sweep
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
			free(gData0);					
			gData0 = malloc(tmpi*sizeof(gdata0_t));
			fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
			preIdx = n1;
			preChr = ichr;
		}		
		if(qe>gData0[0].start){						//sorted by start
			i=0;
			while(i<=tmpi && gData0[i].end<=qs)
				i++;//		
			while(i<tmpi && gData0[i].start<qe){
				if(gData0[i].end>qs)
					hits[gData0[i].idx]++;
				i++;
			}		
		}
		if(n2>n1){									//n2>n1: qs<bd
			int32_t bd = IGD->nbp*(n1+1);			//only keep the first		
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
				tmpi = IGD->nCnt[ichr][j];
				if(tmpi>0){
					if(j!=preIdx || ichr!=preChr){
						fseek(fP, IGD->tIdx[ichr][j], SEEK_SET);			
						free(gData0);					
						gData0 = malloc(tmpi*sizeof(gdata0_t));
						fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
						preIdx = j;
						preChr = ichr;
					}				
					if(qe>gData0[0].start){//qs<maxE	
						i = 0;
						while(i<tmpi && gData0[i].start<bd)
							i++;		//qs<bdi<tmpi && 
						while(i<tmpi && gData0[i].end<=qs)
							i++;		
						while(i<tmpi && gData0[i].start<qe){
							if(gData0[i].end>qs)
								hits[gData0[i].idx]++;
							i++;
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

int64_t getOverlaps0_m0(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, jj, ichr, n1, m=0;	//define boundary!
	int32_t tE, tS, tmpi, bd, qe, qs;
	int64_t nols = 0;
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;
			tmpi = IGD->nCnt[ichr][n1];			
			if(m%1000==0)
				printf("d0_m0 %i\t%i\t%i\n", n1, m, tmpi);		
			if(tmpi>0){							
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData0);					
				gData0 = malloc(tmpi*sizeof(gdata0_t));
				fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
				for(j=0;j<tmpi;j++){
					qe = gData0[j].end;
					qs = gData0[j].start;
					if(qe>gData0[0].start){
						jj = gData0[j].idx;
						if(tmpi<16){
							tE = tmpi-1;
							while(gData0[tE].start>=qe)
								tE--;
						}
						else
							tE = bSearch0(gData0, 0, tmpi-1, qe);	//idx
						tS = 0;
						if(qs<bd)
							while(tS<=tE && gData0[tS].start<bd)tS++;		//exclude 
						for(i=tE; i>=tS; i--){ 
							if(gData0[i].end>qs){
								nols++;
								hitmap[jj][gData0[i].idx]++;
							}
						} 
					}
				}	
			}
			m++;
		}
	}
	//-----------------------------------------------------	
    return nols;
}

int64_t getOverlaps0_m(int64_t **hitmap)
{	//single pass sweep
	int i, j, ii, jj, ichr, n1, mm, m=0;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if(m%1000==0)
				printf("m0 %i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData0);					
				gData0 = malloc(tmpi*sizeof(gdata0_t));
				fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
				//in-bin search: O(n)
				for(j=0;j<tmpi;j++){
					qe = gData0[j].end;
					qs = gData0[j].start;
					jj = gData0[j].idx;			
					if(qs>=bd){
						hitmap[jj][jj]++;
						i=j+1;
					}
					else{//skip duplications
						i=j+1;
						while(i<tmpi && gData0[i].start<bd)
							i++;
					}
					while(i<tmpi && gData0[i].start<qe){
						ii = gData0[i].idx;
						if(ii>jj)hitmap[jj][ii]++;
						else if(jj>ii)hitmap[ii][jj]++;
						else hitmap[ii][jj]+=2;
						i++;
					}
				}
			}			
			m++;
		}
	}
    return nols;	
}

//using AIList: no decomp
int64_t getOverlaps0_m1(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, ii, jj, ichr, n1, mm, m=0;	
	int32_t tE, tS, tmpi, tmpi1, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			tmpi1=tmpi-1;
			if(m%1000==0)
				printf("%i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData0);					
				gData0 = malloc(tmpi*sizeof(gdata0_t));
				fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
				//construct ailist--------------------------------
				maxE = malloc(tmpi*sizeof(int32_t));
				tmax = gData0[0].end;
				for(i=0;i<tmpi;i++){
					if(gData0[i].end>tmax)tmax = gData0[i].end;
					maxE[i]=tmax;
				}
				j=0;
				while(j<tmpi){
					qe = gData0[j].end;
					qs = gData0[j].start;
					//mm = 1;
					if(qe>gData0[0].start){					
						jj = gData0[j].idx;	
						//mm+=j;
						//while(mm<tmpi && gData0[mm].end==qe && gData0[mm].start==qs && gData0[mm].idx==jj)
						//	mm++;
						//mm-=j;	
						tS = j;
						if(qs<bd)
							while(tS<tmpi && gData0[tS].start<bd)tS++;		//exclude 
						i = bSearch0(gData0, tS, tmpi1, qe);	//idx
						while(i>=tS && maxE[i]>qs){
							if(gData0[i].end>qs){
								ii = gData0[i].idx;
								if(ii>jj)hitmap[jj][ii]++;//=mm;
								else if(jj>ii) hitmap[ii][jj]++;//=mm; 
								else hitmap[ii][jj]+=2;//!!!	
							}
							i--;
						} 
					}
					j++;//=mm;
				}
				free(maxE);	
			}			
			m++;
		}
	}
    return nols;
}

//using AIList: decomp
int64_t getOverlaps0_m2(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, ii, jj, k, ichr, n1, m=0, mm;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax, rs, re;
	int64_t nols = 0;
	int32_t *maxE;
	int nc=1, lenC[MAXC], idxC[MAXC];		//components
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if((m++)%1000==0) printf("--%s\t%i\t%i\t%i\n", IGD->cName[ichr], n1, m, tmpi);			
			if(tmpi>0){	//ichr==12 && n1>1000 && 						
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData0);					
				gData0 = malloc(tmpi*sizeof(gdata0_t));
				fread(gData0, sizeof(gdata0_t)*tmpi, 1, fP);
				//construct ailist-----------------------------------------------------------------
				maxE = malloc(tmpi*sizeof(int32_t));				
				construct0(gData0, tmpi, &nc, idxC, lenC, maxE, 20);	
				//printf("%i\t%i\t%i\t %i\n", ichr, n1, tmpi, nc);				
				j = 0;	
				while(j<tmpi){				
					qe = gData0[j].end, qs = gData0[j].start, jj = gData0[j].idx;
					mm=j+1;
					while(mm<tmpi && gData0[mm].end==qe && gData0[mm].start==qs && gData0[mm].idx==jj)
						mm++;
					mm-=j;						
					//for each component-----------------------------------------------------------
					for(k=0;k<nc;k++){
						rs = MAX(idxC[k],j), re = idxC[k]+lenC[k];			
						if(rs<re && qe>gData0[rs].start){		
							if(qs<bd)
								while(rs<re && gData0[rs].start<bd)rs++;		//exclude 
							i = bSearch0(gData0, rs, re-1, qe);				//idx
							while(i>=rs && maxE[i]>qs){
								if(gData0[i].end>qs){
									ii = gData0[i].idx;
									if(ii>jj)hitmap[jj][ii]+=mm;
									else if(jj>ii)hitmap[ii][jj]+=mm; 		//nols++; 
									else hitmap[ii][jj]+=mm+mm;   	
								}
								i--;
							} 
						}
					}
					j+=mm;
				}
				free(maxE);	
			}
		}
	}
    return nols;
}

//=============================================================================
//--------------------------for gdata_t----------------------------------------
int32_t get_overlaps_f(char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{   //FJoin for test
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
		if(qe>gData[0].start){//sorted by start
			i=0;				
			while(i<tmpi && gData[i].end<=qs)
				i++;	// 	
			while(i<tmpi && gData[i].start<qe){
				if(gData[i].end>qs)
					hits[gData[i].idx]++;
				i++;
			}
		}
		if(n2>n1){									//n2>n1
			int32_t bd = IGD->nbp*(n1+1);			//only keep the first		
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
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
					if(qe>gData[0].start){//bSearch!=-1							
						i = 0;
						while(i<tmpi && gData[i].start<bd)
							i++;		//qs<bd
						while(i<tmpi && gData[i].end<=qs)
							i++;//		
						while(i<tmpi && gData[i].start<qe){
							if(gData[i].end>qs)
								hits[gData[i].idx]++;
							i++;
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

/*
int32_t get_overlaps(char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{   //for gdat0_t
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
				if(gData[i].end>qs){
					//nols++;
					hits[gData[i].idx]++;
				} 
			}
		}
		if(n2>n1){									//n2>n1
			int32_t bd = IGD->nbp*(n1+1);			//only keep the first		
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
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
					if(qe>gData[0].start){//bSearch!=-1		
						if(tmpi<16){
							tE = tmpi-1;
							while(gData[tE].start>=qe)tE--;
						}
						else
							tE = bSearch(gData, 0, tmpi-1, qe);	//idx					
						tS = 0;
						while(tS<tmpi && gData[tS].start<bd)tS++;		//qs<bd
						for(i=tE; i>=tS; i--){ 
							if(gData[i].end>qs){		
								//nols++;
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
}*/


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

int64_t getOverlaps_m0(int64_t **hitmap, int32_t v)
{	//load igd tile one by one
	int i, j, jj, ichr, n1, m=0;	//define boundary!
	int32_t tE, tS, tmpi, bd, qe, qs;
	int64_t nols = 0;
	if(v==0)v=-1;
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;
			tmpi = IGD->nCnt[ichr][n1];			
			//if(m%1000==0)
				printf("%i\t%i\t%i\n", n1, m, tmpi);		
			if(tmpi>0){							
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				for(j=0;j<tmpi;j++){
					if(gData[j].value>=v){
						qe = gData[j].end;
						qs = gData[j].start;
						if(qe>gData[0].start){
							jj = gData[j].idx;
							if(tmpi<16){
								tE = tmpi-1;
								while(gData[tE].start>=qe)
									tE--;
							}
							else
								tE = bSearch(gData, 0, tmpi-1, qe);	//idx
							tS = 0;
							if(qs<bd)
								while(tS<=tE && gData[tS].start<bd)tS++;		//exclude 
							for(i=tE; i>=tS; i--){ 
								if(gData[i].end>qs && gData[i].value>=v){
									nols++;
									hitmap[jj][gData[i].idx]++;
								}
							} 
						}
					}
				}	
			}
			m++;
		}
	}
	//-----------------------------------------------------	
    return nols;
}

int64_t getOverlaps_m(int64_t **hitmap)
{
	int i, j, ii, jj, ichr, n1, mm, m=0;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if(m%1000==0)
				printf("m %i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//in-bin search: O(n)
				for(j=0;j<tmpi;j++){
					qe = gData[j].end;
					qs = gData[j].start;
					jj = gData[j].idx;
					if(qs>=bd){
						hitmap[jj][jj]++;
						i=j+1;
					}
					else{//skip duplications
						i=j+1;
						while(i<tmpi && gData[i].start<bd)
							i++;
					}
					while(i<tmpi && gData[i].start<qe){
						ii = gData[i].idx;
						hitmap[jj][ii]++;
						hitmap[ii][jj]++;						
						//if(ii>jj)hitmap[jj][ii]++;
						//else if(jj>ii)hitmap[ii][jj]++;
						//else hitmap[ii][jj]+=2;
						i++;
					}
				}
			}			
			m++;
		}
	}
    return nols;	
}

int64_t getOverlaps_m_v(int64_t **hitmap, int32_t v)
{
	int i, j, ii, jj, ichr, n1, mm, m=0;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if(m%1000==0)
				printf("m_v %i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//in-bin search: O(n)
				for(j=0;j<tmpi;j++){
					if(gData[j].value>=v){
						qe = gData[j].end;
						qs = gData[j].start;
						jj = gData[j].idx;	
						if(qs>=bd){
							hitmap[jj][jj]++;
							i=j+1;
						}
						else{//skip duplications
							i=j+1;
							while(i<tmpi && gData[i].start<bd)
								i++;
						}
						while(i<tmpi && gData[i].start<qe){
							if(gData[i].value>=v){
								ii = gData[i].idx;
								hitmap[jj][ii]++;
								hitmap[ii][jj]++;
							}
							i++;
						}						
					}
				}
			}			
			m++;
		}
	}
    return nols;	
}

int64_t getOverlaps0_m_x(int64_t **hitmap, int32_t x)
{	//extend all regions then self mapping==>subtract m_v 's results
	int i, j, k, j1, j2, ii, tmpi, ix, iy, ichr, n1, jj, m=0;	//define boundary!
	int32_t tE, tS, len1, len2, lenG, lenH, lenT, bd, bd1, rs, re, qe, qs, tmax;
	int64_t nols = 0;
	gdata_t *g1=NULL, *g2=NULL, *gG=NULL, *gH=NULL;		
	for(i=0;i<IGD->nFiles;i++)
		IGD->finfo[i].nr = 1;
	for(ichr=0; ichr<IGD->nCtg; ichr++){
		int32_t nT = IGD->nTile[ichr];
		for(n1=0; n1<nT; n1++){	
			bd = IGD->nbp*n1;
			bd1 = bd+IGD->nbp;
			//1. setup tile data relay
			if(n1>0 && n1<nT-1){		//relay g2->g1, load next as g2
				if(len1>0){
					free(g1), g1=NULL;
				}
				if(len2>0){
					g1 = g2, len1  = len2, g2=NULL;
				}
				else len1 = 0;
				len2 = IGD->nCnt[ichr][n1+1];
				if(len2>0){
					fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
					g2 = malloc(len2*sizeof(gdata0_t));
					fread(g2, sizeof(gdata0_t)*len2, 1, fP);
				}					
			}
			else if(n1==0){				//load first 2
				len1 = IGD->nCnt[ichr][n1];
				if(len1>0){
					fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);								
					g1 = malloc(len1*sizeof(gdata0_t));
					fread(g1, sizeof(gdata0_t)*len1, 1, fP);
				}				
				len2 = IGD->nCnt[ichr][n1+1];
				if(len2>0){
					fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
					g2 = malloc(len2*sizeof(gdata0_t));
					fread(g2, sizeof(gdata0_t)*len2, 1, fP);
				}
				lenH=0, lenG=0;	
			}
			else{						
				if(len1>0){
					free(g1), len1=0, g1=NULL;
				}
				if(len2>0){
					g1 = g2, len1 = len2, g2 = NULL, len2=0;
				}
			}
			
			//2. build gG, gH: relay gHead	
			//printf("m_v_x %i\t%i\t%i\t%i\t%i\n", n1, len1, len2, lenG, lenH);						
			if(gG!=NULL){
				free(gG), lenG=0, gG=NULL;
			}
			if(g1!=NULL){	
				gG = malloc((len1+len2+lenH)*sizeof(gdata0_t));
				lenG = 0;
				if(gH!=NULL){
					memcpy(&gG[0], &gH[0], lenH*sizeof(gdata0_t));
					lenG = lenH;
					free(gH);
					gH=NULL, lenH=0;
				}
				gH = malloc(len1*sizeof(gdata0_t));	
				lenH = 0;	
				k = lenG; j=0;
				for(i=0;i<len1;i++){
					gG[k].start = g1[i].start-x;
					gG[k].end   = g1[i].end+x;
					gG[k].idx   = g1[i].idx;
					if(gG[k].end>bd1 && g1[i].end<=bd1){
						memcpy(&gH[j], &gG[k], sizeof(gdata0_t));					
						j++;
					}
					k++;	
				}	
				lenH = j;	
				lenG = k;//len1+lenH
			}
			else{
				if(gH!=NULL){	//NULL--0
					free(gH), lenH=0, gH=NULL;
				}
			}
			if(gG!=NULL && g2!=NULL){
				i = 0, k = lenG;
				while(i<len2 && g2[i].start-x<bd1){
					if(g2[i].start>=bd1){
						gG[k].start = g2[i].start-x;
						gG[k].end   = g2[i].end+x;
						gG[k].idx   = g2[i].idx;
						k++;
					}
					i++;
				}
				lenG=k;
			}
			
			//3. map gt1	
			/*tmpi = lenG;
			if(m%1000==0)
				printf("m_x %i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){								
				//in-bin search: O(n)
				for(j=0;j<tmpi;j++){
					qe = gG[j].end;
					qs = gG[j].start;
					jj = gG[j].idx;	
					if(qs>=bd){
						hitmap[jj][jj]++;
						i=j+1;
					}
					else{//skip duplications
						i=j+1;
						while(i<tmpi && gG[i].start<bd)
							i++;
					}
					while(i<tmpi && gG[i].start<qe){
						ii = gG[i].idx;
						hitmap[jj][ii]++;
						hitmap[ii][jj]++;
						i++;
					}					
				}
			}*/
			
			if(m%1000==0)
				printf("m_x %i\t%i\t%i\t%i\n", n1, m, len1, lenG);	
			//----------------------------------------------
			//
			if(lenG>0 && len1>0){
				ix=0, iy=0;	//ix: g1; iy: gG
				while(ix<len1 && iy<lenG){
					if(g1[ix].end>gG[iy].start){
						qe = g1[ix].end;
						qs = g1[ix].start;
						jj = g1[ix].idx;
						i = iy;
						while(i<lenG && qe>gG[i].start){
							if(qs<gG[i].end){
								if(gG[i].start>=bd || qs>=bd){	
									ii = gG[i].idx;					
									hitmap[ii][jj]++;
								}
							}
							i++;
						}	
						//move iy
						while(iy<lenG && gG[iy].end<=qs)
							iy++;		
					}
					ix++;
				}
			}	

			//---------------------------------------------
			
			tmpi = len1;
			if(tmpi>0){								
				//in-bin search: O(n)
				for(j=0;j<tmpi;j++){
					qe = g1[j].end;
					qs = g1[j].start;
					jj = g1[j].idx;	
					IGD->finfo[jj].nr++;
					if(qs>=bd){
						hitmap[jj][jj]--;
						i=j+1;
					}
					else{//skip duplications
						i=j+1;
						while(i<tmpi && g1[i].start<bd)
							i++;
					}
					while(i<tmpi && g1[i].start<qe){
						ii = g1[i].idx;
						hitmap[jj][ii]--;
						hitmap[ii][jj]--;
						i++;
					}
				}
			}			
			m++;
		}//for n1
		if(g1!=NULL){
			free(g1), len1=0, g1=NULL;
		}
		if(g2!=NULL){	
			free(g2), len2=0, g2=NULL;
		}
		if(gG!=NULL){
			free(gG), lenG = 0, gG = NULL;
		}
		if(gH!=NULL){
			free(gH), lenH = 0, gH = NULL;
		}
	}//ichr
    return nols;    
}

int64_t getOverlaps_m_v_x(int64_t **hitmap, int32_t v, int32_t x)
{	//extend all regions then self mapping==>subtract m_v 's results
	int i, j, k, ix, iy, j1, j2, ii, tmpi, ichr, n1, jj, m=0;	//define boundary!
	int32_t tE, tS, len1, len2, lenG, lenH, lenT, bd, bd1, rs, re, qe, qs, tmax;
	int64_t nols = 0;
	gdata_t *g1=NULL, *g2=NULL, *gG=NULL, *gH=NULL;		
	for(i=0;i<IGD->nFiles;i++)
		IGD->finfo[i].nr = 1;
	for(ichr=0; ichr<IGD->nCtg; ichr++){
		int32_t nT = IGD->nTile[ichr];
		for(n1=0; n1<nT; n1++){	
			bd = IGD->nbp*n1;
			bd1 = bd+IGD->nbp;
			//1. setup tile data relay
			if(n1>0 && n1<nT-1){		//relay g2->g1, load next as g2
				if(len1>0){
					free(g1), g1=NULL;
				}
				if(len2>0){
					g1 = g2, len1  = len2, g2=NULL;
				}
				else len1 = 0;
				len2 = IGD->nCnt[ichr][n1+1];
				if(len2>0){
					fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
					g2 = malloc(len2*sizeof(gdata_t));
					fread(g2, sizeof(gdata_t)*len2, 1, fP);
				}					
			}
			else if(n1==0){				//load first 2
				len1 = IGD->nCnt[ichr][n1];
				if(len1>0){
					fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);								
					g1 = malloc(len1*sizeof(gdata_t));
					fread(g1, sizeof(gdata_t)*len1, 1, fP);
				}				
				len2 = IGD->nCnt[ichr][n1+1];
				if(len2>0){
					fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
					g2 = malloc(len2*sizeof(gdata_t));
					fread(g2, sizeof(gdata_t)*len2, 1, fP);
				}
				lenH=0, lenG=0;	
			}
			else{						
				if(len1>0){
					free(g1), len1=0, g1=NULL;
				}
				if(len2>0){
					g1 = g2, len1 = len2, g2 = NULL, len2=0;
				}
			}
			
			//2. build gG, gH: relay gHead	
			//printf("m_v_x %i\t%i\t%i\t%i\t%i\n", n1, len1, len2, lenG, lenH);						
			if(gG!=NULL){
				free(gG), lenG=0, gG=NULL;
			}
			if(g1!=NULL){	
				gG = malloc((len1+len2+lenH)*sizeof(gdata_t));
				lenG = 0;
				if(gH!=NULL){
					memcpy(&gG[0], &gH[0], lenH*sizeof(gdata_t));
					lenG = lenH;
					free(gH);
					gH=NULL, lenH=0;
				}
				gH = malloc(len1*sizeof(gdata_t));	
				lenH = 0;	
				k = lenG; j=0;
				for(i=0;i<len1;i++){
					gG[k].start = g1[i].start-x;
					gG[k].end   = g1[i].end+x;
					gG[k].idx   = g1[i].idx;
					gG[k].value = g1[i].value;
					if(gG[k].end>bd1 && g1[i].end<=bd1){
						memcpy(&gH[j], &gG[k], sizeof(gdata_t));					
						j++;
					}
					k++;	
				}	
				lenH = j;	
				lenG = k;//len1+lenH
			}
			else{
				if(gH!=NULL){	//NULL--0
					free(gH), lenH=0, gH=NULL;
				}
			}
			if(gG!=NULL && g2!=NULL){
				i = 0, k = lenG;
				while(i<len2 && g2[i].start-x<bd1){
					if(g2[i].start>=bd1){
						gG[k].start = g2[i].start-x;
						gG[k].end   = g2[i].end+x;
						gG[k].idx   = g2[i].idx;
						gG[k].value = g2[i].value;	
						k++;
					}
					i++;
				}
				lenG=k;
			}
			
			//3. map gG*g1 then subtract g1*g1	
			if(m%1000==0)
				printf("m_v_x %i\t%i\t%i\t%i\n", n1, m, len1, lenG);	
			//----------------------------------------------
			//
			if(lenG>0 && len1>0){
				ix=0, iy=0;	//ix: g1; iy: gG
				while(ix<len1 && iy<lenG){
					if(g1[ix].end>gG[iy].start && g1[ix].value>=v){
						qe = g1[ix].end;
						qs = g1[ix].start;
						jj = g1[ix].idx;
						i = iy;
						while(i<lenG && qe>gG[i].start){
							if(qs<gG[i].end && gG[i].value>=v){
								if(gG[i].start>=bd || qs>=bd){	
									ii = gG[i].idx;					
									hitmap[ii][jj]++;
								}
							}
							i++;
						}	
						//move iy
						while(iy<lenG && (gG[iy].end<=qs || gG[iy].value<v))
							iy++;		
					}
					ix++;
				}
			}	

			//---------------------------------------------
			tmpi = len1;
			if(tmpi>0){								
				//in-bin search: O(n)
				for(j=0;j<tmpi;j++){
					if(g1[j].value>=v){
						qe = g1[j].end;
						qs = g1[j].start;
						jj = g1[j].idx;
						IGD->finfo[jj].nr++;			
						if(qs>=bd){
							hitmap[jj][jj]--;
							i=j+1;
						}
						else{//skip duplications
							i=j+1;
							while(i<tmpi && g1[i].start<bd)
								i++;
						}
						while(i<tmpi && g1[i].start<qe){
							if(g1[i].value>=v){
								ii = g1[i].idx;
								hitmap[jj][ii]--;
								hitmap[ii][jj]--;
							}
							i++;
						}	
					}
				}
			}			
			m++;
		}//for n1
		if(g1!=NULL){
			free(g1), len1=0, g1=NULL;
		}
		if(g2!=NULL){	
			free(g2), len2=0, g2=NULL;
		}
		if(gG!=NULL){
			free(gG), lenG = 0, gG = NULL;
		}
		if(gH!=NULL){
			free(gH), lenH = 0, gH = NULL;
		}
	}//ichr
    return nols;    
}

//using fjoin:
int64_t getOverlaps_m3(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, ii, jj, ichr, n1, mm, done, m=0;	
	int32_t tE, tS, tmpi, tmpi1,bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *Wx, *Wy, nx, ny, ix, iy, ti, tt;//windows	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			tmpi1=tmpi-1;
			if(m%1000==0)
				printf("m3 %i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//construct ailist---------------------------------------------------
				Wx = malloc(tmpi*sizeof(int32_t));
				Wy = malloc(tmpi*sizeof(int32_t));
				nx=0, ny=0;	
				//-------------------------------------------------------------------
				ix=0, iy=0;
				done = 0;
				while(done==0){
					if(ix<tmpi && gData[ix].start<=gData[iy].start){
						//compare ix with Wy
						tt=0; 				//number of removals
						for(j=0;j<ny;j++){
							ti = Wy[j];
							if(gData[ti].end>gData[ix].start){//Wy's start<ix.end
								if(gData[ix].start>=bd || gData[ti].start>=bd)
									hitmap[gData[ti].idx][gData[ix].idx]++;
							}
							else{			//remove Wy[j]: costy; move j+
								Wy[j]=-1;	//marker
								tt++;	
							}
						}
						//re-arrange Wy
						if(tt>0){
							tt=0;
							for(j=0;j<ny;j++){
								if(Wy[j]==-1)
									tt++;	
								else if(tt>0)
									Wy[j-tt]=Wy[j];
							}
							ny-=tt;
						}
						//add ix to Wx if ix intersectable iy
						if(gData[ix].end>gData[iy].start){
							Wx[nx] = ix;
							nx++;
						}
						ix++;
					}
					else if(iy<tmpi){
						//compare iy with Wx
						tt=0; 				//number of removals
						for(j=0;j<nx;j++){
							ti = Wx[j];
							if(gData[ti].end>gData[iy].start){//Wx's start<iy.end
								if(gData[iy].start>=bd || gData[ti].start>=bd)							
									hitmap[gData[iy].idx][gData[ti].idx]++;
							}
							else{			//remove Wy[j]: costy; move j+
								Wx[j]=-1;	//marker
								tt++;	
							}
						}
						//re-arrange Wx
						if(tt>0){
							tt=0;
							for(j=0;j<nx;j++){
								if(Wx[j]==-1)
									tt++;	
								else if(tt>0)
									Wx[j-tt]=Wx[j];
							}
							nx-=tt;
						}
						//add iy to Wy
						if(gData[iy].end>gData[ix].start){
							Wy[ny] = iy;
							ny++;
						}
						iy++;		
					}
					else
						done=1;
				}
				free(Wx);
				free(Wy);	
			}			
			m++;
		}
	}
    return nols;
}

//using FJoin:
int64_t getOverlaps_mf(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, ii, jj, ichr, n1, mm, done, m=0;	
	int32_t tE, tS, tmpi, bd, qe, qs, rs, re, tmax;
	int64_t nols = 0;
	int32_t ix, iy;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if(m%1000==0)
				printf("mf %i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//-------------------------------------------------------------------
				ix=0, iy=0;
				while(ix<tmpi && iy<tmpi){
					//a------y as q, x as r---------------
					//scan ix.. to x.end>qs
					qs = gData[iy].start;
					while(ix<tmpi && gData[ix].end<=qs)ix++;
					//check ix against iy..
					qe = gData[ix].end;
					qs = gData[ix].start;
					j = gData[ix].idx;
					i = iy;
					while(i<tmpi && qe>gData[i].start){
						if(qs<gData[i].end){
							if(gData[i].start>=bd || qs>=bd)							
								hitmap[gData[i].idx][j]++;
						}
						i++;
					}
					//done with ix
					ix++;
					if(ix==tmpi)break;
					
					//b------x as q, y as r---------------
					//scan and move iy.. to y.end>qs
					qs = gData[ix].start;
					while(iy<tmpi && gData[iy].end<=qs)iy++;
					//check iy against ix..
					qe = gData[iy].end;
					qs = gData[iy].start;
					j = gData[iy].idx;
					i = ix;
					while(i<tmpi && qe>gData[i].start){
						if(qs<gData[i].end){
							if(gData[i].start>=bd || qs>=bd)							
								hitmap[j][gData[i].idx]++;
						}
						i++;
					}
					//done with iy
					iy++;
				}
			}			
			m++;
		}
	}
    return nols;
}


//using AIList: no decomp
int64_t getOverlaps_m1(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, ii, jj, ichr, n1, mm, m=0;	
	int32_t tE, tS, tmpi, tmpi1,bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			tmpi1=tmpi-1;
			if(m%1000==0)
				printf("%i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//construct ailist---------------------------------------------------
				maxE = malloc(tmpi*sizeof(int32_t));
				tmax = gData[0].end;
				for(i=0;i<tmpi;i++){
					if(gData[i].end>tmax)tmax = gData[i].end;
					maxE[i]=tmax;
				}
				j=0;
				while(j<tmpi){
					qe = gData[j].end;
					qs = gData[j].start;
					if(qe>gData[0].start){					
						jj = gData[j].idx;	
						tS = j;
						if(qs<bd){
							while(tS<tmpi && gData[tS].start<bd)tS++;	//exclude 
						}
						if(tS==j){
							hitmap[jj][jj]++;	
							tS++;
						}
						i = bSearch(gData, tS, tmpi1, qe);	
						while(i>=tS && maxE[i]>qs){
							if(gData[i].end>qs){
								ii = gData[i].idx;
								hitmap[jj][ii]++;
								hitmap[ii][jj]++; 
							}
							i--;
						} 
					}
					j++;
				}
				free(maxE);	
			}			
			m++;
		}
	}
    return nols;
}

//for comparison
int64_t getOverlaps_m1a(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, ii, jj, ichr, n1, mm, m=0;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if(m%1000==0)
				printf("%i\t%i\t%i\n", n1, m, tmpi);			
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
				j=0;
				while(j<tmpi){
					qe = gData[j].end;
					qs = gData[j].start;
					if(qe>gData[0].start){					
						jj = gData[j].idx;		
						tS = 0;
						if(qs<bd)
							while(tS<tmpi && gData[tS].start<bd)tS++;		//exclude 
						i = bSearch(gData, tS, tmpi-1, qe);	//idx
						while(i>=tS && maxE[i]>qs){
							if(gData[i].end>qs){
								hitmap[jj][gData[i].idx]++;
							}
							i--;
						} 
					}
					j++;
				}
				free(maxE);	
			}			
			m++;
		}
	}
    return nols;
}

//using AIList: decomp
int64_t getOverlaps_m2(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, ii, jj, k, ichr, n1, m=0, mm;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax, rs, re;
	int64_t nols = 0;
	int32_t *maxE;
	int nc=1, lenC[MAXC], idxC[MAXC];		//components
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if((m++)%1000==0) printf("--%s\t%i\t%i\t%i\n", IGD->cName[ichr], n1, m, tmpi);			
			if(tmpi>0){	//ichr==12 && n1>1000 && 						
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//construct ailist-----------------------------------------------------------------
				maxE = malloc(tmpi*sizeof(int32_t));				
				construct(gData, tmpi, &nc, idxC, lenC, maxE, 20);	
				//printf("%i\t%i\t%i\t %i\n", ichr, n1, tmpi, nc);				
				j = 0;	
				while(j<tmpi){				
					qe = gData[j].end, qs = gData[j].start, jj = gData[j].idx;					
					//for each component-----------------------------------------------------------
					for(k=0;k<nc;k++){
						rs = MAX(idxC[k],j), re = idxC[k]+lenC[k];			
						if(rs<re && qe>gData[rs].start){		
							if(qs<bd){
								while(rs<re && gData[rs].start<bd)rs++;		//exclude 
							}
							if(rs==j){
								hitmap[jj][jj]++;
								rs++;
							}
							i = bSearch(gData, rs, re-1, qe);				//idx
							while(i>=rs && maxE[i]>qs){
								if(gData[i].end>qs){
									ii = gData[i].idx;
									hitmap[jj][ii]++;
									hitmap[ii][jj]++;  	
								}
								i--;
							} 
						}
					}
					j++;
				}
				free(maxE);	
			}
		}
	}
    return nols;
}

//using AIList: decomp && simple case of split
int64_t getOverlaps_m2a(int64_t **hitmap)
{	//load igd tile one by one
	int i, j, ii, jj, k, ichr, n1, m=0, mm;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax, rs, re;
	int64_t nols = 0;
	int32_t *maxE;
	int nc=1, lenC[MAXC], idxC[MAXC];		//components
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if((m++)%1000==0) printf("%i\t%i\t%i\n", n1, m, tmpi);			
			if(tmpi>0){							
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//construct ailist--------------------------------
				maxE = malloc(tmpi*sizeof(int32_t));				
				construct(gData, tmpi, &nc, idxC, lenC, maxE, 20);					
				j = 0;	
				while(j<tmpi){				
					qe = gData[j].end, qs = gData[j].start, jj = gData[j].idx;
					mm=j+1;
					while(mm<tmpi && gData[mm].end==qe && gData[mm].start==qs && gData[mm].idx==jj)
						mm++;
					mm-=j;						
					//for each component:
					for(k=0;k<nc;k++){
						rs = MAX(idxC[k],j), re = idxC[k]+lenC[k];			
						if(rs<re && qe>gData[rs].start){		
							if(qs<bd)
								while(rs<re && gData[rs].start<bd)rs++;		//exclude 
							i = bSearch(gData, rs, re-1, qe);				//idx
							while(i>=rs && maxE[i]>qs){
								if(gData[i].end>qs){
									ii = gData[i].idx;
									hitmap[jj][ii]+=mm;
									hitmap[ii][jj]+=mm; 		   	
								}
								i--;
							} 
						}
					}
					j+=mm;
				}
				free(maxE);	
			}
		}
	}
    return nols;
}

//using AIList: no decomp
int64_t getOverlaps_m1_v(int64_t **hitmap, int32_t v)
{	//load igd tile one by one
	//define: >=v
	int i, j, ii, jj, ichr, n1, m=0;	
	int32_t tE, tS, tmpi, tmpi1, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	//for(i=0;i<IGD->nFiles;i++)
	//	IGD->finfo[i].nr = 1;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			tmpi1=tmpi-1;
			if(m%1000==0)
				printf("%i\t%i\t%i\n", n1, m, tmpi);			
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
				j=0;
				while(j<tmpi){				
					if(gData[j].value>=v){
						qe = gData[j].end;
						qs = gData[j].start;						
						jj = gData[j].idx;						
						//IGD->finfo[jj].nr++;
						if(qe>gData[0].start){	
							tS = j;
							if(qs<bd){
								while(tS<tmpi && gData[tS].start<bd)tS++;	//exclude 
							}
							if(tS==j){
								hitmap[jj][jj]++;	
								tS++;
							}
							i = bSearch(gData, tS, tmpi1, qe);	
							while(i>=tS && maxE[i]>qs){
								if(gData[i].end>qs && gData[i].value>=v){
									ii = gData[i].idx;
									hitmap[jj][ii]++;
									hitmap[ii][jj]++;
								}
								i--;
							} 						
						}
					}
					j++;
				}
				free(maxE);	
			}			
			m++;
		}
	}
    return nols;
}

//using AIList: decomp
int64_t getOverlaps_m2_v(int64_t **hitmap, int32_t v)
{	//load igd tile one by one
	int i, j, ii, jj, k, ichr, n1, m=0, mm;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax, rs, re;
	int64_t nols = 0;
	int32_t *maxE;
	int nc=1, lenC[MAXC], idxC[MAXC];		//components
	//for(i=0;i<IGD->nFiles;i++)
	//	IGD->finfo[i].nr = 1;
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
			if((m++)%1000==0) printf("--%s\t%i\t%i\t%i\n", IGD->cName[ichr], n1, m, tmpi);			
			if(tmpi>0){	//ichr==12 && n1>1000 && 						
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				//construct ailist-----------------------------------------------------------------
				maxE = malloc(tmpi*sizeof(int32_t));				
				construct(gData, tmpi, &nc, idxC, lenC, maxE, 20);	
				//printf("%i\t%i\t%i\t %i\n", ichr, n1, tmpi, nc);				
				j = 0;	
				while(j<tmpi){	
					if(gData[j].value>=v){			
						qe = gData[j].end, qs = gData[j].start, jj = gData[j].idx;	
						//IGD->finfo[jj].nr++;				
						//for each component-----------------------------------------------------------
						for(k=0;k<nc;k++){
							rs = MAX(idxC[k],j), re = idxC[k]+lenC[k];			
							if(rs<re && qe>gData[rs].start){		
								if(qs<bd){
									while(rs<re && gData[rs].start<bd)rs++;		//exclude 
								}
								if(rs==j){
									hitmap[jj][jj]++;
									rs++;
								}	
								i = bSearch(gData, rs, re-1, qe);				//idx
								while(i>=rs && maxE[i]>qs){
									if(gData[i].end>qs && gData[i].value>=v){
										ii = gData[i].idx;
										hitmap[jj][ii]++;
										hitmap[ii][jj]++;	   	
									}
									i--;
								} 
							}
						}
					}
					j++;
				}
				free(maxE);	
			}
		}
	}
    return nols;
}

int64_t getOverlaps_m0_x(int64_t **hitmap, int32_t v, int32_t x)
{	//flanking x	
	int i, j, ichr, n1, jj, m=0;	//define boundary!
	int32_t tE, tS, tmpi, bd, qe, qs;
	int64_t nols = 0;
	if(v==0)v=-1;
	for(i=0;i<IGD->nFiles;i++)
		IGD->finfo[i].nr = 1;	
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
				for(j=0;j<tmpi;j++){
				if(gData[j].value>=v){
					jj = gData[j].idx;
					IGD->finfo[jj].nr++;
					//flanking
					//1.Left: qe>bd (case 1: qs>=bd; case 2: qs<bd)
					qe = gData[j].start;
					if(qe>bd){
						qs = MAX(0,qe-x);						
						if(qs<bd && n1>0){//rare case:cross boundary-->load (n1-1)th tile
							int tmpt = IGD->nCnt[ichr][n1-1];
							fseek(fP, IGD->tIdx[ichr][n1-1], SEEK_SET);								
							gdata_t *gtmp = malloc(tmpt*sizeof(gdata_t));
							fread(gtmp, sizeof(gdata_t)*tmpt, 1, fP);
							tE = bSearch(gtmp, 0, tmpt-1, qe);	//idx
							tS = 0;
							for(i=tE; i>=tS;i--){ 
								if(gtmp[i].end>qs && gtmp[i].value>=v){//left boundary
									nols++;
									hitmap[jj][gtmp[i].idx]++;
								}
							} 
							free(gtmp);
						}
						if(qe>gData[0].start){
							tE = bSearch(gData, 0, tmpi-1, qe);	//idx
							tS = 0;
							if(qs<bd)
								while(tS<tmpi && gData[tS].start<bd)tS++;		//exclude 
							for(i=tE; i>=tS;i--){ 
								if(gData[i].end>qs && gData[i].value>=v){// && gData[i].start>bd){//left boundary
									nols++;
									hitmap[jj][gData[i].idx]++;
								}
							} 
						}						

					}//qe<=bd belongs to the previous tile
					
					//2.Right:qs<bd1 (case 1: qe<=bd1; case2--cross bd1: qe>bd1)
					qs = gData[j].end;
					int32_t bd1 = bd+IGD->nbp;
					if(qs<bd1){						
						qe = qs+x;
						if(qe>gData[0].start){
							tE = bSearch(gData, 0, tmpi-1, qe);	//idx
							tS = 0;
							for(i=tE; i>=tS;i--){ 
								if(gData[i].end>qs && gData[i].value>=v){// && gData[i].start>bd){//left boundary
									nols++;
									hitmap[jj][gData[i].idx]++;
								}
							} 
						}		
						if(qe>bd1){//cross boundary: n1+1<nTiles?
							int tmpt = IGD->nCnt[ichr][n1+1];
							fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
							gdata_t *gtmp = malloc(tmpt*sizeof(gdata_t));
							fread(gtmp, sizeof(gdata_t)*tmpt, 1, fP);
							if(tmpt<16){
								tE = tmpt-1;
								while(gtmp[tE].start>=qe)tE--;
							}
							else
								tE = bSearch(gtmp, 0, tmpt-1, qe);	//idx
							tS = 0;
							if(qs<bd1)
								while(tS<tmpt && gtmp[tS].start<bd1)tS++;		//exclude 
							for(i=tE; i>=tS;i--){ 
								if(gtmp[i].end>qs && gtmp[i].value>=v){
									nols++;
									hitmap[jj][gtmp[i].idx]++;
								}
							} 
							free(gtmp);
						}
					}
				}//gData[j].value>v
				}	
			}
		}
	}
    return nols;
}

int64_t getOverlaps_m1_x(int64_t **hitmap, int32_t v, int32_t x)
{	//flanking x	
	int i, j, ichr, n1, jj, m=0;	//define boundary!
	int32_t tE, tS, len1, len2, len3, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *E1, *E2, *E3;
	gdata_t *g1, *g2, *g3;
	for(i=0;i<IGD->nFiles;i++)
		IGD->finfo[i].nr = 1;
	for(ichr=0; ichr<IGD->nCtg; ichr++){
		int32_t nT = IGD->nTile[ichr];
		//--------------------------------------------			
		for(n1=0; n1<nT; n1++){
			//work on 3 consective tiles
			if(n1<nT-1 && n1>0){//general
				if(len1>0){
					free(g1), free(E1), len1=0;
				}
				if(len2>0){
					g1=g2, E1=E2, len1 = len2;
					g2=NULL, E2=NULL, len2=0;
				}
				if(len3>0){
					g2=g3, E2=E3, len2 = len3;
					g3=NULL, E3=NULL, len3=0;
				}
				//load g3, e3
				len3 = IGD->nCnt[ichr][n1+1];
				if(len3>0){
					fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
					g3 = malloc(len3*sizeof(gdata_t));
					fread(g3, sizeof(gdata_t)*len3, 1, fP);
					E3 = malloc(len3*sizeof(int32_t));
					tmax = g3[0].end;
					for(i=0;i<len3;i++){
						if(g3[i].end>tmax)tmax = g3[i].end;
						E3[i]=tmax;
					}
				}				
			}			
			else if(n1==0){//first
				len1=0;
				len2 = IGD->nCnt[ichr][n1];
				len3 = IGD->nCnt[ichr][n1+1];
				if(len2>0){
					fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);								
					g2 = malloc(len2*sizeof(gdata_t));
					fread(g2, sizeof(gdata_t)*len2, 1, fP);
					E2 = malloc(len2*sizeof(int32_t));
					tmax = g2[0].end;
					for(i=0;i<len2;i++){
						if(g2[i].end>tmax)tmax = g2[i].end;
						E2[i]=tmax;
					}
				}	
				if(len3>0){
					fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
					g3 = malloc(len3*sizeof(gdata_t));
					fread(g3, sizeof(gdata_t)*len3, 1, fP);
					E3 = malloc(len3*sizeof(int32_t));
					tmax = g3[0].end;
					for(i=0;i<len3;i++){
						if(g3[i].end>tmax)tmax=g3[i].end;
						E3[i]=tmax;
					}			
				}			
			}			
			else{//n1==nT-1: last
				if(len1>0){
					free(g1), free(E1), len1=0;
				}
				if(len2>0){
					g1=g2, E1=E2, len1 = len2;
					g2=NULL, E2=NULL, len2=0;
				}
				if(len3>0){
					g2=g3, E2=E3, len2 = len3;
					g3=NULL, E3=NULL;
				}
				len3 = 0;			
			}
			//---------------------------------------------------------
			//printf("%i\t%i\t%i\t%i\n", n1, len1,len2,len3);
			bd = IGD->nbp*n1;
			if(m%1000==0)
				printf("%i\n", m);			
			m++;
			if(len2>0){												
				for(j=0;j<len2;j++){
				if(g2[j].value>=v){
					jj = g2[j].idx;
					IGD->finfo[jj].nr++;					
					//flanking
					//1.Left: qe>bd (case 1: qs>=bd; case 2: qs<bd)
					qe = g2[j].start;
					if(qe>bd){
						qs = MAX(0,qe-x);						
						if(qs<bd && len1>0){
							i = bSearch(g1, 0, len1-1, qe);	//idx
							while(i>=0 && E1[i]>qs){ 
								if(g1[i].end>qs && g1[i].value>=v){//left boundary
									nols++;
									hitmap[jj][g1[i].idx]++;
								}
								i--;
							} 
						}
						if(qe>g2[0].start){
							i = bSearch(g2, 0, len2-1, qe);	//idx
							tS = 0;
							if(qs<bd)
								while(g2[tS].start<bd)tS++;		//exclude 
							while(i>=tS && E2[i]>qs){ 
								if(g2[i].end>qs && g2[i].value>=v){
									nols++;
									hitmap[jj][g2[i].idx]++;
								}
								i--;
							} 
						}	
					}
					
					//2.Right:qs<bd1 (case 1: qe<=bd1; case2--cross bd1: qe>bd1)
					qs = g2[j].end;
					int32_t bd1 = bd+IGD->nbp;
					if(qs<bd1){						
						qe = qs+x;
						if(qe>g2[0].start){
							i = bSearch(g2, 0, len2-1, qe);	//idx
							while(i>=0 && E2[i]>qs){ 
								if(g2[i].end>qs && g2[i].value>=v){// && gData[i].start>bd){//left boundary
									nols++;
									hitmap[jj][g2[i].idx]++;
								}
								i--;
							} 
						}		
						if(qe>bd1 && len3>0){
							i = bSearch(g3, 0, len3-1, qe);	//idx
							tS = 0;
							if(qs<bd1)
								while(g3[tS].start<bd1)tS++;		//?tS<tmpt-1
							while(i>=tS && E3[i]>qs){ 
								if(g3[i].end>qs && g3[i].value>=v){// && gData[i].start>bd){//left boundary
									nols++;
									hitmap[jj][g3[i].idx]++;
								}
								i--;
							} 
						}
					}
				}//.value>v
				}//j
			}
		}//n1
		if(len1>0){
			free(g1), free(E1), len1=0;
		}
		if(len2>0){
			free(g2), free(E2), len2=0;
		}
	}//ichr
    return nols;
}

int64_t getOverlaps_m2_x(int64_t **hitmap, int32_t v, int32_t x)
{	//tbd
	int i, j, k, ichr, n1, jj, m=0;	//define boundary!
	int32_t tE, tS, tmpi, bd, qe, qs, rs, re, len1, len2, len3;
	int64_t nols = 0;
	int nc1=1, lenC1[MAXC], idxC1[MAXC];		//components
	int nc2=1, lenC2[MAXC], idxC2[MAXC];		//components
	int nc3=1, lenC3[MAXC], idxC3[MAXC];		//components		
	int32_t *E1, *E2, *E3;
	gdata_t *g1, *g2, *g3;
	if(v==0) v = -1;
	for(i=0;i<IGD->nFiles;i++)
		IGD->finfo[i].nr = 1;
	for(ichr=0; ichr<IGD->nCtg; ichr++){
		int32_t nT = IGD->nTile[ichr];
		//--------------------------------------------			
		for(n1=0; n1<nT; n1++){
			if(n1<nT-1 && n1>0){
				if(len1>0){
					free(g1), free(E1), len1=0, nc1=0;
				}
				if(len2>0){
					g1=g2, E1=E2, len1 = len2, nc1 = nc2;					
					g2=NULL, E2=NULL, len2=0;	nc2=0;
					memcpy(lenC1, lenC2, sizeof(lenC2));
					memcpy(idxC1, idxC2, sizeof(idxC2));
				}
				if(len3>0){
					g2=g3, E2=E3, len2 = len3, nc2=nc3;
					g3=NULL, E3=NULL, len3=0;	nc3=0;
					memcpy(lenC2, lenC3, sizeof(lenC3));
					memcpy(idxC2, idxC3, sizeof(idxC3));
				}				
				//load g3, e3
				len3 = IGD->nCnt[ichr][n1+1];
				if(len3>0){
					fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
					g3 = malloc(len3*sizeof(gdata_t));
					fread(g3, sizeof(gdata_t)*len3, 1, fP);
					E3 = malloc(len3*sizeof(int32_t));
					construct(g3, len3, &nc3, idxC3, lenC3, E3, 20);	
				}				
			}			
			else if(n1==0){
				len1=0, nc1=0, nc2=0, nc3=0;
				len2 = IGD->nCnt[ichr][n1];
				len3 = IGD->nCnt[ichr][n1+1];
				if(len2>0){
					fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);								
					g2 = malloc(len2*sizeof(gdata_t));
					fread(g2, sizeof(gdata_t)*len2, 1, fP);
					E2 = malloc(len2*sizeof(int32_t));
					construct(g2, len2, &nc2, idxC2, lenC2, E2, 20);	
				}		
				if(len3>0){
					fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
					g3 = malloc(len3*sizeof(gdata_t));
					fread(g3, sizeof(gdata_t)*len3, 1, fP);
					E3 = malloc(len3*sizeof(int32_t));
					construct(g3, len3, &nc3, idxC3, lenC3, E3, 20);	
				}				
			}			
			else{//n1==nT-1
				if(len1>0){
					free(g1), free(E1), len1=0, nc1=0;
				}
				if(len2>0){
					g1=g2, E1=E2, len1 = len2, nc1 = nc2;					
					g2=NULL, E2=NULL, len2=0;	nc2=0;
					memcpy(lenC1, lenC2, sizeof(lenC2));
					memcpy(idxC1, idxC2, sizeof(idxC2));
				}
				if(len3>0){
					g2=g3, E2=E3, len2 = len3, nc2=nc3;
					g3=NULL, E3=NULL;
					memcpy(lenC2, lenC3, sizeof(lenC3));
					memcpy(idxC2, idxC3, sizeof(idxC3));
				}									
				//load g3, e3
				len3 = 0, nc3 = 0;		
			}
			//---------------------------------------------------------
			bd = IGD->nbp*n1;
			if(m%1000==0)
				printf("%i\n", m);			
			m++;
			if(len2>0){												
				for(j=0;j<len2;j++){
				if(g2[j].value>=v){
					jj = g2[j].idx;	
					IGD->finfo[jj].nr++;
					//flanking
					//1.Left: qe>bd (case 1: qs>=bd; case 2: qs<bd)
					qe = g2[j].start;
					if(qe>bd){
						qs = MAX(0,qe-x);						
						if(qs<bd && len1>0){
							//for each component:
							for(k=0;k<nc1;k++){
								rs = idxC1[k], re = rs+lenC1[k];		
								if(qe>g1[rs].start){	
									i = bSearch(g1, rs, re-1, qe);	//idx
									while(i>=rs && E1[i]>qs){
										if(g1[i].end>qs && g1[i].value>=v){
											nols++;
											hitmap[jj][g1[i].idx]++;
										}
										i--;
									} 
								}
							}						
						}
						if(qe>g2[0].start){						
							//for each component:
							for(k=0;k<nc2;k++){
								rs = idxC2[k], re = rs+lenC2[k];		
								if(qe>g2[rs].start){						
									tS = rs;
									if(qs<bd)
										while(tS<re && g2[tS].start<bd)tS++;	//exclude 
									i = bSearch(g2, tS, re-1, qe);	//idx
									while(i>=tS && E2[i]>qs){
										if(g2[i].end>qs && g2[i].value>=v){
											nols++;
											hitmap[jj][g2[i].idx]++;
										}
										i--;
									} 
								}
							}						 
						}	
					}//qe<=bd belongs to the previous tile
					
					//2.Right:qs<bd1 (case 1: qe<=bd1; case2--cross bd1: qe>bd1)
					qs = g2[j].end;
					int32_t bd1 = bd+IGD->nbp;
					if(qs<bd1){					
						qe = qs+x;
						if(qe>g2[0].start){	
							for(k=0;k<nc2;k++){
								rs = idxC2[k], re = rs+lenC2[k];
								if(qe>g2[rs].start){
									i = bSearch(g2, rs, re-1, qe);	//idx
									while(i>=rs && E2[i]>qs){
										if(g2[i].end>qs && g2[i].value>=v){
											nols++;
											hitmap[jj][g2[i].idx]++;
										}
										i--;
									} 
								}
							}												 
						}		
						if(qe>bd1 && len3>0){//cross boundary: n1+1<nTiles?
							for(k=0;k<nc3;k++){
								rs = idxC3[k], re = rs+lenC3[k];		
								if(qe>g3[rs].start){						
									tS = rs;
									if(qs<bd1)
										while(tS<re && g3[tS].start<bd1)tS++;		//exclude 
									i = bSearch(g3, tS, re-1, qe);	//idx
									while(i>=tS && E3[i]>qs){
										if(g3[i].end>qs && g3[i].value>=v){
											nols++;
											hitmap[jj][g3[i].idx]++;
										}
										i--;
									} 
								}
							}
						}
					}
				}//.value>v
				}//j
			}
		}				
		if(len1>0){
			free(g1), free(E1), len1=0;
		}
		if(len2>0){
			free(g2), free(E2), len2=0;
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
    int i, j, checking=0, mode=-1, mt=0, ext=0, xlen=0,  mv=0, mx=0, ichr, k;
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
        else if(strcmp(argv[i], "-o")==0){
            if(i+1<argc)
                strcpy(out, argv[i+1]);
        }   
        else if(strcmp(argv[i], "-c")==0){
            checking = 1;
        } 
        else if(strcmp(argv[i], "-m")==0){
            mode = 0;				//default
            if(i+1<argc)
                mt = atoi(argv[i+1]);//map mode type: return 0 if not digit
        }  
        else if(strcmp(argv[i], "-x")==0){
        	mx = 1;
            if(i+1<argc){
                ext = 1; //find cofactors
                xlen = atoi(argv[i+1]);//symmetric: .end+, .start-
            }
        }                               
    }  
     
    //----------------------------------------------------------
	fP = fopen(igdName, "rb");				//share
    if(mode==0){	//mapping: gType==1
    	int64_t **hitmap = malloc(nfiles*sizeof(int64_t*));
    	for(i=0;i<nfiles;i++)
    		hitmap[i] = calloc(nfiles, sizeof(int64_t));
    	double **fmap = malloc(nfiles*sizeof(double*));
    	for(i=0;i<nfiles;i++)
    		fmap[i] = calloc(nfiles, sizeof(double));    		
    	if(IGD->gType==0){
 			if(mt==1)
				getOverlaps0_m1(hitmap);
			else if(mt==2)
				getOverlaps0_m2(hitmap);
			else //if(mt==4)
				getOverlaps0_m(hitmap);				
		 	//else
			//	getOverlaps0_m0(hitmap);   
    	}
    	else{   		
			if(mt==0){	//default			
				if(mv==0 && mx==0)
					getOverlaps_m(hitmap);
				else if(mx==0)
					getOverlaps_m_v(hitmap, v);
				else
					getOverlaps_m_v_x(hitmap, v, xlen);

			}
			else if(mt==1){
				if(mv==0)
					getOverlaps_m1(hitmap);
				else if(mx>0)
					getOverlaps_m1_x(hitmap, v, xlen);
				else
					getOverlaps_m1_v(hitmap, v);
			}
		 	else if(mt==2){
				if(mv==0)
					getOverlaps_m2(hitmap);
				else if(mx>0)
					getOverlaps_m2_x(hitmap, v, xlen);
				else
					getOverlaps_m2_v(hitmap, v);
			} 
			else if(mt==3)
				getOverlaps_m3(hitmap);
			else if(mt==4)
				getOverlaps_mf(hitmap);	
			else{
				if(mx>0)
						getOverlaps_m0_x(hitmap, v, xlen);
				else
					getOverlaps_m0(hitmap, v);			
			}
			//symmetry
			//if(mt<3){
			//	for(j=0;j<nfiles;j++){
			//		for(i=j+1;i<nfiles;i++){
			//			hitmap[i][j]=hitmap[j][i];//uint32_t	
			//		}
			//	}
			//}
    	}
    	//calculate J-index
    	if(mx>0){
			for(j=0;j<nfiles;j++){
				fmap[j][j] = 0.0;
				for(i=j+1;i<nfiles;i++){
					if(hitmap[j][i]>0){//ratio of increased intersections over the number of sites
						fmap[j][i] = (double)hitmap[j][i]/(double)(IGD->finfo[i].nr+IGD->finfo[j].nr);
						fmap[i][j] = fmap[j][i];
					}
				}
			}
    	}
    	else{
			for(j=0;j<nfiles;j++){
				fmap[j][j] = 1.0;
				for(i=j+1;i<nfiles;i++){
					if(hitmap[j][i]>0){
						fmap[j][i] = (double)hitmap[j][i]/(double)(hitmap[j][j]+hitmap[i][i]-hitmap[j][i]);
						fmap[i][j] = fmap[j][i];
					}
				}
			}    	
    	}	   		
    	FILE *fp;
		if(strlen(out)<2)strcpy(out,"Hitmap");
		fp = fopen(out, "w");
	    if(fp==NULL)
	        printf("Can't open file %s\n", out);
	    else{
	        fprintf(fp, "%u\t%u\t%u\n", nfiles, nfiles, v);
	        for(i=0;i<nfiles;i++){
	            for(j=0;j<nfiles;j++)
	                fprintf(fp, "%lld\t", (long long)hitmap[i][j]); 
	            fprintf(fp, "\n");
	        } 
	        fclose(fp);
	    }  
	    //ratio: 
	    char out1[64];
	    strcpy(out1, out);
	    strcat(out1, "_j"); 
	    fp = fopen(out1, "w");
	    if(fp==NULL)
	        printf("Can't open file %s\n", out1);
	    else{
	        fprintf(fp, "%u\t%u\t%u\n", nfiles, nfiles, v);
	        for(i=0;i<nfiles;i++){
	            for(j=0;j<nfiles;j++)
	                fprintf(fp, "%.8f\t", fmap[i][j]); 
	            fprintf(fp, "\n");
	        } 
	        fclose(fp);
	    }    
	
    	for(i=0;i<nfiles;i++){
    		free(hitmap[i]);
    		free(fmap[i]);
    	}
    	free(fmap);
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

