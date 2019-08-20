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
"             -x <extended length on query ends>\n"
"             -m heatmap with igd database itself\n"
"             -c display all intersects\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
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
                if(j1<cLen1) memcpy(&L2[len++], &L0[t], sizeof(gdata_t));
                else memcpy(&L1[k++], &L0[t], sizeof(gdata_t));                 
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

int32_t get_overlaps(char *chrm, int32_t qs, int32_t qe, int32_t *hits)
{   //no need to store every overlaps, only get the number of hits for each file
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
					nols++;
					hits[gData[i].idx]++;
				} 
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
				
				if(qe>gData[0].start){
					if(tmpi<16){
						tE = tmpi-1;
						while(gData[tE].start>=qe)tE--;
					}
					else
						tE = bSearch(gData, 0, tmpi-1, qe);	//idx
					tS = 0;
					while(gData[tS].start<bd)tS++;		//exclude 
					for(i=tE; i>=tS; i--){ 
						if(gData[i].end>qs){// && gData[i].start>bd){//left boundary		
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

//for gData with v
int32_t get_overlaps_v(char *chrm, int32_t qs, int32_t qe, int32_t v, int32_t *hits)
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
				if(gData[i].end>qs && gData[i].value>v){
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
					while(gData[tS].start<bd)tS++;		//exclude <left boundary
					for(i=tE; i>=tS;i--){ 
						if(gData[i].end>qs && gData[i].value>v){
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

int64_t getOverlaps(char *qFile, int32_t *hits)
{	//for gdata_t
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

int64_t getOverlaps_v(char *qFile, int32_t *hits, int32_t v)
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

int64_t getOverlaps_m0(uint32_t **hitmap)
{	//load igd tile one by one
	int i, j, jj, ichr, n1;	//define boundary!
	int32_t tE, tS, tmpi, bd, qe, qs;
	int64_t nols = 0;
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;
			tmpi = IGD->nCnt[ichr][n1];	
			if(tmpi>0){							
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				for(j=0;j<tmpi;j++){
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
							if(gData[i].end>qs){// && gData[i].start>bd){//left boundary
								nols++;
								hitmap[jj][gData[i].idx]++;
							}
						} 
					}
				}	
			}
		}
	}
	//-----------------------------------------------------	
    return nols;
}

//using AIList: no decomp
int64_t getOverlaps_m1(uint32_t **hitmap)
{	//load igd tile one by one
	int i, j, jj, ichr, n1;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax;
	int64_t nols = 0;
	int32_t *maxE;	
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;		
			tmpi = IGD->nCnt[ichr][n1];
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

//using AIList: decomp
int64_t getOverlaps_m2(uint32_t **hitmap)
{	//load igd tile one by one
	int i, j, jj, k, ichr, n1;	
	int32_t tE, tS, tmpi, bd, qe, qs, tmax, rs, re;
	int64_t nols = 0;
	int32_t *maxE;
	int nc=1, lenC[MAXC], idxC[MAXC];		//components
	int32_t m = 0;	
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
				construct(gData, tmpi, &nc, idxC, lenC, maxE, 20);
				/*printf("%i\t%i\t%i\n", ichr, n1, nc);
				for(k=0;k<nc;k++){
					rs = idxC[k], re = rs+lenC[k];
					printf("--%i\t%i\n", rs, re);
					for(j=rs;j<re;j++)
						printf(":%i\t%i\t%u\n",gData[j].start, gData[j].end, maxE[j]);
				}*/
				for(j=0;j<tmpi;j++){
					qe = gData[j].end;
					qs = gData[j].start;
					jj = gData[j].idx;							
					//for each component:
					for(k=0;k<nc;k++){
						rs = idxC[k], re = rs+lenC[k];		
						if(qe>gData[rs].start){						
							tS = rs;
							if(qs<bd)
								while(tS<re && gData[tS].start<bd)tS++;		//exclude 
							i = bSearch(gData, tS, re-1, qe);	//idx
							while(i>=tS && maxE[i]>qs){
								if(gData[i].end>qs){
									//nols++;
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

int64_t getOverlaps_m0_x(uint32_t **hitmap, int32_t x)
{	//flanking x	
	int i, j, ichr, n1;	//define boundary!
	int32_t tE, tS, tmpi, bd, qe, qs;
	int64_t nols = 0;
	for(ichr=0; ichr<IGD->nCtg; ichr++){	
		for(n1=0; n1<IGD->nTile[ichr]; n1++){
			bd = IGD->nbp*n1;
			tmpi = IGD->nCnt[ichr][n1];
			if(tmpi>0){								
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData);					
				gData = malloc(tmpi*sizeof(gdata_t));
				fread(gData, sizeof(gdata_t)*tmpi, 1, fP);
				for(j=0;j<tmpi;j++){
					//flanking
					//1.Left: qe>bd (case 1: qs>=bd; case 2: qs<bd)
					qe = gData[j].start;
					if(qe>bd){
						qs = MAX(0,qe-x);						
						if(qs<bd && n1>0){//rare case:cross boundary-->load (n1-1)th tile
							int tmpt = IGD->nCnt[ichr][n1-1];
							fseek(fP, IGD->tIdx[ichr][n1-1], SEEK_SET);								
							gdata_t *gtmp = malloc(tmpi*sizeof(gdata_t));
							fread(gtmp, sizeof(gdata_t)*tmpt, 1, fP);
							if(tmpt<16){
								tE = tmpt-1;
								while(gtmp[tE].start>=qe)
									tE--;
							}
							else
								tE = bSearch(gtmp, 0, tmpt-1, qe);	//idx
							tS = 0;
							for(i=tE; i>=tS;i--){ 
								if(gtmp[i].end>qs){//left boundary
									nols++;
									hitmap[gData[j].idx][gData[i].idx]++;
								}
							} 
							free(gtmp);
						}
						if(qe>gData[0].start){
							if(tmpi<16){
								tE = tmpi-1;
								while(gData[tE].start>=qe)
									tE--;
							}
							else
								tE = bSearch(gData, 0, tmpi-1, qe);	//idx
							tS = 0;
							if(qs<bd)
								while(gData[tS].start<bd)tS++;		//exclude 
							for(i=tE; i>=tS;i--){ 
								if(gData[i].end>qs){// && gData[i].start>bd){//left boundary
									nols++;
									hitmap[gData[j].idx][gData[i].idx]++;
								}
							} 
						}						

					}//qe<=bd belongs to the previous tile
					
					//2.Right:qs<bd1 (case 1: qe<=bd1; case2--cross bd1: qe>bd1)
					qs = gData[j].end;
					int32_t bd1 = bd+IGD->nbp;
					if(qs<bd1){
						qs = MAX(0,qe-x);						
						qe = qs+x;
						if(qe>gData[0].start){
							if(tmpi<16){
								tE = 0;
								while(gData[tE].start>=qe)
									tE--;
							}
							else
								tE = bSearch(gData, 0, tmpi-1, qe);	//idx
							tS = 0;
							for(i=tE; i>=tS;i--){ 
								if(gData[i].end>qs){// && gData[i].start>bd){//left boundary
									nols++;
									hitmap[gData[j].idx][gData[i].idx]++;
								}
							} 
						}		
						if(qe>bd1){//cross boundary: n1+1<nTiles?
							int tmpt = IGD->nCnt[ichr][n1+1];
							fseek(fP, IGD->tIdx[ichr][n1+1], SEEK_SET);								
							gdata_t *gtmp = malloc(tmpi*sizeof(gdata_t));
							fread(gtmp, sizeof(gdata_t)*tmpt, 1, fP);
							if(tmpt<16){
								tE = tmpt-1;
								while(gtmp[tE].start>=qe)tE--;
							}
							else
								tE = bSearch(gtmp, 0, tmpt-1, qe);	//idx
							tS = 0;//qs<bd1
							while(gtmp[tS].start<bd1)tS++;		//exclude 
							for(i=tE; i>=tS;i--){ 
								if(gtmp[i].end>qs){// && gData[i].start>bd){//left boundary
									nols++;
									hitmap[gData[j].idx][gData[i].idx]++;
								}
							} 
							free(gtmp);
						}
					}
				}	
			}
		}
	}
    return nols;
}

int64_t getOverlaps_m1_x(uint32_t **hitmap, int32_t x)
{	//tbd
	int64_t nols=0;
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
    int checking=0, mode=-1, mt=0, ext=0, xlen=0,  ichr, k;
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
    int32_t *hits = calloc(nfiles, sizeof(int32_t));  
    //-----------------------------------------------------
    for(int i=3; i<argc; i++){
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
        else if(strcmp(argv[i], "-v")==0){
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
            mode = 0;
            if(i+1<argc)
                mt = atoi(argv[i+1]);//mode type
        }  
        else if(strcmp(argv[i], "-x")==0){
            if(i+1<argc){
                ext = 1; //find cofactors
                xlen = atoi(argv[i+1]);//symmetric: .end+, .start-
            }
        }                               
    }  
    //check if fils exist

      
    //----------------------------------------------------------
	fP = fopen(igdName, "rb");				//share
    if(mode==0){	//mapping
    	uint32_t **hitmap = malloc(nfiles*sizeof(uint32_t*));
    	for(int i=0;i<nfiles;i++)
    		hitmap[i] = calloc(nfiles, sizeof(uint32_t));
    	if(mt==0)
    		getOverlaps_m0(hitmap);
    	else if(mt==1)
    		getOverlaps_m1(hitmap);
     	else if(mt==2)
    		getOverlaps_m2(hitmap); 
    	FILE *fp;
		if(strlen(out)<2)strcpy(out,"Hitmap");
		fp = fopen(out, "w");
	    if(fp==NULL)
	        printf("Can't open file %s\n", out);
	    else{
	        fprintf(fp, "%u\t%u\t%u\n", nfiles, nfiles, v);
	        for(int i=0;i<nfiles;i++){
	            for(int j=0;j<nfiles;j++)
	                fprintf(fp, "%u\t", hitmap[i][j]); 
	            fprintf(fp, "\n");
	        } 
	        fclose(fp);
	    }     
	
    	for(int i=0;i<nfiles;i++)
    		free(hitmap[i]);
    	free(hitmap);  		
        //if(ext==0)
        //    search_self(igdName, v, out);
        //else
        //    search_self_ext(igdName, v, out, xlen);
    }
    else if(mode==1){//for a query dataset (file)
    	getOverlaps(qfName, hits);
    	printf("index\t File_name\t number of regions\t number of hits\n"); 
    	int64_t total = 0;       
        for(int i=0;i<nfiles;i++){
        	if(hits[i]>0)
            	printf("%i\t%i\t%i\t%s\n", i, IGD->finfo[i].nr, hits[i], IGD->finfo[i].fileName); 
        	total += hits[i];
        }
        //char *qtype = qfName + strlen(qfName) - 4;    
        //if(strcmp(".bed", qtype)!=0){
        //    printf("%s is not a bed file", qfName);
        //    return EX_OK;
        //}
        //fi = fopen(qfName, "rb");
        //if(!fi){
        //    printf("%s does not exist", qfName);
        //    return EX_OK;
        //}
        //fclose(fi);     
        //search(igdName, qfName, v, out, checking);
        printf("Total: %lld\n", (long long)total);
    }
    else if(mode==2){//mode 2 for a single region
    	ols = get_overlaps(chrm, qs, qe, hits); 
    	printf("index\t File_name\t number of regions\t number of hits\n");        
        for(int i=0;i<nfiles;i++)
            printf("%i\t%i\t%i\t%s\n", i, IGD->finfo[i].nr, hits[i], IGD->finfo[i].fileName);
    }
    else
        return search_help(EX_OK);      

	fclose(fP);
    free(IGD->nTile);
    for(int i=0;i<IGD->nCtg;i++){
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

