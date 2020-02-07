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
	if(mode==1){//for a query dataset (file)  
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

