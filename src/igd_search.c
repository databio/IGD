//=====================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//
//time ./igd_search Test110000.bed /media/john/Extra/ucsc_igd/ucsc.igd
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

int32_t get_overlaps_r(char *chrm, int32_t qs, int32_t qe, int32_t *hits)
{   //no need to store every overlaps, only get the number of hits
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;
	int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tS, tmpi, mTile = IGD->nTile[ichr]-1;
	int32_t nols = 0;
	if(n1>mTile) 
		return 0;
	n2 = MIN(n2, mTile);	
	if(n2==n1){										//most often
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
			if(qs<gData[tmpi-1].end){				//sorted by end
				if(tmpi<16){
					tS = 0;
					while(gData[tS].end<=qs)
						tS++;
				}
				else
					tS = bSearch(gData, tmpi, qs);	//idx
				for(i=tS; i<tmpi; i++){
					if(gData[i].start< qe){
						nols++;
						hits[gData[i].idx]++;
					} 
				}
			}
		}
	}
	else{											//n2>n1
		int32_t bd = IGD->nbp*(n1+1);
		for(j=n1; j<n2; j++){
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
				if(qs<gData[tmpi-1].end){
					if(tmpi<16){
						tS = 0;
						while(gData[tS].end<=qs)
							tS++;
					}
					else
						tS = bSearch(gData, tmpi, qs);	//idx
					for(i=tS; i<tmpi;i++){ 
						if(gData[i].start<qe && gData[i].end <= bd){
							nols++;
							hits[gData[i].idx]++;
						}
					} 
				}
			}
			bd+=IGD->nbp;		
		}	
		tmpi = IGD->nCnt[ichr][n2];
		if(tmpi>0){
			free(gData);					
			gData = malloc(tmpi*sizeof(gdata_t));
			fread(gData, sizeof(gdata_t) * tmpi, 1, fP);
			preIdx=n2;
			if(qs<gData[tmpi-1].end){
				if(tmpi<16){
					tS = 0;
					while(gData[tS].end<=qs)
						tS++;
				}
				else
					tS = bSearch(gData, tmpi, qs);	//idx	
				for(i=tS; i<tmpi; i++){
					if(gData[i].start<qe){
					 	nols++;
					 	hits[gData[i].idx]++;
					}
				} 
			}
		}
	}
	//-----------------------------------------------------	
    return nols;
}

int32_t get_overlaps_r1(char *chrm, int32_t qs, int32_t qe, int32_t *hits)
{   //no need to store every overlaps, only get the number of hits
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;
	int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tS, tmpi, mTile = IGD->nTile[ichr]-1;
	int32_t nols = 0;
	if(n1>mTile) 
		return 0;
	n2 = MIN(n2, mTile);	
	if(n2==n1){										//most often
		tmpi = IGD->nCnt[ichr][n1];
		if(tmpi>0){
			if(n1!=preIdx || ichr!=preChr){
				fseek(fP, IGD->tIdx[ichr][n1], SEEK_SET);			
				free(gData1);					
				gData1 = malloc(tmpi*sizeof(gdata1_t));
				fread(gData1, sizeof(gdata1_t)*tmpi, 1, fP);
				preIdx = n1;
				preChr = ichr;
			}
			if(qs<gData1[tmpi-1].end){				//sorted by end
				if(tmpi<16){
					tS = 0;
					while(gData1[tS].end<=qs)
						tS++;
				}
				else
					tS = bSearch1(gData1, tmpi, qs);	//idx
				for(i=tS; i<tmpi; i++){
					if(gData1[i].start< qe){
						nols++;
						hits[gData1[i].idx]++;
					} 
				}
			}
		}
	}
	else{											//n2>n1
		int32_t bd = IGD->nbp*(n1+1);
		for(j=n1; j<n2; j++){
			tmpi = IGD->nCnt[ichr][j];
			if(tmpi>0){
				if(j!=preIdx || ichr!=preChr){
					fseek(fP, IGD->tIdx[ichr][j], SEEK_SET);			
					free(gData1);					
					gData1 = malloc(tmpi*sizeof(gdata1_t));
					fread(gData1, sizeof(gdata1_t)*tmpi, 1, fP);
					preIdx = j;
					preChr = ichr;
				}
				if(qs<gData1[tmpi-1].end){
					if(tmpi<16){
						tS = 0;
						while(gData1[tS].end<=qs)
							tS++;
					}
					else
						tS = bSearch1(gData1, tmpi, qs);	//idx
					for(i=tS; i<tmpi;i++){ 
						if(gData1[i].start<qe && gData1[i].end <= bd){
							nols++;
							hits[gData1[i].idx]++;
						}
					} 
				}
			}
			bd+=IGD->nbp;		
		}	
		tmpi = IGD->nCnt[ichr][n2];
		if(tmpi>0){
			free(gData1);					
			gData1 = malloc(tmpi*sizeof(gdata1_t));
			fread(gData1, sizeof(gdata1_t) * tmpi, 1, fP);
			preIdx=n2;
			if(qs<gData1[tmpi-1].end){
				if(tmpi<16){
					tS = 0;
					while(gData1[tS].end<=qs)
						tS++;
				}
				else
					tS = bSearch1(gData1, tmpi, qs);	//idx	
				for(i=tS; i<tmpi; i++){
					if(gData[i].start<qe){
					 	nols++;
					 	hits[gData1[i].idx]++;
					}
				} 
			}
		}
	}
	//-----------------------------------------------------	
    return nols;
}

int64_t get_overlaps(char *qFile, int32_t *hits)
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
			nl = get_overlaps_r(chrm, st, en, hits);
			ols += nl;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);	
	return ols;
}

int64_t get_overlaps1(char *qFile, int32_t *hits)
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
			nl = get_overlaps_r1(chrm, st, en, hits);
			ols += nl;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);	
	return ols;
}

int64_t get_overlaps_v(char *qFile, int32_t *hits, int32_t v)
{	//tbd
	int64_t nols=0;
	return nols;
}

int32_t get_overlaps_r2(char *chrm, int32_t qs, int32_t qe, int32_t *hits)
{
	int ichr = get_id(chrm);
	if(ichr<0)
		return 0;
	int i, j, n1 = qs/IGD->nbp, n2 = (qe-1)/IGD->nbp;	//define boundary!
	int32_t tS, tmpi, mTile = IGD->nTile[ichr]-1;
	int32_t nols = 0;
	if(n1>mTile) 
		return 0;
	n2 = MIN(n2, mTile);	
	if(n2==n1){										//most often
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
			if(qs<gData[tmpi-1].end){				//sorted by end
				if(tmpi<16){
					tS = 0;
					while(gData[tS].end<=qs)
						tS++;
				}
				else
					tS = bSearch(gData, tmpi, qs);	//idx
				for(i=tS; i<tmpi; i++){
					if(gData[i].start< qe){
						nols++;
						hits[gData[i].idx]++;
					} 
				}
			}
		}
	}
	else{											//n2>n1
		int32_t bd = IGD->nbp*(n1+1);
		for(j=n1; j<n2; j++){
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
				if(qs<gData[tmpi-1].end){
					if(tmpi<16){
						tS = 0;
						while(gData[tS].end<=qs)
							tS++;
					}
					else
						tS = bSearch(gData, tmpi, qs);	//idx
					for(i=tS; i<tmpi;i++){ 
						if(gData[i].start<qe && gData[i].end <= bd){
							nols++;
							hits[gData[i].idx]++;
						}
					} 
				}
			}
			bd+=IGD->nbp;		
		}	
		tmpi = IGD->nCnt[ichr][n2];
		if(tmpi>0){
			free(gData);					
			gData = malloc(tmpi*sizeof(gdata_t));
			fread(gData, sizeof(gdata_t) * tmpi, 1, fP);
			preIdx=n2;
			if(qs<gData[tmpi-1].end){
				if(tmpi<16){
					tS = 0;
					while(gData[tS].end<=qs)
						tS++;
				}
				else
					tS = bSearch(gData, tmpi, qs);	//idx	
				for(i=tS; i<tmpi; i++){
					if(gData[i].start<qe){
					 	nols++;
					 	hits[gData[i].idx]++;
					}
				} 
			}
		}
	}
	//-----------------------------------------------------	
    return nols;
}

int64_t get_overlaps_m(char *qFile, int32_t *hits, int32_t v)
{	//tbd
	int64_t nols=0;
	return nols;
}

int64_t get_overlaps_m_x(char *qFile, int32_t *hits, int32_t v, int32_t x)
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
    int checking=0, mode=-1, ext=0, xlen=0,  ichr, k;
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
    int32_t *hits = calloc(IGD->nFiles, sizeof(int32_t));    
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
        /*else if(strcmp(argv[i], "-v")==0){
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
        }  
        else if(strcmp(argv[i], "-x")==0){
            if(i+1<argc){
                ext = 1; //find cofactors
                xlen = atoi(argv[i+1]);//symmetric: .end+, .start-
            }
        } */                               
    }  
    //check if fils exist

      
    //----------------------------------------------------------
	fP = fopen(igdName, "rb");				//share
    if(mode==0){
        //if(ext==0)
        //    search_self(igdName, v, out);
        //else
        //    search_self_ext(igdName, v, out, xlen);
    }
    else if(mode==1){
    	if(IGD->gType==0)
    		get_overlaps(qfName, hits);
    	else
    		get_overlaps1(qfName, hits);
    	printf("index\t File_name\t number of regions\t number of hits\n"); 
    	int64_t total = 0;       
        for(int i=0;i<IGD->nFiles;i++){
            printf("%i %s %i %i\n", i, IGD->finfo[i].fileName, IGD->finfo[i].nr, hits[i]); 
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
    	if(IGD->gType==0)
    		ols = get_overlaps_r(chrm, qs, qe, hits); 
    	else
    		ols = get_overlaps_r1(chrm, qs, qe, hits);
    	printf("index\t File_name\t number of regions\t number of hits\n");        
        for(int i=0;i<IGD->nFiles;i++)
            printf("%i %s %i %i\n", i, IGD->finfo[i].fileName, IGD->finfo[i].nr, hits[i]);
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

