//===================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//database intervals sorted by _start: 8/12/2019
//-----------------------------------------------------------------------------------
#include "igd_create.h"
KHASH_MAP_INIT_STR(str, int32_t)
typedef khash_t(str) strhash_t;

int create_help(int exit_code)
{
    fprintf(stderr,
"%s, v%s\n"
"usage:   %s create <input dir> <output dir> <output igd name> [options] \n"
"             -s  <Type of data structure> \n"
"                   0 for [index, start, end]\n"
"                   1 for [index, start, end, value], default\n"
"                   2 single combined BED file, datatype 1\n"
"             -f  (iPath is a file that lists paths of data src files) \n"
"             -b  <Tile size in power of 2 (default 14)> \n"
"             -c  < .BED column as value >=4 (default 4) \n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

//create ucsc igd from gz
void create_igd(char *iPath, char *oPath, char *igdName)
{  	
	//0. Initialize igd
	igd_t *igd = igd_init();
	printf("igd_create 0\n");
	
    //1. Get the files  
    glob_t gResult;
    int rtn = glob(iPath, 0, NULL, &gResult);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath);
        return;
    }
    char** file_ids = gResult.gl_pathv;
    int32_t n_files = gResult.gl_pathc; 
    if(n_files<1)   
        printf("Too few files (add to path /*): %i\n", n_files);    
    int32_t *nr = calloc(n_files, sizeof(int32_t));
    double *avg = calloc(n_files, sizeof(double));
    printf("igd_create 1: %i\n", n_files);        
    //2. Read files
    int nCols=16;
    unsigned char buffer[1024];     
    int32_t i, j, k, ig, i0=0, i1=0, L0=0, L1=1, m, nL; //int64_t? 
    char **splits = malloc((nCols+1)*sizeof(char *));    
    while(i0<n_files){
        //2.1 Start from (i0, L0): read till (i1, L1)
        ig = i0; 
        m = 0;   
		//2.2 Read ~4GB data from files
        while(m==0 && ig<n_files){   	//m>0 defines breaks when reading maxCount      
			//printf("%i, %i, %i, %s\n", i0, ig, nL, file_ids[ig]);
			gzFile fp;
			if ((fp = gzopen(file_ids[ig], "r")) == 0)
				return;                             
		    nL = 0; 
		    if(ig==i0 && L0>0){  		 //pass L0 lines of a big file
		        while(nL<L0 && gzgets(fp, buffer, 1024)!=NULL)
		            nL++;              
		    }      			
			while (m==0 && gzgets(fp, buffer, 1024)!=NULL) {
				str_splits(buffer, &nCols, splits); 
				int32_t  st = atol(splits[1]), en = atol(splits[2]), va = 0;
				if(nCols>4) va = atol(splits[4]);
				igd_add(igd, splits[0], st, en, va, ig);				
				nr[ig]++;
				avg[ig]+=en-st;
				nL++;
                if(igd->total>maxCount){
                    m = 1;
                    i1 = ig;
                    L1 = nL;    		//number of total lines or next line
                }
			}
			gzclose(fp);
			if(m==0) ig++;
		}
        //2.3 Save/append tiles to disc, add cnts tp Cnts 			
        
		igd_saveT(igd, oPath);
        i0 = ig; 
        L0 = L1;
        L1 = 0;        
	} 		
    
	printf("igd_create 2\n");
	
	//3. save _index.tsv: 4 columns--index, filename, nr, avg
    //Also has a header line: 
    char idFile[128];
    char *tchr;   
    sprintf(idFile, "%s%s%s", oPath, igdName, "_index.tsv");    
    FILE *fpi = fopen(idFile, "w");
    if(fpi==NULL)
        printf("Can't open file %s", idFile);     
    fprintf(fpi, "Index\tFile\tNumber of regions\tAvg size\n");    
    for(i=0; i<n_files; i++){
        tchr = strrchr(file_ids[i], '/');
        if(tchr!=NULL)
            tchr += 1;
        else
            tchr = file_ids[i];
        fprintf(fpi, "%i\t%s\t%i\t%f\n", i, tchr, nr[i], avg[i]/nr[i]);     
    }
    fclose(fpi);   
    free(nr);
    free(avg);    
    printf("igd_create 3\n");
    
	//4. Sort tile data and save into single files per ctg
	igd_save(igd, oPath, igdName);	
	globfree(&gResult); 
	printf("igd_create 4\n");	
	free(splits); 
}

//create igd from a file list
void create_igd_f(char *iPath, char *oPath, char *igdName)
{  	
	//0. Initialize igd
	igd_t *igd = igd_init();
	printf("igd_create 0\n");
	
    //1. Get and check the files  
    FILE *fp0 = fopen(iPath, "r");
    if(fp0==NULL)
        printf("Can't open file %s", iPath);
        
	char buf[1024];
    unsigned char buffer[1024]; 	
    int n_files=0;  
    while(fgets(buf, 1024, fp0)!=NULL)
		n_files++;    
    if(n_files<1){   
        printf("Too few files (add to path /*): %i\n", n_files);  
        return;
	}      
    char** file_ids = malloc(n_files*sizeof(char *));
    fseek(fp0, 0, SEEK_SET);   
    //int nCols=16; 
    char *ctg;  
    int32_t st=0, en=0; 
    gzFile fp;
    int ix=0;    
    while(fgets(buf, 1024, fp0)!=NULL){	
    	buf[strcspn(buf, "\n")] = 0;  
        //check if openable and valid
		if ((fp = gzopen(buf, "r"))!=NULL){
			gzgets(fp, buffer, 1024);
			//str_splits(buffer, &nCols, splits);
			ctg = parse_bed(buffer, &st, &en);
			if(ctg){     
				file_ids[ix] = strdup(buf); 
				ix++;
			}
			gzclose(fp);
		}
    }        
    fclose(fp0);
    n_files = ix;

    int32_t *nr = calloc(n_files, sizeof(int32_t));
    double *avg = calloc(n_files, sizeof(double));
    printf("igd_create 1: %i\n", n_files); 
           
    //2. Read files
    int32_t va, i, j, k, ig, i0=0, i1=0, L0=0, L1=1, m, nL; //int64_t?        		
    while(i0<n_files){
        //2.1 Start from (i0, L0): read till (i1, L1)
        ig = i0; 
        m = 0;    
		//2.2 Read ~4GB data from files
        while(m==0 && ig<n_files){   	//m>0 defines breaks when reading maxCount       
			printf("%i, %i, %i, %s\n", i0, ig, nL, file_ids[ig]);
			fp = gzopen(file_ids[ig], "r");                           
		    nL = 0; 
		    if(ig==i0 && L0>0){  		 //pass L0 lines of a big file
		        while(nL<L0 && gzgets(fp, buffer, 1024)!=NULL)
		            nL++;              
		    }      			
			while (m==0 && gzgets(fp, buffer, 1024)!=NULL) {
				ctg = parse_bed(buffer, &st, &en);				
				if(ctg && st>=0 && en<321000000){				
					igd_add(igd, ctg, st, en, va, ig);				
					nr[ig]++;
					avg[ig]+=en-st;
				}
				nL++;
	            if(igd->total>maxCount){
	                m = 1;
	                i1 = ig;
	                L1 = nL;    		//number of total lines or next line
	            }
			}
			gzclose(fp);
			if(m==0) ig++;
		}
        //2.3 Save/append tiles to disc, add cnts tp Cnts 
		//printf("---------------%i\t %i\n", i0, ig);        				         
		igd_saveT(igd, oPath);
        i0 = ig; 
        L0 = L1;
        L1 = 0;        
	} 
   
	printf("igd_create 2\n");
	
	//3. save _index.tsv: 4 columns--index, filename, nr, avg
    //Also has a header line: 
    char idFile[1024];
    char *tchr;   
    sprintf(idFile, "%s%s%s", oPath, igdName, "_index.tsv");    
    FILE *fpi = fopen(idFile, "w");
    if(fpi==NULL)
        printf("Can't open file %s", idFile);     
    fprintf(fpi, "Index\tFile\tNumber of regions\tAvg size\n");    
    for(i=0; i<n_files; i++){
        tchr = strrchr(file_ids[i], '/');
        if(tchr!=NULL)
            tchr += 1;
        else
            tchr = file_ids[i];
        fprintf(fpi, "%i\t%s\t%i\t%f\n", i, tchr, nr[i], avg[i]/nr[i]);     
    }
    fclose(fpi);   
    free(nr);
    free(avg);    
    printf("igd_create 3\n");
    
	//4. Sort tile data and save into single files per ctg
	igd_save(igd, oPath, igdName);	
	printf("igd_create 4\n");
	free(file_ids);	
	//free(splits); 
}

//create ucsc igd from gz
void create_igd0(char *iPath, char *oPath, char *igdName)
{  	
	//0. Initialize igd
	igd0_t *igd = igd0_init();
	printf("igd_create 0\n");
	
    //1. Get the files  
    glob_t gResult;
    int rtn = glob(iPath, 0, NULL, &gResult);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath);
        return;
    }
    char** file_ids = gResult.gl_pathv;
    int32_t n_files = gResult.gl_pathc; 
    if(n_files<1)   
        printf("Too few files (add to path /*): %i\n", n_files);    
    int32_t *nr = calloc(n_files, sizeof(int32_t));
    double *avg = calloc(n_files, sizeof(double));
    printf("igd_create 1: %i\n", n_files);        
    //2. Read files
    int nCols=16;
    unsigned char buffer[256];     
    int32_t i, j, k, ig, i0=0, i1=0, L0=0, L1=1, m, nL; //int64_t? 
	char **splits = malloc((nCols+1)*sizeof(char *));     
    while(i0<n_files){
        //2.1 Start from (i0, L0): read till (i1, L1)
        ig = i0; 
        m = 0;    
               
		//2.2 Read ~4GB data from files
        while(m==0 && ig<n_files){   	//m>0 defines breaks when reading maxCount      
			//printf("%i, %i, %i, %s\n", i0, ig, nL, file_ids[ig]);
			gzFile fp;
			if ((fp = gzopen(file_ids[ig], "r")) == 0)
				return;                             
		    nL = 0; 
		    if(ig==i0 && L0>0){  		 //pass L0 lines of a big file
		        while(nL<L0 && gzgets(fp, buffer, 256)!=NULL)
		            nL++;              
		    }      			
			while (m==0 && gzgets(fp, buffer, 256)!=NULL) {
				str_splits(buffer, &nCols, splits); 
				int32_t  st = atol(splits[1]), en = atol(splits[2]);
				igd0_add(igd, splits[0], st, en, ig);				
				nr[ig]++;
				avg[ig]+=en-st;
				nL++;
                if(igd->total>maxCount){
                    m = 1;
                    i1 = ig;
                    L1 = nL;    		//number of total lines or next line
                }
			}

			gzclose(fp);
			if(m==0) ig++;
		}
        //2.3 Save/append tiles to disc, add cnts tp Cnts 			
		igd0_saveT(igd, oPath);
        i0 = ig; 
        L0 = L1;
        L1 = 0;        
	}  		
	printf("igd_create 2\n");
	
	//3. save _index.tsv: 4 columns--index, filename, nr, avg
    //Also has a header line: 
    char idFile[1024];
    char *tchr;   
    sprintf(idFile, "%s%s%s", oPath, igdName, "_index.tsv");    
    FILE *fpi = fopen(idFile, "w");
    if(fpi==NULL)
        printf("Can't open file %s", idFile);     
    fprintf(fpi, "Index\tFile\tNumber of regions\tAvg size\n");    
    for(i=0; i<n_files; i++){
        tchr = strrchr(file_ids[i], '/');
        if(tchr!=NULL)
            tchr += 1;
        else
            tchr = file_ids[i];
        fprintf(fpi, "%i\t%s\t%i\t%f\n", i, tchr, nr[i], avg[i]/nr[i]);     
    }
    fclose(fpi);   
    free(nr);
    free(avg);    
    printf("igd_create 3\n");
    
	//4. Sort tile data and save into single files per ctg
	igd0_save(igd, oPath, igdName);	
	globfree(&gResult); 
	printf("igd_create 4\n");	
	free(splits);  
}

//create igd from a single .bed.gz Jaspar file with dataset index at 4th column
void create_igd_bed4(char *iPath, char *oPath, char *igdName)
{   //Process line by line: ipath is the file name
    //1. Get file_ids, n_files  
    int nCols=32;
    unsigned char buffer[1024];        
    char **splits = malloc((nCols+1)*sizeof(char *));     
    //-------------------------------------------------------------------------
    printf("igd_create 1\n");    
    igd_t *igd = igd_init(); 
    igd->gType = 1;
    int maxF = 360000;
    int32_t *nr = calloc(maxF, sizeof(int32_t));
    double *avg = calloc(maxF, sizeof(double));
    char** file_ids = malloc(maxF*sizeof(char *));
    //--------------------------------------------
    void *fHash = kh_init(str);
    khint_t k;
    //--------------------------------------------         
    //2. Read file  
    int64_t nL = 0;
    int32_t st, en, n_files = 0; 
    int i, j=0, idx, absent;                          
    gzFile zfile = gzopen(iPath, "r");              
    while(gzgets(zfile, buffer, 1024)!=NULL){ 
        str_splits(buffer, &nCols, splits);	      
        //Set up hash table/find file id ---------
        strhash_t *fH = (strhash_t*)fHash;
        k = kh_put(str, fH, splits[3], &absent);
        if(absent){
        	kh_key(fH, k) = strdup(splits[3]);
        	kh_val(fH, k) = n_files;
        	file_ids[n_files] = strdup(splits[3]);
        	n_files++;
        }
        idx = kh_val(fH, k);
        //----------------------------------------
        st = atol(splits[1]), en = atol(splits[2]); 	
		igd_add(igd, splits[0], st, en, atol(splits[4]), idx);			
		nr[idx]++;
		avg[idx]+=en-st;
		nL++;	
        //Save/append tiles to disc, add cnts tp Cnts 
        if(igd->total>=maxCount){ 
        	j++;
        	printf("--igd_saveT1--%i, %lld\n", j, (long long)igd->total);
			igd_saveT(igd, oPath);
		    nL = 0;        
        }
	}	
	gzclose(zfile);     
	kh_destroy(str, fHash);
	if(nL>0)
		igd_saveT(igd, oPath);  
	printf("igd_create 2\n");
	
	//3. save _index.tsv: 4 columns--index, filename, nr, avg
    //Also has a header line: 
    char idFile[1024];
    char *tchr;   
    sprintf(idFile, "%s%s%s", oPath, igdName, "_index.tsv");    
    FILE *fpi = fopen(idFile, "w");
    if(fpi==NULL)
        printf("Can't open file %s", idFile);     
    fprintf(fpi, "Index\tFile\tNumber of regions\tAvg size\n");    
    for(i=0; i<n_files; i++){  	
        tchr = strrchr(file_ids[i], '/');
        if(tchr!=NULL)
            tchr += 1;
        else
            tchr = file_ids[i];
        fprintf(fpi, "%i\t%s\t%i\t%f\n", i, tchr, nr[i], avg[i]/nr[i]);     
    }

    fclose(fpi);   
    free(nr);
    free(avg); 
    printf("igd_create 3\n");
	//4. Sort tile data and save into a single file
	igd_save(igd, oPath, igdName);	
	printf("igd_create 4\n");   
	free(file_ids); 
	free(splits);       
}

//-------------------------------------------------------------------------------------
int igd_create(int argc, char **argv)
{
    if (argc < 5) 
        return create_help(EX_OK);       

    int32_t i, j, n;
    char ipath[1024];
    char opath[1024];
    char *s1, *s2;
    strcpy(ipath, argv[2]);
    strcpy(opath, argv[3]);
    char *dbname = argv[4]; 
    int dtype = 1, ftype = 0;	//file type 1: list of bed file
    
    tile_size = 16384;								//b14 default
    for(i=5; i<argc; i++){
        if(strcmp(argv[i], "-s")==0 && i+1<argc)	//data structure type
            dtype = atoi(argv[i+1]); 
        if(strcmp(argv[i], "-b")==0 && i+1<argc){
            n = atoi(argv[i+1]);
            if(n>10 && n<20)
            	tile_size = pow(2, n);  
        } 
        if(strcmp(argv[i], "-f")==0) 
        	ftype = 1;                                
    }
    
	if(opath[strlen(opath)-1]!='/'){
	    strcat(opath, "/");
	}          
		  
	if(ftype==0 && dtype!=2){                            
	    if(ipath[strlen(ipath)-1]=='/'){
	        strcat(ipath, "*");
	    }
	    else if(ipath[strlen(ipath)-1]!='*'){
	        strcat(ipath, "/*");
	    }
	}  
    
    //check if the subfolders exist:    
    char ftmp[1024];      
    struct stat st = {0};  
    
    sprintf(ftmp, "%s%s%s", opath, dbname, ".igd");
    if(stat(ftmp, &st) == 0)
        printf("The igd database file %s exists!\n", ftmp);  
    else{
        if (stat(opath, &st) == -1){
            mkdir(opath, 0777);    
        }
        sprintf(ftmp, "%s%s", opath, "data0");
        if (stat(ftmp, &st) == -1)
            mkdir(ftmp, 0777);
        if(dtype==0)
        	create_igd0(ipath, opath, dbname);
        else if(dtype==2)
         	create_igd_bed4(ipath, opath, dbname);
        else if(ftype==1)
        	create_igd_f(ipath, opath, dbname);
        else //default
        	create_igd(ipath, opath, dbname);         
    } 

    return EX_OK;
}

