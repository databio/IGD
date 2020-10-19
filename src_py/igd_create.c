//===================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//database intervals sorted by _start: 8/12/2019
//-----------------------------------------------------------------------------------
#include "igd_create.h"

int create_help(int exit_code)
{
    fprintf(stderr,
"%s, v%s\n"
"usage:   %s create <input dir> <output dir> <output igd name> [options] \n"
"             -b  <Tile size in power of 2 (default 14)> \n"
"             -c  < .BED column as value >=4 (default 4) \n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

void create_iGD(iGD_t *iGD, char *iPath, char *oPath, char *igdName, int tile_size)
{  	
    //*. Check if the subfolders exist: 
    if(oPath[strlen(oPath)-1]!='/'){
        strcat(oPath, "/");
    }          
                                 
    if(iPath[strlen(iPath)-1]=='/'){
        strcat(iPath, "*");
    }
    else if(iPath[strlen(iPath)-1]!='*'){
        strcat(iPath, "/*");
    }   
       
    char ftmp[128];      
    struct stat st = {0};  
    
    sprintf(ftmp, "%s%s%s", oPath, igdName, ".igd");
    if(stat(ftmp, &st) == 0){
        printf("The igd database file %s exists!\n", ftmp); 
        return EX_OK;
    } 
    else{
        if (stat(oPath, &st) == -1){
            mkdir(oPath, 0777);    
        }
        sprintf(ftmp, "%s%s", oPath, "data0");
        if (stat(ftmp, &st) == -1)
            mkdir(ftmp, 0777);        
    } 

	//0. Initialize igd
	igd_t *igd = igd_init(tile_size);
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
    while(i0<n_files){
        //2.1 Start from (i0, L0): read till (i1, L1)
        ig = i0; 
        m = 0;    
        char **splits = malloc((nCols+1)*sizeof(char *));                
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
		free(splits);          
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
	//5. Form iGD
	sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");
	open_iGD(iGD, idFile);
}
