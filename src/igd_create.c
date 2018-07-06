//=====================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <err.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <glob.h>
#include <zlib.h>
#include <errno.h>
#include <sysexits.h>
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
//-------------------------------------------------------------------------------------
struct igd_data
{
    uint32_t i_idx;        			//genomic object--data set index
    uint32_t r_start;      			//region start
    uint32_t r_end;        			//region end
    uint32_t g_val;        			//signal level
};

struct query_data
{
    uint32_t q_idx;        			//tile index--locus position
    uint32_t r_start;      			//region start
    uint32_t r_end;        			//region end
    uint32_t g_val;        			//signal level
};

struct igd_mix
{    
    uint32_t m_idx;                 //tile index--chr 
    struct igd_data igd;
};

struct igd_info
{
    char* fileName;
    uint32_t nd;
    double md;
};

char** str_split( char* str, char delim, int *nmax);
char** str_split_t( char* str, int nItems);
struct query_data* get_igdlist(char *qfName, uint32_t *nblocks, uint32_t *nRegions, double *mRegion);	
struct igd_info* get_igdinfo(char *ifName, uint32_t *nFiles);
void create_igd(char *path, char* igdName);
void append_igd(struct igd_mix *mdata, uint32_t *counts, struct igd_info *fInfo, int nFiles, char* igdName);
void store_igd(char*igdName);
void str_splits( char* str, int *nmax, char **splits);

char *fileBase = "_b14_";         	//14 bits block
uint32_t nmax[] = {15940, 15580, 12760, 12240, 11670, 10990, 10260, 9370, 8860, 8610, 8710, 
        8580, 7300, 6840, 6510, 5830, 5370, 5160, 3820, 4150, 2980, 3240, 9880, 3510};
char *folder[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
     	 "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
uint32_t gstart[25] = {0, 15940, 31520, 44280, 56520, 68190, 79180, 89440, 98810, 107670, 
        116280, 124990, 133570, 140870, 147710, 154220, 160050, 165420, 170580, 174400, 
        178550, 181530, 184770, 194650, 198160}; 
uint32_t nTiles = 198160;
uint32_t nbp = 16384, bgz_buf = 1024;
uint64_t maxCount = 536870912;//134217728;//maxCount*16=2GB,67108864,33554432, 67108864;//
uint32_t *g2ichr;

//-------------------------------------------------------------------------------------
int compare_iidx(const void *a, const void *b)
{
    struct igd_data *pa = (struct igd_data *) a;
    struct igd_data *pb = (struct igd_data *) b;
    return pa->i_idx - pb->i_idx;
}

int compare_rend(const void *a, const void *b)
{
    struct igd_data *pa = (struct igd_data *) a;
    struct igd_data *pb = (struct igd_data *) b;
    return pa->r_end - pb->r_end;
}

int compare_qidx(const void *a, const void *b)
{
    struct query_data *pa = (struct query_data *) a;
    struct query_data *pb = (struct query_data *) b;
    return pa->q_idx - pb->q_idx;
}

int compare_rstart(const void *a, const void *b)
{
    struct query_data *pa = (struct query_data *) a;
    struct query_data *pb = (struct query_data *) b;
    return pa->r_start - pb->r_start;
}

int compare_midx(const void *a, const void *b)
{
    struct igd_mix *pa = (struct igd_mix *) a;
    struct igd_mix *pb = (struct igd_mix *) b;
    return pa->m_idx - pb->m_idx;
}

char** str_split_t( char* str, int nItems)
{
    char **splits;
    char *tmp;
    int i;
    if (str == NULL)
        return NULL;
    else {    
        splits = malloc((nItems+1)*sizeof(*splits)); 
        i=0;
        tmp = strtok(str, "\t");
        while(tmp!=NULL && i<nItems){
            splits[i] = tmp;
            tmp = strtok(NULL, "\t");
            i++;
        }
    }
    //printf("%s %s %s \n", splits[0], splits[1], splits[2]);
    return splits;
}

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

char** str_split( char* str, char delim, int *nmax)
{   //slightly faster than _t
    char** splits;
    char* ch;    
    int ns;
    if (str == NULL || delim=='\0')
        return NULL;
    else {
        splits = malloc((*nmax+1) * sizeof(*splits));
        splits[*nmax] = NULL;
        ch = str;
        ns = 1;
        splits[0] = str;
        do {
            if (*ch == delim)
            {
                splits[ns++] = &ch[1];
                *ch = '\0';
            }
            ch++;
        } while (*ch != '\0' && ns < *nmax+1);
    }
    *nmax = ns;
    return splits;
}

void append_igd(struct igd_mix *mdata, uint32_t *counts, struct igd_info *fInfo, int nFiles, char* igdName)
{   //binary: single folder (path without\)
    //both head file and data files (tile index)
    //Head file: path/igdName_index.tsv; datafiles: path/igdName_fileBase_tileNo.igd
    char idFile[128];
    sprintf(idFile, "%s%s%s", "igdata/", igdName, "_index.tsv");    
    FILE *fp = fopen(idFile, "a");
    if(fp==NULL)
        printf("Can't open file %s", idFile);
    uint32_t i, j;
    for(i=0; i<nFiles; i++){
        fprintf(fp, "%s %u %f\n", fInfo[i].fileName, fInfo[i].nd, fInfo[i].md);     
    }
    fclose(fp);
    
    //---append igd data: sorted already
    uint64_t nI = 0;
    //struct igd_data *td;
    for(i=0;i<nTiles;i++){
        if(counts[i]>0){
            sprintf(idFile, "%s%s%s%u%s", "igdata/data0/", igdName, fileBase, i, ".igd");
            fp = fopen(idFile, "ab");
            if(fp==NULL)
                printf("Can't open file %s", idFile);           
            for(j=nI;j<nI+counts[i];j++)
                fwrite(&mdata[j].igd, sizeof(struct igd_data), 1, fp);
            fclose(fp);           
            nI += counts[i];
        }
    }
}

//reload tile igd files, sort them and save them into a single file
void store_igd(char *igdName)
{   //first data, last counts: more efficient!
    char idFile[128], iname[128];
    sprintf(idFile, "%s%s%s%s", "igdata/data1/", igdName, fileBase, ".igd");    
    FILE *fp1 = fopen(idFile, "ab"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t i, nrec;  
    FILE *fp0; 
    struct igd_data *gdata;
    uint32_t *counts = calloc(nTiles, sizeof(uint32_t));
    for(i=0;i<nTiles;i++){
        sprintf(iname, "%s%s%s%u%s", "igdata/data0/", igdName, fileBase, i, ".igd");
        fp0 = fopen(iname, "rb");
        if(fp0!=NULL){   
            fseek(fp0, 0, SEEK_END);
            nrec = ftell(fp0)/sizeof(struct igd_data);
            fseek(fp0, 0, SEEK_SET);
            if(nrec>0){
                gdata = malloc(nrec*sizeof(struct igd_data));
                fread(gdata, sizeof(struct igd_data), nrec, fp0);
                qsort(gdata, nrec, sizeof(struct igd_data), compare_rend);
                //append the data to fp1
                fwrite(gdata, sizeof(struct igd_data), nrec, fp1);
                counts[i] = nrec;
                free(gdata);
            }      
            fclose(fp0);
        }
    }
    fwrite(counts, sizeof(uint32_t), nTiles, fp1);
    fclose(fp1);
    free(counts);
}

//create ucsc igd from gz
void create_igd_gz(char *iPath, char *oPath, char *igdName)
{   //Process line by line
    //1. Get the files  
    glob_t gResult;
    //strcat(iPath, "*");
    int rtn = glob(iPath, 0, NULL, &gResult);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath);
        return;
    }
    
    char** file_ids = gResult.gl_pathv;
    uint32_t n_files = gResult.gl_pathc; 
    if(n_files<6)   
        printf("Too few files (add to path /*): %u\n", n_files); 
     
    uint32_t *nd = calloc(n_files, sizeof(uint32_t));
    float *md = calloc(n_files, sizeof(float)); 
    //2. Read region data
    uint32_t i, ii, j, k, t, df1, df2, df4, ichr, n1, n2, ti, nR;
    uint32_t *counts = calloc(nTiles, sizeof(uint32_t));    //134217728*16=2G
    uint32_t *Counts = calloc(nTiles, sizeof(uint32_t));    //total
    double delta;    
    char idFile[256];  
    clock_t start, end;
    start = clock();  
    end = clock();    
    printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);   
    
    //-------------------------------------------------------------------------
    char *ftype;
    int nCols=5, err;                    
    int tlen, bytes_read;
    unsigned char buffer[bgz_buf]; 
    uint64_t cTotal, nL;
    //1. find record lines (from one part of big file or many samll files) of size ~maxCount
    //2. procese and save (n0: starting file index, n1 ending file' n0 start line, n2 ending line)
    //(i0, L0)-->(i1, L1): i1 can be the same as i0 if filesize>maxCount
    uint32_t i0=0, i1=0, L0=0, L1=1, m;  
    struct igd_data **gData; 
    char **splits = malloc((nCols+1)*sizeof(char *)); 
    while(i0<n_files){
        //1. start from (i0, L0): find (i1, L1)
        cTotal = 0;
        i = i0; 
        m = 0;              
        memset(counts, 0, nTiles*sizeof(uint32_t));
        printf("start: %u\n", i0);
        while(i<n_files && m==0){   //n>0 defines breaks when reading a big file         
            ftype = file_ids[i] + strlen(file_ids[i]) - 7;
            if(strcmp(".bed.gz", ftype)==0 || strcmp(".txt.gz", ftype)==0){        
                //a. Prepare: get the counts               
                //printf("%s", file_ids[i]);
                gzFile zfile;
                zfile = gzopen (file_ids[i], "r");
                if (!zfile) {
                    fprintf (stderr, "gzopen of '%s' failed: %s.\n", file_ids[i],
                             strerror (errno));
                    exit (EXIT_FAILURE);
                }             
                
                printf("%u %u\t", i, (uint32_t)cTotal);
                
                nL = 0; 
                if(i==i0 && L0>0){   //pass n0 lines of a big file
                    while(nL<L0 && gzgets(zfile, buffer, bgz_buf)!=NULL)
                        nL++;              
                }                 
                while(m==0 && gzgets(zfile, buffer, bgz_buf)!=NULL){
                    nL++;
                    //splits = str_split_t(buffer, nCols); 
                    str_splits(buffer, &nCols, splits);  
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0],"chrY")==0)
                            ichr = 23;
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    if(ichr>=0 && ichr<=23){
                        df1 = (uint32_t)atoi(splits[1]);               
                        df2 = (uint32_t)atoi(splits[2]); 
                        n1 = df1/nbp;
                        n2 = df2/nbp-n1;       
                        for(j=0;j<=n2;j++){
                            if(n1+j<nmax[ichr]){
                                ti = n1+j+gstart[ichr];
                                counts[ti]++; 
                                cTotal++;                          
                            }       
                        }
                        if(cTotal>maxCount){
                            m = 1;
                            i1 = i;
                            L1 = nL;    //number of total lines or next line
                        }
                    } 
                    //free(splits);               
                }   //while getLine 
                gzclose(zfile);      
            }   //if bed.gz
            if(m==0)
                i++;
        }   //while i<n_files
        //---------------------------------------------------------------------
        //2. process files i0--i           
        gData = malloc(nTiles*sizeof(struct igd_data*));
        for(j=0; j<nTiles; j++){
            if(counts[j]>0)
                gData[j] = malloc(counts[j]*sizeof(struct igd_data));  
        }  
        if(gData==NULL)
            printf("Error: memory allocation with igd_data");        
        memset(counts, 0, nTiles*sizeof(uint32_t));
        printf("(%u, %u) (%u, %u) \n", i0, L0, i, L1);           
        for(ii=i0; ii<i; ii++){   //n>0 defines breaks when reading a big file
            ftype = file_ids[ii] + strlen(file_ids[ii]) - 7;
            if(strcmp(".bed.gz", ftype)==0 || strcmp(".txt.gz", ftype)==0){        
                gzFile zfile = gzopen (file_ids[ii], "r"); 
                nL = 0; 
                if(ii==i0 && L0>0){   //pass n0 lines of a big file
                    while(nL<L0 && gzgets(zfile, buffer, bgz_buf)!=NULL)
                        nL++;              
                }                
                while(gzgets(zfile, buffer, bgz_buf)!=NULL){//tbd: read 1Mb each time
                    //splits = str_split(buffer,'\t', &nCols);  
                    str_splits(buffer, &nCols, splits);
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0],"chrY")==0)
                            ichr = 23;
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    if(ichr>=0 && ichr<=23){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
                        if(nCols>4)
                            df4 = (uint32_t)atoi(splits[4]);
                        else
                            df4 = 100;
                        n1 = df1/nbp;
                        n2 = df2/nbp-n1;  
                        //-----------------------------------------------------
                        for(j=0;j<=n2;j++){
                            t = n1+j;
                            if(t<nmax[ichr]){
                                ti = t+gstart[ichr];
                                k = counts[ti]; 
                                gData[ti][k].i_idx = ii;
                                gData[ti][k].r_start = df1;
                                gData[ti][k].r_end = df2;
                                gData[ti][k].g_val = df4;
                                counts[ti]++;  
                            }
                        }
                    }  
                    //free(splits);
                }   //while gzgets           
                gzclose (zfile);               
            }   //gz file
        }   //ii
        //--------------------------------------------------------------------- 
        if(m>0){
            ii=i;  //m>0 defines breaks when reading a big file
            ftype = file_ids[ii] + strlen(file_ids[ii]) - 7;
            if(strcmp(".bed.gz", ftype)==0 || strcmp(".txt.gz", ftype)==0){        
                gzFile zfile = gzopen (file_ids[ii], "r"); 
                nL = 0;
                if(ii==i0 && L0>0){
                    while(nL<L0 && gzgets(zfile, buffer, bgz_buf)!=NULL)
                        nL++; 
                }                            
                while(nL<L1 && gzgets(zfile, buffer, bgz_buf)!=NULL){//tbd: read 1Mb each time
                    //splits = str_split(buffer,'\t', &nCols);  
                    str_splits(buffer, &nCols, splits);
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0],"chrY")==0)
                            ichr = 23;
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    if(ichr>=0 && ichr<=23){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
                        if(nCols>4)
                            df4 = (uint32_t)atoi(splits[4]);
                        else
                            df4 = 100;
                        n1 = df1/nbp;
                        n2 = df2/nbp-n1;  
                        //-----------------------------------------------------
                        for(j=0;j<=n2;j++){
                            t = n1+j;
                            if(t<nmax[ichr]){
                                ti = t+gstart[ichr];
                                k = counts[ti]; 
                                gData[ti][k].i_idx = ii;
                                gData[ti][k].r_start = df1;
                                gData[ti][k].r_end = df2;
                                gData[ti][k].g_val = df4;
                                counts[ti]++;  
                            }
                        }
                    }  
                    //free(splits);
                    nL++;
                }            
                gzclose (zfile);               
            }   //gz file  
        }    
        //---------------------------------------------------------------------           
        end = clock();    
        printf("File %u processing time: %f \n", i, ((double)(end-start))/CLOCKS_PER_SEC);        
        //save gData
        for(j=0;j<nTiles;j++){
            //if(j%1000==0)
            //    printf("%u %u\n", j, counts[j]);
            if(counts[j]>0){    
                Counts[j] += counts[j];
                k = g2ichr[j];                      
                sprintf(idFile, "%s%s%s/%s_%u%s", oPath, "data0/", folder[k], igdName, j-gstart[k], ".igd");
                FILE *fp = fopen(idFile, "ab");
                if(fp==NULL)
                    printf("Can't open file %s", idFile);
                fwrite(gData[j], sizeof(struct igd_data), counts[j], fp);
                fclose(fp); 
                //free(gData[j]); 
            }
        }   
        for(j=0;j<nTiles;j++){;
            if(counts[j]>0)
                free(gData[j]);   
        }
        free(gData);
        end = clock();    
        printf("Saving time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
        gData = NULL;                                      
        //---------------------------------------------------------------------
        i0 = i; 
        L0 = L1;
        L1 = 0;
    }  
    //save _index.tsv: 4 columns--index, filename, nd, md
    //Also has a header line: 
    char *tchr;   
    sprintf(idFile, "%s%s%s", oPath, igdName, "_index.tsv");    
    FILE *fpi = fopen(idFile, "w");
    if(fpi==NULL)
        printf("Can't open file %s", idFile);
        
    fprintf(fpi, "Index\tFile\tNumber of regions\tAvg size\n");    
    for(i=0; i<n_files; i++){
        tchr = strrchr(file_ids[i], '/');
        if(tchr!=NULL)
            tchr+=1;
        else
            tchr = file_ids[i];
        fprintf(fpi, "%u\t%s\t%u\t%f\n", i, tchr, nd[i], md[i]/nd[i]);     
    }
    fclose(fpi);   
    free(nd);
    free(md);    
    end = clock();    
    printf("TSV saved: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
    //-------------------------------------------------------------------------
    //Reload tile data, sort and save them into a single file
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");    
    FILE *fp1 = fopen(idFile, "wb"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t nrec;  
    FILE *fp0;    
    fwrite(Counts, sizeof(uint32_t), nTiles, fp1);
    char iname[256]; 
    struct igd_data *gdata;            
    for(i=0;i<nTiles;i++){
        if(Counts[i]>0){
            k = g2ichr[i];
            sprintf(iname, "%s%s%s/%s_%u%s", oPath, "data0/", folder[k], igdName, i-gstart[k], ".igd");
            fp0 = fopen(iname, "rb");
            if(fp0!=NULL){   
                nrec = Counts[i];
                gdata = malloc(nrec*sizeof(struct igd_data));
                fread(gdata, sizeof(struct igd_data), nrec, fp0);
                qsort(gdata, nrec, sizeof(struct igd_data), compare_rend);
                //append the data to fp1
                fwrite(gdata, sizeof(struct igd_data), nrec, fp1);
                free(gdata);
            }      
            fclose(fp0);
        }
    }
    fclose(fp1); 
    end = clock();    
    printf("igd_w finished: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);    
    //------------------------------------------------------------------------- 
    free(splits);  
    free(counts);
    free(Counts);
    globfree(&gResult); 
}

//create ucsc igd from gz
void create_igd_gz0(char *iPath, char *oPath, char *igdName)
{   //Process line by line
    //1. Get the files  
    glob_t gResult;
    int rtn = glob(iPath, 0, NULL, &gResult);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath);
        return;
    }

    char** file_ids = gResult.gl_pathv;
    uint32_t n_files = gResult.gl_pathc;
     
    //2. Read region data
    uint32_t i, j, k, t, df1, df2, df4, ichr, n1, n2, ti, nR;
    uint32_t *counts = calloc(nTiles, sizeof(uint32_t));    //134217728*16=2G
    struct igd_data **gData;
    struct igd_info *fi = malloc(1*sizeof(struct igd_info));
    double delta;    
    
    printf("nfiles:%i\n", (int)gResult.gl_pathc);        
    //-------------------------------------------------------------------------
    char *ftype;
    char idFile[256];
    char **splits;
    int nCols=5, err;                    
    int tlen, bytes_read;
    unsigned char buffer[bgz_buf];  
    for(i=0; i<n_files; i++){
        ftype = file_ids[i] + strlen(file_ids[i]) - 7;
        if(strcmp(".bed.gz", ftype)==0 || strcmp(".txt.gz", ftype)==0){        
            //a. Prepare: get the counts               
            //printf("%s", file_ids[i]);
            memset(counts, 0, nTiles*sizeof(uint32_t));
            gzFile zfile;
            zfile = gzopen (file_ids[i], "r");
            if (!zfile) {
                fprintf (stderr, "gzopen of '%s' failed: %s.\n", file_ids[i],
                         strerror (errno));
                exit (EXIT_FAILURE);
            }                   
            while(gzgets(zfile, buffer, bgz_buf)!=NULL){
                splits = str_split(buffer,'\t', &nCols);  
                ichr = -1;
                tlen = strlen(splits[0]);
                if(tlen<6 && tlen>3){
                    if(strcmp(splits[0], "chrX")==0)
                        ichr = 22;
                    else if(strcmp(splits[0],"chrY")==0)
                        ichr = 23;
                    else{
                        rtn = atoi(&splits[0][3]);
                        if(rtn!=0)
                            ichr = (uint32_t)(rtn-1);
                    }
                }
                if(ichr>=0 && ichr<=23){
                    df1 = (uint32_t)atoi(splits[1]);               
                    df2 = (uint32_t)atoi(splits[2]); 
                    n1 = df1/nbp;
                    n2 = df2/nbp-n1;       
                    for(j=0;j<=n2;j++){
                        if(n1+j<nmax[ichr]){
                            ti = n1+j+gstart[ichr];
                            counts[ti]++; 
                        }       
                    }
                }                
            }
                                  
            //b. Get the igd
            gData = malloc(nTiles*sizeof(struct igd_data*));
            for(j=0; j<nTiles; j++){
                if(counts[j]>0)
                    gData[j] = malloc(counts[j]*sizeof(struct igd_data));  
            }          
            memset(counts, 0, nTiles*sizeof(uint32_t));
            nR = 0;
            delta = 0.0;
            gzseek(zfile, 0, SEEK_SET);                  
            while(gzgets(zfile, buffer, bgz_buf)!=NULL){
                splits = str_split(buffer,'\t', &nCols);  
                ichr = -1;
                tlen = strlen(splits[0]);
                if(tlen<6 && tlen>3){
                    if(strcmp(splits[0], "chrX")==0)
                        ichr = 22;
                    else if(strcmp(splits[0],"chrY")==0)
                        ichr = 23;
                    else{
                        rtn = atoi(&splits[0][3]);
                        if(rtn!=0)
                            ichr = (uint32_t)(rtn-1);
                    }
                }
                if(ichr>=0 && ichr<=23){
                    df1 = (uint32_t)atoi(splits[1]);
                    df2 = (uint32_t)atoi(splits[2]);
                    if(nCols>4)
                        df4 = (uint32_t)atoi(splits[4]);
                    else
                        df4 = 100;
                    delta += df2-df1;
                    n1 = df1/nbp;
                    n2 = df2/nbp-n1;  
                    //---------------------------------------------------------
                    for(j=0;j<=n2;j++){
                        t = n1+j;
                        if(t<nmax[ichr]){
                            ti = t+gstart[ichr];
                            k = counts[ti]; 
                            gData[ti][k].i_idx = i;
                            gData[ti][k].r_start = df1;
                            gData[ti][k].r_end = df2;
                            gData[ti][k].g_val = df4;
                            counts[ti]++;  
                        }
                    }
                    nR++;
                }                
            }
            gzclose (zfile);     
            //-----------------------------------------------------------------
            //save gData
            for(j=0;j<nTiles;j++){
                if(counts[j]>0){             
                    sprintf(idFile, "%s%u%s", igdName, j, ".igd");
                    //if(j==2)
                    //    printf("%s\n", idFile);
                    //FILE *fp = fopen(idFile, "ab");
                    //fwrite(gData[j], sizeof(struct igd_data), counts[j], fp);
                    //fclose(fp);  
                }
            }
            //qsort(gData, nt, sizeof(struct igd_mix), compare_midx);
            printf("%i %u \n", i, nR);         
            //fi[0].fileName = file_ids[i];
            //fi[0].nd = nR;
            //fi[0].md = delta/nR;
            //append_igd(gData, counts, fi, 1, igdName);
            for(j = 0; j < nTiles; j++){
                if(counts[j]>0)
                    free(gData[j]);
            }
            free(gData);
        }
        else if(strcmp(".vcf.gz", ftype)==0){
        
        }
        else{
            printf("File type not supported: %s\n", ftype);
        }
        //printf("%i %s %llu\n", i, file_ids[i], (unsigned long long)f_size[i]);
    }
    //qsort(query, *nblocks, sizeof(struct query_data), compare_qidx); 
    //free(fi->fileName);
    //free(fi);    
    free(counts);
    globfree(&gResult); 
}

//create ucsc igd from bed files
void create_igd(char *path, char *igdName)
{   //For large data sources (>2GB): save to individual-tile files in append
    //   mode (igd mode 0); when done: reload each tile file and sort it and save
    //   to a single disk file (igd mode 1) for better search performance.
    //   In case of a single source data > 2GB(mem limit), deal with it chr by chr. 
    //1. Get the files  
    glob_t gResult;
    int rtn = glob(path, 0, NULL, &gResult);     
    if(rtn!=0){
        printf("wrong dir path: %s", path);
        return;
    }

    char** file_ids = gResult.gl_pathv;
    uint32_t n_files = gResult.gl_pathc;
     
    //2. Read region data
    uint32_t i, j, k, t, df1, df2, df4, ichr, n1, n2, ti, nR, nt;
    uint32_t *counts = calloc(nTiles, sizeof(uint32_t));    //134217728*16=2G
    struct igd_mix *gData;
    struct igd_info *fi = malloc(1*sizeof(struct igd_info));
    double delta;    
    
    printf("nfiles:%i\n", (int)gResult.gl_pathc);        
    //----------------------------------------------------------------------
    char *ftype;
    char **splits;
    int nCols=5, err;                    
    int bytes_read;
    unsigned char buffer[bgz_buf];  
    for(i=0; i<n_files; i++){
        ftype = file_ids[i] + strlen(file_ids[i]) - 7;
        memset(counts, 0, nTiles*sizeof(uint32_t));
        nt = 0;
        nR = 0;
        delta = 0.0;
        if(strcmp(".bed.gz", ftype)==0 || strcmp(".txt.gz", ftype)==0){
            //printf("%s\n", file_ids[i]);
            gzFile zfile;
            zfile = gzopen (file_ids[i], "r");
            if (! zfile) {
                fprintf (stderr, "gzopen of '%s' failed: %s.\n", file_ids[i],
                         strerror (errno));
                exit (EXIT_FAILURE);
            }                   
            while(gzgets(zfile, buffer, bgz_buf)!=NULL){
                splits = str_split(buffer,'\t', &nCols);       
                if(strlen(splits[0])<6){
                    if(strcmp(splits[0], "chrX")==0)
                        ichr = 22;
                    else if(strcmp(splits[0],"chrY")==0)
                        ichr = 23;
                    else
                        ichr = (uint32_t)(atoi(&splits[0][3])-1);
                    df1 = (uint32_t)atoi(splits[1]);
                    df2 = (uint32_t)atoi(splits[2]);
                    if(nCols>4)
                        df4 = (uint32_t)atoi(splits[4]);
                    else
                        df4 = 600;
                    delta += df2-df1;
                    n1 = df1/nbp;
                    n2 = df2/nbp-n1;  
                    //---------------------------------------------------------
                    for(j=0;j<=n2;j++){
                        ti = n1+j+gstart[ichr];
                        gData[nt].m_idx = ti;
                        counts[ti]++;
                        gData[nt].igd.i_idx = i;
                        gData[nt].igd.r_start = df1;
                        gData[nt].igd.r_end = df2;
                        gData[nt].igd.g_val = df4;
                        nt++;
                    }
                    nR++;
                }                
            }
            gzclose (zfile);     
            //-----------------------------------------------------------------
            qsort(gData, nt, sizeof(struct igd_mix), compare_midx);         
            fi[0].fileName = file_ids[i];
            fi[0].nd = nR;
            fi[0].md = delta/nR;
            append_igd(gData, counts, fi, 1, igdName);
        }
        else if(strcmp(".vcf.gz", ftype)==0){
        
        }
        else{
            printf("File type not supported: %s\n", ftype);
        }
        //printf("%i %s %llu\n", i, file_ids[i], (unsigned long long)f_size[i]);
    }
    //qsort(query, *nblocks, sizeof(struct query_data), compare_qidx); 
    free(fi->fileName);
    free(fi);     
    free(gData);
    free(counts);
    globfree(&gResult); 
}

struct query_data* get_igdlist(char *qfName, uint32_t *nblocks, uint32_t *nRegions, double *mRegion)
{  
    FILE *fp = fopen(qfName, "r");
    if(fp==NULL)
        return NULL;
    char buf[1024], ch;
    int nregions = 0;
    while((ch=fgetc(fp))!=EOF){
        if(ch=='\n')
        	    nregions++;
    }
    uint8_t *df0 = malloc(nregions*sizeof(uint8_t));
    uint32_t *df1 = malloc(nregions*sizeof(uint32_t));
    uint32_t *df2 = malloc(nregions*sizeof(uint32_t)); 
    uint32_t *n1 = malloc(nregions*sizeof(uint32_t));
    uint8_t *n2 = malloc(nregions*sizeof(uint8_t)); 

    char **splits;
    int nCols = 5;  
    fseek(fp, 0, SEEK_SET);
    int k, nm, i=0, nextras=0;
    double delta=0.0;//number of regions that cross the tile boundary
    while(fgets(buf, 1024, fp)!=NULL){	
        //splits = str_split_t(buf, nItems);
        splits = str_split(buf,'\t', &nCols);       
        if(strlen(splits[0])<6){
    	   if(strcmp(splits[0], "chrX")==0)
            	df0[i] = 22;
    	   else if(strcmp(splits[0],"chrY")==0)
            	df0[i] = 23;
    	   else
        	   df0[i] = (uint8_t)(atoi(&splits[0][3])-1);
    	   df1[i] = (uint32_t)atoi(splits[1]);
    	   df2[i] = (uint32_t)atoi(splits[2]);
    	   delta += df2[i]-df1[i];
    	   n1[i] = df1[i]/nbp;
    	   n2[i] = df2[i]/nbp-n1[i];  
    	   nextras += n2[i];     
    	   i++;
        }
    }
    fclose(fp);
    nm = i;
    
    *nblocks = nm+nextras;
    *nRegions = nm;
    *mRegion = delta/(double)nm;
    struct query_data *query = malloc(*nblocks*sizeof(struct query_data));
    uint32_t j=0;
    for(i=0;i<nm;i++){
        query[j].q_idx   = n1[i] + gstart[df0[i]];
        query[j].r_start = df1[i];
        query[j].r_end   = df2[i];
        j++;
        if(n2[i]>0){
            for(k=1; k<=n2[i]; k++){
            	query[j].q_idx   = k + n1[i] + gstart[df0[i]];
            	query[j].r_start = df1[i];
            	query[j].r_end   = df2[i];	
            	j++;		
            }
        }
    }
    qsort(query, *nblocks, sizeof(struct query_data), compare_rstart);
    free(df0);
    free(df1);
    free(df2);
    free(n1);
    free(n2);
    return query;   
}

struct igd_info* get_igdinfo(char *ifName, uint32_t *nFiles)
{   //read head file __index.tsv to get info
    FILE *fp = fopen(ifName, "r");
    if(fp==NULL){
        printf("file not found:%s\n", ifName);
        return NULL;
    }
    char buf[1024], ch;
    uint32_t i, nfiles = -1;
    while((ch=fgetc(fp))!=EOF){
        if(ch=='\n')
        	    nfiles++;
    }
    struct igd_info *fi = (struct igd_info*)malloc(nfiles*sizeof(struct igd_info));
    //char** fNames = malloc(nfiles*sizeof(char*));
    //uint32_t *nd = malloc(nfiles*sizeof(uint32_t));
    //double *md = malloc(nfiles*sizeof(double));
    //------------
    char **splits;
    int ncols = 4;  
    fseek(fp, 0, SEEK_SET);
    i=0;
    fgets(buf, 1024, fp);   //header
    while(fgets(buf, 1024, fp)!=NULL){	
        splits = str_split(buf,'\t', &ncols);        
        fi[i].fileName = (char *)calloc(strlen(splits[1]) + 1, sizeof(char));
        strcpy(fi[i].fileName, splits[1]);
        fi[i].nd = (uint32_t)atoi(splits[3]);
        fi[i].md = (double)atoi(splits[2]);   
        i++;
    }        
    //printf("%s %i %f \n", fNames[6], nd[6], md[6]);  
    //printf("Total file: %i\n", i);
    *nFiles = nfiles;
    fclose(fp);
    return fi;
}

//-------------------------------------------------------------------------------------
int igd_create(int argc, char **argv)
{
    if (argc < 4) 
        errx(1, "usage:\t%s <input path><ouput path><dbName>\n", argv[0]);   
    int usingw = 0;
    if (argc==3)
        usingw = 1;     
    //convert block index to chr index for convenience
    g2ichr = malloc(nTiles*sizeof(uint32_t));
    uint32_t i, j;
    for(i=0; i<24; i++){  
        	for(j=gstart[i]; j<gstart[i+1]; j++)      
        	    g2ichr[j] = i;
    }
    //---------------------------------------------------------------------------------  
    //char *ipath = "/media/john/CE30F6EE30F6DC81/ucsc_data/parsed_tracks_sorted/*";//ucscweb_sort/*";//ucsc_data/database/*";//ucscweb_sort/*";//roadmap_sort/*";
    //char *opath = "/media/john/CE30F6EE30F6DC81/ucsc_igd_b14/";  
    //char *ipath = "/scratch/jf4x/rme/*";//ucscweb_sort/*";//ucsc_data/database/*";//ucscweb_sort/*";//roadmap_sort/*";    
    //char *opath = "/media/john/CE30F6EE30F6DC81/rme_igd_b14/";       
    //char *opath = "/scratch/jf4x/rme_igd/";//ucscweb_sort/*";//ucsc_data/database/*";//ucscweb_sort/*";//roadmap_sort/*"; 
    char *ipath = argv[1];
    char *opath = argv[2];
    char *dbname = argv[3];
    create_igd_gz(ipath, opath, dbname);    
    //if(usingw)
    //    search(qfName, igdName);      
        //testMain(fname, igdName);

    free(g2ichr);
    return EX_OK;
}

