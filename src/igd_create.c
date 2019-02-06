//=====================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include "igd_base.h"
#include <zlib.h>
#include <sys/stat.h>
#include <sys/types.h>

//-------------------------------------------------------------------------------------
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
void create_igd_gz(char *iPath, char *oPath, char *igdName, int mode)
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
    if(n_files<1)   
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
    //printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);   
    
    //-------------------------------------------------------------------------
    char *ftype;
    int nCols=16, err;                    
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
        //printf("start: %u\n", i0);
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
                
                //printf("%u %u\t", i, (uint32_t)cTotal);
                
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
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
        //printf("(%u, %u) (%u, %u) \n", i0, L0, i, L1);           
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){ 
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
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
        //printf("File %u processing time: %f \n", i, ((double)(end-start))/CLOCKS_PER_SEC);        
        //save gData
        //printf("igdName %s", igdName);
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
        //printf("Saving time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
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
    //printf("TSV saved: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
    //-------------------------------------------------------------------------
    //Reload tile data, sort and save them into a single file
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");    
    FILE *fp1 = fopen(idFile, "wb"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t nrec;  
    FILE *fp0;       
    i=1122330000;   //type 0 data
    fwrite(&i, sizeof(uint32_t), 1, fp1);      
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
            if(mode==0)
                remove(iname);              
        }
    }
    fclose(fp1); 
    end = clock();    
    //printf("igd_w finished: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);       
    free(splits);  
    free(counts);
    free(Counts);
    globfree(&gResult); 
}

//create ucsc igd from gz
void create_igd_gz1(char *iPath, char *oPath, char *igdName, int mode)
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
    if(n_files<1)   
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
    //printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);   
    
    //-------------------------------------------------------------------------
    char *ftype;
    int nCols=16, err;                    
    int tlen, bytes_read;
    unsigned char buffer[bgz_buf]; 
    uint64_t cTotal, nL;
    //1. find record lines (from one part of big file or many samll files) of size ~maxCount
    //2. procese and save (n0: starting file index, n1 ending file' n0 start line, n2 ending line)
    //(i0, L0)-->(i1, L1): i1 can be the same as i0 if filesize>maxCount
    uint32_t i0=0, i1=0, L0=0, L1=1, m;  
    struct igd_data1 **gData; 
    char **splits = malloc((nCols+1)*sizeof(char *)); 
    while(i0<n_files){
        //1. start from (i0, L0): find (i1, L1)
        cTotal = 0;
        i = i0; 
        m = 0;              
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("start: %u\n", i0);
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
                
                //printf("%u %u\t", i, (uint32_t)cTotal);
                
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
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
        gData = malloc(nTiles*sizeof(struct igd_data1*));
        for(j=0; j<nTiles; j++){
            if(counts[j]>0)
                gData[j] = malloc(counts[j]*sizeof(struct igd_data1));  
        }  
        if(gData==NULL)
            printf("Error: memory allocation with igd_data");        
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("(%u, %u) (%u, %u) \n", i0, L0, i, L1);           
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
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
        //printf("File %u processing time: %f \n", i, ((double)(end-start))/CLOCKS_PER_SEC);        
        //save gData
        //printf("igdName %s", igdName);
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
                fwrite(gData[j], sizeof(struct igd_data1), counts[j], fp);
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
        //printf("Saving time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
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
    //printf("TSV saved: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
    //-------------------------------------------------------------------------
    //Reload tile data, sort and save them into a single file
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");    
    FILE *fp1 = fopen(idFile, "wb"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t nrec;  
    FILE *fp0;            
    i=1122331111;   //type 3 data
    fwrite(&i, sizeof(uint32_t), 1, fp1);        
    fwrite(Counts, sizeof(uint32_t), nTiles, fp1);
    char iname[256]; 
    struct igd_data1 *gdata;        
    for(i=0;i<nTiles;i++){
        if(Counts[i]>0){
            k = g2ichr[i];
            sprintf(iname, "%s%s%s/%s_%u%s", oPath, "data0/", folder[k], igdName, i-gstart[k], ".igd");
            fp0 = fopen(iname, "rb");
            if(fp0!=NULL){   
                nrec = Counts[i];
                gdata = malloc(nrec*sizeof(struct igd_data1));
                fread(gdata, sizeof(struct igd_data1), nrec, fp0);
                qsort(gdata, nrec, sizeof(struct igd_data1), compare_rend1);
                //append the data to fp1
                fwrite(gdata, sizeof(struct igd_data1), nrec, fp1);
                free(gdata);
            }      
            fclose(fp0);
            if(mode==0)
                remove(iname);              
        }
    }
    fclose(fp1); 
    end = clock();    
    //printf("igd_w finished: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);    
    //-------------------------------------------------------------------------     
    free(splits);  
    free(counts);
    free(Counts);
    globfree(&gResult); 
}

void constructAIList(struct igd_data2* B, int nB, struct igd_data2* aiL)
{   //return the ailist with head: [nSub, lSub], cLen=20; fixed head size: 4*8 Bytes==2 elements of igd_data2
    //aiL size>=2 if exist; nB>0  
    int lenT, len, iter, i, i1, j, j1, k, k0, t, tS, tE; 
    uint32_t maxE, tt; 
    int numC, idxC[7], lenC[7];      
    struct igd_data2* aiT;
    struct igd_data2 tData; 
    int maxC=7, minL = 64, cLen = 30, cLen1=10;  //max number of components, minL>cLen      
        
    for(i=0;i<7;i++){
        lenC[i]=0;
        idxC[i]=0;
    }
    
    qsort(B, nB, sizeof(struct igd_data2), compare_rstart2);              
    if(nB<=minL){                     
        memcpy(&aiL[2], B, nB*sizeof(struct igd_data2)); 
        //printf("nB in %i\n", nB);           
        numC = 1;
        lenC[0] = nB;         
    }
    else{ 
        aiT = malloc(nB*sizeof(struct igd_data2));    
        iter = 0;
        lenT = nB;
        k = 0;
        k0 = 0;
        while(iter<maxC && lenT>minL){   
            len = 0;            
            for(t=0; t<lenT-cLen; t++){
                tt = B[t].r_end;
                j=1;    j1=1;
                while(j<cLen && j1<cLen1){
                    if(B[j+t].r_end>=tt)  
                        j1++;
                    j++;
                }
                if(j1<cLen1){//number of items not satisfying
                    memcpy(&aiT[len], &B[t], sizeof(struct igd_data2));
                    len++;
                }
                else{
                    memcpy(&aiL[k+2], &B[t], sizeof(struct igd_data2));
                    k++;
                }                    
            } 
            memcpy(&aiL[k+2], &B[lenT-cLen], cLen*sizeof(struct igd_data2));
               
            k += cLen;
            lenT = len;                
            idxC[iter] = k0;
            lenC[iter] = k-k0;
            k0 = k;
            if(lenT<=minL || iter==maxC-2){//exit: add aiT to the end
                iter++;
                if(lenT>0){
                    memcpy(&aiL[k+2], aiT, lenT*sizeof(struct igd_data2));
                    idxC[iter] = k;
                    lenC[iter] = lenT;
                    iter++;
                    numC = iter;
                }
                else
                    numC = iter;                   
            }
            else{
                memcpy(B, aiT, lenT*sizeof(struct igd_data2));
                iter++;
            }
        }
        free(aiT);     
    }  
    
    //for(i=0;i<7;i++)
    //    printf("lenC 2 %i\n", lenC[i]);
    
    //augment  
    aiL[0].i_idx   = numC;
    aiL[0].r_start = lenC[0];    
    aiL[0].r_end   = lenC[1];
    aiL[0].r_max   = lenC[2];
    aiL[1].i_idx   = lenC[3];
    aiL[1].r_start = lenC[4];    
    aiL[1].r_end   = lenC[5];
    aiL[1].r_max   = lenC[6]; 
           
    for(j=0; j<numC; j++){ 
        k0 = idxC[j];
        k = k0 + lenC[j];
        maxE = aiL[k0+2].r_end;
        for(t=k0+2; t<k+2; t++){
            if(aiL[t].r_end > maxE)
                maxE = aiL[t].r_end;
            aiL[t].r_max = maxE;  
        }             
    }           
}

//create ucsc igd from gz: ailist in-bin
void create_igd_gz2(char *iPath, char *oPath, char *igdName, int mode)
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
    if(n_files<1)   
        printf("No file (add to path /*): %u\n", n_files); 
     
    uint32_t *nd = calloc(n_files, sizeof(uint32_t));
    float *md = calloc(n_files, sizeof(float)); 

    //2. Read region data
    uint32_t i, ii, j, k, t, df1, df2, df4, ichr, n1, n2, ti, nR;
    uint32_t *counts = calloc(nTiles, sizeof(uint32_t));    //134217728*16=2G
    uint32_t *Counts = calloc(nTiles, sizeof(uint32_t));    //total
    double delta;    
    char idFile[256];  
    //clock_t start, end;
    //start = clock();  
    //end = clock();    
    //printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);   
    
    //-------------------------------------------------------------------------
    char *ftype;
    int nCols=16, err;                    
    int tlen, bytes_read;
    unsigned char buffer[bgz_buf]; 
    uint64_t cTotal, nL;
    //1. find record lines (from one part of big file or many samll files) of size ~maxCount
    //2. procese and save (n0: starting file index, n1 ending file' n0 start line, n2 ending line)
    //(i0, L0)-->(i1, L1): i1 can be the same as i0 if filesize>maxCount
    uint32_t i0=0, i1=0, L0=0, L1=1, m;  
    struct igd_data2 **gData; 
    char **splits = malloc((nCols+1)*sizeof(char *)); 
    while(i0<n_files){
        //1. start from (i0, L0): find (i1, L1)
        cTotal = 0;
        i = i0; 
        m = 0;              
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("start: %u\n", i0);
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
                
                //printf("%u %u\t", i, (uint32_t)cTotal);
                
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
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
        gData = malloc(nTiles*sizeof(struct igd_data2*));//easy to copy
        for(j=0; j<nTiles; j++){
            if(counts[j]>0)
                gData[j] = malloc(counts[j]*sizeof(struct igd_data2));  
        }  
        if(gData==NULL)
            printf("Error: memory allocation with igd_data");        
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("(%u, %u) (%u, %u) \n", i0, L0, i, L1);           
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
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
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
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
        //end = clock();    
        //printf("File %u processing time: %f \n", i, ((double)(end-start))/CLOCKS_PER_SEC);        
        //save gData
        //printf("igdName %s", igdName);
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
                fwrite(gData[j], sizeof(struct igd_data2), counts[j], fp);
                fclose(fp); 
                //free(gData[j]); 
            }
        }   
        for(j=0;j<nTiles;j++){;
            if(counts[j]>0)
                free(gData[j]);   
        }
        free(gData);
        //end = clock();    
        //printf("Saving time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
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
    //end = clock();    
    //printf("TSV saved: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
    //-------------------------------------------------------------------------
    //Reload tile data, sort and save them into a single file
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");    
    FILE *fp1 = fopen(idFile, "wb"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t nrec;  
    FILE *fp0;        
    i=1122332222;   //type 6 data
    fwrite(&i, sizeof(uint32_t), 1, fp1);   
    fwrite(Counts, sizeof(uint32_t), nTiles, fp1);  //
    char iname[256]; 
    struct igd_data2 *gdata, *gdata0;      
    for(i=0;i<nTiles;i++){
        if(Counts[i]>0){
            k = g2ichr[i];
            sprintf(iname, "%s%s%s/%s_%u%s", oPath, "data0/", folder[k], igdName, i-gstart[k], ".igd");
            fp0 = fopen(iname, "rb");
            if(fp0!=NULL){   
                nrec = Counts[i];
                gdata0 = malloc(nrec*sizeof(struct igd_data2));
                gdata  = malloc((nrec+2)*sizeof(struct igd_data2));                
                fread(gdata0, sizeof(struct igd_data2), nrec, fp0);
                //-----construct gdata6
                constructAIList(gdata0, nrec, gdata);
                //----------------------------------------------------------
                //append the data to fp1
                fwrite(gdata, sizeof(struct igd_data2), nrec+2, fp1);
                free(gdata);
                free(gdata0);
            }      
            fclose(fp0);
            if(mode==0)
                remove(iname);              
        }
    }
    fclose(fp1); 
    //end = clock();    
    //printf("igd_w finished: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);    
    //-------------------------------------------------------------------------       
    free(splits);  
    free(counts);
    free(Counts);
    globfree(&gResult); 
}

//create ucsc igd from plain text files
void create_igd(char *iPath, char *oPath, char *igdName, int mode)
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
    if(n_files<1){   
        printf("No data file: %u\n", n_files); 
        return;
    }
     
    uint32_t *nd = calloc(n_files, sizeof(uint32_t));   //number of datasets
    float *md = calloc(n_files, sizeof(float));         //mean 

    //2. Read region data
    uint32_t i, ii, j, k, t, df1, df2, df4, ichr, n1, n2, ti, nR;
    uint32_t *counts = calloc(nTiles, sizeof(uint32_t));    //134217728*16=2G
    uint32_t *Counts = calloc(nTiles, sizeof(uint32_t));    //total
    double delta;    
    char idFile[256];  
    clock_t start, end;
    start = clock();  
    end = clock();    
    //printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);   
    
    //-------------------------------------------------------------------------
    char *ftype;
    int nCols=16, err;                    
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
        //printf("start: %u\n", i0);
        while(i<n_files && m==0){   //m>0 defines breaks when reading a big file         
            ftype = file_ids[i] + strlen(file_ids[i]) - 4;
            if(strcmp(".bed", ftype)==0){        
                //a. Prepare: get the counts               
                //printf("%s", file_ids[i]);
                FILE *fp;
                fp = fopen(file_ids[i], "r");
                if (!fp) {
                    fprintf (stderr, "open of '%s' failed: %s.\n", file_ids[i],
                             strerror (errno));
                    exit (EXIT_FAILURE);
                }             
                
                //printf("%u %u\t", i, (uint32_t)cTotal);
                
                nL = 0; 
                if(i==i0 && L0>0){   //pass n0 lines of a big file
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++;              
                }                 
                while(m==0 && fgets(buffer, bgz_buf, fp)!=NULL){
                    nL++;
                    //splits = str_split_t(buffer, nCols); 
                    str_splits(buffer, &nCols, splits);  
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);               
                        df2 = (uint32_t)atoi(splits[2]); 
                        n1 = df1/nbp;
                        n2 = df2/nbp-n1;    //tile start at boundary   
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
                fclose(fp);      
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
        if(gData==NULL){
            printf("Error: memory allocation with igd_data"); 
            return;
        }       
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("(%u, %u) (%u, %u) \n", i0, L0, i, L1);           
        for(ii=i0; ii<i; ii++){   //n>0 defines breaks when reading a big file
            ftype = file_ids[ii] + strlen(file_ids[ii]) - 4;
            if(strcmp(".bed", ftype)==0){        
                FILE *fp = fopen (file_ids[ii], "r"); 
                nL = 0; 
                if(ii==i0 && L0>0){   //skip n0 lines of a big file
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++;              
                }                
                while(fgets(buffer, bgz_buf, fp)!=NULL){//tbd: read 1Mb each time
                    //splits = str_split(buffer,'\t', &nCols);  
                    str_splits(buffer, &nCols, splits);
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
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
                fclose(fp);               
            }   //.bed file
        }   //ii
        //--------------------------------------------------------------------- 
        if(m>0){
            ii=i;  //m>0 defines breaks when reading a big file
            //ftype = file_ids[ii] + strlen(file_ids[ii]) - 7;
            //if(strcmp(".bed.gz", ftype)==0 || strcmp(".txt.gz", ftype)==0){ 
            ftype = file_ids[i] + strlen(file_ids[i]) - 4;
            if(strcmp(".bed", ftype)==0){    
                FILE *fp = fopen (file_ids[ii], "r"); 
                nL = 0;
                if(ii==i0 && L0>0){
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++; 
                }                            
                while(nL<L1 && fgets(buffer, bgz_buf, fp)!=NULL){//tbd: read 1Mb each time
                    //splits = str_split(buffer,'\t', &nCols);  
                    str_splits(buffer, &nCols, splits);
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
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
                fclose(fp);               
            }   //bed file  
        }    
        //---------------------------------------------------------------------           
        end = clock();    
        //printf("File %u processing time: %f \n", i, ((double)(end-start))/CLOCKS_PER_SEC);        
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
        //printf("Saving time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
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
    if(fpi==NULL){
        printf("Can't open file %s", idFile);
        return;
    }
        
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
    //printf("TSV saved: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
    //-------------------------------------------------------------------------
    //Reload tile data, sort and save them into a single file
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");    
    FILE *fp1 = fopen(idFile, "wb"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t nrec;  
    FILE *fp0;    
    i=1122330000;   //type 0 data
    fwrite(&i, sizeof(uint32_t), 1, fp1);  
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
            if(mode==0)
                remove(iname);              
        }
    }
    fclose(fp1); 
    end = clock();    
    //printf("igd_w finished: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);    
    //-------------------------------------------------------------------------    
    free(splits);  
    free(counts);
    free(Counts);
    globfree(&gResult);
}

//create ucsc igd from plain text files
void create_igd1(char *iPath, char *oPath, char *igdName, int mode)
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
    if(n_files<1){   
        printf("Too few files (add to path /*): %u\n", n_files); 
        return;
    }
     
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
    //printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);   
    
    //-------------------------------------------------------------------------
    char *ftype;
    int nCols=16, err;                    
    int tlen, bytes_read;
    unsigned char buffer[bgz_buf]; 
    uint64_t cTotal, nL;
    //1. find record lines (from one part of big file or many samll files) of size ~maxCount
    //2. procese and save (n0: starting file index, n1 ending file' n0 start line, n2 ending line)
    //(i0, L0)-->(i1, L1): i1 can be the same as i0 if filesize>maxCount
    uint32_t i0=0, i1=0, L0=0, L1=1, m;  
    struct igd_data1 **gData; 
    char **splits = malloc((nCols+1)*sizeof(char *)); 
    while(i0<n_files){
        //1. start from (i0, L0): find (i1, L1)
        cTotal = 0;
        i = i0; 
        m = 0;        
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("start: %u\n", i0);
        while(i<n_files && m==0){   //n>0 defines breaks when reading a big file         
            ftype = file_ids[i] + strlen(file_ids[i]) - 4;
            if(strcmp(".bed", ftype)==0){        
                //a. Prepare: get the counts               
                //printf("%s", file_ids[i]);
                FILE *fp;
                fp = fopen(file_ids[i], "r");
                if (!fp) {
                    fprintf (stderr, "open of '%s' failed: %s.\n", file_ids[i],
                             strerror (errno));
                    exit (EXIT_FAILURE);
                }             
                
                //printf("%u %u\t", i, (uint32_t)cTotal);
                
                nL = 0; 
                if(i==i0 && L0>0){   //pass n0 lines of a big file
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++;              
                }                 
                while(m==0 && fgets(buffer, bgz_buf, fp)!=NULL){
                    nL++;
                    //splits = str_split_t(buffer, nCols); 
                    str_splits(buffer, &nCols, splits);  
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
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
                fclose(fp);      
            }   //if bed.gz
            if(m==0)
                i++;
        }   //while i<n_files
        //---------------------------------------------------------------------
        //2. process files i0--i           
        gData = malloc(nTiles*sizeof(struct igd_data1*));
        for(j=0; j<nTiles; j++){
            if(counts[j]>0)
                gData[j] = malloc(counts[j]*sizeof(struct igd_data1));  
        }  
        if(gData==NULL){
            printf("Error: memory allocation with igd_data"); 
            return;
        }       
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("(%u, %u) (%u, %u) \n", i0, L0, i, L1);           
        for(ii=i0; ii<i; ii++){   //n>0 defines breaks when reading a big file
            ftype = file_ids[ii] + strlen(file_ids[ii]) - 4;
            if(strcmp(".bed", ftype)==0){        
                FILE *fp = fopen (file_ids[ii], "r"); 
                nL = 0; 
                if(ii==i0 && L0>0){   //pass n0 lines of a big file
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++;              
                }                
                while(fgets(buffer, bgz_buf, fp)!=NULL){//tbd: read 1Mb each time
                    //splits = str_split(buffer,'\t', &nCols);  
                    str_splits(buffer, &nCols, splits);
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
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
                                counts[ti]++;  
                            }
                        }
                    }  
                    //free(splits);
                }   //while gzgets           
                fclose(fp);               
            }   //.bed file
        }   //ii
        //--------------------------------------------------------------------- 
        if(m>0){
            ii=i;  //m>0 defines breaks when reading a big file
            //ftype = file_ids[ii] + strlen(file_ids[ii]) - 7;
            //if(strcmp(".bed.gz", ftype)==0 || strcmp(".txt.gz", ftype)==0){ 
            ftype = file_ids[i] + strlen(file_ids[i]) - 4;
            if(strcmp(".bed", ftype)==0){    
                FILE *fp = fopen (file_ids[ii], "r"); 
                nL = 0;
                if(ii==i0 && L0>0){
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++; 
                }                            
                while(nL<L1 && fgets(buffer, bgz_buf, fp)!=NULL){//tbd: read 1Mb each time
                    //splits = str_split(buffer,'\t', &nCols);  
                    str_splits(buffer, &nCols, splits);
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
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
                                counts[ti]++;  
                            }
                        }
                    }  
                    //free(splits);
                    nL++;
                }            
                fclose(fp);               
            }   //bed file  
        }    
        //---------------------------------------------------------------------           
        end = clock();    
        //printf("File %u processing time: %f \n", i, ((double)(end-start))/CLOCKS_PER_SEC);        
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
                fwrite(gData[j], sizeof(struct igd_data1), counts[j], fp);
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
        //printf("Saving time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
        gData = NULL;                                      
        //---------------------------------------------------------------------
        i0 = i; 
        L0 = L1;
        L1 = 0;
    }  
    //save _index.tsv: 4 columns--index, filename, nd, md
    //Also has a header line: 
    //printf("Line 1803\n");
    
    char *tchr;   
    sprintf(idFile, "%s%s%s", oPath, igdName, "_index.tsv");    
    FILE *fpi = fopen(idFile, "w");
    if(fpi==NULL){
        printf("Can't open file %s", idFile);
        return;
    }
        
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
    //printf("TSV saved: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
    //-------------------------------------------------------------------------
    //Reload tile data, sort and save them into a single file
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");    
    FILE *fp1 = fopen(idFile, "wb"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t nrec;  
    FILE *fp0;  
    i=1122331111;   //type 3 data
    fwrite(&i, sizeof(uint32_t), 1, fp1);           
    fwrite(Counts, sizeof(uint32_t), nTiles, fp1);
    char iname[256]; 
    struct igd_data1 *gdata;  
              
    for(i=0;i<nTiles;i++){
        if(Counts[i]>0){
            k = g2ichr[i];
            sprintf(iname, "%s%s%s/%s_%u%s", oPath, "data0/", folder[k], igdName, i-gstart[k], ".igd");
            fp0 = fopen(iname, "rb");
            if(fp0!=NULL){   
                nrec = Counts[i];
                gdata = malloc(nrec*sizeof(struct igd_data1));
                fread(gdata, sizeof(struct igd_data), nrec, fp0);
                qsort(gdata, nrec, sizeof(struct igd_data1), compare_rend1);
                //append the data to fp1
                fwrite(gdata, sizeof(struct igd_data1), nrec, fp1);
                free(gdata);
            }      
            fclose(fp0);
            if(mode==0)
                remove(iname);              
        }
    }  
    
    fclose(fp1); 
    end = clock();    
    //printf("igd_w finished: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);    
    //-------------------------------------------------------------------------    
    free(splits);  
    free(counts);
    free(Counts);
    globfree(&gResult);
}

//create ucsc igd from plain text files
void create_igd2(char *iPath, char *oPath, char *igdName, int mode)
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
    if(n_files<1){   
        printf("Too few files (add to path /*): %u\n", n_files); 
        return;
    }
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
    //printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);   
    
    //-------------------------------------------------------------------------
    char *ftype;
    int nCols=16, err;                    
    int tlen, bytes_read;
    unsigned char buffer[bgz_buf]; 
    uint64_t cTotal, nL;
    //1. find record lines (from one part of big file or many samll files) of size ~maxCount
    //2. procese and save (n0: starting file index, n1 ending file' n0 start line, n2 ending line)
    //(i0, L0)-->(i1, L1): i1 can be the same as i0 if filesize>maxCount
    uint32_t i0=0, i1=0, L0=0, L1=1, m;  
    struct igd_data2 **gData; 
    char **splits = malloc((nCols+1)*sizeof(char *)); 
    while(i0<n_files){
        //1. start from (i0, L0): find (i1, L1)               
        //printf("%s", file_ids[i0]);
        cTotal = 0;
        i = i0; 
        m = 0;              
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("start: %u\n", i0);
        while(i<n_files && m==0){   //n>0 defines breaks when reading a big file         
            ftype = file_ids[i] + strlen(file_ids[i]) - 4;
            if(strcmp(".bed", ftype)==0){        
                //a. Prepare: get the counts               

                FILE *fp;
                fp = fopen(file_ids[i], "r");
                if (!fp) {
                    fprintf (stderr, "open of '%s' failed: %s.\n", file_ids[i],
                             strerror (errno));
                    exit (EXIT_FAILURE);
                }             
                
                //printf("%u %u\t", i, (uint32_t)cTotal);           
                nL = 0; 
                if(i==i0 && L0>0){   //skip n0 lines of a big file
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++;              
                }                 
                while(m==0 && fgets(buffer, bgz_buf, fp)!=NULL){
                    nL++;
                    //splits = str_split_t(buffer, nCols); 
                    str_splits(buffer, &nCols, splits);  
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
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
                fclose(fp);      
            }   //if bed.gz
            else{
                printf("File %s type not .bed\n", file_ids[i]);
            }
            if(m==0)
                i++;
        }   //while i<n_files
        //---------------------------------------------------------------------
        //2. process files i0--i           
        gData = malloc(nTiles*sizeof(struct igd_data2*));
        for(j=0; j<nTiles; j++){
            if(counts[j]>0){
                gData[j] = malloc(counts[j]*sizeof(struct igd_data2));
            }  
        }  
        if(gData==NULL){
            printf("Error: memory allocation with igd_data"); 
            return;
        }       
        memset(counts, 0, nTiles*sizeof(uint32_t));
        //printf("(%u, %u) (%u, %u) \n", i0, L0, i, L1);           
        for(ii=i0; ii<i; ii++){   //n>0 defines breaks when reading a big file
            ftype = file_ids[ii] + strlen(file_ids[ii]) - 4;
            if(strcmp(".bed", ftype)==0){        
                FILE *fp = fopen (file_ids[ii], "r"); 
                nL = 0; 
                if(ii==i0 && L0>0){   //pass n0 lines of a big file
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++;              
                }                
                while(fgets(buffer, bgz_buf, fp)!=NULL){//tbd: read 1Mb each time
                    //splits = str_split(buffer,'\t', &nCols);  
                    str_splits(buffer, &nCols, splits);
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
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
                                counts[ti]++;  
                            }
                        }
                    }  
                    //free(splits);
                }   //while gzgets           
                fclose(fp);               
            }   //.bed file
        }   //ii
        //--------------------------------------------------------------------- 
        if(m>0){
            ii=i;  //m>0 defines breaks when reading a big file
            //ftype = file_ids[ii] + strlen(file_ids[ii]) - 7;
            //if(strcmp(".bed.gz", ftype)==0 || strcmp(".txt.gz", ftype)==0){ 
            ftype = file_ids[i] + strlen(file_ids[i]) - 4;
            if(strcmp(".bed", ftype)==0){    
                FILE *fp = fopen (file_ids[ii], "r"); 
                nL = 0;
                if(ii==i0 && L0>0){
                    while(nL<L0 && fgets(buffer, bgz_buf, fp)!=NULL)
                        nL++; 
                }                            
                while(nL<L1 && fgets(buffer, bgz_buf, fp)!=NULL){//tbd: read 1Mb each time
                    //splits = str_split(buffer,'\t', &nCols);  
                    str_splits(buffer, &nCols, splits);
                    ichr = -1;
                    tlen = strlen(splits[0]);
                    if(tlen<6 && tlen>3){
                        if(strcmp(splits[0], "chrX")==0)
                            ichr = 22;
                        else if(strcmp(splits[0], "chrY")==0)
                            ichr = 23;
                        else if(strcmp(splits[0], "chrM")==0)
                            ichr = 24;                           
                        else{
                            rtn = atoi(&splits[0][3]);
                            if(rtn!=0)
                                ichr = (uint32_t)(rtn-1);
                        }
                    }
                    else{
                        k=25;
                        while(k<94 && strcmp(splits[0], folder[k])!=0)
                            k++;
                        if(k<94)
                            ichr = k;
                    }
                    if(ichr>=0){
                        df1 = (uint32_t)atoi(splits[1]);
                        df2 = (uint32_t)atoi(splits[2]);
                        nd[ii]++;                        
                        md[ii] += df2-df1;
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
                                counts[ti]++;  
                            }
                        }
                    }  
                    //free(splits);
                    nL++;
                }            
                fclose(fp);               
            }   //bed file  
        }    
        //---------------------------------------------------------------------           
        end = clock();    
        //printf("File %u processing time: %f \n", i, ((double)(end-start))/CLOCKS_PER_SEC);        
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
                fwrite(gData[j], sizeof(struct igd_data2), counts[j], fp);
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
        //printf("Saving time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
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
    if(fpi==NULL){
        printf("Can't open file %s", idFile);
        return;
    }
        
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
    //printf("TSV saved: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC); 
    //-------------------------------------------------------------------------
    //Reload tile data, sort and save them into a single file
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");    
    FILE *fp1 = fopen(idFile, "wb"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t nrec;  
    FILE *fp0;   
    i=1122332222;   //type 6 data
    fwrite(&i, sizeof(uint32_t), 1, fp1);          
    fwrite(Counts, sizeof(uint32_t), nTiles, fp1);
    char iname[256]; 
    struct igd_data2 *gdata, *gdata0; 
    for(i=0;i<nTiles;i++){
        if(Counts[i]>0){

            k = g2ichr[i];
            sprintf(iname, "%s%s%s/%s_%u%s", oPath, "data0/", folder[k], igdName, i-gstart[k], ".igd");
            fp0 = fopen(iname, "rb");
            if(fp0!=NULL){   
                nrec = Counts[i];
                gdata0 = malloc(nrec*sizeof(struct igd_data2));
                gdata  = malloc((nrec+2)*sizeof(struct igd_data2));                
                fread(gdata0, sizeof(struct igd_data2), nrec, fp0);
                //-----construct gdata6
                constructAIList(gdata0, nrec, gdata);
                //----------------------------------------------------------
                //append the data to fp1
                fwrite(gdata, sizeof(struct igd_data2), nrec+2, fp1);
                free(gdata);
                free(gdata0);                
            } 
            fclose(fp0);            
            if(mode==0)
                remove(iname);     
        }
    }
    fclose(fp1); 
    end = clock();    
    //printf("igd_w finished: time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);        
    free(splits);  
    free(counts);
    free(Counts);
    globfree(&gResult);
}

//-------------------------------------------------------------------------------------
int igd_create(int argc, char **argv)
{
    if (argc < 5) 
        errx(1, "usage:\t%s <input path><ouput path><dbName><option>\n", argv[0]);       
    //convert block index to chr index for convenience
    g2ichr = malloc(nTiles*sizeof(uint32_t));
    uint32_t i, j;
    for(i=0; i<24; i++){  
        for(j=gstart[i]; j<gstart[i+1]; j++)      
	    g2ichr[j] = i;
    }
    char ipath[128];
    char opath[128];
    strcpy(ipath, argv[2]);
    strcpy(opath, argv[3]);
    char *dbname = argv[4]; 
    int mode = 0;//mode 1: save tile data
    int dtype = 0;//default:g_data 4 items, -1 old, 0 [default], 1 [3 items], 2[ailist]
    int gztxt = 0;//input type .gz or .bed /txt
    
    for(i=5; i<argc; i++){
        if(strcmp(argv[i], "-t")==0)
            gztxt = 1; 
        else if(strcmp(argv[i], "-m")==0)
            mode = 1; 
        else if(strcmp(argv[i], "-s")==0 && i+1<argc)//data structure type
            dtype = atoi(argv[i+1]);              
    }
                                   
    if(opath[strlen(opath)-1]!='/'){
        strcat(opath, "/");
    }
    if(ipath[strlen(ipath)-1]=='/'){
        strcat(ipath, "*");
    }
    else if(ipath[strlen(ipath)-1]!='*'){
        strcat(ipath, "/*");
    }
    
    //check if the subfolders exist:    
    char ftmp[128];      
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
        for (i=0;i<24;i++){	    
            sprintf(ftmp, "%s%s%s", opath, "data0/", folder[i]);
            if (stat(ftmp, &st) == -1)   
                mkdir(ftmp, 0777);
        }
        if(dtype==0){
            if(gztxt==1)
                create_igd(ipath, opath, dbname, mode); //test file
            else
                create_igd_gz(ipath, opath, dbname, mode); 
        }
        else if(dtype==1){
            if(gztxt==1)
                create_igd1(ipath, opath, dbname, mode); //test file
            else
                create_igd_gz1(ipath, opath, dbname, mode); 
        } 
        else if(dtype==2){
            if(gztxt==1)
                create_igd2(ipath, opath, dbname, mode); //test file
            else
                create_igd_gz2(ipath, opath, dbname, mode); 
        }               
    } 

    free(g2ichr);
    return EX_OK;
}

