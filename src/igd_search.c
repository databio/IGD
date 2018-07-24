//=====================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//
//time ./igd_search Test110000.bed /media/john/Extra/ucsc_igd/ucsc.igd
//-------------------------------------------------------------------------------------
#include "igd_search.h"
//-------------------------------------------------------------------------------------

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

struct igd_mix* get_overlaps(struct query_data *query, uint32_t nblocks, uint32_t nmax, uint32_t* nOL)
{
    uint32_t i, j, i1, i2, bk, ichr, k, nrec, nols=0, q1, q2, m;
    char iname[128];
    struct igd_mix *overlaps = malloc(nmax*sizeof(struct igd_mix));
    struct igd_data *gdata;
    m = 0;
    while (m<nblocks) {    
        i1 = m;
        bk = query[i1].q_idx;
        while(m<nblocks-1 && query[m+1].q_idx==bk)
            m++;
        i2 = m;
        //-----------------------------------------------------------------------------
        ichr = (uint8_t)g2ichr[bk];
        k = bk - gstart[ichr];
        sprintf(iname, "%s%s%s%i%s", "igdata/",folder[ichr], "/b14_", k, ".igd");
        FILE *fp = fopen(iname, "rb");

        if(fp!=NULL){
            fseek(fp, 0, SEEK_END);
            nrec = ftell(fp)/sizeof(struct igd_data);
            fseek(fp, 0, SEEK_SET);
            if(nrec>0){
                gdata = malloc(nrec*sizeof(struct igd_data));
                fread(gdata, sizeof(struct igd_data), nrec, fp);
         	    //printf("%i %i %i \n", nrec, gdata[0].r_start, gdata[0].r_end);               
                for(i=i1; i<=i2; i++){
                    q1 = query[i].r_start;
                    q2 = query[i].r_end;
                    for(j=0;j<nrec;j++){
                        if(!(q1>gdata[j].r_end || q2<gdata[j].r_start)){
                    			overlaps[nols].igd = gdata[j];
                    			overlaps[nols].m_idx   = bk;
                    			nols++;
                        }
                    }	    
                }
                free(gdata);
            }
            fclose(fp);
        }
        m++;
    }
    *nOL = nols;
    return overlaps;
}

//using single-file igd_data
struct igd_mix* get_overlaps_w(struct query_data *query, uint32_t nblocks, char *igdName, uint32_t nmax, uint32_t* nOL)
{   //add alignment: padding 
    //assume in-tile igdata is sorted by region end 
    FILE* fp = fopen(igdName, "rb");
    if(!fp)
        return NULL;    
        
    uint32_t i, j, t1, t2, i1, i2, bk, ichr, nols=0, q1, q2, m;
    struct igd_mix *overlaps = malloc(nmax*sizeof(struct igd_mix));
    struct igd_data *gdata;
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);//number of struct
    uint64_t *mloc = malloc((nTiles+1)*sizeof(uint64_t));
    //fseek(fp, 0, SEEK_SET);
    fread(counts, sizeof(uint32_t), nTiles, fp);   
    for(i=0; i<nTiles; i++)
        mloc[i+1] = counts[i]*16;
   
    mloc[0]=len0;
    for(i=1; i<nTiles; i++)
        mloc[i] += mloc[i-1];    
    m = 0;
    while (m<nblocks) {    
        i1 = m;
        bk = query[i1].q_idx;
        while(m<nblocks-1 && query[m+1].q_idx==bk)
            m++;
        i2 = m;
        if(bk<nTiles && counts[bk]>0){      
            ichr = (uint8_t)g2ichr[bk];
            fseek(fp, mloc[bk], SEEK_SET);
            gdata = malloc(counts[bk]*16);
            fread(gdata, 16, counts[bk], fp);
            //update the in-tile db start j0: not really faster
            for(i=i1; i<=i2; i++){
        	    q1 = query[i].r_start;
        	    q2 = query[i].r_end;             
            	for(j=0;j<counts[bk];j++){
          		    t2 = gdata[j].r_end;
        		    if(q1<=t2){
            	       t1 = gdata[j].r_start;
        		       if(q2>=t1){    		          		    
                    		overlaps[nols].igd = gdata[j];
                    		overlaps[nols].m_idx   = bk;
                			nols++;
            			}
        		    }
        		}	    
    	    }
    	    free(gdata);
        }
        m++;
    }  
    //check for everlaps: only needed for duplicated queries       
    fclose(fp);
    *nOL = nols;
    free(counts);
    free(mloc);
    return overlaps;
}

//using single-file igd_data
uint64_t get_overlaps_n(char *qfName, char *igdName, uint32_t *nregions, double *mean_size, uint32_t *hits)
{   //no need to store every overlaps, only get the number of hits
    //assume in-tile igdata is sorted by region end 
    //.Reuse igddata if current query is in the same tile as the previous
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;    
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
     
    int ichr;      
    uint32_t i, j, t1, t2, q1, q2, m;
    uint32_t n1, n2, idx, idx0, nRegions=0, nCols = 5, bd;  
    struct igd_data *gdata = malloc(1*sizeof(struct igd_data));
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);//number of struct
    uint64_t *mloc = malloc((nTiles+1)*sizeof(uint64_t));
    //fseek(fp, 0, SEEK_SET);
    fread(counts, sizeof(uint32_t), nTiles, fi);   
    for(i=0; i<nTiles; i++)
        mloc[i+1] = counts[i]*16;
    mloc[0]=len0;
    for(i=1; i<nTiles; i++)
        mloc[i] += mloc[i-1];  
   
   //----read qf for a region--------
    char buf[1024], ch;
    char **splits;
    double delta=0.0;
    uint64_t nols=0;
    idx0 = nTiles+10;
    while(fgets(buf, 1024, fq)!=NULL){	
	//printf("%u %s",(uint32_t)nols, buf); 
        splits = str_split(buf,'\t', &nCols); 
        if(strlen(splits[0])>5)
            ichr = -1;  
        else if(strcmp(splits[0], "chrX")==0)
            ichr = 22;
        else if(strcmp(splits[0], "chrY")==0)
            ichr = 23;
        else{
            ichr = (int)(atoi(&splits[0][3])-1);
        }           
        if(ichr>=0){
            q1  = (uint32_t)atoi(splits[1]);
            q2  = (uint32_t)atoi(splits[2]);
            delta += q2 - q1;
            nRegions++;
            n1 = q1/nbp;
            n2 = q2/nbp-n1;   
            idx = n1 + gstart[ichr];
            //find overlaps with this region  
            if(n2==0){
                if(counts[idx]>0){
                    if(idx!=idx0){
                        free(gdata);
                        fseek(fi, mloc[idx], SEEK_SET);
                        gdata = malloc(counts[idx]*16);
                        fread(gdata, 16, counts[idx], fi);
                        idx0 = idx;
                    }
                    //update the in-tile db start j0: not really faster           
                    for(j=0;j<counts[idx0];j++){
                        t2 = gdata[j].r_end;
                        if(q1<=t2){
                            t1 = gdata[j].r_start;
                            if(q2>=t1){    		          		    
                                hits[gdata[j].i_idx]++;
                                nols++;
                            }
                        }
                    } 
                }          
            }
            else{ 
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(idx+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    if(counts[m]>0){
                        if(m!=idx0){
                            free(gdata);
                            fseek(fi, mloc[m], SEEK_SET);
                            gdata = malloc(counts[m]*16);
                            fread(gdata, 16, counts[m], fi);
                            idx0 = m;
                        }
                        //update the in-tile db start j0: not really faster           
                        for(j=0;j<counts[m];j++){
                            t2 = gdata[j].r_end;
                            if(q1<=t2 && t2<bd){
                                t1 = gdata[j].r_start;
                                if(q2>=t1){ 
                                    hits[gdata[j].i_idx]++; 
                                    nols++;
                                }
                            }
                        }
                    }                   
                    bd += nbp;
                }	
                //--last tile: normal------------------------------
                m=idx+n2;
                if(counts[m]>0){
                    if(m!=idx0){
                        free(gdata);
                        fseek(fi, mloc[m], SEEK_SET);
                        gdata = malloc(counts[m]*16);
                        fread(gdata, 16, counts[m], fi);
                        idx0 = m;
                    }
                    //update the in-tile db start j0: not really faster           
                    for(j=0;j<counts[m];j++){
                        t2 = gdata[j].r_end;
                        if(q1<=t2){
                            t1 = gdata[j].r_start;
                            if(q2>=t1){    		          		    
                                hits[gdata[j].i_idx]++; 
                                nols++;
                            }
                        }
                    }
                }
            }   
        }   //if
        free(splits);     
    }   //while  
    free(gdata);
    *mean_size = delta/nRegions;
    *nregions = nRegions;
    fclose(fq);   
    fclose(fi);
    free(counts);
    free(mloc);
    return nols;
}

//using single-file igd_data: dynamic
uint64_t get_overlaps_v(char *qfName, char *igdName, uint32_t v, uint32_t *nregions, double *mean_size, uint32_t *hits)
{   //no need to store every overlaps, only get the number of hits
    //assume in-tile igdata is sorted by region end 
    //.Reuse igddata if current query is in the same tile as the previous
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;    
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
     
    int ichr;      
    uint32_t i, j, t1, t2, q1, q2, m;
    uint32_t n1, n2, idx, idx0, nRegions=0, nCols = 5, bd;  
    struct igd_data *gdata = malloc(1*sizeof(struct igd_data));
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);//number of struct
    uint64_t *mloc = malloc((nTiles+1)*sizeof(uint64_t));
    //fseek(fp, 0, SEEK_SET);
    fread(counts, sizeof(uint32_t), nTiles, fi);   
    for(i=0; i<nTiles; i++)
        mloc[i+1] = counts[i]*16;
    mloc[0]=len0;
    for(i=1; i<nTiles; i++)
        mloc[i] += mloc[i-1];  
   
   //----read qf for a region--------
    char buf[1024], ch;
    char **splits;
    double delta=0.0;
    uint64_t nols=0;
    idx0 = nTiles+10;
    while(fgets(buf, 1024, fq)!=NULL){	
	//printf("%u %s",(uint32_t)nols, buf); 
        splits = str_split(buf,'\t', &nCols); 
        if(strlen(splits[0])>5)
            ichr = -1;  
        else if(strcmp(splits[0], "chrX")==0)
            ichr = 22;
        else if(strcmp(splits[0], "chrY")==0)
            ichr = 23;
        else{
            ichr = (int)(atoi(&splits[0][3])-1);
        }           
        if(ichr>=0){
            q1  = (uint32_t)atoi(splits[1]);
            q2  = (uint32_t)atoi(splits[2]);
            delta += q2 - q1;
            nRegions++;
            n1 = q1/nbp;
            n2 = q2/nbp-n1;   
            idx = n1 + gstart[ichr];
            //find overlaps with this region  
            if(n2==0){
                if(counts[idx]>0){
                    if(idx!=idx0){
                        free(gdata);
                        fseek(fi, mloc[idx], SEEK_SET);
                        gdata = malloc(counts[idx]*16);
                        fread(gdata, 16, counts[idx], fi);
                        idx0 = idx;
                    }
                    //update the in-tile db start j0: not really faster           
                    for(j=0;j<counts[idx0];j++){
                        t2 = gdata[j].r_end;
                        if(gdata[j].g_val>v && q1<=t2){
                            t1 = gdata[j].r_start;
                            if(q2>=t1){    		          		    
                                hits[gdata[j].i_idx]++;
                                nols++;
                            }
                        }
                    } 
                }          
            }
            else{ 
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(idx+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    if(counts[m]>0){
                        if(m!=idx0){
                            free(gdata);
                            fseek(fi, mloc[m], SEEK_SET);
                            gdata = malloc(counts[m]*16);
                            fread(gdata, 16, counts[m], fi);
                            idx0 = m;
                        }
                        //update the in-tile db start j0: not really faster           
                        for(j=0;j<counts[m];j++){
                            t2 = gdata[j].r_end;
                            if(gdata[j].g_val>v && q1<=t2 && t2<bd){
                                t1 = gdata[j].r_start;
                                if(q2>=t1){ 
                                    hits[gdata[j].i_idx]++; 
                                    nols++;
                                }
                            }
                        }
                    }                   
                    bd += nbp;
                }	
                //--last tile: normal------------------------------
                m=idx+n2;
                if(counts[m]>0){
                    if(m!=idx0){
                        free(gdata);
                        fseek(fi, mloc[m], SEEK_SET);
                        gdata = malloc(counts[m]*16);
                        fread(gdata, 16, counts[m], fi);
                        idx0 = m;
                    }
                    //update the in-tile db start j0: not really faster           
                    for(j=0;j<counts[m];j++){
                        t2 = gdata[j].r_end;
                        if(gdata[j].g_val>v && q1<=t2){
                            t1 = gdata[j].r_start;
                            if(q2>=t1){    		          		    
                                hits[gdata[j].i_idx]++; 
                                nols++;
                            }
                        }
                    }
                }
            }   
        }   //if
        free(splits);     
    }   //while  
    free(gdata);
    *mean_size = delta/nRegions;
    *nregions = nRegions;
    fclose(fq);   
    fclose(fi);
    free(counts);
    free(mloc);
    return nols;
}

void search(char* qfName, char* igdName, uint32_t v)
{   //name standard: igd_file(dbname.igd), index_file(dbname_index.tsv)
    uint32_t i, nq=1, nFiles, nCols=2, genome_size=3095677412;
    double mq = 1.0;
    char tmp[128];
    strcpy(tmp, igdName);
    tmp[strrchr(tmp, '.')-tmp] = '\0';
    char *idFile = tmp;//str_split(tmp, '.', &nCols)[0];
    strcat(idFile, "_index.tsv");     
 
    struct igd_info *fi = get_igdinfo(idFile, &nFiles); 
    printf("nfiles: %u\n", nFiles);  
    
    clock_t start, end;
    start = clock();
    uint32_t *hits = calloc(nFiles, sizeof(uint32_t));

    if(v>0)
        uint64_t nOL = get_overlaps_v(qfName, igdName, v, &nq, &mq, hits);
    else   
        uint64_t nOL = get_overlaps_n(qfName, igdName, &nq, &mq, hits);   
    end = clock();   
    
    printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);
    printf("%u %u %f \n", (uint32_t)nOL, nq, mq);
    printf("index\t File_name\t number of regions\t mean-region-size \t number of hits\n");
    for(i=0;i<10;i++)
        printf("%i %s %u %u %u\n", i, fi[i].fileName, fi[i].nd, (uint32_t)fi[i].md, hits[i]);
    //---------------------------------------------------------------------------------
    free(fi->fileName);
    free(fi);
    free(hits);
}

//-------------------------------------------------------------------------------------
int igd_search(int argc, char **argv)
{   //igd[0] search[1] query100.bed[2] home/john/iGD/rme_igd/roadmap.igd[3]
    uint32_t v = 0;   
    //convert block index to chr index for convenience
    
    g2ichr = malloc(nTiles*sizeof(uint32_t));
    uint32_t i, j;
    for(i=0; i<24; i++){  
        	for(j=gstart[i]; j<gstart[i+1]; j++)      
        	    g2ichr[j] = i;
    }
    char *qfName = argv[2];
    char *igdName = argv[3];
    if(argc>4)
        v = argv[4];
        
    search(qfName, igdName, v);      

    free(g2ichr);
    return EX_OK;
}

