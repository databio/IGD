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
    int i, nfiles = -1;//first row is title
    while((ch=fgetc(fp))!=EOF){
        if(ch=='\n')
        	    nfiles++;
    }  
    struct igd_info *fi = (struct igd_info*)malloc(nfiles*sizeof(struct igd_info));
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
        fi[i].nd = (uint32_t)atoi(splits[2]);
        fi[i].md = (double)atoi(splits[3]);   
        i++;
    }        
    //printf("%s %i %f \n", fNames[6], nd[6], md[6]);  
    //printf("Total file: %i\n", i);
    *nFiles = (uint32_t)nfiles;
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
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;       
    int ichr, tc, tL, tR, tM, tS, tlen, rtn;      
    uint32_t i, j, k, t1, t2, q1, q2, m;
    uint32_t n1, n2, idx, idx0, nRegions=0, nCols = 16, bd;  
    struct igd_data *gdata = malloc(1*sizeof(struct igd_data));
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);//number of struct
    uint64_t *mloc = malloc((nTiles+1)*sizeof(uint64_t));
    //fseek(fp, 0, SEEK_SET);
    fread(counts, sizeof(uint32_t), nTiles, fi);   
    for(i=0; i<nTiles; i++){
        mloc[i+1] = counts[i]*16;//bytes
    }
    mloc[0]=len0;
    for(i=1; i<nTiles; i++)
        mloc[i] += mloc[i-1];  
   
   //----read qf for a region--------
    char buf[1024], ch;
    char **splits;
    double delta=0.0;
    uint64_t nols=0;
    idx0 = nTiles+10;
    
    i=0;
    while(fgets(buf, 1024, fq)!=NULL){	
	//printf("%u %s",(uint32_t)nols, buf); 
        ichr = -1;
        splits = str_split(buf,'\t', &nCols); 
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
            while(k<93 && strcmp(splits[0], folder[k])!=0)
                k++;
            if(k<93)
                ichr = k;
        }          
        if(ichr>=0){
            q1  = (uint32_t)atoi(splits[1]);
            q2  = (uint32_t)atoi(splits[2]);
            delta += q2 - q1;
            nRegions++;
            n1 = q1/nbp;
            n2 = q2/nbp-n1;   
            idx = n1 + gstart[ichr];
            
            if(n2==0){
                tc = counts[idx];
                if(tc>0){
                    if(idx!=idx0){
                        free(gdata);
                        fseek(fi, mloc[idx], SEEK_SET);
                        gdata = malloc(tc*16);
                        fread(gdata, 16, tc, fi);
                        idx0 = idx;
                    }
                    //optimize the search: b-search 
                    if(q1>=gdata[tc-1].r_end){
                        //no overlap:do nothing tL>nc;tR<0
                    } 
                    else if(tc<32){
                        tS=0;
                        while(gdata[tS].r_end<=q1)
                            tS++;
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                            }
                        }                     
                    }
                    else{//search tS: the 1st, from Left, t_end satisfies [tS].end>q1
                        tL=0;   tR=tc-1;  
                        tS = -1;    //no exclusion; tL<nc-1
                        while(tL<tR-1){
                            tM = (tL+tR)/2; 
                            if(gdata[tM].r_end>q1)
                                tR = tM;
                            else
                                tL = tM+1;
                        }
                        if(gdata[tL].r_end>q1)
                            tS = tL;
                        else if(gdata[tR].r_end>q1)
                            tS = tR;
                        //------------------------------
                        //if(tS>0){should be
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                            }
                        } 
                        //}                                      
                    } 
                }          
            }
            else{ //n2!=0: tile start at lower boundary, end below upper boundary (open region)
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(n1+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    tc = counts[m];
                    if(tc>0){
                        if(m!=idx0){
                            free(gdata);
                            fseek(fi, mloc[m], SEEK_SET);
                            gdata = malloc(tc*16);
                            fread(gdata, 16, tc, fi);
                            idx0 = m;
                        }
                        //update the in-tile db start j0: not really faster 
                        if(q1>=gdata[tc-1].r_end){  //!!!0204
                            //printf("Line 560 gdata[tc-1] %u m %u \n", gdata[tc-1].r_end, q1);
                            //no overlap:do nothing tL>nc;tR<0
                        } 
                        else{
                            tS=0;
                            while(gdata[tS].r_end<=q1) 
                                tS++;                                 
                            for(j=tS;j<tc;j++){
                                t2 = gdata[j].r_end;
                                if(t2<bd && q2>gdata[j].r_start){
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
                tc = counts[m];
                if(tc>0){
                    if(m!=idx0){
                        free(gdata);
                        fseek(fi, mloc[m], SEEK_SET);
                        gdata = malloc(tc*16);
                        fread(gdata, 16, tc, fi);
                        idx0 = m;
                    }
                    if(tc<32){
                        tS=0;
                        //while(gdata[tS].r_end<=q1)
                        //    tS++;
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                            }
                        }                                          
                    }
                    else{//half dual-binary search 
                        tS = 0;
                        if(gdata[tL].r_end>q1)
                            tS = tL;
                        else if(gdata[tR].r_end>q1)
                            tS = tR;
                        //------------------------------
                        //if(tS>0){should be
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                            }
                        }                                                            
                    } //else
                }//if tc>0
            }   //else n2>0            
        }   //if n2=0
        free(splits);
        i++;  
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

//using single-file igd_data
uint64_t get_overlaps_n0(char *qfName, char *igdName, uint32_t *nregions, double *mean_size, uint32_t *hits)
{   //no need to store every overlaps, only get the number of hits
    //assume in-tile igdata is sorted by region end 
    //.Reuse igddata if current query is in the same tile as the previous  
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;       
    int ichr, tc, tL, tR, tM, tS, tlen, rtn;      
    uint32_t i, j, k, t1, t2, q1, q2, m;
    uint32_t n1, n2, idx, idx0, nRegions=0, nCols = 16, bd;  
    struct igd_data *gdata = malloc(1*sizeof(struct igd_data));
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);//number of struct
    uint64_t *mloc = malloc((nTiles+1)*sizeof(uint64_t));
    //fseek(fp, 0, SEEK_SET);
    fread(&i, sizeof(uint32_t), 1, fi);
    fread(counts, sizeof(uint32_t), nTiles, fi);   
    for(i=0; i<nTiles; i++){
        mloc[i+1] = counts[i]*16;
    }
    mloc[0]=len0 + 4;
    for(i=1; i<nTiles; i++)
        mloc[i] += mloc[i-1];  
   
   //----read qf for a region--------
    char buf[1024], ch;
    char **splits;
    double delta=0.0;
    uint64_t nols=0;
    idx0 = nTiles+10;
    
    i=0;
    while(fgets(buf, 1024, fq)!=NULL){	
	//printf("%u %s",(uint32_t)nols, buf); 
        //tHits = 0;
        ichr = -1;
        splits = str_split(buf,'\t', &nCols); 
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
            while(k<93 && strcmp(splits[0], folder[k])!=0)
                k++;
            if(k<93)
                ichr = k;
        }                   
        if(ichr>=0){
            q1  = (uint32_t)atoi(splits[1]);
            q2  = (uint32_t)atoi(splits[2]);
            delta += q2 - q1;
            nRegions++;
            n1 = q1/nbp;
            n2 = q2/nbp-n1;   
            idx = n1 + gstart[ichr];
            
            if(n2==0){
                tc = counts[idx];
                if(tc>0){
                    if(idx!=idx0){
                        free(gdata);
                        fseek(fi, mloc[idx], SEEK_SET);
                        gdata = malloc(tc*16);
                        fread(gdata, 16, tc, fi);
                        idx0 = idx;
                    }
                    //optimize the search: b-search 
                    if(q1>=gdata[tc-1].r_end){
                        //no overlap:do nothing tL>nc;tR<0
                    } 
                    else if(tc<32){
                        tS=0;
                        while(gdata[tS].r_end<=q1)
                            tS++;
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                            }
                        }                     
                    }
                    else{//search tS: the 1st, from Left, t_end satisfies [tS].end>q1
                        tL=0;   tR=tc-1;  
                        tS = -1;    //no exclusion; tL<nc-1
                        while(tL<tR-1){
                            tM = (tL+tR)/2; 
                            if(gdata[tM].r_end>q1)
                                tR = tM;
                            else
                                tL = tM+1;
                        }
                        if(gdata[tL].r_end>q1)
                            tS = tL;
                        else if(gdata[tR].r_end>q1)
                            tS = tR;
                        //------------------------------
                        //if(tS>0){should be
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                            }
                        } 
                        //}                                      
                    } 
                }          
            }
            else{ //n2!=0: tile start at lower boundary, end below upper boundary (open region)
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(n1+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    tc = counts[m];
                    if(tc>0){
                        if(m!=idx0){
                            free(gdata);
                            fseek(fi, mloc[m], SEEK_SET);
                            gdata = malloc(tc*16);
                            fread(gdata, 16, tc, fi);
                            idx0 = m;
                        }
                        //update the in-tile db start j0: not really faster 
                        if(q1>=gdata[tc-1].r_end){  //!!!0204
                            //printf("Line 560 gdata[tc-1] %u m %u \n", gdata[tc-1].r_end, q1);
                            //no overlap:do nothing tL>nc;tR<0
                        } 
                        else{
                            tS=0;
                            while(gdata[tS].r_end<=q1) 
                                tS++;                                 
                            for(j=tS;j<tc;j++){
                                t2 = gdata[j].r_end;
                                if(t2<bd && q2>gdata[j].r_start){
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
                tc = counts[m];
                if(tc>0){
                    if(m!=idx0){
                        free(gdata);
                        fseek(fi, mloc[m], SEEK_SET);
                        gdata = malloc(tc*16);
                        fread(gdata, 16, tc, fi);
                        idx0 = m;
                    }
                    if(tc<32){
                        tS=0;
                        //while(gdata[tS].r_end<=q1)
                        //    tS++;
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                                //tHits++;
                            }
                        }                                          
                    }
                    else{
                        tS = 0;//q1<bd
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                            }
                        }                                                            
                    } //else
                }//if tc>0
            }   //else n2>0            
        }   //if n2=0
        free(splits);
        i++;  
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

//using single-file igd_data
uint64_t get_overlaps_n1(char *qfName, char *igdName, uint32_t *nregions, double *mean_size, uint32_t *hits)
{   //no need to store every overlaps, only get the number of hits
    //assume in-tile igdata is sorted by region end 
    //.Reuse igddata if current query is in the same tile as the previous  
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;       
    int ichr, tc, tL, tR, tM, tS, tlen, rtn;      
    uint32_t i, j, k, t1, t2, q1, q2, m;
    uint32_t n1, n2, idx, idx0, nRegions=0, nCols = 16, bd;  
    struct igd_data1 *gdata = malloc(1*sizeof(struct igd_data1));
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);//number of struct
    uint64_t *mloc = malloc((nTiles+1)*sizeof(uint64_t));
    //fseek(fp, 0, SEEK_SET);
    fread(&i, sizeof(uint32_t), 1, fi);
    fread(counts, sizeof(uint32_t), nTiles, fi);   
    for(i=0; i<nTiles; i++){
        mloc[i+1] = counts[i]*12;//bytes
        //if(counts[i]>0)
        //    printf("i %u, counts %u\n", i, counts[i]);
    }
    mloc[0]=len0 + 4;
    for(i=1; i<nTiles; i++)
        mloc[i] += mloc[i-1];  
   
   //----read qf for a region--------
    char buf[1024], ch;
    char **splits;
    double delta=0.0;
    uint64_t nols=0;
    idx0 = nTiles+10;
    
    //int tHits;
    i=0;
    while(fgets(buf, 1024, fq)!=NULL){	
	//printf("%u %s",(uint32_t)nols, buf); 
        //tHits = 0;
        ichr = -1;
        splits = str_split(buf,'\t', &nCols); 
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
            while(k<93 && strcmp(splits[0], folder[k])!=0)
                k++;
            if(k<93)
                ichr = k;
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
            //printf("%i, %i, %u, %u,  %u, %u\n", i, ichr, n1, n2, idx, counts[idx]);
            
            if(n2==0){
                tc = counts[idx];
                if(tc>0){
                    if(idx!=idx0){
                        free(gdata);
                        fseek(fi, mloc[idx], SEEK_SET);
                        gdata = malloc(tc*12);
                        fread(gdata, 12, tc, fi);
                        idx0 = idx;
                    }
                    //optimize the search: b-search 
                    if(q1>=gdata[tc-1].r_end){
                        //no overlap:do nothing tL>nc;tR<0
                    } 
                    else if(tc<32){
                        tS=0;
                        while(gdata[tS].r_end<=q1)
                            tS++;
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                                //tHits++;
                            }
                        }                     
                    }
                    else{//search tS: the 1st, from Left, t_end satisfies [tS].end>q1
                        tL=0;   tR=tc-1;  
                        tS = -1;    //no exclusion; tL<nc-1
                        while(tL<tR-1){
                            tM = (tL+tR)/2; 
                            if(gdata[tM].r_end>q1)
                                tR = tM;
                            else
                                tL = tM+1;
                        }
                        if(gdata[tL].r_end>q1)
                            tS = tL;
                        else if(gdata[tR].r_end>q1)
                            tS = tR;
                        //------------------------------
                        //if(tS>0){should be
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                                //tHits++;
                            }
                        } 
                        //}                                      
                    } 
                }          
            }
            else{ //n2!=0: tile start at lower boundary, end below upper boundary (open region)
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(n1+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    tc = counts[m];
                    if(tc>0){
                        if(m!=idx0){
                            free(gdata);
                            fseek(fi, mloc[m], SEEK_SET);
                            gdata = malloc(tc*12);
                            fread(gdata, 12, tc, fi);
                            idx0 = m;
                        }
                        //update the in-tile db start j0: not really faster 
                        if(q1>=gdata[tc-1].r_end){  //!!!0204
                            //printf("Line 560 gdata[tc-1] %u m %u \n", gdata[tc-1].r_end, q1);
                            //no overlap:do nothing tL>nc;tR<0
                        } 
                        else{
                            tS=0;
                            while(gdata[tS].r_end<=q1) 
                                tS++;                                 
                            for(j=tS;j<tc;j++){
                                t2 = gdata[j].r_end;
                                if(t2<bd && q2>gdata[j].r_start){
                                    hits[gdata[j].i_idx]++;                               
                                    nols++;
                                    //tHits++;
                                }
                            }
                        }
                    }                   
                    bd += nbp;
                }	
                //--last tile: normal------------------------------
                m=idx+n2;
                tc = counts[m];
                if(tc>0){
                    if(m!=idx0){
                        free(gdata);
                        fseek(fi, mloc[m], SEEK_SET);
                        gdata = malloc(tc*12);
                        fread(gdata, 12, tc, fi);
                        idx0 = m;
                    }
                    //update the in-tile db start j0: not really faster 
                    //optimize the search: b-search 
                    //if(q1>=gdata[tc-1].r_end){
                        //no overlap:do nothing tL>nc;tR<0
                    //} 
                    //else 
                    if(tc<32){
                        tS=0;
                        //while(gdata[tS].r_end<=q1)
                        //    tS++;
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                                //tHits++;
                            }
                        }                                          
                    }
                    else{//half dual-binary search 
                        tS=0;
                        //------------------------------
                        //if(tS>0){should be
                        for(j=tS;j<tc;j++){
                            if(q2>gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;                               
                                nols++;
                                //tHits++;
                            }
                        }                                                            
                    } //else
                }//if tc>0
            }   //else n2>0            
        }   //if n2=0
        free(splits);
        //if(tHits>0)  
        //    printf("%i, %i, %u, %u, %i\n", i, ichr, q1, q2, tHits);
        i++;  
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

//assume .end not inclusive:component
int bSearch(struct igd_data2* As, int idxS, int idxE, uint32_t qe)
{   //find the index of the first .s<qe start from right
    int tL=idxS, tR=idxE-1, tM, tE=-1;
    if(As[tR].r_start < qe)
        return tR;
    else if(As[tL].r_start >= qe)
        return -1;
    while(tL<tR-1){
        tM = (tL+tR)/2; 
        if(As[tM].r_start >= qe)
            tR = tM-1;
        else
            tL = tM;
    }
    if(As[tR].r_start < qe)
        tE = tR;
    else if(As[tL].r_start < qe)
        tE = tL;       
    return tE;   //tE: index of the first item satisfying .s<qe from right
}

//using single-file igd_data: in-bin ailist search
uint64_t get_overlaps_n2(char *qfName, char *igdName, uint32_t *nregions, double *mean_size, uint32_t *hits)
{   //no need to store every overlaps, only get the number of hits
    //assume in-tile igdata is sorted by region end 
    //.Reuse igddata if current query is in the same tile as the previous  
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;       
    int ichr, t, tc, tL, tR, tM, tS, tlen, rtn;      
    uint32_t i, j, k, t1, t2, q1, q2, m;
    uint32_t n1, n2, idx, idx0, nRegions=0, nCols = 5, bd;  
    struct igd_data2 *gdata = malloc(1*sizeof(struct igd_data2));
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);
    uint64_t *mloc = malloc((nTiles+1)*sizeof(uint64_t));
    //fseek(fp, 0, SEEK_SET);
    fread(&i, sizeof(uint32_t), 1, fi);
    //printf("i = %u\n", i);
    fread(counts, sizeof(uint32_t), nTiles, fi);   
    for(i=0; i<nTiles; i++){
        if(counts[i]>0)
            mloc[i+1] = counts[i]*16 + 32;
        else
            mloc[i+1] = 0;
    }
    mloc[0] = len0 + 4;//first 4 bytes
    for(i=1; i<nTiles; i++)
        mloc[i] += mloc[i-1];  
   
   //----read qf for a region--------
    char buf[1024], ch;
    char **splits;
    double delta=0.0;
    uint32_t header[8], rs, re;
    uint64_t nols=0;
    idx0 = nTiles+10;
    while(fgets(buf, 1024, fq)!=NULL){	
	//printf("%u %s",(uint32_t)nols, buf); 
        ichr = -1;
        splits = str_split(buf,'\t', &nCols); 
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
            while(k<93 && strcmp(splits[0], folder[k])!=0)
                k++;
            if(k<93)
                ichr = k;

        }          
        if(ichr>=0){
            q1  = (uint32_t)atoi(splits[1]);
            q2  = (uint32_t)atoi(splits[2]);
            delta += q2 - q1;
            nRegions++;
            n1 = q1/nbp;
            n2 = q2/nbp-n1;   
            idx = n1 + gstart[ichr];            
            //if(idx > 200660)
            //    printf("str %s, ichr %i, idx %i, cnt %u\n", splits[0], ichr, idx, counts[idx]);
            //find overlaps with this region
            if(n2==0){
                tc = counts[idx];
                if(tc>0){
                    if(idx!=idx0){
                        free(gdata);
                        fseek(fi, mloc[idx], SEEK_SET);
                        fread(&header, 4, 8, fi);
                        gdata = malloc(tc*16);
                        fread(gdata, 16, tc, fi);
                        idx0 = idx;
                    }
                    //optimize the search: b-search 
                    if(q2 < gdata[0].r_start){
                        //no overlap:do nothing tL>nc;tR<0
                    } 
                    else{
                        //BSearch takes very little time!!! 2M bSearch(50M)~0.5s!!!  
                        rs=0;   //sublist index range in aiList
                        for(k=0; k<header[0]; k++){
                            re = rs+header[k+1];
                            t = bSearch(gdata, rs, re, q2); //inline not better 
                            while(t >= rs && gdata[t].r_max > q1){
                                if(gdata[t].r_end>q1){
                                    hits[gdata[t].i_idx]++;  
                                    //printf("%s %s %s %u %u %u\n", splits[0], splits[1], splits[2], gdata[t].r_start, gdata[t].r_end, gdata[t].i_idx);                             
                                    nols++;              
                                } 
                                t--;
                            }                       
                            rs = re;    
                        }
                    }
                }          
            }
            else{ 
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(n1+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    tc = counts[m];
                    if(tc>0){
                        if(m!=idx0){
                            free(gdata);
                            fseek(fi, mloc[m], SEEK_SET);
                            fread(&header, 4, 8, fi);
                            gdata = malloc(tc*16);
                            fread(gdata, 16, tc, fi);
                            idx0 = m;
                        }
                        rs=0;   //q2>all start
                        for(k=0; k<header[0]; k++){
                            re = rs+header[k+1];
                            t = bSearch(gdata, rs, re, q2); //upper limit
                            while(t >= rs && gdata[t].r_max > q1){
                                if(gdata[t].r_end < bd && gdata[t].r_end>q1){
                                    //printf("%s %s %s %u %u %u\n", splits[0], splits[1], splits[2], gdata[t].r_start, gdata[t].r_end, gdata[t].i_idx);                                    
                                    hits[gdata[t].i_idx]++;                               
                                    nols++;              
                                } 
                                t--;
                            }                       
                            rs = re;    
                        }                        
                    }                   
                    bd += nbp;
                }	
                //--last tile: normal------------------------------
                m=idx+n2;
                tc = counts[m];
                if(tc>0){
                    if(m!=idx0){
                        free(gdata);
                        fseek(fi, mloc[m], SEEK_SET);
                        gdata = malloc(tc*16);
                        fread(&header, 4, 8, fi);
                        fread(gdata, 16, tc, fi);
                        idx0 = m;
                    }   
                    rs=0;   
                    for(k=0; k<header[0]; k++){
                        re = rs+header[k+1];
                        t = bSearch(gdata, rs, re, q2); //Ie 
                        while(t >= rs && gdata[t].r_max > q1){
                            if(gdata[t].r_end>q1){
                                //printf("%s %s %s %u %u %u\n", splits[0], splits[1], splits[2], gdata[t].r_start, gdata[t].r_end, gdata[t].i_idx);                                
                                hits[gdata[t].i_idx]++;                               
                                nols++;              
                            } 
                            t--;
                        }                       
                        rs = re;    
                    }             
                }//if tc>0
            }   //else n2>0            
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
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;       
    int ichr, tc, tL, tR, tM, tS, tlen, rtn;   
    uint32_t i, j, k, t1, t2, q1, q2, m;
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
        ichr = -1;
        splits = str_split(buf,'\t', &nCols); 
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
            while(k<93 && strcmp(splits[0], folder[k])!=0)
                k++;
            if(k<93)
                ichr = k;
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
                tc = counts[idx];
                if(tc>0){
                    if(idx!=idx0){
                        free(gdata);
                        fseek(fi, mloc[idx], SEEK_SET);
                        gdata = malloc(tc*16);
                        fread(gdata, 16, tc, fi);
                        idx0 = idx;
                    }
                    if(q1>gdata[tc-1].r_end){
                        //no overlap:do nothing tL>nc;tR<0
                    } 
                    else if(tc<64){
                        tS=0;
                        while(gdata[tS].r_end<q1)
                            tS++;                                          
                        for(j=tS;j<tc;j++){
                            if(gdata[j].g_val>v && q2>=gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;
                                nols++;
                            }
                        }                     
                    }
                    else{//half dual-binary search 
                        //search tS: the 1st t_end on the left of q_start
                        tL=0;   tR=tc-1;  
                        tS = -1;    //no exclusion; tL<nc-1
                        while(tL<tR-1){
                            tM = (tL+tR)/2; 
                            if(gdata[tM].r_end<q1)
                                tL = tM;
                            else
                                tR = tM - 1;
                        }
                        if(gdata[tR].r_end<q1)
                            tS = tR;
                        else if(gdata[tL].r_end<q1)
                            tS = tL;
                        tS++; 
                        for(j=tS;j<tc;j++){
                            if(gdata[j].g_val>v && q2>=gdata[j].r_start){                          		          		    
                                hits[gdata[j].i_idx]++;
                                nols++;
                            }
                        }                                       
                    } //else
                }          
            }
            else{ 
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(n1+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    tc = counts[m];
                    if(tc>0){
                        if(m!=idx0){
                            free(gdata);
                            fseek(fi, mloc[m], SEEK_SET);
                            gdata = malloc(tc*16);
                            fread(gdata, 16, tc, fi);
                            idx0 = m;
                        }
                        tS=0;
                        while(gdata[tS].r_end<q1)
                            tS++;            
                        for(j=tS;j<tc;j++){
                            t2 = gdata[j].r_end;
                            if(gdata[j].g_val>v && t2<bd){
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
                tc = counts[m];
                if(tc>0){
                    if(m!=idx0){
                        free(gdata);
                        fseek(fi, mloc[m], SEEK_SET);
                        gdata = malloc(tc*16);
                        fread(gdata, 16, tc, fi);
                        idx0 = m;
                    }
                    //update the in-tile db start j0: not really faster
                    if(q1>gdata[tc-1].r_end){
                        //no overlap:do nothing tL>nc;tR<0
                    } 
                    else if(tc<64){
                        tS=0;
                        while(gdata[tS].r_end<q1)
                            tS++;                                          
                        for(j=tS;j<tc;j++){
                            if(gdata[j].g_val>v && q2>=gdata[j].r_start){    		          		    
                                hits[gdata[j].i_idx]++;
                                nols++;
                            }
                        }                     
                    }
                    else{//half dual-binary search 
                        //search tS: the 1st t_end on the left of q_start
                        tL=0;   tR=tc-1;  
                        tS = -1;    //no exclusion; tL<nc-1
                        while(tL<tR-1){
                            tM = (tL+tR)/2; 
                            if(gdata[tM].r_end<q1)
                                tL = tM;
                            else
                                tR = tM - 1;
                        }
                        if(gdata[tR].r_end<q1)
                            tS = tR;
                        else if(gdata[tL].r_end<q1)
                            tS = tL;
                        tS++; 
                        for(j=tS;j<tc;j++){
                            if(gdata[j].g_val>v && q2>=gdata[j].r_start){                          		          		    
                                hits[gdata[j].i_idx]++;
                                nols++;
                            }
                        }                                       
                    } //else                   
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

void search(char* qfName, char* igdName, uint32_t v, char *out)
{   //name standard: igd_file(dbname.igd), index_file(dbname_index.tsv)
    uint32_t i, nq=1, nFiles, nCols=2, genome_size=3095677412;
    double mq = 1.0;
    char tmp[128];
    strcpy(tmp, igdName);
    tmp[strrchr(tmp, '.')-tmp] = '\0';
    char *idFile = tmp;//str_split(tmp, '.', &nCols)[0];
    strcat(idFile, "_index.tsv");     
 
    struct igd_info *fi = get_igdinfo(idFile, &nFiles); 
    //printf("nfiles: %u\n", nFiles);  
    //determine igd_data type: data structure length
    int type = 0;
    FILE* fp = fopen(igdName, "rb");
    if(!fp)
        return;  
    else{
        fread(&i, 4, 1, fp);

        if(i==1122331111)
            type = 1;
        else if(i==1122332222)
            type = 2;
        else if(i==1122330000)
            type = 0;
        else
            type = -1;//old style without type
        fclose(fp);
    }        
    //printf("igdData type %i, %u\n", type, i);      
    clock_t start, end;
    start = clock();
    uint32_t *hits = calloc(nFiles, sizeof(uint32_t));
    uint64_t nOL;
    if(v>0){
        nOL = get_overlaps_v(qfName, igdName, v, &nq, &mq, hits);       
    }
    else{  
        if(type==0)
            nOL = get_overlaps_n0(qfName, igdName, &nq, &mq, hits);
        else if(type==1)
            nOL = get_overlaps_n1(qfName, igdName, &nq, &mq, hits);   
        else if(type==2)
            nOL = get_overlaps_n2(qfName, igdName, &nq, &mq, hits); 
        else if(type==-1) 
            nOL = get_overlaps_n(qfName, igdName, &nq, &mq, hits);                        
    } 
    end = clock();  
 
    if(strlen(out)>1){
        FILE *fp = fopen(out, "w");
        if(fp==NULL)
            printf("Can't open file %s\n", idFile);
        else{
            fprintf(fp, "Number of overlaps: %u, number of query regions: %u, mean query size: %f \n", (uint32_t)nOL, nq, mq);     
            fprintf(fp, "index\t File_name\t number of regions\t mean-region-size \t number of hits\n");
            for(i=0;i<nFiles;i++)
                fprintf(fp, "%i %s %u %u %u\n", i, fi[i].fileName, fi[i].nd, (uint32_t)fi[i].md, hits[i]);  
            fclose(fp);
        }     
    }
    else{
        //printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);
        //printf("Number of overlaps: %u, number of query regions: %u, mean query size: %f \n", (uint32_t)nOL, nq, mq);  
        printf("index\t File_name\t number of regions\t mean-region-size \t number of hits\n");        
        for(i=0;i<nFiles;i++)
            printf("%i %s %u %u %u\n", i, fi[i].fileName, fi[i].nd, (uint32_t)fi[i].md, hits[i]);         
    }   
    
    //printf("out=%s\n", out);
    //printf("v=%u\n", v); 
    //---------------------------------------------------------------------------------
    free(fi->fileName);
    free(fi);
    free(hits);
}

//-------------------------------------------------------------------------------------
int igd_search(int argc, char **argv)
{   //igd[0] search[1] query100.bed[2] home/john/iGD/rme_igd/roadmap.igd[3]
    uint32_t v = 0;
    char out[64]="";   
    //convert block index to chr index for convenience
    g2ichr = malloc(nTiles*sizeof(uint32_t));
    uint32_t i, j;
    for(i=0; i<93; i++){  
    	for(j=gstart[i]; j<gstart[i+1]; j++)      
    	    g2ichr[j] = i;
    }
    char *qfName = argv[2];
    char *igdName = argv[3];
    
    //check if fils exist
    char *ftype = igdName + strlen(igdName) - 4;
    char *qtype = qfName + strlen(qfName) - 4;    
    if(strcmp(".igd", ftype)!=0){
        printf("%s is not an igd database", igdName);
        return EX_OK;
    }
    if(strcmp(".bed", qtype)!=0){
        printf("%s is not a bed file", qfName);
        return EX_OK;
    }
            
    FILE* fi = fopen(igdName, "rb");
    if(!fi){
        printf("%s does not exist", igdName);
        return EX_OK;
    }
    fclose(fi); 
    fi = fopen(qfName, "rb");
    if(!fi){
        printf("%s does not exist", qfName);
        return EX_OK;
    }
    fclose(fi); 
       
    if(argc>=6){
        if(strcmp(argv[4], "-v")==0)
            v = atoi(argv[5]);
        else if(strcmp(argv[4], "-o")==0)
            strcpy(out, argv[5]);
    }
    if(argc>=8){   
        if(strcmp(argv[6], "-v")==0)
            v = atoi(argv[7]);
        else if(strcmp(argv[6], "-o")==0)
            strcpy(out, argv[7]);
    }
    search(qfName, igdName, v, out);      

    free(g2ichr);
    return EX_OK;
}

