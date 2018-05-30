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

struct igd_overlap
{
    uint32_t i_idx;        			//genomic object
    uint32_t r_start;      			//region start
    uint32_t r_end;        			//region end
    uint32_t g_val;        			//signal level
    uint8_t r_chr;                  //chr:0-23
};

struct igd_info
{
    char* fileName;
    uint32_t nd;
    double md;
};

char** str_split( char* str, char delim, int nmax);
char** str_split_t( char* str, int nItems);
struct query_data* get_igdlist(char *fname, uint32_t *nblocks, uint32_t *nRegions, double *mRegion);	
struct igd_overlap* get_overlaps(struct query_data *query, uint32_t nblocks, uint32_t nmax, uint32_t *nOL);
struct igd_overlap* get_overlaps_w(struct query_data *query, uint32_t nblocks, char *igdName, uint32_t nmax, uint32_t *nOL);
struct igd_info* get_igdinfo(char *fname, uint32_t *nFiles);
uint64_t get_overlaps_n(char *qfName, char *igdName, uint32_t *nregions, double *mean_size, uint32_t *hits);
char *fileBase = "b14";         	//14 bits block
uint32_t nmax[] = {15940, 15580, 12760, 12240, 11670, 10990, 10260, 9370, 8860, 8610, 8710, 
        8580, 7300, 6840, 6510, 5830, 5370, 5160, 3820, 4150, 2980, 3240, 9880, 3510};
char *folder[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
     	"chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
uint32_t gstart[25] = {0, 15940, 31520, 44280, 56520, 68190, 79180, 89440, 98810, 107670, 
        116280, 124990, 133570, 140870, 147710, 154220, 160050, 165420, 170580, 174400, 
        178550, 181530, 184770, 194650, 198160}; 
uint32_t nTiles = 198160;
uint32_t nbp = 16384;
uint32_t *g2ichr;
//-------------------------------------------------------------------------------------
int compare_iidx(const void *a, const void *b)
{
    struct igd_data *pa = (struct igd_data *) a;
    struct igd_data *pb = (struct igd_data *) b;
    return pa->i_idx - pb->i_idx;
}

int compare_qidx(const void *a, const void *b)
{
    struct query_data *pa = (struct query_data *) a;
    struct query_data *pb = (struct query_data *) b;
    return pa->q_idx - pb->q_idx;
}

int compare_oidx(const void *a, const void *b)
{
    struct igd_overlap *pa = (struct igd_overlap *) a;
    struct igd_overlap *pb = (struct igd_overlap *) b;
    return pa->i_idx - pb->i_idx;
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


char** str_split( char* str, char delim, int nmax)
{   //slightly faster than _t
    char** splits;
    char* ch;    
    int ns;
    if (str == NULL || delim=='\0')
        return NULL;
    else {
        splits = malloc((nmax+1) * sizeof(*splits));
        splits[nmax] = NULL;
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
        } while (*ch != '\0' && ns < nmax+1);
    }
    return splits;
}

struct query_data* get_igdlist(char *fname, uint32_t *nblocks, uint32_t *nRegions, double *mRegion)
{   //---------------------------------------------------------------------------------
    FILE *fp = fopen(fname, "r");
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
    int nItems = 5;  
    fseek(fp, 0, SEEK_SET);
    int nm, i=0, nextras=0;
    double delta=0.0;//number of regions that cross the tile boundary
    while(fgets(buf, 1024, fp)!=NULL){	
        //splits = str_split_t(buf, nItems);
        splits = str_split(buf,'\t', nItems);       
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
            for(int k=1; k<=n2[i]; k++){
            	query[j].q_idx   = k + n1[i] + gstart[df0[i]];
            	query[j].r_start = df1[i];
            	query[j].r_end   = df2[i];	
            	j++;		
            }
        }
    }
    qsort(query, *nblocks, sizeof(struct query_data), compare_qidx);
    free(df0);
    free(df1);
    free(df2);
    free(n1);
    free(n2);
    return query;   
}

struct igd_info* get_igdinfo(char *fname, uint32_t *nFiles)
{
    FILE *fp = fopen(fname, "r");
    if(fp==NULL)
        return NULL;
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
    double delta=0.0;//number of regions that cross the tile boundary
    fgets(buf, 1024, fp);   //header
    while(fgets(buf, 1024, fp)!=NULL){	
        splits = str_split(buf,'\t', ncols);        
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


struct igd_overlap* get_overlaps(struct query_data *query, uint32_t nblocks, uint32_t nmax, uint32_t* nOL)
{
    uint32_t i, j, i1, i2, bk, ichr, k, nrec, nols=0, q1, q2, m;
    char iname[128];
    struct igd_overlap *overlaps = malloc(nmax*sizeof(struct igd_overlap));
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
                			overlaps[nols].i_idx   = gdata[j].i_idx;
                			overlaps[nols].r_start = gdata[j].r_start;
                			overlaps[nols].r_end   = gdata[j].r_end;
                			//overlaps[nols].g_val   = gdata[j].g_val;
                			overlaps[nols].r_chr   = ichr;
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
struct igd_overlap* get_overlaps_w(struct query_data *query, uint32_t nblocks, char *igdName, uint32_t nmax, uint32_t* nOL)
{   //add alignment: padding 
    //assume in-tile igdata is sorted by region end 
    FILE* fp = fopen(igdName, "rb");
    if(!fp)
        return NULL;    
        
    uint32_t i, j, t1, t2, i1, i2, bk, ichr, nols=0, q1, q2, m;
    struct igd_overlap *overlaps = malloc(nmax*sizeof(struct igd_overlap));
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
                			overlaps[nols].i_idx   = gdata[j].i_idx;
                			overlaps[nols].r_start = t1;
                			overlaps[nols].r_end   = t2;
                			//overlaps[nols].g_val   = gdata[j].g_val;
                			overlaps[nols].r_chr   = ichr;
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
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;    
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
           
    uint32_t i, j, t1, t2, ichr, q1, q2, m;
    uint32_t n1, n2, idx, nRegions=1, nItems = 5, bd;  
    struct igd_data *gdata;
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
    while(fgets(buf, 1024, fq)!=NULL){	
        splits = str_split(buf,'\t', nItems);       
        if(strlen(splits[0])<6){
            if(strcmp(splits[0], "chrX")==0)
                ichr = 22;
            else if(strcmp(splits[0], "chrY")==0)
                ichr = 23;
            else
                ichr = (uint32_t)(atoi(&splits[0][3])-1);
            q1  = (uint32_t)atoi(splits[1]);
            q2  = (uint32_t)atoi(splits[2]);
            delta += q2 - q1;
            nRegions++;
            n1 = q1/nbp;
            n2 = q2/nbp-n1;   
            idx = n1 + gstart[ichr];
            //find overlaps with this region  
            if(n2==0){
                fseek(fi, mloc[idx], SEEK_SET);
                gdata = malloc(counts[idx]*16);
                fread(gdata, 16, counts[idx], fi);
                //update the in-tile db start j0: not really faster           
                for(j=0;j<counts[idx];j++){
                    t2 = gdata[j].r_end;
                    if(q1<=t2){
                        t1 = gdata[j].r_start;
                        if(q2>=t1){    		          		    
                            hits[gdata[j].i_idx]++;
                            nols++;
                        }
                    }
                } 
                free(gdata);          
            }
            else{ 
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(idx+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    fseek(fi, mloc[m], SEEK_SET);
                    gdata = malloc(counts[m]*16);
                    fread(gdata, 16, counts[m], fi);
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
                    bd += nbp;
                    free(gdata);
                }	
                //--last tile: normal------------------------------
                m=idx+n2;
                fseek(fi, mloc[m], SEEK_SET);
                gdata = malloc(counts[m]*16);
                fread(gdata, 16, counts[m], fi);
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
                free(gdata);
            }   
        }   //if
    }   //while  
    *mean_size = delta/nRegions;
    *nregions = nRegions;
    fclose(fq);   
    fclose(fi);
    free(counts);
    free(mloc);
    return nols;
}

//-------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 2) 
        errx(1, "usage:\t%s <input path>\n", argv[0]);   
    int usingw = 0;
    if (argc==3)
        usingw = 1;
        
    clock_t start0, end0;
    start0 = clock();       
    //convert block index to chr index for convenience
    g2ichr = malloc(nTiles*sizeof(uint32_t));
    uint32_t i, j;
    for(i=0; i<24; i++){  
        	for(j=gstart[i]; j<gstart[i+1]; j++)      
        	    g2ichr[j] = i;
    }
    //---------------------------------------------------------------------------------
    char *fname = argv[1];
    char *igdName = argv[2];
    uint32_t nblocks=0, nq=1, nFiles;
    double mq = 1.0;
    struct query_data *query = get_igdlist(fname, &nblocks, &nq, &mq);
    //-----open idfile-----
    char tmp[128];
    strcpy(tmp, igdName);
    char *idFile = str_split(tmp, '_', 2)[0];
    strcat(idFile, "_index.tsv");  
    struct igd_info *fi = get_igdinfo(idFile, &nFiles);
    printf("%s %i %f\n", fi[10].fileName, fi[10].nd, fi[10].md);     
    //---------------------
    end0 = clock(); 
    
    clock_t start, end;
    start = clock();
    uint32_t nOL=0, nmax=12000000;
    struct igd_overlap *overlaps;   
    if(usingw>0)
        overlaps = get_overlaps_w(query, nblocks, igdName, nmax, &nOL);
    else
        overlaps = get_overlaps(query, nblocks, nmax, &nOL);
    
    end = clock();
    printf("time: %f\n", ((double)(end0-start0))/CLOCKS_PER_SEC);    
    printf("time: %f\n", ((double)(end-start))/CLOCKS_PER_SEC);
    printf("nblock: %i %i %f \n", nblocks, nq, mq);
    printf("number of overlaps:%i\n", nOL);
    for(i=nOL-1000;i<nOL;i++)
     	printf("%i %i %i %i\n", overlaps[i].i_idx, overlaps[i].r_start, overlaps[i].r_end, overlaps[i].r_chr);
    //---------------------------------------------------------------------------------
    free(fi->fileName);
    free(fi);
    free(query);
    free(overlaps);
    free(g2ichr);
}
