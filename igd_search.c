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
struct igd_mix* get_overlaps(struct query_data *query, uint32_t nblocks, uint32_t nmax, uint32_t *nOL);
struct igd_mix* get_overlaps_w(struct query_data *query, uint32_t nblocks, char *igdName, uint32_t nmax, uint32_t *nOL);
struct igd_info* get_igdinfo(char *ifName, uint32_t *nFiles);
uint64_t get_overlaps_n(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);
uint64_t get_overlaps_n1(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);
void search(char* qfName, char* igdName);
void str_splits(char* str, int *nmax, char **splits);

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
uint32_t nbp = 16384;
uint32_t maxCount = 67108864, bgz_buf = 1024;//maxCount*16=2GB134217728,
uint32_t *g2ichr;

//-------------------------------------------------------------------------------------
//This section is taken from giggle package
/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 * kfunc.h
 *  Created on: May 1, 2015
 *      Author: nek3d
 */
long double _lbinom(long long n, long long k);
long double _hypergeo(long long n11, long long n1_, long long n_1, long long n);

typedef struct {
    long long n11, n1_, n_1, n;
    long double p;
} _hgacc_t;

// incremental version of hypergenometric distribution
long double _hypergeo_acc(long long n11, long long n1_, long long n_1, long long n, _hgacc_t *aux);
long double _kt_fisher_exact(long long n11, long long n12, long long n21, long long n22, long double *_left, long double *_right, long double *two);
 
// log\binom{n}{k}
long double _lbinom(long long n, long long k)
{
    if (k == 0 || n == k) return 0;
    return lgammal(n+1) - lgammal(k+1) - lgammal(n-k+1);
}

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n

// hypergeometric distribution
long double _hypergeo(long long n11, long long n1_, long long n_1, long long n)
{
    //***DEBUG***
    return expl(_lbinom(n1_, n11) + _lbinom(n-n1_, n_1-n11) - _lbinom(n, n_1));
}

// incremental version of hypergenometric distribution
long double _hypergeo_acc(long long n11, long long n1_, long long n_1, long long n, _hgacc_t *aux)
{
    if (n1_ || n_1 || n) {
        aux->n11 = n11; aux->n1_ = n1_; aux->n_1 = n_1; aux->n = n;
    } else { // then only n11 changed; the rest fixed
        if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
            if (n11 == aux->n11 + 1) { // incremental
                aux->p *= (long double)(aux->n1_ - aux->n11) / n11
                    * (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
                aux->n11 = n11;
                return aux->p;
            }
            if (n11 == aux->n11 - 1) { // incremental
                aux->p *= (long double)aux->n11 / (aux->n1_ - n11)
                    * (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
                aux->n11 = n11;
                return aux->p;
            }
        }
        aux->n11 = n11;
    }
    aux->p = _hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);

    return aux->p;
}

long double _kt_fisher_exact(long long n11,
                             long long n12,
                             long long n21,
                             long long n22,
                             long double *_left,
                             long double *_right,
                             long double *two)
{
    long long i, j, max, min;
    long double p, q, left, right;
    _hgacc_t aux;
    long long n1_, n_1, n;

    n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n

    max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
    min = n1_ + n_1 - n;    // not sure why n11-n22 is used instead of min(n_1,n1_)
    if (min < 0) min = 0; // min n11, for left tail
    *two = *_left = *_right = 1.;

    if (min == max) return 1.; // no need to do test


    q = _hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
    if (q < 1e-200) q = 1e-200;

    // left tail
    p = _hypergeo_acc(min, 0, 0, 0, &aux);
    for (left = 0., i = min + 1; p < 0.99999999 * q && i<=max; ++i) // loop until underflow
        left += p, p = _hypergeo_acc(i, 0, 0, 0, &aux);
    --i;
    if (p < 1.00000001 * q) left += p;
    else --i;
    // right tail
    p = _hypergeo_acc(max, 0, 0, 0, &aux);
    for (right = 0., j = max - 1; p < 0.99999999 * q && j>=0; --j) // loop until underflow
        right += p, p = _hypergeo_acc(j, 0, 0, 0, &aux);
    ++j;
    if (p < 1.00000001 * q) right += p;
    else ++j;
    // two-tail
    *two = left + right;
    if (*two > 1.) *two = 1.;
    // adjust left and right
    if (labs((long) (i - n11)) < labs((long) (j - n11)) && q != 0.0) right = 1. - left + q;
    else left = 1.0 - right + q;
    *_left = left; *_right = right;
    return q;
}

//The following 2 are from giggle
//{{{double log2fc(double ratio)
double log2fc(double ratio)
{
    if (fabs(ratio) < 0.0001)
        return 0.0;

    if (ratio < 1) {
        ratio = 1.0/ratio;
        return -1.0 * log2(ratio);
    }

    return log2(ratio);
}
//}}}

//{{{double neglog10p(double sig)
long double neglog10p(long double sig)
{
    if (fabsl(sig) < -DBL_MAX)
        return 10.0;
    return -1.0 * log10l(sig);
}
//}}}

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

void str_splits(char* str, int *nmax, char **splits)
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

    int nCols = 5;    
    char **splits = malloc((nCols+1)*sizeof(char *)); 
  
    fseek(fp, 0, SEEK_SET);
    int k, nm, i=0, nextras=0;
    double delta=0.0;//number of regions that cross the tile boundary
    while(fgets(buf, 1024, fp)!=NULL){	
        //splits = str_split_t(buf, nItems);
        str_splits(buf, &nCols, splits);       
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
    qsort(query, *nblocks, sizeof(struct query_data), compare_qidx);
    free(splits);
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
    //-----------------------------------------------------
    int ncols = 4;    
    char **splits = malloc((ncols+1)*sizeof(char *));     
    
    fseek(fp, 0, SEEK_SET);
    i=0;
    fgets(buf, 1024, fp);   //header
    while(fgets(buf, 1024, fp)!=NULL){	
        str_splits(buf, &ncols, splits);        
        fi[i].fileName = (char *)calloc(strlen(splits[1]) + 1, sizeof(char));
        strcpy(fi[i].fileName, splits[1]);
        fi[i].nd = (uint32_t)atoi(splits[3]);
        fi[i].md = (double)atoi(splits[2]);   
        i++;
    }        
    //printf("%s %i %f \n", fNames[6], nd[6], md[6]);  
    //printf("Total file: %i\n", i);
    free(splits);
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
uint64_t get_overlaps_n1(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits)
{   //no need to store every overlaps, only get the number of hits
    //faster if in-tile igdata is sorted by region end 
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;    
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
    int ichr;           
    uint32_t i, j, t1, t2, q1, q2, m;
    uint32_t n1, n2, idx, nRegions=0, nCols = 5, bd;  
    struct igd_data *gdata;
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);//number of struct
    uint64_t *mloc = calloc((nTiles+1),sizeof(uint64_t));
    //fseek(fp, 0, SEEK_SET);
    fread(counts, sizeof(uint32_t), nTiles, fi);   
    for(i=0; i<nTiles; i++)
        mloc[i+1] = counts[i]*16;
    mloc[0]=len0;
    for(i=1; i<nTiles; i++)
        mloc[i] += mloc[i-1];  
   
   //----read qf for a region--------
    char buf[1024];
    char **splits;
    double delta=0.0;
    uint64_t nols=0;
    while(fgets(buf, 1024, fq)!=NULL){	
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
                //printf("%u\t %u\t %u\t", counts[idx], q1, q2); 
                if(counts[idx]>0){          
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
                                //if(gdata[j].i_idx==12)
                                //    printf("%i %i %i %i\n", q1, q2, t1, t2);
                                nols++;
                            }
                        }
                    } 
                    free(gdata);  
                }        
            }
            else{ 
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(n1+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    if(counts[m]>0){                
                        //printf("--%u\t %u\t %u\t", counts[m], q1, q2); //cause core dumped! 
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
                                    //if(gdata[j].i_idx==12)
                                    //    printf("--%i %i %i %i\n", q1, q2, t1, t2);
                                    nols++;
                                }
                            }
                        }
                        free(gdata);
                    }                        
                    bd += nbp;
                }	
                //--last tile: normal------------------------------
                m=idx+n2;
                if(counts[m]>0){
                    //printf("--%u\t %u\t %u\t", counts[m], q1, q2); //cause core dumped!               
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
                                //if(gdata[j].i_idx==12)
                                //    printf("----%i %i %i %i\n", q1, q2, t1, t2); 
                                nols++;
                            }
                        }
                    }
                    free(gdata);
                }
            }   
        }   //if
    }   //while  
    *mq = delta/nRegions;
    *nq = nRegions;
    fclose(fq);   
    fclose(fi);
    free(counts);
    free(mloc);
    return nols;
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
    }   //while  
    if(splits!=NULL)
        free(splits);           
    free(gdata);
    *mean_size = delta/nRegions;
    *nregions = nRegions;
    fclose(fq);   
    fclose(fi);
    free(counts);
    free(mloc);
    return nols;
}

void search(char* qfName, char* igdName)
{   //name standard: igd_file(dbname.igd), index_file(dbname_index.tsv)
    uint32_t i, nq=1, nFiles, nCols=2, genome_size=3095677412;
    double mq = 1.0;
    char tmp[128];
    strcpy(tmp, igdName);
    tmp[strrchr(tmp, '.')-tmp] = '\0';
    char *idFile = tmp;//str_split(tmp, '.', &nCols)[0];
    strcat(idFile, "_index.tsv");     
 
    struct igd_info *fi = get_igdinfo(idFile, &nFiles); 
    printf("nfiles: %u", nFiles);  
    
    clock_t start, end;
    start = clock();
    uint32_t *hits = calloc(nFiles, sizeof(uint32_t));

    uint64_t nOL = get_overlaps_n(qfName, igdName, &nq, &mq, hits);   
    end = clock();   
    
    printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);
    printf("%u %u %u %f \n", nFiles, (int)nOL, nq, mq);
    printf("index\t File_name\t number of regions\t mean-region-size \t number of hits\n");
    for(i=0;i<10;i++)
        printf("%i %s %u %u %u\n", i, fi[i].fileName, fi[i].nd, (uint32_t)fi[i].md, hits[i]);
    //---------------------------------------------------------------------------------
    free(fi->fileName);
    free(fi);
    free(hits);
}

//-------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 2) 
        errx(1, "usage:\t%s <input path>\n", argv[0]);   
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
    char *qfName = argv[1];
    char *igdName = argv[2];
    
    if(usingw)
        search(qfName, igdName);      

    free(g2ichr);
}

//-------------------------------------------------------------------------------------
//notes section
/*
.Initialized on may 10, 2018
.ulimit -s: soft limit of system-->core dumped?
    $ulimit -s 16384 # 16MB stack
    $ulimit -s unlimited # "unlimited" stack
.giggle -s option for significance analysis with query file:
    use without -s for pure search comparing 
.https://codeforwin.org/2015/05/list-of-all-format-specifiers-in-c-programming.html
.Inlcude <float.h> for DBL_MAX definition
.Using -lm  option to compile linked to math lib
    gcc -o get_overlaps_b14 get_overlaps_b14.c -lm
    use .lz to link zlib
.calloc() is similar to malloc() but the former does initialize memory to 0
.may not free gdata but slower
.chrM: core dumped
.glob--path for files in a path with path/*
.!/bin/sh may not work =>/bin/bash --version
.if python parse.py not work =>try python2 parse.py 
.Error info (giggle index): Could not open...=> ulimit -Sn 16384

r or rb
    Open file for reading.
w or wb
    Truncate to zero length or create file for writing.
a or ab
    Append; open or create file for writing at end-of-file.
r+ or rb+ or r+b
    Open file for update (reading and writing).
w+ or wb+ or w+b
    Truncate to zero length or create file for update.
a+ or ab+ or a+b
    Append; open or create file for update, writing at end-of-file. 
john@JCloud:/media/john/CE30F6EE30F6DC81/ucsc_data$ time giggle  index   
   -i "parsed_tracks_sorted/*gz" -o parsed_tracks_sorted_b -s -f
Indexed 7456805968 intervals.
real	909m59.405s
user	367m22.618s
sys	29m51.091s
  rme_data/split_sort$ time ls *.bed.gz | gargs "tabix {}"

real	4m34.347s
user	2m8.513s
sys	0m20.360s
rme_data$ time giggle index -i "split_sort/*gz" -o split_sort_b -f -s
Indexed 55605005 intervals.

real	1m59.529s
user	1m31.234s
sys	0m7.898s
  ucsc_data$ time giggle search -i parsed_tracks_sorted_b -q ucsc_r100.bed.gz >/dev/null 
real	1m44.226s
user	0m0.196s
sys	0m1.983s
ucsc_data$ time giggle search -i parsed_tracks_sorted_b -q ucsc_r10000.bed.gz >/dev/null 
real	33m59.485s
user	0m12.011s
sys	0m34.640s 
 ucsc_data$ time giggle search -i parsed_tracks_sorted_b -q ucsc_r100000.bed.gz >/dev/null 
real	54m26.230s
user	1m22.248s
sys	0m52.168s
ucsc_data$ time giggle search -i parsed_tracks_sorted_b -q ucsc_r1000.bed.gz >/dev/null 
real	8m42.814s
user	0m1.373s
sys	0m7.795s
ucsc_data$ time giggle search -i parsed_tracks_sorted_b -q ucsc_r1000000.bed.gz >/dev/null 
Killed
real	114m56.288s
user	2m5.064s
sys	1m24.374s

/rme_data$     for Q_SIZE in $Q_SIZES; do
>         speed_test.sh \
>             rme_r$Q_SIZE.bed.gz \
>             split_sort \
>             rme.human.hg19.genome
>     done \
>     > $RESULTS
0.09 real	0.01 user	0.03 sys
29.25 real	28.29 user	0.69 sys
33.83 real	20.90 user	8.49 sys

0.13 real	0.01 user	0.04 sys
33.09 real	32.17 user	0.68 sys
55.95 real	41.66 user	9.74 sys

0.25 real	0.05 user	0.08 sys
35.32 real	34.38 user	0.68 sys
275.11 real	258.66 user	11.72 sys

1.30 real	0.64 user	0.31 sys
40.56 real	39.53 user	0.83 sys
5845.14 real	2457.11 user	29.63 sys

189.47 real	14.39 user	3.30 sys
285.91 real	101.11 user	2.02 sys
25046.15 real	24023.17 user	173.08 sys

112.97 real	65.53 user	3.74 sys
477.85 real	428.48 user	2.13 sys



   fprintf(fp, "%s %s %s %d", "We", "are", "in", 2012);

  


printf outputs to the standard output stream (stdout)

fprintf goes to a file handle (FILE*)

sprintf goes to a buffer you allocated. (char*)
   
*/
