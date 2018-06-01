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

struct igd_overlap
{
    uint32_t i_idx;        			//genomic object
    uint32_t r_start;      			//region start
    uint32_t r_end;        			//region end
    uint32_t g_val;        			//signal level
    uint8_t r_chr;                 //chr:0-23
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
uint64_t get_overlaps_n(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);
uint64_t get_overlaps_n1(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);
void test1(char* fname, char* igdName);
void testMain(char* fname, char* igdName);
void search(char* fname, char* igdName);
void create_roadmap();

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
//This section is used here for comparison with giggle etc
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

//create roadmap igd from gz
void create_igd(char *path)
{   
    //0. Get the files   
    glob_t gresult;
    int rtn = glob(path, 0, NULL, &gresult);     
    if(rtn!=0){
        printf("wrong dir path: %s", path);
        return;
    }
    char** file_ids = gresult.gl_pathv;
    uint32_t n_files = gresult.gl_pathc;
    
    //1. Read head info: assume a .bed or .vcf file 
    char* tmp = file_ids[0] + strlen(file_ids[0]) - 4;
    int fileType = 0, colC, colS, colE, colV;
    if(strcmp(".bed", tmp)==0){
        fileType = 1;
        colC = 0;
        colS = 1;
        colE = 2;
        colV = 4;
    }
    else if(strcmp(".vcf", tmp)==0){
        filetype = 2; 
        colC = 0;
        colS = 1;
        colE   
    }
    else{
        printf("File type not supported");
        return;   
    }
    regionData = pd.read_csv(file_path+file_ids[0], delimiter='\t', header=None) 
    nrows, ncols = regionData.shape
    print("nFiles, rows[0], cols[0], ", n_files, nrows, ncols, "\n") 
        
    //2. Read region data: read int64 default--int32 should be better
    nRegions = np.zeros(n_files, dtype='uint32')
    mRegion = np.zeros(n_files, dtype='uint32')
    count = np.zeros(nTiles, dtype=np.uint32)    
    data = np.empty(nTiles, dtype=object)        #bytearray        
    for i, id_ in tqdm.tqdm(enumerate(file_ids)):
        file = file_path + id_
        regionData = pd.read_csv(file, delimiter='\t', header=None)       
        df = regionData.sort_values(by=[0, 1])   #first by str, then by start
        n1 = df[1].values//nbp
        n2 = df[2].values//nbp-n1 
        itmp = len(n1)
        nRegions[i] = itmp
        tmp = np.sum(df[2]-df[1])/itmp    #64-bit
        mRegion[i]= int(tmp)
        rchr, ridx, rcnt = np.unique(df[0].values, return_index=True, return_counts=True)        
        //if a record crosses the block boundary, list it under both blocks (duplicates)
        //the start and end values are kept for fast processing (np): serialization and deserial..
        rc1 = df[1].values #.astype('uint32')
        rc2 = df[2].values #.astype('uint32')
        if ncols<5:
            rc3 = np.ones(len(df[1]), dtype='uint32')
        else:
            rc3 = df[4].values #.astype('uint32')
        #rec_bytes = np.array(rc1, dtype=int)
        for m in range(0, len(rchr)):
            ichr = -1
            if rchr[m] == 'chrX':
                ichr = 22
            elif rchr[m] == 'chrY':
                ichr = 23
            else:
                tmps = rchr[m][3:]
                if tmps.isdigit():
                    ichr = int(tmps)-1
            if ichr<24 and ichr>=0:
                for k in range(0, rcnt[m]):
                    idx0 = k+ridx[m]
                    idx = n1[idx0]+gstart[ichr]
                    #16 bytes for fast pack/unpack
                    rec = struct.pack('IIII', i, rc1[idx0], rc2[idx0], rc3[idx0])          
                    for j in range(0,n2[idx0]+1):
                        if idx+j<nTiles:
                            if data[idx+j]==None:
                                data[idx+j] = rec
                            else:
                                data[idx+j] += rec 

    //save all in a single file
    headInfo = {'File-id':file_ids, 'number_of_regions':nRegions, 'mean_region_size':mRegion}
    headInfo = pd.DataFrame(headInfo)
    headInfo.to_csv('igdata/roadmap_index.tsv', sep='\t')    
    //headInfo.to_csv('igdata/roadmap_index.tsv', columns=['File_id', 'number_of_regions', 'mean_region_size'], sep='\t')    
    //---------------------------------------------------------
    t0 = time.time()
    file = open('igdata/roadmap_'+fileBase+'.igd', 'wb')
    //Write header info: (1587200, count*14)
    for m in range(nTiles):
        if data[m]!=None:
            count[m]=len(data[m])/16 # number of struct
        else:
            count[m]=0
       
    file.write(count.tostring())
    tmpd = []
    for m in range(nTiles):
        if count[m]>0:
            //to list then sort and then repack
            tmpd = list(struct.iter_unpack('IIII', data[m]))
            tmpd.sort(key=itemgetter(2))
            tmp = bytearray()
            for i in range(count[m]):
                tmp += struct.pack('IIII', *tmpd[i])           
            file.write(tmp)       
    file.close()
    print('t_save=', time.time()-t0)   
    
    globfree(&gresult); 
}

struct query_data* get_igdlist(char *fname, uint32_t *nblocks, uint32_t *nRegions, double *mRegion)
{  
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
    uint32_t n1, n2, idx, nRegions=0, nItems = 5, bd;  
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
        splits = str_split(buf,'\t', nItems);   
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
            if(nRegions>100000)
                printf("%u\t %u\t %u\t \n", counts[idx], q1, q2);          
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
    FILE* fi = fopen(igdName, "rb");
    if(!fi)
        return 0;    
    FILE* fq = fopen(qfName, "r");
    if(!fq)
        return 0;    
     
    int ichr;      
    uint32_t i, j, t1, t2, q1, q2, m;
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
            }
            else{ 
                //deal with duplicates: find the unique list before or after 
                bd = nbp*(idx+1); 
                //in tiles (m=0, m=n2-1): t2<bd(m)-----------------        
                for(m=idx; m<idx+n2; m++){
                    if(counts[m]>0){
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
                        free(gdata);
                    }                   
                    bd += nbp;
                }	
                //--last tile: normal------------------------------
                m=idx+n2;
                if(counts[m]>0){
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

void test1(char* fname, char* igdName)
{
    uint32_t nblocks=0, nq=1, nFiles;
    double mq = 1.0;
    struct query_data *query = get_igdlist(fname, &nblocks, &nq, &mq);
    
    //-----open idfile-----
    char tmp[128];
    strcpy(tmp, igdName);
    char *idFile = str_split(tmp, '_', 2)[0];
    strcat(idFile, "_index.tsv");  
    struct igd_info *fi = get_igdinfo(idFile, &nFiles);
    printf("%s %u %f\n", fi[10].fileName, fi[10].nd, fi[10].md);     
    
    clock_t start, end;
    start = clock();
    uint32_t nOL=0, nmax=15000000;
    struct igd_overlap *overlaps;   
    overlaps = get_overlaps_w(query, nblocks, igdName, nmax, &nOL);
    
    end = clock();   
    printf("time: %f\n", ((double)(end-start))/CLOCKS_PER_SEC);
    printf("nblock: %u %u %f \n", nblocks, nq, mq);
    printf("number of overlaps:%i\n", nOL);
    //---------------------------------------------------------------------------------
    free(fi->fileName);
    free(fi);
    free(query);
    free(overlaps);
}

void testMain(char* fname, char* igdName)
{
    uint32_t nq=1, nFiles, genome_size=3095677412;
    double mq = 1.0;
    char tmp[128];
    strcpy(tmp, igdName);
    char *idFile = str_split(tmp, '_', 2)[0];
    strcat(idFile, "_index.tsv");  
    struct igd_info *fi = get_igdinfo(idFile, &nFiles);   
    
    clock_t start, end;
    start = clock();
    uint32_t *hits = calloc(sizeof(uint32_t), nFiles);
    uint64_t nOL = get_overlaps_n(fname, igdName, &nq, &mq, hits);   
    end = clock();   
    
    printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);
    printf("%u %u %u %f \n", nFiles, (int)nOL, nq, mq);
    long long n11, n12, n21, n22, n22_full, n3;
    double comp_mean, ratio;
    long double left, right, two, r;
    for(int i=0;i<nFiles;i++){
        //printf("%i %s %i %i %i\n", i, fi[i].fileName, fi[i].nd, (int)fi[i].md, hits[i]);
        n11 = (long long)(hits[i]);
        n12 = (long long)(MAX(0,nq-hits[i]));
        n21 = (long long)(MAX(0,fi[i].nd-hits[i]));
        comp_mean = fi[i].md + mq;
        n3 = n11 + n12 + n21;
        n22_full = (long long)MAX(n3, genome_size/comp_mean);
        n22 = MAX(1, n22_full - n3);
        left, right, two;
        //printf("%i %u %f %lli %lli %lli %lli\n", i, hits[i], comp_mean, n11, n12, n21, n22);
        r = _kt_fisher_exact(n11,n12,n21,n22,&left,&right,&two);
        ratio = (((double)(n11 + 1)/(double)MAX(1,n12))/((double)(n21 + 1)/(double)(n22 + 1)));
        printf("%s\t %u\t %u\t %.17g\t %.17Lg\t %.17Lg\t %.17Lg\t %.17Lg\t \n",
            fi[i].fileName, fi[i].nd, hits[i], ratio, two, left, right, 
            log2fc(ratio) * neglog10p(two));
    }
    //---------------------------------------------------------------------------------
    free(fi->fileName);
    free(fi);
    free(hits);
}

void search(char* fname, char* igdName)
{
    uint32_t nq=1, nFiles, genome_size=3095677412;
    double mq = 1.0;
    char tmp[128];
    strcpy(tmp, igdName);
    char *idFile = str_split(tmp, '_', 2)[0];
    strcat(idFile, "_index.tsv");  
    struct igd_info *fi = get_igdinfo(idFile, &nFiles);   
    
    clock_t start, end;
    start = clock();
    uint32_t *hits = calloc(nFiles, sizeof(uint32_t));
    uint64_t nOL = get_overlaps_n1(fname, igdName, &nq, &mq, hits);   
    end = clock();   
    
    printf("time: %f \n", ((double)(end-start))/CLOCKS_PER_SEC);
    printf("%u %u %u %f \n", nFiles, (int)nOL, nq, mq);
    printf("index\t File_name\t number of regions\t mean-region-size \t number of hits\n");
    for(int i=0;i<nFiles/20;i++)
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
    //---------------------------------------------------------------------------------
    
    //char* path = "/media/john/CE30F6EE30F6DC81/roadmap_sort/";
    
    char *fname = argv[1];
    char *igdName = argv[2];
    if(usingw)
        search(fname, igdName);      
        //testMain(fname, igdName);

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
.calloc() is similar to malloc() but the former does initialize memory to 0
.may not free gdata but slower
.chrM: core dumped
*/