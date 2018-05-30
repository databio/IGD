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
    uint8_t r_chr;                    //chr:0-23
};

char** str_split( char* str, char delim, int nmax);
struct query_data* get_igdlist(char *fname, uint32_t *nblocks);	
struct igd_overlap* get_overlaps(struct query_data *query, uint32_t nblocks, uint32_t nmax, uint32_t *nOL);
struct igd_overlap* get_overlaps_w(struct query_data *query, uint32_t nblocks, uint32_t nmax, uint32_t *nOL);

char *fileBase = "bb12";         		//14 bits block
uint32_t nbp = 4096;
uint32_t nmax[] = {63748, 62292, 51012, 48924, 46648, 43932, 41008, 37444, 35400, 
	34432, 34820, 34304, 29164, 27324, 26004, 23304, 21456, 20608, 15256, 
	16576, 11908, 12932, 39504, 14024};
char *folder[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
 	"chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
uint32_t gstart[25] = {0, 63748, 126040, 177052, 225976, 272624, 316556, 357564, 
	395008, 430408, 464840, 499660, 533964,563128, 590452, 616456, 639760, 
	661216, 681824, 697080, 713656, 725564, 738496, 778000, 792024};
uint32_t nTiles = 792024;
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

char** str_split( char* str, char delim, int nmax)
{
    char** splits;
    char* c;    
    int ns;
    if (str == NULL || delim=='\0')
        return NULL;
    else {
        splits = malloc((nmax+1) * sizeof(*splits));
        splits[nmax] = NULL;
        c = str;
        ns = 1;
        splits[0] = str;
        do {
            if (*c == delim)
            {
                splits[ns++] = &c[1];
                *c = '\0';
            }
            c++;
        } while (*c != '\0' && ns < nmax+1);
    }
    return splits;
}

struct query_data* get_igdlist(char *fname, uint32_t *nblocks)
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
    int nmax = 5;  
    fseek(fp, 0, SEEK_SET);
    int nm, i=0, nextras = 0;//number of regions that cross the tile boundary
    while(fgets(buf, 1024, fp)!=NULL){	
        	splits = str_split(buf,'\t', nmax);
        	if(strlen(splits[0])<6){
        	    if(strcmp(splits[0], "chrX")==0)
            		df0[i] = 22;
        	    else if(strcmp(splits[0],"chrY")==0)
            	 	df0[i] = 23;
        	    else
        	      	df0[i] = (uint8_t)(atoi(&splits[0][3])-1);
        	    df1[i] = (uint32_t)atoi(splits[1]);
        	    df2[i] = (uint32_t)atoi(splits[2]); 
        	    n1[i] = df1[i]/nbp;
        	    n2[i] = df2[i]/nbp-n1[i];  
        	    nextras += n2[i];     
        	    i++;
        	}
    }
    fclose(fp);
    nm = i;
    
    *nblocks = nm+nextras;
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
    	sprintf(iname, "%s%s%s%i%s", "igdata/",folder[ichr], "/bb12_", k, ".igd");
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
struct igd_overlap* get_overlaps_w(struct query_data *query, uint32_t nblocks, uint32_t nmax, uint32_t* nOL)
{
    char *iname = "igdata/bb12.igd";   
    FILE* fp = fopen(iname, "rb");
    if(!fp)
        return NULL;    
        
    uint32_t i, j, i1, i2, bk, ichr, k, nrec, nols=0, q1, q2, m;
    struct igd_overlap *overlaps = malloc(nmax*sizeof(struct igd_overlap));
    struct igd_data *gdata;
    uint32_t len0 = nTiles*sizeof(uint32_t);
    uint32_t *counts = malloc(len0);//number of bytes
    uint32_t *mloc = malloc(len0 + sizeof(uint32_t));
    //fseek(fp, 0, SEEK_SET);
    fread(counts, sizeof(uint32_t), nTiles, fp);   
    for(i=0; i<nTiles; i++){
        mloc[i+1] = counts[i];
        counts[i]/=16;
    }    
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
        	    gdata = malloc(counts[bk]*sizeof(struct igd_data));
        	    fread(gdata, sizeof(struct igd_data), counts[bk], fp);
            //printf("%i %i %i \n", counts[bk], gdata[0].r_start, gdata[0].r_end);  
        	    for(i=i1; i<=i2; i++){
                q1 = query[i].r_start;
                q2 = query[i].r_end;
                for(j=0;j<counts[bk];j++){
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
        m++;
    }  
             
    fclose(fp);
    *nOL = nols;
    free(counts);
    free(mloc);
    return overlaps;
}

//-------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 2) 
        errx(1, "usage:\t%s <input path>\n", argv[0]);   
    int usingw = 0;
    if (argc==3)
        usingw = atoi(argv[2]);
    //convert block index to chr index for convenience
    g2ichr = malloc(nTiles*sizeof(uint32_t));
    uint32_t i, j;
    for(i=0; i<24; i++){  
        	for(j=gstart[i]; j<gstart[i+1]; j++)      
        	    g2ichr[j] = i;
    }
    //---------------------------------------------------------------------------------
    char *fname = argv[1];
    uint32_t nblocks=0;
    struct query_data *query = get_igdlist(fname, &nblocks);
    clock_t start, end;
    start = clock();
    uint32_t nOL=0, nmax=6000000;
    struct igd_overlap *overlaps;   
    if(usingw>0)
        overlaps = get_overlaps_w(query, nblocks, nmax, &nOL);
    else
        overlaps = get_overlaps(query, nblocks, nmax, &nOL);
   
    end = clock();
    printf("time: %f\n", ((double)(end-start))/CLOCKS_PER_SEC);
    printf("nblock: %i\n", nblocks);
    printf("number of overlaps:%i\n", nOL);
    for(i=0;i<10;i++)
     	printf("%i %i %i\n", overlaps[i].i_idx, overlaps[i].r_start, overlaps[i].r_end);
    //---------------------------------------------------------------------------------
    free(query);
    free(overlaps);
    free(g2ichr);
}
