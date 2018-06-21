//=====================================================================================
//Read igd region data and query data, and then find all overlaps 
//by Jianglin Feng  05/12/2018
//igd_tiles2w: put all tiles data into a single file 
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
#include <errno.h>

//-------------------------------------------------------------------------------------
struct igd_data
{
    uint32_t i_idx;        			//genomic object--data set index
    uint32_t r_start;      			//region start
    uint32_t r_end;        			//region end
    uint32_t g_val;        			//signal level
};

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
uint64_t maxCount = 134217728;//maxCount*16=2GB,67108864,33554432, 67108864;//
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

//reload tile igd files, sort them and save them into a single file
void tiles2w(char *tPath, char* oPath, char *igdName)
{   //tPath: for igd tile data files
    //oPath: for the single igd file 
    int i, k;
    char idFile[256];
    sprintf(idFile, "%s%s%s", oPath, igdName, ".igd");    
    FILE *fp1 = fopen(idFile, "wb"); 
    if(fp1==NULL)
        printf("Can't open file %s", idFile);     
    uint32_t nrec;  
    FILE *fp0; 
    
    clock_t start, end;
    start = clock();        
    //struct igd_data *gdata;
    uint32_t *counts = calloc(nTiles, sizeof(uint32_t)); 
    fwrite(counts, sizeof(uint32_t), nTiles, fp1);

    char iname[256]; 
    struct igd_data *gdata; 
    printf("Starting ...\n");         
    for(i=0;i<nTiles;i++){
        k = g2ichr[i];
        sprintf(iname, "%s%s%s/%s_%u%s", tPath, "data0/", folder[k], igdName, i-gstart[k], ".igd");
        fp0 = fopen(iname, "rb");
        if(fp0!=NULL){   
            fseek(fp0, 0, SEEK_END);
            nrec = ftell(fp0)/sizeof(struct igd_data);
            fseek(fp0, 0, SEEK_SET);
            if(nrec>0){
                gdata = malloc(nrec*sizeof(struct igd_data));
                fread(gdata, sizeof(struct igd_data), nrec, fp0);
                fclose(fp0);
                qsort(gdata, nrec, sizeof(struct igd_data), compare_rend);
                //append the data to fp1
                fwrite(gdata, sizeof(struct igd_data), nrec, fp1);
                counts[i] = nrec;
                free(gdata);
                //gdata = NULL;
            }      
            //fclose(fp0);
        }
        if(i%1000==0)
            printf("%u\t%u\t%f \n", i, nrec, (double)(clock()-start)/CLOCKS_PER_SEC);
    }
    fseek(fp1, 0, SEEK_SET);
    fwrite(counts, sizeof(uint32_t), nTiles, fp1);
    fclose(fp1); 
    end = clock();    
    free(counts);   
    printf("igd_w finished: time: %f \n", (double)(end-start)/CLOCKS_PER_SEC);    
}

//-------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 4) 
        errx(1, "usage:\t%s <tile-file path><ouput path><dbName>\n", argv[0]);   
    //convert block index to chr index for convenience
    g2ichr = malloc(nTiles*sizeof(uint32_t));
    uint32_t i, j;
    for(i=0; i<24; i++){  
        	for(j=gstart[i]; j<gstart[i+1]; j++)      
        	    g2ichr[j] = i;
    }
    char *tpath = argv[1];
    char *opath = argv[2];
    char *dbname = argv[3];
    tiles2w(tpath, opath, dbname);    

    free(g2ichr);
}

/*
echo cfq > /sys/block/hda/queue/scheduler

*/