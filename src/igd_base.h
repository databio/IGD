//=================================================================================
//Common structs, parameters, functions
//by Jianglin Feng  05/12/2018
//re-designed 7/1/2019
//database intervals sorted by _start: 8/12/2019
//---------------------------------------------------------------------------------
#ifndef __IGD_BASE_H__
#define __IGD_BASE_H__

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
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include "khash.h"
#include "kseq.h"

#define PROGRAM_NAME  "igd"
#define MAJOR_VERSION "0"
#define MINOR_VERSION "1"
#define REVISION_VERSION "1"
#define BUILD_VERSION "0"
#define VERSION MAJOR_VERSION "." MINOR_VERSION "." REVISION_VERSION
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define maxCount 268435456	//16* = 4GB memory
#define MAXC 10							//max number of components
//---------------------------------------------------------------------------------
typedef struct{							//default 
    int32_t idx;        				//genomic object--data set index
    int32_t start;      				//region start
    int32_t end;        				//region end
    int32_t value;
} gdata_t;

typedef struct{							//default 
    int32_t idx;        				//genomic object--data set index
    int32_t start;      				//region start
    int32_t end;        				//region end
} gdata0_t;

typedef struct{
    char* fileName;						//dataset file
    int32_t nr;							//number regions/dataset
    double md;    						//average width of the regions
} info_t;

typedef struct{
	int32_t ncnts, nCnts, mcnts;		//batch counts, total, max
	gdata_t *gList;
} tile_t;

typedef struct{
	int32_t ncnts, nCnts, mcnts;		//batch counts, total, max
	gdata0_t *gList;
} tile0_t;

typedef struct{
	char *name;    						//name of the contig
	int32_t mTiles;						//determined by the interval start and end 
	tile_t *gTile;                  	//tile data
} ctg_t;

typedef struct{
	char *name;    						//name of the contig
	int32_t mTiles;						//determined by the interval start and end 
	tile0_t *gTile;                  	//tile data
} ctg0_t;

typedef struct{		
	int32_t nbp, gType, nctg, mctg;		// number of base pairs, data type: 0, 1, 2 etc; size differs	
	int64_t total;						//total region in each ctg
	ctg_t *ctg;        					//list of contigs (of size _n_ctg_) 
} igd_t;

typedef struct{		
	int32_t nbp, gType, nctg, mctg;		// number of base pairs, data type: 0, 1, 2 etc; size differs	
	int64_t total;						//total region in each ctg
	ctg0_t *ctg;        					//list of contigs (of size _n_ctg_) 
} igd0_t;

typedef struct{							//for retrieving from disk file
	int32_t nFiles;
	info_t *finfo;
	char fname[64];
	int32_t nbp, gType, nCtg;			//data type: 0, 1, 2 etc; size differs	
	char **cName;						//name of ctgs
	int32_t *nTile;						//num of tiles in each ctg
	int32_t **nCnt;						//num of counts in each tile
	int64_t **tIdx;  					//tile index *sizeof -> location in .igd file
} iGD_t;

//---------Globals-----------------------------------------------------------------
extern void *hc;						//dict for converting contig names to int
extern iGD_t *IGD;
extern gdata_t *gData;
extern gdata0_t *gData0;
extern int32_t preIdx, preChr, tile_size;
extern FILE *fP;
//---------------------------------------------------------------------------------
//Parse a line of BED file
void str_splits( char* str, int *nmax, char **splits);
char *parse_bed(char *s, int32_t *st_, int32_t *en_);

//Binary search
int32_t bSearch(gdata_t *gdata, int32_t t0, int32_t tc, int32_t qe);
int32_t bSearch0(gdata0_t *gdata, int32_t t0, int32_t tc, int32_t qe);

//Add an interval
void igd_add(igd_t *igd, const char *chrm, int32_t s, int32_t e, int32_t v, int32_t idx);
void igd0_add(igd0_t *igd, const char *chrm, int32_t s, int32_t e, int32_t idx);

//Get id from igd dict
int32_t get_id(const char *chrm);

//Get file info from .tsv
info_t *get_fileinfo(char *ifName, int32_t *nFiles);

//Get igd info from .igd
iGD_t *get_igdinfo(char *igdFile);

//Initialize igd_t
igd_t *igd_init(void);
igd0_t *igd0_init(void);

//Save tile data
void igd_saveT(igd_t *igd, char *oPath);
void igd0_saveT(igd0_t *igd, char *oPath);

//Sort and save igd
void igd_save(igd_t *igd, char *oPath, char *igdName);
void igd0_save(igd0_t *igd, char *oPath, char *igdName);

//Free ailist data
void igd_destroy(igd_t *igd);
void igd0_destroy(igd0_t *igd);
//---------------------------------------------------------------------------------
//The following section taken from Dr Heng Li's cgranges
// (https://github.com/lh3/cgranges)

KSTREAM_INIT(gzFile, gzread, 0x10000)
/**************
 * Radix sort *
 **************/
#define RS_MIN_SIZE 64
#define RS_MAX_BITS 8

#define KRADIX_SORT_INIT(name, rstype_t, rskey, sizeof_key) \
	typedef struct { \
		rstype_t *b, *e; \
	} rsbucket_##name##_t; \
	void rs_insertsort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		rstype_t *i; \
		for (i = beg + 1; i < end; ++i) \
			if (rskey(*i) < rskey(*(i - 1))) { \
				rstype_t *j, tmp = *i; \
				for (j = i; j > beg && rskey(tmp) < rskey(*(j-1)); --j) \
					*j = *(j - 1); \
				*j = tmp; \
			} \
	} \
	void rs_sort_##name(rstype_t *beg, rstype_t *end, int n_bits, int s) \
	{ \
		rstype_t *i; \
		int size = 1<<n_bits, m = size - 1; \
		rsbucket_##name##_t *k, b[1<<RS_MAX_BITS], *be = b + size; \
		assert(n_bits <= RS_MAX_BITS); \
		for (k = b; k != be; ++k) k->b = k->e = beg; \
		for (i = beg; i != end; ++i) ++b[rskey(*i)>>s&m].e; \
		for (k = b + 1; k != be; ++k) \
			k->e += (k-1)->e - beg, k->b = (k-1)->e; \
		for (k = b; k != be;) { \
			if (k->b != k->e) { \
				rsbucket_##name##_t *l; \
				if ((l = b + (rskey(*k->b)>>s&m)) != k) { \
					rstype_t tmp = *k->b, swap; \
					do { \
						swap = tmp; tmp = *l->b; *l->b++ = swap; \
						l = b + (rskey(tmp)>>s&m); \
					} while (l != k); \
					*k->b++ = tmp; \
				} else ++k->b; \
			} else ++k; \
		} \
		for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k-1)->e; \
		if (s) { \
			s = s > n_bits? s - n_bits : 0; \
			for (k = b; k != be; ++k) \
				if (k->e - k->b > RS_MIN_SIZE) rs_sort_##name(k->b, k->e, n_bits, s); \
				else if (k->e - k->b > 1) rs_insertsort_##name(k->b, k->e); \
		} \
	} \
	void radix_sort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		if (end - beg <= RS_MIN_SIZE) rs_insertsort_##name(beg, end); \
		else rs_sort_##name(beg, end, RS_MAX_BITS, (sizeof_key - 1) * RS_MAX_BITS); \
	}

/*********************
 * Convenient macros *
 *********************/

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#define EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((a), (m)); \
	}while (0) 

#endif

