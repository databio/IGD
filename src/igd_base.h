//=====================================================================================
//Common structs, parameters, functions
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
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
#include <errno.h>
#include <sysexits.h>

#define PROGRAM_NAME  "igd"
#define MAJOR_VERSION "0"
#define MINOR_VERSION "1"
#define REVISION_VERSION "1"
#define BUILD_VERSION "0"
#define VERSION MAJOR_VERSION "." MINOR_VERSION "." REVISION_VERSION
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

extern char *fileBase;   			//14 bits block
extern char gfile[];   			//gfile:genome-sizes file
extern uint32_t nmax0[], *nmax;
extern char *folder0[], **folder;	//genome-chr name
extern uint32_t gstart0[], *gstart;
extern uint32_t nTiles;
extern uint32_t nbp, nChr;			//number of chr groups
extern uint32_t bgz_buf;
extern uint64_t maxCount;
extern uint32_t *g2ichr;

//-------------------------------------------------------------------------------------
struct igd_data
{   //default data
    uint32_t i_idx;        			//genomic object--data set index
    uint32_t r_start;      			//region start
    uint32_t r_end;        			//region end
    uint32_t g_val;        			//signal level
};

struct igd_data1
{   //region only data
    uint32_t i_idx;        			//genomic object--data set index
    uint32_t r_start;      			//region start
    uint32_t r_end;        			//region end
};

struct igd_data2                   //ailist data structure
{   //region only data
    uint32_t i_idx;        			//genomic object--data set index
    uint32_t r_start;      			//region start
    uint32_t r_end;        			//region end
    uint32_t r_max;               //augment  
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
void str_splits( char* str, int *nmax, char **splits);
int setup_igd(char* g_file);			//setup basic parameters
//-------------------------------------------------------------------------------------
int compare_iidx(const void *a, const void *b);

int compare_rend(const void *a, const void *b);

int compare_iidx1(const void *a, const void *b);

int compare_rend1(const void *a, const void *b);

int compare_qidx(const void *a, const void *b);

int compare_rstart(const void *a, const void *b);

int compare_rstart2(const void *a, const void *b);

int compare_midx(const void *a, const void *b);

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
long double _hypergeo(long long nA, long long nAC, long long nAB, long long n);

typedef struct {
    long long nA, nAC, nAB, n;
    long double p;
} _hgacc_t;

// incremental version of hypergenometric distribution
long double _hypergeo_acc(long long nA, long long nAC, long long nAB, long long n, _hgacc_t *aux);
long double _kt_fisher_exact(long long nA, long long nC, long long nB, long long nD, long double *_left, long double *_right, long double *two);
double log2fc(double ratio);
long double neglog10p(long double sig);

#endif

