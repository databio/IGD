#ifndef __IGD_SEARCH_H__
#define __IGD_SEARCH_H__

//=====================================================================================
//Search the igd database for all overlaps with queries
//by Jianglin Feng  05/12/2018
//
//time ./igd_search Test110000.bed /media/john/Extra/ucsc_igd/ucsc.igd
//-------------------------------------------------------------------------------------
#include "igd_base.h"
//-------------------------------------------------------------------------------------
struct query_data* get_igdlist(char *qfName, uint32_t *nblocks, uint32_t *nRegions, double *mRegion);	

struct igd_info* get_igdinfo(char *ifName, uint32_t *nFiles);

struct igd_mix* get_overlaps(struct query_data *query, uint32_t nblocks, uint32_t nmax, uint32_t *nOL);

struct igd_mix* get_overlaps_w(struct query_data *query, uint32_t nblocks, char *igdName, uint32_t nmax, uint32_t *nOL);

uint64_t get_overlaps_self(char *igdName, uint32_t nFiles, uint32_t **hitmap);

uint64_t get_overlaps_n(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);//obsolete?

uint64_t get_overlaps_n0(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);

uint64_t get_overlaps_n1(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);

uint64_t get_overlaps_n2(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);

uint64_t get_overlaps_n0_c(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);

uint64_t get_overlaps_n1_c(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);

uint64_t get_overlaps_n2_c(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);

uint64_t get_overlaps_n0_r(char *igdName, int ichr, uint32_t qs, uint32_t qe);

uint64_t get_overlaps_n1_r(char *igdName, int ichr, uint32_t qs, uint32_t qe);

uint64_t get_overlaps_n2_r(char *igdName, int ichr, uint32_t qs, uint32_t qe);

int bSearch(struct igd_data2* As, int idxS, int idxE, uint32_t qe);

uint64_t get_overlaps_v(char *qfName, char *igdName, uint32_t v, uint32_t *nq, double *mq, uint32_t *hits);
uint64_t get_overlaps_self_v(char *igdName, uint32_t nFiles, uint32_t v, uint32_t **hitmap);

uint64_t get_overlaps_self_v_x(char *igdName, uint32_t nFiles, uint32_t v, int xlen, uint32_t *countf, uint32_t **hitmap);
uint64_t get_overlaps_self_x(char *igdName, uint32_t nFiles, int xlen, uint32_t **hitmap);

void search(char* igdName, char* qfName, uint32_t v, char* out, int checking);

void search_r(char* igdName, int ichr, uint32_t qs, uint32_t qe);

int igd_search(int argc, char **argv);

#endif

