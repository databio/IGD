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

uint64_t get_overlaps_n(char *qfName, char *igdName, uint32_t *nq, double *mq, uint32_t *hits);

uint64_t get_overlaps_v(char *qfName, char *igdName, uint32_t v, uint32_t *nq, double *mq, uint32_t *hits);

void search(char* qfName, char* igdName, uint32_t v);

int igd_search(int argc, char **argv);

#endif

