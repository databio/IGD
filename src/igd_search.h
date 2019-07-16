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
//Single query
int32_t get_overlaps_r(char *chrm, int32_t qs, int32_t qe, int32_t *hits);
int32_t get_overlaps_r1(char *chrm, int32_t qs, int32_t qe, int32_t *hits);
int32_t get_overlaps_r2(char *chrm, int32_t qs, int32_t qe, int32_t *hits);

//query file: call _r
int64_t get_overlaps(char *qFile, int32_t *hits);
int64_t get_overlaps1(char *qFile, int32_t *hits);

//search with value, mapping
int64_t get_overlaps_v(char *qFile, int32_t *hits, int32_t v);		//dtype 1
int64_t get_overlaps_m(char *qFile, int32_t *hits, int32_t v); 		//map, dtype
int64_t get_overlaps_m_x(char *qFile, int32_t *hits, int32_t v, int32_t x); //map, q extended

//search main
int igd_search(int argc, char **argv);

#endif

