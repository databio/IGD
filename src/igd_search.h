#ifndef __IGD_SEARCH_H__
#define __IGD_SEARCH_H__

//=====================================================================================
//Search the igd database for all overlaps with queries
//by Jianglin Feng  05/12/2018
//
//time ./igd_search Test110000.bed /media/john/Extra/ucsc_igd/ucsc.igd
//database intervals sorted by _start: 8/12/2019
//-------------------------------------------------------------------------------------
#include "igd_base.h"
//-------------------------------------------------------------------------------------
//Single query
int32_t get_overlaps(char *chrm, int32_t qs, int32_t qe, int32_t *hits);
int32_t get_overlaps_v(char *chrm, int32_t qs, int32_t qe, int32_t v, int32_t *hits);

//query file: call _r
int64_t getOverlaps(char *qFile, int32_t *hits);
int64_t getOverlaps_v(char *qFile, int32_t *hits, int32_t v);	

//construct for mapping
void construct(gdata_t *glist, int32_t nr, int32_t *nc, int32_t *idxC, int32_t *lenC, int32_t *maxE, int cLen);

//search with value, mapping
int64_t getOverlaps_m0(uint32_t **hitmap, int32_t v); 	//map, dtype
int64_t getOverlaps_m1(uint32_t **hitmap); 	//map, dtype, test
int64_t getOverlaps_m2(uint32_t **hitmap); 	//map, dtype, test
int64_t getOverlaps_m1_v(uint32_t **hitmap, int32_t v); 	//map, dtype, test
int64_t getOverlaps_m2_v(uint32_t **hitmap, int32_t v); 	//map, dtype, test
int64_t getOverlaps_m0_x(uint32_t **hitmap, int32_t v, int32_t x); //map, q extended
int64_t getOverlaps_m1_x(uint32_t **hitmap, int32_t v, int32_t x); //map, q extended
int64_t getOverlaps_m2_x(uint32_t **hitmap, int32_t v, int32_t x); //map, q extended

//search main
int igd_search(int argc, char **argv);

#endif

