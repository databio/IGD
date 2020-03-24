#ifndef __IGD_SEARCH_H__
#define __IGD_SEARCH_H__

//=====================================================================================
//Search the igd database for all overlaps with queries
//by Jianglin Feng  05/12/2018
//
//time ./igd_search Test110000.bed /media/john/Extra/ucsc_igd/ucsc.igd
//database intervals sorted by _start: 8/12/2019
//Reorganized: 11/06/2019
//-------------------------------------------------------------------------------------
#include "igd_base.h"
//-------------------------------------------------------------------------------------
//for gdata0_t only-----
int32_t get_overlaps0(char *chrm, int32_t qs, int32_t qe, int64_t *hits);
int64_t getOverlaps0(char *qFile, int64_t *hits);

//for seqpare
void seq_overlaps(char *chrm, int32_t qs, int32_t qe, overlap_t *hits, uint32_t *n, uint32_t *m);
void seqOverlaps(char *qFile, overlap_t *hits);

//Single query
int32_t get_overlaps(char *chrm, int32_t qs, int32_t qe, int64_t *hits);
int32_t get_overlaps_v(char *chrm, int32_t qs, int32_t qe, int32_t v, int64_t *hits);

//query file: call _r
int64_t getOverlaps(char *qFile, int64_t *hits);
int64_t getOverlaps_v(char *qFile, int64_t *hits, int32_t v);	

//mapping
int64_t getMap(uint32_t **hitmap);
int64_t getMap_v(uint32_t **hitmap, int32_t v);

//search main
int igd_search(int argc, char **argv);

#endif

