#ifndef __IGD_SEARCH_H__
#define __IGD_SEARCH_H__

//=====================================================================================
//Search the igd database for all overlaps with queries
//by Jianglin Feng  05/12/2018
//
//database intervals sorted by _start: 8/12/2019
//Reorganized: 11/06/2019
//-------------------------------------------------------------------------------------
#include "igd_base.h"
//-------------------------------------------------------------------------------------
//Single query
void get_overlaps(iGD_t *iGD, char *chrm, int32_t qs, int32_t qe, int64_t *hits);

//query file: call _r
void getOverlaps(char **igdFile, char **qFile, int64_t *hits);

#endif

