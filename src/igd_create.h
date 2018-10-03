#ifndef __IGD_CREATE_H__
#define __IGD_CREATE_H__

//=====================================================================================
//Create igd database 
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include "igd_base.h"

//-------------------------------------------------------------------------------------
void append_igd(struct igd_mix *mdata, uint32_t *counts, struct igd_info *fInfo, int nFiles, char* igdName);

//reload tile igd files, sort them and save them into a single file
void store_igd(char *igdName);

//create ucsc igd from gz
void create_igd_gz3(char *iPath, char *oPath, char *igdName, int mode);

//create ucsc igd plain bed files
void create_igd3(char *iPath, char *oPath, char *igdName, int mode);

//main
int igd_create(int argc, char **argv);

#endif
