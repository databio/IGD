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

//create igd from .bed.gz files
void create_igd_gz(char *iPath, char *oPath, char *igdName, int mode);

//create igd from .bed.gz
void create_igd_gz1(char *iPath, char *oPath, char *igdName, int mode);

//create igd from .bed.gz: ailist data structure for in-bin search
void create_igd_gz2(char *iPath, char *oPath, char *igdName, int mode);

//create igd from a single .bed.gz file with the name/index at 4th column: Jaspar2018
void create_igd_gz3(char *iPath, char *oPath, char *igdName, int mode);

//create igd from .bed text files
void create_igd(char *iPath, char *oPath, char *igdName, int mode);

void create_igd1(char *iPath, char *oPath, char *igdName, int mode);

void create_igd2(char *iPath, char *oPath, char *igdName, int mode);

void constructAIList(struct igd_data2* B, int nB, struct igd_data2* aiL);

//main
int igd_create(int argc, char **argv);

#endif
