#ifndef __IGD_CREATE_H__
#define __IGD_CREATE_H__

//=====================================================================================
//Create igd database 
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include "igd_base.h"

//create igd from .bed.gz files
void create_igd(char *iPath, char *oPath, char *igdName);
void create_igd0(char *iPath, char *oPath, char *igdName);
void create_igd_f(char *iPath, char *oPath, char *igdName);//create from a file list
//create igd from a single .bed.gz file with dataset index at 4th column
void create_igd_Jaspar(char *iPath, char *oPath, char *igdName);

//main
int igd_create(int argc, char **argv);

#endif
