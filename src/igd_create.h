#ifndef __IGD_CREATE_H__
#define __IGD_CREATE_H__

//=====================================================================================
//Create igd database 
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include "igd_base.h"

//create igd from .bed.gz files
void create_igd(char *iPath, char *oPath, char *igdName, int gtype);

//create igd from a single .bed.gz file with dataset index at 4th column
void create_igd1(char *iPath, char *oPath, char *igdName);

//main
int igd_create(int argc, char **argv);

#endif
