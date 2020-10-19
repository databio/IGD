#ifndef __IGD_CREATE_H__
#define __IGD_CREATE_H__

//=====================================================================================
//Create igd database 
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include "igd_base.h"

//create igd from .bed.gz files
void create_iGD(iGD_t *iGD, char *iPath, char *oPath, char *igdName, int tile_size);

#endif
