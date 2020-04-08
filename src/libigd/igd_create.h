#ifndef __IGD_CREATE_H__
#define __IGD_CREATE_H__

//=====================================================================================
//Create igd database
//by Jianglin Feng  05/12/2018
//-------------------------------------------------------------------------------------
#include "igd_base.h"

//create igd from .bed.gz files


/**
 * @brief Create an IGD database object (*.igd) from a folder
 *
 * @param i_path Input directory
 * @param o_path Output directory
 * @param igd_name Filename for output *.igd file
 * @param tile_size Tile size, must be power of 2 (default is 14)
 * @return A vector of overlap counts for each database file.
 */
void create_IGD(char **i_path, char **o_path, char **igd_name, int *tile_size);

int create_IGD_from_task(CreateTask_t* cTask);


#endif
