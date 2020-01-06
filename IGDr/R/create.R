#==========================iGD in R========================================================
#----------------- Copyright (C) 2019 Jianglin Feng  --------------------------------------
#
#   This file is a part of the package IGDr. This function creates a igd
#   database from a folder of .bed files
#
#   The IGDr package is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   any later version.
#
#   The IGDr package is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#-----------------------------------------------------------------------------------------

#' Function to create an IGD database from a folder of .bed or .bed.gz files,
#' or a list of such folders
#'
#' @param iPath	folder where your input files are stored
#' @param oPath	the folder that the created IGD database will be stored
#' @param igdName the name you give to the IGD database (.igd will be added to it)
#' @param binsize the size in basepairs for the bin (block) used in the database:
#' usually 8192, 16384, 32768, ... as a power of 2; default 16384
#' @return an igd database will be created in the specified folder
#' @export
#' @examples
#' library("IGDr")
#' IGDr::createIGD("data/rme3", "data/testigd", "roadmap_b14")
createIGD <- function(iPath, oPath, igdName, binsize=16384) {
    .C("create_iGD", as.character(iPath), as.character(oPath), as.character(igdName), as.integer(binsize))
}

#' Function to create an IGD database from a list of source files (.bed or .bed.gz)
#'
#' @param iPath	path to a txt file that lists the paths of all the source files
#' @param oPath	the folder that the created IGD database will be stored
#' @param igdName the name you give to the IGD database (.igd will be added to it)
#' @param binsize the size in basepairs for the bin (block) used in the database:
#' usually 8192, 16384, 32768, ... as a power of 2
#' @return an igd database will be created in the specified folder
#' @export
createIGD_f <- function(iPath, oPath, igdName, binsize=16384) {
    .C("create_iGD_f", as.character(iPath), as.character(oPath), as.character(igdName), as.integer(binsize))
}
