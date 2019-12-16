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

createIGD <- function(iPath, oPath, igdName, binsize=16384) {
    .C("create_iGD", as.character(iPath), as.character(oPath), as.character(igdName), as.integer(binsize))
}
