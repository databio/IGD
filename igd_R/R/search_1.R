#==========================iGD in R========================================================
#----------------- Copyright (C) 2019 Jianglin Feng  --------------------------------------
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
search_1 <- function(gdata, nd, ichr, qs, qe, hits) {
  #.C("search_1", as.integer(gdata), as.integer(nd), as.integer(ichr), as.integer(qs), as.integer(qe), hits=as.integer64(hits))$hits
}
