from libc.stdint cimport int64_t, int32_t, uint32_t
from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np

cdef extern from "igd_base.h":
	ctypedef struct iGD_t:
		pass
	iGD_t *iGD_init()
	int32_t get_nFiles(iGD_t *iGD)
	void open_iGD(iGD_t *iGD, char *igdFile) 
	void close_iGD(iGD_t *iGD)
    
cdef extern from "igd_create.h":
	void create_iGD(iGD_t *iGD, char *iPath, char *oPath, char *igdName, int tile_size)
	
cdef extern from "igd_search.h":
	void get_overlaps(iGD_t *iGD, char *chrm, int32_t qs, int32_t qe, int64_t *hits)	  
	int64_t getOverlaps(iGD_t *iGD, char *qFile, int64_t *hits)
    
cdef class igd_py:
	cdef iGD_t *iGD
	
	def __cinit__(self):
		self.iGD = iGD_init()

	def __dealloc__(self):
		close_iGD(self.iGD)
		
	def get_nFiles(self):
		return get_nFiles(self.iGD)
	
	def create(self, iPath, oPath, igdName, bin_size):
		create_iGD(self.iGD, str.encode(iPath), str.encode(oPath), str.encode(igdName), int(bin_size))
	
	def open(self, igdFile):
		open_iGD(self.iGD, str.encode(igdFile))
		
	def search_1(self, chrm, np.int32_t qs, np.int32_t qe, np.ndarray[np.int64_t, ndim=1, mode="c"] hits not None):
		get_overlaps(self.iGD, str.encode(chrm), qs, qe, &hits[0])
				
	def search_n(self, qFile, np.ndarray[np.int64_t, ndim=1, mode="c"] hits not None):
		cdef int64_t nols = getOverlaps(self.iGD, str.encode(qFile), &hits[0])
		return nols
    
