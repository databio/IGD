#!/usr/bin/env python

import igd_py as iGD
from timeit import default_timer as timer
import sys
import numpy as np
import pandas as pd

def main(argv):
	if len(argv) < 4:
		print("To create: igd_test.py create <path to source folder> <path to output folder> <name for igd> \n \
			To search: igd_test.py search <path to igd file> <query file>")
		sys.exit(1)	

	igd = iGD.igd_py()
	if argv[1]=="create" and len(argv)>=5:
		igd.create(argv[2], argv[3], argv[4], 16384)	
	
	elif argv[1]=="search" and len(argv)>=4:
		igd.open(argv[2])
		nFiles = igd.get_nFiles()
		hits = np.zeros(nFiles, dtype='int64')
		total = igd.search_n(argv[3], hits)		
		print("Total: ", total, "\n")
	
	print("nFiles: ", igd.get_nFiles(), "\n")
		
if __name__ == "__main__":
	main(sys.argv)
