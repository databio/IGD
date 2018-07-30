---
title: "Using iGD"
author: "John Feng"
date: "06/18/2018"
---

This vignette shows how to create an iGD database from roadmap dataset and then search the iGD database for overlaps with a query file. Make sure your computer memory is at least 4GB.

## Create roadmap iGD database

First download the entire `rme` folder of the [roadmap data source](http://big.databio.org/igd/data) to the current working directory. The `rme` folder contains 1905 `.bed.gz` data files.

Then:  
```
mkdir rme_igd
igd create "rme/*" "rme_igd/" "roadmap"
```

This will generate the following in the output folder `rme_igd`:
1. a subfolder `data0` and its subfolders `chr1,..., chr22, chrX, chrY` with ~180,000 igd bin files (mode 0); 
2. a single igd database file (mode 1) `roadmap.igd` and dataset index file `roadmap_index.tsv`.


## Search iGD for overlaps

Download a sample query file [query100.bed](http://big.databio.org/igd/data) to the same directory as above.
 
Then:
```
igd search query100.bed rme_igd/roadmap.igd
```

