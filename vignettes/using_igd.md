---
title: "Using iGD"
author: "John Feng"
date: "06/18/2018"
updated: "5/8/2020"
---

This vignette shows how to create an iGD database from roadmap dataset and then search the iGD database for overlaps with a query file. Make sure your computer memory is at least 4GB.

## Create roadmap iGD database

First download the roadmap data source [rme.tgz](http://big.databio.org/igd/data/rme.tgz) and extract it:
```
wget http://big.databio.org/igd/data/rme.tgz
tar -xzf rme.tgz
```
The `rme` folder contains 1905 `.bed.gz` data files.

Then: 
```
mkdir rme_igd
igd create "rme/*" "rme_igd/" "roadmap"
```

This will generate the following in the output folder `rme_igd`:
a single igd database file (mode 1) `roadmap.igd` and dataset index file `roadmap_index.tsv`.


## Search iGD for overlaps

Download a sample query file [query100.bed](http://big.databio.org/igd/data/query100.bed) to the same directory as above.
 
Then:
```
wget http://big.databio.org/igd/data/query100.bed
igd search rme_igd/roadmap.igd -q query100.bed 
```

