---
title: "benchmarks"
author: "John Feng"
date: "07/18/2018"
---

This vignette shows benchmarks of iGD for searching igd databases `roadmap.igd` and `ucsc.igd` with `igd search`. The corressponding benchmarks of [giggle](https://github.com/ryanlayer/giggle) index file folders `rme_index` and `ucsc_index` are also presented. 

## Benchmarks for roadmap.igd

Database `roadmap.igd` is created from Roadmap data source by using iGD, see vignettes `using_igd.Rmd` for detailed information. The compressed file `rme_igd.tgz` can be downloaded [here](http://big.databio.org/igd). Unzip it after download(`~1.3GB` with folder `rme_igd`).

A series query files [queries.tgz](http://big.databio.org/igd) of various number of region sets are used for testing:

```
for i in 10 100 1000 10000 100000 1000000
do
time igd search queries/r$i.bed rme_igd/roadmap.igd >/dev/null 
done
```
The following speed results are from a laptop of 2.5GHz cpu:
0.009s, 0.009s, 0.014s, 0.051s, 0.338s, 2.085s.

For giggle: download the index `rme_index.tgz` and queries `queries_gz.tgz` from http://big.databio.org/igd. The size of unziped rme_index is ~2.3GB.

```
for i in 10 100 1000 10000 100000 1000000
do
time giggle search -i rme_index -q queries_gz/ucsc_r$i.bed.gz >/dev/null 
done
```
The corresponding times are:
0.013s, 0.017s, 0.064s, 0.582s, 5.241s, 45.439s. 


## Benchmarks for ucsc.igd

Database `ucsc.igd` is created from ucsc data source and it can be downloaded from [databio.org/igd](http://big.databio.org/igd). The compressed file `ucsc.igd.tgz` is 22GB and unzipped database file ucsc.igd is about ~123GB in folder `ucsc_igd`. 

Using the same query series as above for testing:

```
for i in 10 100 1000 10000 100000 1000000
do
time igd search queries/r$i.bed ucsc_igd/ucsc.igd >/dev/null
done
```
The times are:
0.206s, 1.671s, 16.239s, 2m24.222s, 13m26.413s, 17m6.142s. 

For giggle: download the index `ucsc_index.tgz` at [here](http://big.databio.org/igd). The size of unziped ucsc_index is ~620GB.

```
for i in 10 100 1000 10000 100000 1000000
do
time giggle search -i ucsc_index -q queries_gz/ucsc_r$i.bed.gz >/dev/null 
done
```
The results are:
17.169s, 1m18.631s, 5m31.372s, 18m27.102s, 73m50.718s, 138m20.332(killed for unknown reason).

