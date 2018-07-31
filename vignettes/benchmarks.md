---
title: "benchmarks"
author: "John Feng"
date: "07/18/2018"
---

This vignette shows benchmarks of iGD for searching igd databases `roadmap.igd` and `ucsc.igd` with `igd search`. The corressponding benchmarks of [giggle](https://github.com/ryanlayer/giggle) index files `rme_giggle` and `ucsc_giggle` are also presented. 

## Benchmarks for roadmap data source

### iGD:
We have pre-built the roadmap.igd database using Roadmap data. Download and unzip the roadmap database (1.3GB) like this:
```
wget http://big.databio.org/igd/rme_igd.tgz
tar -xf rme_igd.tgz
```
A series query files of various number of region sets are used for testing:
```
wget http://big.databio.org/igd/queries_igd.tgz
tar -xf queries_igd.tgz
```
To run the speed test:
```
for i in 10 100 1000 10000 100000 1000000
do
time igd search queries_igd/r$i.bed rme_igd/roadmap.igd >/dev/null 
done
```
The following results are from a laptop of 2.5GHz cpu:
```
0.009s, 0.009s, 0.014s, 0.051s, 0.338s, 2.085s.
```
### giggle:
Download and unzip the pre-indexed roadmap giggle database (~2.3GB) and queries:
```
wget http://big.databio.org/igd/rme_giggle.tgz
tar -xf rme_giggle.tgz
wget http://big.databio.org/igd/queries_giggle.tgz
tar -xf queries_giggle.tgz
```
Run the speed test:
```
for i in 10 100 1000 10000 100000 1000000
do
time giggle search -i rme_giggle -q queries_giggle/ucsc_r$i.bed.gz >/dev/null 
done
```
The corresponding times are:
```
0.013s, 0.017s, 0.064s, 0.582s, 5.241s, 45.439s.
``` 

## Benchmarks for ucsc data source

### iGD:
Download and unzip the pre-built ucsc_igd database (~123GB):
```
wget http://big.databio.org/igd/ucsc_igd.tgz
tar -xf ucsc_igd.tgz
```
Run the speed test:
```
for i in 10 100 1000 10000 100000 1000000
do
time igd search queries_igd/r$i.bed ucsc_igd/ucsc.igd >/dev/null 
done
```
The times are:
```
0.206s, 1.671s, 16.239s, 2m24.222s, 13m26.413s, 17m6.142s. 
```

### giggle:
Download and unzip the pre-indexed ucsc giggle database (~620GB):
```
wget http://big.databio.org/igd/ucsc_giggle.tgz
tar -xf ucsc_giggle.tgz
```
Run the speed test:
```
for i in 10 100 1000 10000 100000 1000000
do
time giggle search -i ucsc_giggle -q queries_giggle/ucsc_r$i.bed.gz >/dev/null 
done
```
The results are:
```
17.169s, 1m18.631s, 5m31.372s, 18m27.102s, 73m50.718s, 138m20.332(killed for unknown reason).
```

