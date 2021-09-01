# IGD: A high-performance search engine for large-scale genomic interval datasets

## Summary

Databases of large-scale genome projects now contain thousands of genomic interval datasets. These data are a critical resource for understanding the function of DNA. However, our ability to examine and integrate interval data of this scale is limited. Here, we introduce the integrated genome database (IGD), a method and tool for searching genome interval datasets more than three orders of magnitude faster than existing approaches, while using only one hundredth of the memory. IGD uses a novel linear binning method that allows us to scale analysis to billions of genomic regions. 

## Citation

If you use IGD in your research, please cite:


Jianglin Feng, Nathan C Sheffield. IGD: high-performance search for large-scale genomic interval datasets. *Bioinformatics*, Volume 37, Issue 1, 1 January 2021, Pages 118–120, https://doi.org/10.1093/bioinformatics/btaa1062

Preprint: https://www.biorxiv.org/content/10.1101/2020.06.08.139758v1


## How to build iGD

If zlib is not already installed, install it:
```
sudo apt-get install libpng12-0
```
Then:
```
git clone https://github.com/databio/iGD.git
cd iGD
make
```
the executable `igd` is in the subfolder `bin`. And then copy it to /usr/local/bin.

## How to run iGD

### 1. Create iGD database
 
#### 1.1 Create iGD database from a genome data source folder
```
igd create "/path/to/data_source_folder/*" "/path/to/igd_folder/" "databaseName" [option]

where:

- "path/to/data_source_folder/" is the path of the folder that contains .bed.gz or .bed data files.

- "path/to/igd_folder/" is the path to the output igd folder;

- "databaseName" is the name you give to the database, for eaxmple, "roadmap"

option:

-b: bin-size (power of 2; default 14, which is 16384 bp)
```
#### 1.2 Create iGD database from a list of source files
 
```
igd create "/path/to/source-list file" "/path/to/igd_folder/" "databaseName" -f [option]

where:

- "/path/to/source-list file" is the path to the file that lists the source files

- "path/to/igd_folder/" is the path to the output igd folder;

- "databaseName" is the name you give to the database, for eaxmple, "roadmap"

option:

-b: bin-size (power of 2; default 14, which is 16384 bp)
```


### 2. Search iGD for overlaps
```
igd search "path/to/igd_data_file" -q "path/to/query_file"

where:

- path/to/igd_data_file is the path to the igd data

- path/to/query_file is the path to the query file (.bed or .bed.gz)

other options:

-r <chrN start end> (a single query)

-v <signal value 0-1000> (signal value > v)

-o <output file Name>

-s (output Seqpare similarity)

-f (output full overlaps, for -q and -r only)

-m (hitsmap of igd datasets)

```

For a detailed example, please check out the `vignettes`.

## R-wrapper of IGD

### 1. Create iGD database 
#### 1.1  from a genome data source
```
> library(IGDr)
> createIGD("/path/to/data_source_folder/*" "/path/to/igd_folder/" "databaseName" [option]

where:

- "path/to/data_source_folder/" is the path of the folder that contains .bed.gz or .bed data files.

- "path/to/igd_folder/" is the path to the output igd folder;

- "databaseName" is the name you give to the database, for eaxmple, "roadmap"

options:

-b: bin size in bp (default 16384)
```
#### 1.2  from a file that contains the list of genome data source files 
```
> library(IGDr)
> createIGD_f("/path/to/source-list file" "/path/to/igd_folder/" "databaseName" [option]

where:

- "path/to/the list file/" is the path to the file that contains the .bed.gz or .bed data files.

- "path/to/igd_folder/" is the path to the output igd folder;

- "databaseName" is the name you give to the database, for eaxmple, "roadmap"

options:

-b: bin size in bp (default 16384)
```

### 2. search the igd database in R (an example for a created igd file)

Search the igd database with a single query:
```
> igd_file = "igdr_b14/roadmap.igd"
> library(IGDr)
> igd <- IGDr::IGDr(igd_file)
> hits <- search_1r(igd, "chr6", 1000000, 10000000)
> hits
```
Search the igd database with n queries:
```
> igd_file = "igdr_b14/roadmap.igd"
> library(IGDr)
> igd <- IGDr::IGDr(igd_file)
> chrms = c("chr6", "chr1", "chr2")
> starts = c(10000, 100000, 1000000)
> ends = (100000, 1000000, 10000000)
> hits <- search_nr(igd, 3, chrms, starts, ends)
> hits
```
Search a whole query file chainRn4.bed
```
> igd_file = "igdr_b14/roadmap.igd"
> query_file = "r10000.bed"
> library(bit64)
> library(IGDr)
> fi = IGDr::getFInfo(igd_file)
> hits = integer64(fi$nFiles)
> ret = IGDr::search_all(igd_file, query_file, hits)
> for(i in 1:fi$nFiles){
  cat(i, "\t", toString(ret[i]), "\t", toString(fi$fInfo[i,2]), "\n")
  }
>
```
