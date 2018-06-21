# iGD: Reshape and integrate large-scale data sources for highly efficient genome analysis
## Summary:
Genome analysis usually requires comparing one set of genomic loci (region set) to many other annotated region sets that may come from one or more large-scale data sources. A wide variety of computational tools, such as LOLA, BART, Bedtools, Tabix, GenomeRunner and Giggle, have been developed for this purpose. Among these tools, Giggle claims to be over three orders of magnitude faster than others in searching millions to billions of genomic intervals. iGD, an integrated genome database, takes a database approach that directly integrates large-scale data sources into an easily searchable database by reshaping the data records. As a database, iGD takes not only the genomic regions, but also the signal levels and/or strand info, which makes the parameterized query (dynamic searching) an easy task. The searching speed of iGD is one to two orders of magnitude faster than Giggle while the database size is several times smaller than that of the Giggle index files. 

## Introduction:
Genome data sources, such as Roadmap, UCSC, and Cistrome, contain thousand of annotated data sets and each dataset has thousands to millions of genomic regions. 
The enrichment analysis is to find overlaps between regions of all data sets in the data source and all regions in the query data set. This becomes complicated when the size of data source regions and the size of query regions are large. For example, UCSC human data source contains ~7 billion of regions or intervals, if the query set has 100,000 interavls, then a total of 100 trillion comparisons may be involved. To effectively search the overlaps, different indexing techniques have been developed. For example, Giggle employs a B+ tree to build a bunch of pure indexing files and then a qurey is carried out on these indexing files instead of the original data. With the indexing, the searching space is significantly reduced. 


 

An example for running iGD:
  1. Download iGD_b14.ipynb to your local computer
  2. Download roadmap_ex from data/ and mkdir (roadmap_igd)
  3. Open a terminal (linux) and cd to iGD_b14.ipynb folder 
  4. Type: jupyter notebook iGD_b14.ipynb
  5. In the last cell: Change the ifilePath to the roadmap_ex folder and ofilePath to the folder of roadmap_igd 
  6. Run all cells--an igd database roadmap_b14.igd and roadmap_index.tsv will be created in the roadmap_igd folder
