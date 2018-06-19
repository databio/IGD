# iGD: Reshape and integrate large-scale data sources for highly efficient genome analysis

Genome analysis usually requires comparing one set of genomic loci (region set) to many other annotated region sets that may come from one or more large-scale data sources. A wide variety of computational tools, such as LOLA, BART, Bedtools, Tabix, GenomeRunner and Giggle, have been developed for this purpose. Among these tools, Giggle claims to be over three orders of magnitude faster than others in searching millions to billions of genomic intervals. iGD takes a database approach that directly integrates large-scale data sources into an easily searchable database by reshaping the data records. As a database, iGD takes not only the genomic regions, but also the signal levels and/or strand info, which makes the parameterized query (dynamic searching) possible. The searching speed of iGD is one to two orders of magnitude faster than Giggle while the database size is serval times smaller than that of the Giggle index files. 

#Basic idea

An example for running iGD:
  1. Download iGD_b14.ipynb to your local computer
  2. Download roadmap_ex from data/ and mkdir (roadmap_igd)
  3. Open a terminal (linux) and cd to iGD_b14.ipynb folder 
  4. Type: jupyter notebook iGD_b14.ipynb
  5. In the last cell: Change the ifilePath to the roadmap_ex folder and ofilePath to the folder of roadmap_igd 
  6. Run all cells--an igd database roadmap_b14.igd and roadmap_index.tsv will be created in the roadmap_igd folder
