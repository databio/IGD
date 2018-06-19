# iGD: Reshape and integrate large-scale data sources for highly efficient genome analysis

Genome analysis usually requires comparing one set of genomic loci (region set) to many other annotated region sets that may come from one or more large-scale data sources. A wide variety of computational tools, such as LOLA, BART, Bedtools, Tabix, GenomeRunner and Giggle, have been developed for this purpose.  

This repository contains ...

An example for running iGD:
  1. Download iGD_b14.ipynb to your local computer
  2. Download roadmap_ex from data/ and mkdir (roadmap_igd)
  3. Open a terminal (linux) and cd to iGD_b14.ipynb folder 
  4. Type: jupyter notebook iGD_b14.ipynb
  5. In the last cell: Change the ifilePath to the roadmap_ex folder and ofilePath to the folder of roadmap_igd 
  6. Run all cells--an igd database roadmap_b14.igd and roadmap_index.tsv will be created in the roadmap_igd folder
