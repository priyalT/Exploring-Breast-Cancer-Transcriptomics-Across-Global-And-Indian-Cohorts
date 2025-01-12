Description
============

This project is divided into various steps:

1. [Preprocessing](https://github.com/priyalT/Multicluster-and-Survival-Analysis/new/main?filename=README.md#preprocessing)
   
We went through each step thoroughly and the files have been arranged in the same way.

Processing Data
==============
* Curated BRCA data from TCGA using TCGABiolinks and it's many functions. Different forms of the same cohort's data was obtained, including gene expression, molecular subtypes, and clinical characteristics.
* Obtained Indian breast cancer data from a study done by IISC. We included gene expression data as well as patient characteristics data. (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40206)
* Indian data was already normalized. To prepare TCGA's genetic expression data for analysis, I had to log-normalize and then scale it. Transposing the dataframe was also required. All of this was achieved with TCGABiolinks and SummaryExperiments.
* Then, I mapped the EMBL ID's in TCGA's dataset to HGNC symbols, which were also seen in the Indian dataset. This prepared the TCGA dataset for further comparative analysis.
* Using intersect(), I subsetted both TCGA and Indian dataset to include common genes. This solidified both our datasets for comparison. 
