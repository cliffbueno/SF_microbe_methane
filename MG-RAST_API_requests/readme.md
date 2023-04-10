## MG-RAST API requests code

This directory contains scripts for interfacing with the MG-RAST API:  http://api.metagenomics.anl.gov/api.html
In particular, scripts here allow one to obtain from a list of mgIDs:

1) Table of gene abundances per sample (here as KO), using API to KO tables script.  Samples are obtained individually with the requests library, and JSON files parsed to obtain single tables per file.
These files are then combined with reduce into a single table.  Here, I rename the samples in the table using an external mapping of my own sample names to mgIDs.

2) The API tables to KO script also produces a list of md5_nr IDs mapped to each KO.  These IDs are unique gene (KO) - species combinations, and their mapping is used by API get Taxonomy.

3) The API get taxonomy script first gets a mapping of md5_nr IDs to GenBank IDs, then full taxonomy of these IDs from MG-RAST's database.

4) The API get taxonomy then merges md5_ID taxonomic information with the md5_ID to KO mapping produced in 2).  This final file contains the abundance of each taxonomic group WITHIN each KO ID for each sample.
However, here I have limited this output to only specific genes of interest.

Several functions of general utility in working with MG-RAST's API are contained in these scripts, but these are not yet a completely modular pipeline.  Some aspects of the scripts are specific to my data, including samples of interest and culled lists of genes focal to my analyses.
Caveat Emptor: Any user of these script is cautioned that intermediate files are quite large (4.5 GB for 130 samples) and numerous.  I wrote many intermediate files to disk, especially API outputs, in order to not need to run these steps again.  Queing data calculated by MG-RAST took nearly one day, and about as much time elapsed in downloading these products for my samples.