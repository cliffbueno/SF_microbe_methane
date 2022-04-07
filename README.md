# Functions and analysis of microbes & CH4 with wetland salinity 

Analysis of 16S rRNA taxa & functional guilds in wetland metagenomes along a salinity gradient in the SF Bay & Delta.  Variation of interest in taxa and functions with soil CH4 fluxes, salinity, and wetland restoration.


## Overview

**Purpose**

This repo is project handoff primarily containing complex plotting / subsetting functions. Complex modeling functions are not yet shared here, and additional python functions / notebooks also are not yet shared.  This is absolutely gradware: "clean enough" to be followed, not much more.  May have bugs, caveat emptor.

Intended as a "pipeline" for passing 16S OTU tables to complex plots for publication.  Though it is an "open" pipeline to permit further customization & exploration on handoff.
Here specifically for addressing recalculation of 16S taxa using an updated Taxonomic assignment pipeline, in order to generate revised figures/tables.  

Original taxonomic data is included, and loading / variable assignment is left in notebooks, commented out.

**Use**
- R modules encapsulate preprocessing and correlation functions.
- Use of these functions for making plots & tables shown in .ipynb files.  
- Some additional features of functions shown for handoff for user exploration & hypothesis generation.

**Requirements**
- Code was tested under R 4.0.2, refactored where developed in earlier versions.  
- See `R_version_package_info.ipynb` for testing environment/package version details (not all packages required)


## Steps in Analysis
- Although precalculated results are shown, the following ordering of notebook runs is required to produce new analyses.
- Further details given in sections below.


### 1) Preprocessing OTU table
- SilvaOTUs are new OTU data.  Preprocessing of this data to fit with older pipeline is shown in `SilvaOTUs/Silva_OTUs_preprocessing.ipynb`
- Preprocessing reformats taxonomic strings, imposes a minimum count filter on OTUs, and imposes a variance stabilized DESeq2 transformation.


### 2) Taxon analysis
**a) Get 16S correlations** 
- req's OTU table.  Produces correlation by taxon rank tables  
- should be updated to log2 transformed data. 

**b) Taxonomic composition and patterns** 
- produces color palette required downstream

**c) NMDS plots of taxonomic based site clusters**
- requires gapfilled data, site colors, 2x grouping variables
- or change to a different function in module for < 2 group vars

**d) 16S Taxa correlation heatmaps** 
- requires precalculated correlations from 2a, taxa colors from 2b


### 3) Guild analysis 
**a) Guilds calculations and barplots** 
- requires OTU table, correct OTU subsetting module, taxa colors, guild colors

**b) Compare 16S & TreeSAPP guilds** 
- requires guild calculations & precalculated TreeSAPP guild counts

**c) Compare 16S guilds from iTagger & SILVA** (for data update) 
- requires guild calculations for each

**d) Guild correlation heatmaps** 
- requires guild calculations; OTU correlation tables from 2a



## Important notes to user

**Guild calculations & versions** 
- calculation of microbial "functional guilds" is performed by `modules\OTU_subsetting...`
- critically v.3.3 contains original taxa -> function mapping (iTagger)
- "new" SILVA taxa -> function mapping in v.3.4 (to be updated)

**Data / data versions**
- All data besides SILVA OTU preprocessing stored in '\data`
- Includes metadata, color palette definitions, other data here not generated programatically (e.g. sites)
- Taxonomic color palette generated in Taxonomic_compositions_and_patterns.ipynb
- 3 versions of metadata given (3.5 is complete; 3 is truncated for publication; 3b_gap_fill_MDS is used for NMDS with gaps filled)
- Post-preprocessed OTU table used has at minimum filled missing taxa values (here w. finest classification), "Consensus.lineage" string for sorting OTUs.
- Figures are NOT SAVED by script, use ggsave(ggsave('FILENAME.pdf', width=6, height = 4)
- Correlation calculations / saves commented out (slow, especially at OTU level)

- two modules to join / combine OTU & metadata & perform common calculations (though not consistently used here):
   "modules/Import_iTagger_OTU_data4plots.R" &  "modules/Import_SILVA_OTU_data4plots.R"

- also two modules to pre-process OTUs: "1_OTU_pre-process_module_0.2.r" & "1_OTU_preprocessing.R"
   latter was used for SILVA preprocessing, but perhaps should have been former (original)... 
   v0.2 includes calculation of "Taxonomy" (Phylum + Proteobact Classes)
   v0.2 also includes transformation of count data to numeric, though not gracefully - may be behind some downstream issues encountered.



## Modules use & dependencies
- In case of future modifications to modules

**1) OTU_preprocessing**
    - used in silvaOTUs/Silva_OTUs_preprocessing.ipynb

**2) OTU_table_to_DESeq2_and_VST_cpm**  
    - used in silvaOTUs/Silva_OTUs_preprocessing.ipynb

**3) OTU_subsetting_modules**
    - used by modules/7_Corr_heatmap_module.R
    - used in guild_analysis/Guilds_calculations_and_barplots.ipynb
    - used in guild_analysis/Guild_correlation_heatmaps.ipynb

**4) OTU_plotting_module_NMDS**
    - used in taxon_analysis/NMDS_plot_of_taxonomic_clusters.ipynb

**5) OTU_barplots_module**
    - used by modules/7_Corr_heatmap_module.R
    - used in taxon_analysis/Taxonomic_composition_and_patterns.ipynb
    - used in guild_analysis/Guilds_calculations_and_barplots.ipynb

**6) Corr_ranks_module**
    - used in taxon_analysis/Get_16S_correlations.ipynb

**7) Corr_heatmap_module**
    - used in taxon_analysis/16S_Taxa_correlation_heatmaps.ipynb
    - used in guild_analysis/Guild_correlation_heatmaps.ipynb
    - requires modules/3_OTU_subsetting_modules_v0.x...R
    - requires modules/5_OTU_barplots_module...R

**8) Guild_corr_heatmaps**
    - used in guild_analysis/Guild_correlation_heatmaps.ipynb

**Not numbered**
    * modules/Import_Silva_OTU_data4plots_v0.1.R" & 
    * modules/Import_iTagger_OTU_data4plots_v0.1.R"
    - used in guild_analysis/Guild_correlation_heatmaps.ipynb
    - load & process otu_V, metaDB, Meta_iTag, site colors (seen in other notebooks)
    - could be used in many other notebooks to streamline OTU + metadata manipulations

