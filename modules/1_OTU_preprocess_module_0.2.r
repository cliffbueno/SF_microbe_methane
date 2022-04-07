
#library(phyloseq)
#library(ggplot2)
#library(RColorBrewer)
library(stringr)
library(reshape2)

#suppressMessages(library(vegan))
#suppressMessages(library(DESeq2))
#suppressMessages(library(gplots))
#suppressMessages(library(Hmisc))

### IMPORT Sample mapping
# metaDB <-read.table("SF_sal_meta_FIX3.txt", sep="\t", header=TRUE)                            # import Mapping    # 

# try keeping all params...
# row.names(metaDB) <- metaDB$Sample               # Row names are samples for phyloseq             #head(map_iTag)
# metaDB = metaDB[,-1]                             # Drop only old index, keep everything else                         

## IMPORT OTU TABLE                              # Raw (w. taxon)
# otu_raw <- read.table("SF_Salinity_gradient_OTU_table.txt", sep='\t', header=TRUE, row.names=1)      # add to fxn below?
#otu_tax <-data.frame(OTU= row.names(otu_raw), Taxonomy = otu_raw$Consensus.lineage)      

# New plaintext OTU ids needed for several transposes below                                    
# row.names(otu_raw)<-paste0("otu_",row.names(otu_raw))                   # added to matrix (for colnames)
#row.names(otu_raw)<-paste0("otu_",OTU_VST_CPM_tax$OTU)                 # added to Taxonomy mtx (for taxon recovery) 

filter_n_otu = function(otu_t, min_reads) {                   # Filter OTU table by minimum count across all samples
                                                              # Assumes last col is Consensus.lineage, OTUs row index
    otu_t$OTU_tot = rowSums(otu_t[,1:(ncol(otu_t)-1)])        # Get row totals
    otu_F = otu_t[(otu_t$OTU_tot > min_reads),]               # Filter by total, min_reads
    otu_F = otu_F[,1:(ncol(otu_F)-1)]                         # drop OTU total col      
    return(otu_F)                                             # return filtered OTU table
}

# Apply function as test
#otu_F = filter_n_otu(otu_raw, 500) 
# dim(otu_F); head(otu_F)

sort_samples_by_meta = function(otu_F, MetaDB, Sample, index_col) { 
# assumes column names in OTU table correspond to "Sample" column in MetaDB

    OTU_samps <- data.frame(Sample=colnames(otu_F))                        # Get Sample names from OTU column names          # I#OTU_samps
    Meta_iTag <- merge(metaDB, OTU_samps, by=Sample)                       # Merge site order and Samples, drops non-matching                  
    rownames(Meta_iTag) <- Meta_iTag$Sample                                # Set row names as Sample 

    # Get site sorting data matrix
    sortr <- c(Sample,index_col)                                           # site, index name list
    Site_sort <-Meta_iTag[,sortr]                                          # slice DF by list 
    colnames(Site_sort) <- c('Sample','Site_index')                        # rename columns   -- needs to be generalized #Site_sort
    dim(Site_sort)

    ### Reorder Samples in OTU table      
    otu_Ft <- t(otu_F)                                                      # transpose OTU table
    otu_FT <- data.frame(Sample=row.names(otu_Ft), otu_Ft)                  # add Sample col

    otu_F_metaT <- merge(otu_FT, Site_sort, by="Sample", all=TRUE)          # Merge with site sorter, keep all
    otu_F_metaT <- otu_F_metaT[order(otu_F_metaT$Site_index),]              # Sort samples by meta index
    otu_F_metaT <- otu_F_metaT[,1:(ncol(otu_F_metaT)-1)]                    # Drop Site index column
    row.names(otu_F_metaT) <- otu_F_metaT$Sample                            # add back rownames
    otu_F_reorder <- data.frame(t(otu_F_metaT[,-1]))                        # transpose while drop new samp column

    return(otu_F_reorder)
}
#dim(otu_F); dim(otu_F_reorder); head(otu_F_reorder); tail(otu_F_reorder)
# head(Meta_iTag)

# otu_S = sort_samples_by_meta(otu_F, MetaDB, "Sample", 'EWsiteHyd_index')
# dim(otu_S); head(otu_S)

clean_OTU_taxon_6ranks = function(OTUin) {
    # Function takes an OTU table with Consensus.lineage, returns OTU table with addtl 6 taxonomic levels for Phyloseq
    # fxn also drops poorly ID'd taxa, cleans some misc delimiters
    # importantly, fxn fills in missing (unknown) taxonomy at each level with next highest taxon, e.g. Nitrospira_CL
    
    ######################################################
    ## 1) Split strings, clean delimiters, misc text ## 

    # Split Taxonomy strings
    TaxSplit <- str_split_fixed(as.character(OTUin$Consensus.lineage),';', 6)
    colnames(TaxSplit)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")       # previously made phyloseq b4 cleaning # TaxTable<-tax_table(NameNumber)   # Make Phyloseq table # head(TaxTable) # attributes(TaxSplit)

    # Clean string delimiters
    TaxSplit <- data.frame(TaxSplit)
    TaxSplit<-as.data.frame(sapply(TaxSplit, gsub,pattern="k__",replace=""))
    TaxSplit<-as.data.frame(sapply(TaxSplit, gsub,pattern="p__",replace=""))
    TaxSplit<-as.data.frame(sapply(TaxSplit, gsub,pattern="c__",replace=""))
    TaxSplit<-as.data.frame(sapply(TaxSplit, gsub,pattern="o__",replace=""))
    TaxSplit<-as.data.frame(sapply(TaxSplit, gsub,pattern="f__",replace=""))
    TaxSplit<-as.data.frame(sapply(TaxSplit, gsub,pattern="g__",replace=""))             # TaxSplit

    # Clean misc. text: 
    TaxSplit<-as.data.frame(sapply(TaxSplit, gsub,pattern="\\s\\(class\\)",replace=""))  # Drop " (class)", note double esc.
    TaxSplit<-as.data.frame(sapply(TaxSplit, gsub,pattern="Candidatus",replace="Cand.")) # "Canditatus" genera -> "Cand."
    row.names(TaxSplit) <-row.names(OTUin)                                               # Add back original OTU numbers                  # head(TaxSplit) #    # TaxSplit

    ######################################################
    ## 2) Drop unwanted groups, poorly ID seq's ## 

    # Drop chloroplast & mitochondrial reads
    TaxSplit2 <-data.frame(TaxSplit)                                                     # not needed anymore, but leave for numbering # levels(TaxSplit2$Kingdom)  # TaxSplit2$Genus
    TaxSplit3 <- TaxSplit2[!grepl("chloroplast", TaxSplit2$Kingdom),]                    # drop chloropl.   
    TaxSplit4 <- TaxSplit3[!grepl("mitochondria", TaxSplit3$Kingdom),]                   # drop mitochon.             # dim(TaxSplit2) # dim(TaxSplit4)

    # Drop reads with no phylum
    TaxSplit5 <- TaxSplit4[!grepl("BacteriaKI", TaxSplit4$Phylum),]                      # assigned KI only 
    TaxSplit5 <- TaxSplit4[!grepl("ArchaeaKI", TaxSplit4$Phylum),]                       # assigned KI only 
    levels(TaxSplit5$Phylum)[levels(TaxSplit5$Phylum)==""] <- NA_character_              # assign NA to Empty Phyla #  Magic, see is.na fxn doc
    TaxClean0 <- TaxSplit5[!is.na(TaxSplit5$Phylum),]                                    # drop NA rows             #  dim(TaxClean0); dim(TaxSplit2)#levels(TaxClean$Phylum)  #  # TaxClean0

    ######################################################
    ## 3) Fill missing taxonomic levels, rebuild Tax Table ## 

    # Define Tax rank cols (needed below, avoiding attach method)
    OTU <- row.names(TaxClean0)
    Kingdom <- TaxClean0$Kingdom
    Phylum <- TaxClean0$Phylum
    class <- TaxClean0$Class
    order <- TaxClean0$Order
    family <- TaxClean0$Family
    genus <- TaxClean0$Genus

    # Fill missing, note capitalization for new vars, fill is sequential using new if missing        # Phylum<-ifelse(TaxClean0$phylum=="",paste(Kingdom,"KI",sep=""),paste(phylum))      # already deleted these
    Class<-ifelse(TaxClean0$Class=="",paste(Phylum,"PH",sep=""),paste(class))
    Order<-ifelse(TaxClean0$Order=="",paste(Class,"CL",sep=""),paste(order))
    Family<-ifelse(TaxClean0$Family=="",paste(Order,"OR",sep=""),paste(family))
    Genus<-ifelse(TaxClean0$Genus=="",paste(Family,"FA",sep=""),paste(genus))

    # Rejoin Taxonomy DS elements
    TaxRedo <-data.frame(OTU, Kingdom,Phylum,Class,Order,Family,Genus)                               # length(OTU);length(Kingdom);length(Phylum);length(Class);length(Order);length(Family);length(Genus)

    # Cleanup multiple sequential tags
    TaxRedo <-as.data.frame(sapply(TaxRedo, gsub,pattern="PHCLORFA",replace="PH"))
    TaxRedo <-as.data.frame(sapply(TaxRedo, gsub,pattern="PHCLOR",replace="PH"))
    TaxRedo <-as.data.frame(sapply(TaxRedo, gsub,pattern="PHCL",replace="PH"))
    TaxRedo <-as.data.frame(sapply(TaxRedo, gsub,pattern="CLORFA",replace="CL"))
    TaxRedo <-as.data.frame(sapply(TaxRedo, gsub,pattern="CLOR",replace="CL"))
    TaxRedo <-as.data.frame(sapply(TaxRedo, gsub,pattern="ORFA",replace="OR"))

    # add OTU as row names
    row.names(TaxRedo)<-TaxRedo$OTU                                                                  # TaxRedo$Genus # dim(TaxRedo)  # TaxRedo

    ######################################################
    ## 4) Prune OTU table by cleaned taxonomy ##
    OTUin$OTU <-row.names(OTUin)                                          # Get origninal OTU table (OTUin), add OTU column       #head(OTUin)
    KeepOTU <- data.frame(OTU=TaxRedo$OTU)                                # Gather OTU numbers from cleaned taxonomy
    row.names(KeepOTU)<-KeepOTU$OTU                                       # Add OTU row names for later.                # KeepOTU

    cleanOTU <- merge(OTUin,TaxRedo, by="OTU")                            # Merge KeepOTUs with original OTU table      # keeping taxonomy for check...  # cleanOTU
    row.names(cleanOTU) <- cleanOTU$OTU                                   # Make OTUs row.names after merge
    return(cleanOTU)
}

# Test taxonomy cleaning
# otu_SFC = clean_OTU_taxon_6ranks(otu_S)
# dim(otu_SFC); head(otu_SFC); tail(otu_SFC)

# New function to do all three steps as one

otu_t_preproc = function(otu_t, min_otu, metaDB, Sample, index_col) {
    # Filters otu_t to min_otu count; then reorders Samples by index_col in metaDB; 
    # then fixes taxonomy and outputs 6 rank + lineage table                                 
    row.names(otu_t)<-paste0("otu_",row.names(otu_t))                         # New plaintext OTU ids needed for several transposes below          # added to matrix (for colnames)

    otu_F = filter_n_otu(otu_t, min_otu)                                      # filter min OTU
    otu_SF = sort_samples_by_meta(otu_F, MetaDB, Sample, index_col)           # Sort samples by meta   
    otu_SFC = clean_OTU_taxon_6ranks(otu_SF)                                  # Clean taxonomy
    
    # Create new "Taxonomy" category
    otu_SFC["Taxonomy"] <- ifelse(otu_P['Phylum']=="Proteobacteria", paste(otu_P[,'Class']), paste(otu_P[,'Phylum']))
    
    # Make counts numeric
    otu_counts <- data.frame(data.matrix(otu_SFC[, 2:(ncol(otu_SFC)-8)]))  # get OTU count data, make numeric
    tax <- otu_SFC[, (ncol(otu_SFC)-7):ncol(otu_SFC)]                      # get OTU taxonomy data
    otu_PP_out <- data.frame(otu_counts, tax)                              # join counts, taxonomy 
    return(otu_PP_out)
}

# Test function
# otu_PP = otu_t_preproc(otu_raw, 500, metaDB, "Sample", "EWsiteHyd_index")
#  dim(otu_PP); head(otu_PP)

#head(otu_SFC)

# Test whether DESeq2 works downstream

# DESeq2 normalize data 
# Import DESeq2 MODULE
#source('2_OTU_table_to_DESeq2_and_VST_cpm_0.6.r')

# Merge site order and Samples
# OTU_samps <- data.frame('Sample'=colnames(otu_PP))                #OTU_samps
# Meta_iTag <- merge(metaDB, OTU_samps, by='Sample')               #colnames(metaDB)
# rownames(Meta_iTag) <- Meta_iTag$Sample


# Calculate DESeq2 normalized VST CPM otu table
#OTU_vstCPM = calc_DESeq2_CPM(otu_PP, metaDB)
# OTU_vstCPM 

# OR
#physeq <- make_Phyloseq_data(otu_PP, Meta_iTag, sample = "Sample")
#physeq
# counts <- sapply(otu_PP, is.numeric)                                       # Select numeric columns
# otu_d <- data.matrix(otu_PP[,counts])                                      # Keep only numeric cols, matrix for # phyloseq  #dim(otu_d); head(otu_d)
# OTU <- otu_table(otu_d, taxa_are_rows=TRUE, errorIfNULL = TRUE)           # Make phyloseq OTU table     # OTU
# OTU# counts
#otu_PP = otu_t_preproc(otu_raw, 500, metaDB, "Sample", "EWsiteHyd_index")
#dim(otu_PP); head(otu_PP)

#otu_S = sort_samples_by_meta(otu_F, MetaDB, "Sample", 'EWsiteHyd_index')
#dim(otu_S); head(otu_S)
