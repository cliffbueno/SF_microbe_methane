
##############################################################################
# OTU table to DESeq2 normalized table dev v.0.6

# returns OTU table as VST CPM OR DESeq object for experimental treatment testing

### Input: 
## - 1) OTU table as text      (otu_P)
#    - best if preprocessed using module to correct taxonomy
    
## - 2) Metadata table         (metadata)
#    - assumes "Sample" (to match OTU cols) is first col, plus one additional column min.
#    - can rename sample = "Sample"

### Output(s): 
## - 1) DESeq VST normalized OTU table, as CPM (Counts Per Million) 
## - 2) DESeq object, if treatment design.  
## - 3) Add VST and CPM to 2) above for plotting contrasts
##############################################################################

##############################################################################

### USER Workflows:

# 1) Get OTU DESeq2 VST CPM table from input OTU table 
#   - OTU_vstCPM = calc_DESeq2_CPM(otu_P, metadata)
#     Makes VST table, requires (MetaDB to have a "Sample" first col and at least one more col)

# 2)  Get OTU DESeq2 object for contrasts, Then get VST CPM table w. taxonomy
#    - DEQSeq_treat = get_DESeq2_treat_effs(otu_P, metadata, treat="LU", ref_level="Ref")
#    - OTU_vstCPM = DESeq2_vstCPM(DEQSeq_treat)

### OUTLINE OF FUNCTIONS (U* are composite user functions)

# 1)  Make phyloseq data
#       physeq = Make_Phyloseq_Data(otu_P, metadata, sample = "Sample")

# 2)  Get DESeq2 object w no treatments
#        OTU_phy2des = get_DESeq2(physeq)

# 3U) Get DESeq2 object w. treatments, rawOTU 
#        DEQSeq_treat = get_DESeq2_treat_effs(otu_P, metaDB, treat="LU", ref_level="Ref")

# 4U)  Get DESeq2 VST CPM data from DSq object
#        OTU_vstCPM = DESeq2_vstCPM(DEQSeq_treat)

# 5U) Get DESeq2 VST CPM data from original OTU, metadata
#        OTU_vstCPM = calc_DESeq2_CPM(otu_P, metaDB)

##############################################################################

# Import libraries
library(phyloseq)
suppressMessages(library(DESeq2))



##############################################################################
### Load TESTING data

### IMPORT Sample mapping
# metaDB <-read.table("SF_sal_meta_FIX3.txt", sep="\t", header=TRUE)           # Import Metadata, keep all
#row.names(metaDB) <- metaDB$Sample                                          # Row names are samples for phyloseq             
#metaDB = metaDB[,-1]                                                        # Drop only old index, keep everything else                                  

### import pre-processed table
# otu_P = read.table("SF_Sal_OTU_table_PP.txt", sep='\t', header=T, row.names=1)
# dim(otu_P);head(otu_P)

# Alternately, pre-process table manually 

# Import module from R source:
# source("OTU_preprocessing.R")  

## IMPORT raw OTU TABLE 
# Raw (w. taxon)
# otu_raw <- read.table("SF_Salinity_gradient_OTU_table.txt", sep='\t', header=TRUE, row.names=1) 

## Run OTU table pre-processing step: Filter reads, sort samples, clean taxnonomy
#otu_PP = otu_t_preproc(otu_raw, 500, metaDB, "Sample", "EWsiteHyd_index")              # dim(otu_PP); head(otu_PP)

# Write PP OTU table
#write.table(otu_PP, "SF_Sal_OTU_table_PP.txt", sep='\t', row.names=FALSE, col.names = TRUE, quote=F)

##############################################################################
### 1) Prepare Phyloseq data

make_Phyloseq_data = function(otu_P, metadata, sample = "Sample"){
    
    # OTU table to phyloseq OTU tab
    counts <- sapply(otu_P, is.numeric)                                       # Select numeric columns
    otu_d <- data.matrix(otu_P[,counts])                                      # Keep only numeric cols, matrix for phyloseq  #dim(otu_d); head(otu_d)
    OTU <- otu_table(otu_d, taxa_are_rows=TRUE, errorIfNULL = TRUE)           # Make phyloseq OTU table     # OTU
    
    # Make Phyloseq Taxonomy table, while dropping OTU col
    tax <- sapply(otu_P, is.factor)                                           # Select factor columns
    FullTax <-otu_P[,tax]                                                     # Keep only factor columns 
    TaxTable<-tax_table(FullTax)                                              # make Phyloseq taxonomy object
    colnames(TaxTable) <- colnames(FullTax)                                   # Add back colnames
    row.names(TaxTable) <- row.names(FullTax)                                 # Add back rownames     
    
    # Get Metadata (sample data) into phyloseq
    row.names(metadata) <- metadata[,sample]                                   # Row names are samples for phyloseq             
    metadata = metadata[,-1]                                                  # Drop only old index, keep everything else                                  
    Map <- sample_data(metadata)                                              # Phyloseq import metadata
    
    ### Make PHYLOSEQ object                                                  # suppressMessages() -- not working
    physeq = phyloseq(OTU,Map,TaxTable)  
    print(physeq)
    return(physeq)
}

# Test Function
# physeq = make_Phyloseq_data(otu_P, metaDB)

# Omitted verbose:
# print(physeq)
# print(sample_variables(physeq))                # Check sample variables  

##############################################################################
### 2) Calculate DESeq2 Phyloseq object
# here using NULL design, only intercept for generating table...  
# NOT for hypothesis testing, but for correlations... ?DESeqDataSet

get_DESeq2 = function(physeq, design = ~ 1){
    design = design                                                            # Design = mean (intercept) 
    OTU_phy2des <- phyloseq_to_deseq2(physeq, design)                         # phyloseq out TO DESeq2 
    OTU_phy2des <- DESeq(OTU_phy2des, test="Wald", fitType="parametric")      # Generate dispersions, VST 

    print('Deseq2 finished computing')
    return(OTU_phy2des)
}

# Test function
#OTU_phy2des = get_DESeq2(physeq)
#OTU_phy2des

##############################################################################
### 3) Calculate DESeq2 Phyloseq object
# for block design with reference treatment

get_DESeq2_treat_effs = function(otu_P, metadata, treat, ref_level){
    physeq = make_Phyloseq_data(otu_P, metaDB)
    sample_data(physeq)$Treat <- relevel(as.factor(get_variable(physeq, treat)), ref=ref_level)
    OTU_phy2des = get_DESeq2(physeq, design = ~ Treat)
    return(OTU_phy2des)
}

# Test function on effects with treatment "LU" and reference level "Ref" -- like a control
# DEQSeq_treat = get_DESeq2_treat_effs(otu_P, metaDB, treat="LU", ref_level="Ref")
# summary(DEQSeq_treat)
# DEQSeq_treat

##############################################################################
### 4) Calculate DESeq2 VST CPM with Taxonomy FROM DESeq2 object

DESeq2_vstCPM = function(OTU_phy2des, physeq){    
    ## REPLACE VST data with COUNTS                                           # validated results same as recalc VST 
    OTU_vs_counts <- counts(OTU_phy2des, normalized=TRUE)                     # get counts                                  
    OTU_table_VST <- data.frame(OTU_vs_counts)                                # make df, shorten name

    # Scale OTU VST table counts to CPM (Counts per million)
    sampTots<-colSums(OTU_table_VST)                                          # Get sample totals                               
    OTU_VST_CPM<-sweep(OTU_table_VST*1000000,2,sampTots,'/')                  # Sweep matrix by:[div. by samp total * 1e+06]    
    OTU_VST_CPM$OTU <- row.names(OTU_VST_CPM)                                 # Add OTU as column for later merges              
    # colSums(OTU_VST_CPM[,1:(ncol(OTU_VST_CPM)-1)],)                         # check sums                                      
    
    ### Add Taxonomy data                                        # source: #head(tax_table(physeq))
    Taxonomy <- data.frame(tax_table(physeq))                    # make DF from phyloseq object 
    Taxonomy$OTU <- row.names(Taxonomy)                          # Add OTU number for merge                        # head(Taxonomy)
    OTU_VST_CPM_tax <- merge(OTU_VST_CPM, Taxonomy, by="OTU")
    
    #print('Deseq2 VST CPM finished computing')
    return(OTU_VST_CPM_tax)
 }    

# Test function
# OTU_vstCPM = DESeq2_vstCPM(DEQSeq_treat)
# OTU_vstCPM = DESeq2_vstCPM(OTU_phy2des)
# head(OTU_vstCPM)

##############################################################################
### 5) Combined: Calculate DESeq2 VST CPM from raw OTU, metadata

calc_DESeq2_CPM = function(otu_P, metadata, sample = "Sample", design = ~ 1){                                  
    physeq = make_Phyloseq_data(otu_P, metadata, sample)                              # to phyloseq format
    OTU_phy2des = get_DESeq2(physeq, design)
    OTU_vstCPM = DESeq2_vstCPM(OTU_phy2des, physeq)                                           # DESeq2 VST CPM calculation
    return(OTU_vstCPM)
}

# test function
# OTU_vstCPM = calc_DESeq2_CPM(otu_P, metaDB)
# head(OTU_vstCPM)
