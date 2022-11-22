# SEM - Lasso semi-auto for SF OTU guilds v1.12
#-  runs under r > 3.6 (glmnet for LASSO)

#- v1.12 clean new functions & messes in v1.11
#- v1.11 testing models with primitive workflows
#- v1.11 some functions for branched models, can always be more
#- v1.11 need to mute branched testing in OLD sections 14, is like 15 but messier
#- v1.11 /some cleanup needed (eg. wft is section II) 
#- v1.11 MOVED test data for composite functions, now in SEM_models/fxn_example_tables/
#  - v1.10 goal: clean fxns and examples to actual fleshed out testing

########################################################################################
# 0) Import libraries, modules

# Graphics packagaes
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

suppressMessages(library(Hmisc))
suppressMessages(library(cowplot))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))
theme_set(theme_grey())             

#suppressMessages(library("ggdendro")) -- note these were moved to bottom as don't work at top, work fine at bottom
#suppressMessages(library('dendextend')) -- problem in dendextend branch clustering k = n, not working, assume overwritten

# utilities
library(parallel)
library(stringr) 
library(reshape2)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(vegan))

# load SEM libraries and functions
suppressMessages(library("lavaan"))
suppressMessages(library(psychTools))
library(AICcmodavg)

# library("semPlot")
# library("DiagrammeR") - seems ok.

suppressMessages(library(devtools))

# install_github('SachaEpskamp/semPlot',  dependencies = F)  - no, dep = T may work or break lib...WAIT
# library(semPlot)

#library("semPlot")
#library("DiagrammeR") 

# load caret and glmnet (r > 3.6)
suppressMessages(library(glmnet)) #install.packages('glmnet')
suppressMessages(library(caret))

# Required packages (move up to top)
#library(ape)
suppressMessages(library("ggdendro"))
suppressMessages(library('dendextend'))

# Import own modules for OTU table processing into guilds

source("Import_SalOTU_dat_Plot_test_v0.1.R")                             # load data and metadata
source("OTU_subsetting_modules_v.0.2_strip.R")                           # gets guilds, contains colors 
suppressMessages(source("Corr_ranks_module_v0.3.2_strip.R"))             # contains simple aggregation function

########################################################################################
### 1) Import data sets

## a) Import metadata -- Read data fom import module for testing convenience here
# source("Import_SalOTU_dat_Plot_test_v0.1.R")
# head(Meta_iTag)

## b) Imported OTU table
# head(otu_V)

## c) Separate Delta Sites (oligo and FW)                                 # levels(Meta_iTag$SALTgroup)
Meta_iTag_FW <- Meta_iTag[Meta_iTag$SALTgroup =="FW",]
Meta_iTag_Oligo <- Meta_iTag[Meta_iTag$SALTgroup =="Oligo",]

Meta_iTag_Delta <-rbind(Meta_iTag_FW, Meta_iTag_Oligo)                    # Meta_iTag_Delta
Delta_sites <-data.frame(Meta_iTag_Delta[,"Sample"])                      # DF of Delta sites
colnames(Delta_sites) [1] <- "Sample"                                     # Rename column "Sample"
#Delta_sites["Samp_Index"] <- seq(1:nrow(Delta_sites))                    # Make sample index for reordering post-merges

# Delta_sites


########################################################################################
### 2) Process guild data and metadata

## 0) import guild counting from OTU table function from module 
# source("OTU_subsetting_modules_v.0.2_strip.R")
# options(repr.plot.width=1.5, repr.plot.height=6) 

## a) get guild counts/sample data from OTU table
Guild_OTUs <- Get_16S_Guilds(otu_V)              # use Get_16S_Guilds to get guilds  
otu_V["OTU"] <- row.names(otu_V)                 # Make OTU number column
# dim(Guild_OTUs); head(Guild_OTUs)

# merge otu table and Guilds 
OTU_guilds <- merge(Guild_OTUs, otu_V, by="OTU")#, all.y=TRUE)
head(OTU_guilds)

## b) Get C guilds from Phyla (based on shotgun data for C)
Firmic <- otu_V[subset(otu_V["Phylum"]=="Firmicutes"),]     # Firmicutes
Firmic['Guild'] <-"Firmic" 

Actino <- otu_V[subset(otu_V["Phylum"]=="Actinobacteria"),] 
Actino['Guild'] <-"Actino" 

Chlorf <- otu_V[subset(otu_V["Phylum"]=="Chloroflexi"),] 
Chlorf['Guild'] <-"Chlorf" 

C_guilds <- rbind(Firmic, Actino, Chlorf)

## c) Process guild data and metadata 
# Get agg_by_cat aggregatopm function from corr module
# source("Corr_ranks_module_v0.3.2_strip.R")

# Aggregate OTU guilds, C guilds & combine  
Guild_agg <- agg_by_cat(OTU_guilds, "Guild")     # aggregate by Guild
C_guild_agg <- agg_by_cat(C_guilds, "Guild")
Guild_agg2 <-rbind(Guild_agg, C_guild_agg)  # Combine OTU and C guild data

## d) calculate new guild ratios, log transform guild data
# prepare data
Guild_agg2 <- data.frame(t(Guild_agg2[,-1]))           # Transpose for new calcs
Guild_agg2 <- replace(Guild_agg2, Guild_agg2==0, 2)    # replace 0s with pseudo counts 

# Add guild ratios 
Guild_agg2$AO_NOB <- (Guild_agg2$AOA + Guild_agg2$AOB) / Guild_agg2$NOB
Guild_agg2$NOB_AO <- Guild_agg2$NOB / (Guild_agg2$AOA + Guild_agg2$AOB)
Guild_agg2$mcr_pmo <- ((Guild_agg2$CH4_ac + Guild_agg2$CH4_H2) / 
                            (Guild_agg2$MOB_I + Guild_agg2$MOB_II + Guild_agg2$MOB_IIa))
Guild_agg2$pmo_mcr <- ((Guild_agg2$MOB_I + Guild_agg2$MOB_II + Guild_agg2$MOB_IIa) / 
                       (Guild_agg2$CH4_ac + Guild_agg2$CH4_H2))

# Log 2 transform data
Guild_agg_L2 <- log2(Guild_agg2)                       # log2
Guild_agg_L2$Sample <- row.names(Guild_agg_L2)
#head(Guild_agg_L2)

#colMeans(Guild_aggT)
#Guild_aggT1000 = Guild_aggT/1000   # Divide OTU counts by 1000 for lavaan
#colMeans(Guild_aggT1000)

## e) get guild colors
guild_colors <- read.table("Guild_color_palette.txt", sep ='\t', header = T)
names(guild_colors)[1] <- 'var'
#guild_colors

## f) Process CH4 and other soil metadata
### Reimport metadata instead of using module loaded metadata
metaDB <- read.table("SF_sal_meta_FIX3b_gap_fill_MDS.txt", sep='\t', header=T)
# metaDB <- read.table("SF_sal_metaLOG_FIX2.txt", sep='\t', header=T)
# head(metaDB); names(metaDB)

# Get only numeric data, LOG 10 transform soil chemistry
metaCHEM <- metaDB[,15:ncol(metaDB)]                                             # should use lapply
metaCHEM_log <-log10(metaCHEM)                                                   # log10 chem
CH4_logn1 <- log10(metaDB[,"CH4_ug_m2_h"]- (min(metaDB[,"CH4_ug_m2_h"])*1.05))   # CH4 logn1
metaCHEM_log["CH4_logn1"] <- CH4_logn1                                           # to metachemLOg       #metaCHEM_log

# reattach non-numeric
metaDB = data.frame(metaDB[,1:14], metaCHEM_log)                                  # metaDB

# Merge site order and Samples
Meta_iTag <- merge(metaDB, OTU_samps, by='Sample')                                # colnames(metaDB)
rownames(Meta_iTag) <- Meta_iTag$Sample

# Reorder location factor
Meta_iTag$Location <-factor(Meta_iTag$Location, levels=c("Sandmound","WestPond","Mayberry","Browns","RushRanch","Joice","Goodyear","WhiteSlough","Tolay","ChinaCamp","Muzzi"))  #head(Meta_iTag)
Meta_iTag$Pl_Sp <-factor(Meta_iTag$Pl_Sp, levels=c("Cattail","Tule","ThreeSq","CattailNL","Phrag","PW","Cord"))

# Resort meta itag by index
indexer = 'EWsiteHyd_index'
Meta_iTag <- Meta_iTag[order(Meta_iTag[indexer]),]
colnames(Meta_iTag) #head(Meta_iTag) #max(Meta_iTag$CH4_ug_m2_h) #plot(x=Meta_iTag$Sample, y=Meta_iTag$CH4_ug_m2_h)

# reduce metadata factors
CH4_samp <- c("Sample", "CH4_ug_m2_h","CH4_logn1",'CO2_mg_m2_h','CO2_soilC_mg_g_d','CH4_CO2', "Bulk_dens", 'H2O_FPS',
               'pH', 'C','N','P','CN','NP', 'NO3_N', 'NH4_N', 'Olsen_P', 'NO3_NH4','NP_ext',
               'Salinity.x', 'Cl', 'SO4', 'SO4_pw', 'Fe', 'Fe_pw','DOC_mg_L')
              
              # 'NO2_pw','NO3_pw', 'NH3_pw', 
              # "C_g_m2")  # not for log data, already been logn1 transf: "CH4_logn1", 
CH4 <- Meta_iTag[,CH4_samp]  
# CH4

## g Merge metadata with guild data 
Guild_CH4 <- merge(CH4, Guild_agg_L2)
# dim(Guild_CH4); head(Guild_CH4)

# Clean NA methane rows from this 
Guild_CH4 <- Guild_CH4[!is.na(Guild_CH4$CH4_ug_m2_h),]           

# Clean NA data bad metadata from this
guild_names <- names(Guild_CH4)
#drop <-c("NO2_pw", "NO3_pw", "NH3_pw", "NH4_N.1", "NO3_N.1")
#keep <- guild_names[!guild_names %in% drop]                      #keep
#Guild_CH4 <- Guild_CH4[keep]

# dim(Guild_CH4); head(Guild_CH4); names(Guild_CH4)

# ggplot(Guild_CH4, aes(x=CH4_ug_m2_h, y=CH4_logn1)) + geom_point()

##  h) Get delta only sites
Guild_CH4_d <- merge(Guild_CH4, Delta_sites)
dim(Guild_CH4_d); head(Guild_CH4_d); colnames(Guild_CH4_d)

## i) Write out convenience tables for later use 
# dir.create("SEM_data")

# write.table(Guild_CH4_d, "SEM_data/SEM_base_log2_10_delta_guild_soil_data.txt")
# write.table(Guild_CH4, "SEM_data/SEM_base_log2_10_all_guild_soil_data.txt")

########################################################################################
### 3) SEM functions for extracting fits, MI covariates

# "Simple model" - 
meth_mod <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h 
             
             CH4_H2 ~ SRB_syn + CO2_mg_m2_h
             MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac
             MOB_I ~ CH4_H2 +CH4_ac'

# Calculate SEM
meth_mod.sem <- sem(meth_mod, data=Guild_CH4_d, fixed.x=FALSE, estimator ="mlm")    # Generates observed var > 1000 x warning, divide otu vas try again

sem_fit_params <- c( "pvalue", "chisq", "df", "npar","aic", "bic", "gfi", "cfi", "rni", "rmsea", "srmr")

# How to get features from Lavaan output, used in functions
# parameterEstimates(meth_mod.sem, standardized=TRUE)  # gets P value for each Regression in summary
# inspect(meth_mod.sem)                                # gets various SEM model params 
# inspect(meth_mod.sem, 'r2')                          # gets R2 for each fitted var

# fitMeasures(meth_mod.sem)                            # gets all model fit measures
# fitMeasures(meth_mod.sem, c("aic", "bic"))           # gets selected model fit measures

# check for SEM variable significance

get_SEM_nonsig_vars = function(sem_mod, p_cut) {

    # get significant only 
    sem_params <- parameterEstimates(sem_mod)                          # Get estimates
    sem_params_bad <- sem_params[sem_params$pvalue > p_cut,]           # cut out significant

    # clean up
    keep <- c('lhs', 'op','rhs','pvalue')                              # params to keep
    sem_params_bad <- sem_params_bad[keep]                             # keep
    sem_params_bad$pvalue <- round(sem_params_bad$pvalue, 2)           #  round p value
    sem_params_bad <- sem_params_bad[complete.cases(sem_params_bad),]  # drop all NA values (comes up in latent mods)
    sem_params_bad <- sem_params_bad[!sem_params_bad$op == "~~",]       # drop covariate fails 
    
    return(sem_params_bad)
    
}

# Demonstrate function
get_SEM_nonsig_vars(meth_mod.sem, p_cut =0.05)

get_SEM_R2s = function(sem_mod){
    
    sem_name <- paste(as.character(deparse(substitute(sem_mod))))
    model <- data.frame(model = sem_name)
    
    R2s <- inspect(sem_mod, 'r2')
    R2s <- data.frame(round(R2s, 3))
    out <- t(R2s)
    
    out <- data.frame(model, out)
    row.names(out) <-NULL
    return(out)
}

# demonstrate function
get_SEM_R2s(meth_mod.sem)


# Define list of fit params
#sem_fit_params <- c("npar", "chisq", "df", "pvalue", "aic", "bic", "gfi", "cfi", "rni", "rmsea", "srmr")
sem_fit_params <- c( "pvalue", "chisq", "df", "npar","aic", "bic", "gfi", "cfi", "rni", "rmsea", "srmr")

# Make function
sem_fit_sum = function(sem_mod, sem_fit_params) {
    
    sem_name <- paste(as.character(deparse(substitute(sem_mod))))
    model <- data.frame(model = sem_name)
    
    fits <- data.frame(fitMeasures(sem_mod, sem_fit_params))
    fits <- round(fits, digits = 3)
    out <- t(fits)
    
    out <- data.frame(model, out)
    row.names(out) <-NULL
   # names(out)[1] <- sem_name
    return(out)
    
    
}

# Test function
sem_fit_sum(meth_mod.sem, sem_fit_params)

# models = c(meth_mod.sem, meth_mod.sem)
meth_mod.sem2 <- meth_mod.sem
models2 = c("meth_mod.sem","meth_mod.sem2")

compare_sem_fits = function(models, sem_fit_params){
    
    out_list <- {}
    
    for (i in seq_along(models)) {

        curr_mod <- get(models[[i]])
        out_list[[i]] <- sem_fit_sum(curr_mod, sem_fit_params)

        }

    out <- do.call("rbind", out_list)
    out$model <- models
    return(out)
    
}

compare_sem_fits(models2, sem_fit_params)

compare_sem_R2s = function(models){
    
    out_list <- {}
    
    for (i in seq_along(models)) {

        curr_mod <- get(models[[i]])
        R2c <- get_SEM_R2s(curr_mod)
        R2c$model <- models[i]
        out_list[[i]] <- R2c
        #out_list[[i]] <- get_SEM_R2s(curr_mod)
        #names(out_list[i]) <- models[i]
        }
    
    out <- Reduce(function(...) merge(..., all=T), out_list)
    
    #out <- Reduce(merge, out_list)
    #out <- do.call("rbind", out_list)
    #out$model <- models
    #out_list
    out
}

# Demonstrate function
#models2 <- c('ch4_mod0_f','ch4_mod0a_f', 'ch4_mod0b_f')
compare_sem_R2s(models2)  # note this doesn't look right with two of the same models names

compare_sem_results = function(models, sem_fit_params){
    
        fits <- compare_sem_fits(models, sem_fit_params)
        R2s <- compare_sem_R2s(models)
    
        out <- data.frame(merge(fits, R2s, by="model"))
        return(out)

}

models2
#meth_mod.sem

# Demonstrate function
compare_sem_results(models2, sem_fit_params)
# testing issues way downstream -- non-matching outputs with running auto composites

man_mod_7 <- c('sem.ch4L0C.107','sem.ch4L0C.113','sem.ch4L0C.116','sem.ch4L0C.118',
              'sem.ch4L0C.121','sem.ch4L0C.69','sem.ch4L0C.79')

auto_mod_7m <- c('sem.107','sem.113','sem.116','sem.118','sem.121','sem.69','sem.79')
auto_mod_7 <- c('sem.8','sem.10','sem.16','sem.34','sem.107','sem.113','sem.116','sem.118','sem.121','sem.69','sem.75','sem.79')

compare_sem_results(man_mod_7, sem_fit_params)
compare_sem_results(auto_mod_7m, sem_fit_params)
compare_sem_results(auto_mod_7, sem_fit_params)


sem_mi_table <- function(sem_mod, mi_lower = 4, mi_sort = "T", head = "T"){
    
    mi_tab <- data.frame(modindices(sem_mod))  # get mod indices table
    mi_tab <- mi_tab[mi_tab$mi > mi_lower,]         # impose cutoff on lower MI 
    
    # sort by mi?  or defaults to covariances first
    if (mi_sort == "T") {out_tab <- mi_tab[order(-mi_tab$mi),]}
        else {out_tab <- mi_tab}
    
    # crop to top 5 most important effects?
    if (head == "T") {out <- head(out_tab)}
        else {out <- out_tab}
    
    return(out)
        
}

# demonstrate function
sem_mi_table(meth_mod.sem)
#sem_mi_table(meth_mod.sem, mi_lower = 4, mi_sort = "F", head = "T")
# set output n results
n_res = 4

# Get covariance matrix, melt
resid <- residuals(meth_mod.sem)$cov
resid[lower.tri(resid, diag = T)] <-0         # get upper triangle only, lower as 0
resid_M <- melt(resid)                        # melt DF

# Clean, sort table
resid_M <- resid_M[!resid_M$value == 0,]      # Drop 0/lower
resid_M$abs <- abs(resid_M$value)             # Get abs for sorting effect size 
resid_M <- resid_M[order(-resid_M$abs),]      # sort by abs
resid_M <- resid_M[,1:3]                      # drop abs col 
resid_M <- resid_M[1:n_res,]                  # get first n results
resid_M# for fun, how many unique variables in there?
cov_uniq <-  unique(c((paste(resid_M$Var1)), c(paste(resid_M$Var2))))
cov_uniq
get_sem_covariances = function(sem_mod, n_res =5){  # output_factors = F - consider output as list not table?
    
    # Get covariance matrix, melt
    resid <- residuals(sem_mod)$cov
    resid[lower.tri(resid, diag = T)] <-0         # get upper triangle only, lower as 0
    resid_M <- melt(resid)                        # melt DF

    # Clean, sort table
    resid_M <- resid_M[!resid_M$value == 0,]      # Drop 0/lower
    resid_M$abs <- abs(resid_M$value)             # Get abs for sorting effect size 
    resid_M <- resid_M[order(-resid_M$abs),]      # sort by abs
    resid_M <- resid_M[,1:3]                      # drop abs col 
    resid_M <- resid_M[1:n_res,]                  # get first n results
    resid_M
              
} 
    


# Demonstrate function
get_sem_covariances(meth_mod.sem, n_res =5)

plot_matrix <- function(matrix_toplot){
corrplot(matrix_toplot, is.corr = FALSE, 
               type = 'lower', 
               order = "original", 
               tl.col='black', tl.cex=.75)
}


options(repr.plot.width=4, repr.plot.height=4) 
plot_matrix(residuals(meth_mod.sem)$cov)

# residuals(meth_mod.sem)$cov
# residuals(meth_mod.sem, type = "cor")$cor   # "keep an eye out for residual corr > 0.1" from http://www.understandingdata.net/2017/03/22/cfa-in-lavaan/
## E) Suggest model updates with external covariates (based on MI and model covariates) 
# prepare data matrix from input
data <- Guild_CH4_d[,-1]
data <- data[!is.na(data$CH4_ug_m2_h),]
data <-data.matrix(data)

# function to get covariates from MI pair
get_mi_covars_pair = function(data, r_cut = 0.5, keep_pair){
    
    corr_m <- cor(data)                                             # get corr
    corr_m[abs(corr_m) < r_cut] <- NA                               # get abs value over threshold only, or NA

    select_cols <- data.frame(corr_m[,keep_pair])                   # Get only cols from MI inputs
    # select_cols <- data.frame(corr_m[,keep_vars])                 -- too many vars at once not helping
    common_data <- t(select_cols[complete.cases(select_cols),])        # get only complete cases -- no NA
    common_data <- round(common_data,2)
    #select_cols
    common_data
}

# Note behavior change with many variables, only common feats kept
# keep_pair = c("CH4_H2", "SRB_syn")
keep_pair = c("CH4_H2", "AOA")
# keep_pair = cov_uniq[1:4]   

# Demonstrate function
get_mi_covars_pair(data, r_cut = 0.6, keep_pair)

# function to get covariates from MI several: drop all NA only
get_mi_covars_many = function(data, r_cut = 0.5, keep_vars){
    
    corr_m <- cor(data)                                             # get corr
    corr_m[abs(corr_m) < r_cut] <- NA                               # get abs value over threshold only, or NA

    # select_cols <- data.frame(corr_m[,keep_pair])                 # Get only cols from MI inputs
    select_cols <- data.frame(corr_m[,keep_vars])                   #  too many vars at once not helping
    # common_data <- select_cols[complete.cases(select_cols),]      # get only complete cases -- no NA
    common_data <- t(select_cols[rowSums(is.na(select_cols)) != ncol(select_cols), ])
    common_data <- round(common_data,2)
    common_data
}

keep_vars <- c("CH4_H2", "MOB_I", "AOA", "NOB", "AO_NOB")
#keep_vars <-cov_uniq

# Demonstrate function
get_mi_covars_many(data, r_cut = 0.6, keep_vars)
# Define simple model" - 
ch4_mod0 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h 
             
             CH4_H2 ~ SRB_syn + CO2_mg_m2_h
             MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac
             MOB_I ~ CH4_H2 +CH4_ac'
# try adding SEM response vars with string defs
CH4 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h'
CH4_H2 <- 'CH4_H2~ SRB_syn + CO2_mg_m2_h'
MOB_IIa <- 'MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac'
MOB_I <- 'MOB_I ~ CH4_H2 +CH4_ac'

# Get vector of response vars as strings
Other <- ""
response_vars <- c("CH4", "CH4_H2", "MOB_I","MOB_IIa", "Other")

define_new_sem_table <- function(mod_name, response_vars, output_prefix = "sem."){

    # Build data frame of sub-formulas from names
    feat_list <- {}
    
    for (i in seq_along(response_vars)) {

        feat_list[[i]] <- get(response_vars[i])

    }
    
    feats <- data.frame(feat_list)
    row.names(feats) <- response_vars
    feats <-data.frame(t(feats))
              
    # Get model formula from feats
    formula <- do.call(paste, c(feats, sep = '\n'))
    feats$formula <- formula
     
    # add model info
    model_out <- paste0(output_prefix, mod_name)
    mod_info <- data.frame(model_out = model_out, model = mod_name, base = mod_name, notes = "base")
    
    out <- data.frame(mod_info, feats)
    row.names(out) <- NULL
    return(out)
        
}

# Demonstrate function
model_test <- define_new_sem_table(mod_name = "ch4_mod0", response_vars)
model_test

# Inspect formula
# cat(model_test$formula)    #model_test$formula

# check formula extraction
formula_test <- model_test$formula
# cat(formula_test)

# Confirm this works in SEM model
ch4_mod0_build_f <- sem(formula_test, data=Guild_CH4_d, fixed.x=FALSE, estimator ="mlm")    # Generates observed var > 1000 x warning, divide otu vas try again
sem_fit_sum(ch4_mod0_build_f, sem_fit_params)
# and paste with newline
# ch4_mod0_build <- paste(CH4, CH4_H2, MOB_IIa, MOB_I, sep = '\n')
# cat(ch4_mod0_build)  - great to visualize / inspect, fails lavaan# seems to work with component list
# response_vars <- c(CH4, CH4_H2, MOB_IIa, MOB_I)
# ch4_mod0_build <- paste(response_vars, sep = '\n')
# cat(ch4_mod0_build)# Confirm this works
# ch4_mod0_build_f <- sem(ch4_mod0_build, data=Guild_CH4_d, fixed.x=FALSE, estimator ="mlm")    # Generates observed var > 1000 x warning, divide otu vas try again
# sem_fit_sum(ch4_mod0_build_f, sem_fit_params)
update_sem_model_table = function(mod_table, base_model, new_mod_name, notes, formula_edits, 
                                  output_prefix = "sem.", 
                                  non_feat_cols = c("model_out","model", "base", "notes", "formula")){
    
     # copy base model & update descriptions
     new_model <- mod_table[mod_table$model == base_model,]
     new_model$base <- base_model
     new_model$model <- new_mod_name
     new_model$model_out <- paste0(output_prefix, new_mod_name) 
     new_model$notes <- notes
    
     # edit formula components
     form_feat <- formula_edits[1]
     form_edit <- formula_edits[2]
    
     new_model[form_feat] <- form_edit
    
     # rebuild full formula
     cols <- names(new_model)
     feat_cols <- cols[!cols %in% non_feat_cols]
     feats <- data.frame(new_model[feat_cols])
                         
     formula <- do.call(paste, c(feats, sep = '\n'))
     new_model$formula <- formula
    
     # append to original table 
     mod_table_updated <- rbind(mod_table, new_model)
                 
     mod_table_updated
     #new_model
}

# Test function
model_test2 <- update_sem_model_table(model_test, base_model = "ch4_mod0", new_mod_name = "ch4_mod0a.1", 
                      notes = "drop CH4~MOB_I",
                      formula_edits = c("CH4", "CH4_ug_m2_h ~ CH4_ac + CH4_H2 + MOB_IIa + CO2_mg_m2_h"))
                      # output_prefix = "sem."#non_feat_cols = c("model", "base", "notes", "formula"))

model_test2
# Test while dropping entire internal feature
model_test3 <- update_sem_model_table(model_test2, base_model = "ch4_mod0a.1", new_mod_name = "ch4_mod0a.2", 
                      notes = "drop MOB_I ~ .",
                      formula_edits = c("MOB_I", ""))

model_test3
dim(model_test2)
#nrow(sem_model_table)
# sem_model_table$formula

sem_table_run_assign = function(sem_mod_table, data, estimator = 'mlm', prefix = "sem."){

        iters <- nrow(sem_mod_table)                                                  # seq_along misbehaving, use iters for loop
    
        for(i in seq(1:iters)) {
         
            # run SEM model            
            curr_model <- sem_mod_table$formula[i]                                    # get ith sem formula
            sem <- sem(curr_model, data=data, fixed.x=FALSE, estimator = estimator, orthogonal =TRUE)   # run SEM; fixed.x not working in fxn
            
            # assign mod_out names to models (GLOBAL ENV)
            sem_name <- paste(sem_mod_table$model_out[i])                             # get model name from model_out
            assign(sem_name, sem, envir=.GlobalEnv)                                   # assign to global env 
        }
        
        # Gather models run and print list to console
        models_run <- paste(sem_mod_table$model_out, collapse =', ')
        print(paste0("ran SEM to create models: ", models_run))
        
        models_run_out_list <- c(paste(sem_mod_table$model_out))
        models_run_out_list
    
    }
    
# Note, could have run through output as list, but not necessary as only assigning vars here  # mods <- {}

# test run and assign models from model table
sem_multi_test <- sem_table_run_assign(model_test2, data=Guild_CH4_d, estimator ="mlm")
sem_multi_test
#cat(x)
#data.frame(x)

# test model summary
compare_sem_results(sem_multi_test, sem_fit_params)


run_compare_sem_models = function(mod_table, data, estimator = "mlm", 
                                  keep_mod_descr = c("model_out","model", "base", "notes"),
                                  sem_fit_params = c("pvalue", "chisq", "df", "npar","aic", "bic", 
                                                      "gfi", "cfi", "rni", "rmsea", "srmr")){
    
    # note mod_table requires columns: model, model_out, formula cols 
    
    # get model table features to include in summary
    mod_descr <- data.frame(mod_table[keep_mod_descr])
    mod_descr
    
    # run assign SEM models in table 
    sem_out_list <- sem_table_run_assign(mod_table, data, estimator)
    #sem_out_list

    # Get SEM results summary as previously
    fit_summary <- data.frame(compare_sem_results(sem_out_list, sem_fit_params))
    
    
    # Clean data summary
    names(fit_summary)[1] <- "model_out"
    fit_summary$chisq <- round(fit_summary$chisq, digits = 1)
    fit_summary$aic <- round(fit_summary$aic, digits = 1)
    fit_summary$bic <- round(fit_summary$bic, digits = 1)
    #fit_summary
    
    # combine model descriptions and fits
    out <- merge(mod_descr, fit_summary, by = "model_out")
    out
    #mod_descr
    #fit_summary
}

# Demonstrate function
out <- run_compare_sem_models(model_test2, data, estimator = "mlm")
out



clean_SEM_models_for_output = function(models, models_out){
    
    ### Clean up models ####  
    models <- models[, -ncol(models)]                           # drop formula column  
    models <- models[, -1]                                      # drop model out columns

    # Hide predictors in formulas
    models$CH4 <- gsub("  ", " ", models$CH4)                   # Note there are double spaces in some, likely result of replace ~.
    models$CH4 <- gsub("CH4_ug_m2_h ~ ", "", models$CH4)        # hide predictor

    models$CH4_H2 <- gsub("  ", " ", models$CH4_H2)             # remove double spaces
    models$CH4_H2 <- gsub("CH4_H2~ ", "", models$CH4_H2)        # hide predictor     

    models$MOB_I <- gsub("MOB_I ~ ", "", models$MOB_I)        

    models$MOB_IIa <- gsub("  ", " ", models$MOB_IIa)           # remove double spaces
    models$MOB_IIa <- gsub("MOB_IIa ~ ", " ", models$MOB_IIa)   # hide predictor          
    # return(models)
    
    # Merge with results
    results_tab <- merge(models_out, models, by = c('model', 'base', 'notes'))
    return(results_tab)
}


# Demonstrate function
# Exper1_results <- clean_SEM_models_for_output(Ex1_sem_model_table, Ex1_sem_out)


########################################################################################
### 4) Simple SEM model testing -- example workflow


# Define simple model" - 
ch4_mod0 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h 
             
             CH4_H2 ~ SRB_syn + CO2_mg_m2_h
             MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac
             MOB_I ~ CH4_H2 +CH4_ac'
# Instead of full function, add SEM response vars with string defs
CH4 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h'
CH4_H2 <- 'CH4_H2~ SRB_syn + CO2_mg_m2_h'
MOB_IIa <- 'MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac'
MOB_I <- 'MOB_I ~ CH4_H2 +CH4_ac'

# Get vector of response vars as strings
Other <- ""
response_vars <- c("CH4", "CH4_H2", "MOB_I","MOB_IIa", "Other")

# Build model
sem_ch4_0_base_table <- define_new_sem_table(mod_name = "ch4_mod0", response_vars)
sem_ch4_0_base_table

ch4_mod0 <- sem_ch4_0_base_table$formula

# run simple model
ch4_mod0_f <- sem(ch4_mod0, data=Guild_CH4_d, fixed.x=FALSE, estimator ="mlm")

#compare_sem_results(ch4_mod0_f, sem_fit_params)

# Get summary of fit
sem_fit_sum(ch4_mod0_f, sem_fit_params)

# Get R2: 
get_SEM_R2s(ch4_mod0_f)

# Get non-significant vars
get_SEM_nonsig_vars(ch4_mod0_f, p_cut =0.05)



########################################################################################
### 5) ML feature selection functions for extracting fits, MI covariates

# SCALE DATA, careful will scale y
# Guild_CH4_dS <- scale(Guild_CH4_d[,-1])

# Remove sample names
Guild_CH4_d0 <- Guild_CH4_d
row.names(Guild_CH4_d0) <- Guild_CH4_d0$Sample
Guild_CH4_d0 <- Guild_CH4_d0[,-1]

# names(Guild_CH4_d)

# prepare data subsets

# get vars
xvars0 <-names(Guild_CH4_d0)

# vars to exclude (maybe)
drop_gas <- c('CH4_ug_m2_h', 'CH4_logn1','CH4_CO2', 'CO2_soilC_mg_g_d') # keep this: ,'CO2_mg_m2_h',
drop_ratios <- c('CN','CP','NP','NP_ext')
drop_pw <- c('SO4_pw', 'Fe_pw','Fe','DOC_mg_L')
drop_guilds <- c('SOxB','MeOB','FeOB','AO_NOB','NOB_AO','mcr_pmo', 'pmo_mcr','Anamx')
drop_taxa <- c('Actino','Chlorf','Firmic')

# Get combined drop lists
drop_vars1 <- c(drop_gas, drop_pw, drop_guilds, drop_taxa, drop_ratios)

# Select guild & chem data:                # Remove drop vars from list
xvars1 <- xvars0[!xvars0 %in% drop_vars1]

# alternately, may want only guild or chem data lists, including ALL

# make glmnet data
data <- Guild_CH4_d0
y_var <- 'CH4_logn1'
x_vars <- xvars1

# SCALE ALL, needed for some functions (SEM screening)
Guild_CH4_dS <- scale(Guild_CH4_d[,-1])              #  careful will scale y 
Guild_CH4_dS_df <- data.frame(Guild_CH4_dS)          #  make df
Guild_CH4_dS_df$CH4_logn1 <- Guild_CH4_d0$CH4_logn1  # fix CH4 scale to original

# scale data, get separate x & y -- not needed, moved into FUNC
# Get SCALED X vars
# CH4_pred1 <- Guild_CH4_d0[xvars1]
# CH4_predS1 <- scale(CH4_pred1)
# names(CH4_pred1)# ; head(CH4_predS1)

# get Y var (CH4_flux)
# CH4_flux <- Guild_CH4_d0$CH4_logn1

# Revised function with force in
glmnet_lassoCV = function(y_var, x_vars, data, force_in = '', max_pred = 13){  #, max_pred = 
    
    # Get x & y data
    y <- data[,y_var]
    x <- data.matrix(data[x_vars])
    xS <- scale(x)                                       # lasso scales X by default?  results differ,
    
    # Get force in vars
    force_in_test <- x_vars %in% force_in
    force_in_vec <- ifelse(force_in_test == FALSE, 1, 0)

    # cross val for best lambda
    set.seed(123) 
    cv <- cv.glmnet(xS, y, alpha = 1, penalty.factor = force_in_vec, dfmax = max_pred)
    print(paste0("best_lambda: ", round(cv$lambda.min, 4)))
    
    # lasso w best lambda
    pfit <- glmnet(xS, y, alpha = 1, lambda = cv$lambda.min, dfmax = max_pred,
                  penalty.factor = force_in_vec) # , exclude = force_out_vec)
    
    # get predictions / fit scores (print)
    predictions <- pfit %>% predict(xS)
    fit <- data.frame(R2(predictions, y), RMSE(predictions, y))
    names(fit) <- c('R2','RMSE')
    fit <-round(fit, 3)
    print(fit)
    
    return(pfit)
    }

# demonstrate function
CH4_lass0 <- glmnet_lassoCV(y_var, x_vars, Guild_CH4_d0)
#CH4_lass0.2 <- glmnet_lassoCV3(y_var, x_vars, Guild_CH4_d0, force_in = c("C","Salinity.x","SRB"))

# dev screen advanced
# CH4_lass0 <- glmnet_lassoCV2(y_var, x_vars, Guild_CH40)
# glmnet_coefs(CH4_lass0)            
### dev advanced check -- coeffs
#glmnet_coefs(CH4_lass00)
# glmnet_coefs(CH4_lass0.2)

# dev advanced check -- plots
CH4_lass0_cplot <- plot_lasso_coefs(CH4_lass0, guild_colors)

options(repr.plot.width=4, repr.plot.height=4)
 CH4_lass0_cplot # + labs(title="Plot of length \n by dose")
# extract lasso coeffs
glmnet_coefs = function(model, output = 'coefs', print_ns = FALSE){
    
    # Extract model coeffs
    coef_m <- data.frame(data.matrix(coef(model)))   # extract coeffs
    names(coef_m) <- 'coef'                          # label coef  
    coef_m$var <- row.names(coef_m)                  # get feats
    
    # get only model selected features
    sig_ests <- coef_m[!coef_m$coef==0,]             # drop NS feats  
    sig_ests <- sig_ests[-1,]                        # drop intercept
    
    sig_feats <- sig_ests$var                        # features only
    
    # non-sig feats
    all_feats <- coef_m$var
    ns_int <- all_feats[!all_feats %in% sig_feats]
    ns <- unlist(paste(ns_int[-1], collapse = ', ' ))
    ns_out <- paste('non-sig: ', ns)
    
    if (print_ns == TRUE){
    print(ns_out)}
    
    if (output == 'vars') {
        data_out <- sig_feats
        } else {data_out <- sig_ests}
      
    return(data_out)
}

# demonstrate function
# glmnet_coefs(CH4_lass0)                             # return table of coefs (e.g. for plotting)
vars <- glmnet_coefs(CH4_lass0, output = 'vars')      # return only list of signif vars
vars

#'CO2_mg_m2_h''Bulk_dens''H2O_FPS''pH''N''SO4''CH4_ac''MOB_II''MOB_IIa''AOB'

# rename col to var
names(guild_colors)[1] <- "var"
#colors <- guild_colors

#model <- CH4_lass0

plot_lasso_coefs = function(model_out, colors, nafill = '#708090'){

    # dev lasso effects plot
    coefs <- glmnet_coefs(model_out)              # Get significant effects

    # merge w colors
    coef_color <- data.frame(merge(coefs, colors, all.x=T, by = 'var'))     # coef_color
    colors2 <- as.character(coef_color$color)                               # only levels of var used  
    colors2[is.na(colors2)] <- '#708090'                                    # replace NA color 
    coef_color$color <- as.factor(colors2)                                  # replace colors (levels)
                         
    # reorder data / factors by coef
    coef_color <- coef_color[order(coef_color$coef),]                       # sort by coef
    coef_color$var <- reorder(coef_color$var, coef_color$coef)              # reorder var by coef
    coef_color$color <- reorder(coef_color$color, coef_color$coef)          # reorder color by coef

    # Make ggplot 
    plot_colors <- c(paste(coef_color$color))                               # Get colors for ggplot

    plot <- ggplot(data=coef_color, aes(x=var, y=coef, fill = var)) +       # base plot
      geom_bar(stat="identity") + scale_fill_manual(values= plot_colors) +  # fill by 
      geom_hline(yintercept=0) +                                            # add y = 0 line
      coord_flip()+ theme_minimal() + theme(legend.position = "none")       # flip; clean plot

    return(plot)
}

# Demonstrate function
names(guild_colors)[1] <- 'var'

CH4_lass0_cplot <- plot_lasso_coefs(CH4_lass0, guild_colors)

options(repr.plot.width=6, repr.plot.height=4)
CH4_lass0_cplot

# to determine models with _all significant features_, not done by leaps, meifly, bestglm packages
# modified from: https://rpubs.com/kaz_yos/exhaustive with coeffs & significance tests added
# could further get R2 if run LM in parallel, not sure its needed.
all_subsets_whh = function(outcome, predictors, dataset, p_cut, sort = "AIC"){
    
    ## Create list of models
    list.of.models <- lapply(seq_along((predictors)), function(n) {

        left.hand.side  <- outcome
        right.hand.side <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")

        paste(left.hand.side, right.hand.side, sep = "  ~  ")
    })

    ## Convert to a vector
    vector.of.models <- unlist(list.of.models)

    ## Fit glm to all models
    list.of.fits <- lapply(vector.of.models, function(x) {

        formula    <- as.formula(x)
        fit        <- glm(formula, data = dataset)
        result.AIC <- extractAIC(fit)

        # all sig filter place
        coeffsM <- data.frame(coef(summary(fit))[,1:4])            # Get coeffs & significance
        coeffsM <- coeffsM[-1,]                                    # drops intercept from NS test

        # test for significance of all feats
        bad_model <- ifelse(any(coeffsM$Pr...t.. > p_cut) == TRUE, "no", "yes")
               
        # get BIC
        BIC <- BIC(fit)
        AICc <- AICc(fit)
        
        ## combine data
        data.frame(num.predictors = result.AIC[1],
               AIC            = result.AIC[2],
               BIC            = BIC,
               AICc           = AICc,
               model          = x,
               all_sig            = bad_model)
    })

    ## Collapse to a data frame
    result <- data.frame(do.call(rbind, list.of.fits))
    result <- result[order(result[,sort]),]
    result <- result[result$all_sig == 'yes',]
    return(result)
}# use function, 10s? / 1023 permutes
subs <- all_subsets_whh('CH4_logn1', vars, Guild_CH4_dS_df, p_cut = 0.05, sort = "AIC")
#dim(subs); subs[1:20,]  # show output

# MOB_mods <- subs[grep("MOB", subs$model),]
# head(MOB_mods)

# library(parallel)

all_subsets_whhP = function(outcome, predictors, dataset, p_cut, sort = "AIC", ncores = 20 ){
    
    ## Create list of models
    list.of.models <- lapply(seq_along((predictors)), function(n) {

        left.hand.side  <- outcome
        right.hand.side <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")

        paste(left.hand.side, right.hand.side, sep = "  ~  ")
    })

    ## Convert to a vector
    vector.of.models <- unlist(list.of.models)

    ## Fit glm to all models
    list.of.fits <- mclapply(vector.of.models, function(x) {

        formula    <- as.formula(x)
        fit        <- glm(formula, data = dataset)
        result.AIC <- extractAIC(fit)

        # all sig filter place
        coeffsM <- data.frame(coef(summary(fit))[,1:4])            # Get coeffs & significance
        coeffsM <- coeffsM[-1,]                                    # drops intercept from NS test

        # test for significance of all feats
        bad_model <- ifelse(any(coeffsM$Pr...t.. > p_cut) == TRUE, "no", "yes")
               
        # get BIC
        BIC <- BIC(fit)
        AICc <- AICc(fit)
        
        ## combine data
        data.frame(n_pred = result.AIC[1],
               AIC            = result.AIC[2],
               BIC            = BIC,
               AICc           = AICc,
               model           = x,
               all_sig            = bad_model)
    }, mc.cores = ncores)

    ## Collapse to a data frame
    #result <- data.frame(do.call(rbind, list.of.fits))
    result <- data.frame(rbindlist(list.of.fits))        # Much faster than do.call!
    result <- result[order(result[,sort]),]
    result <- result[result$all_sig == 'yes',]
    return(result)
}

# use function, 10s? / 1023 permutes single core; 1.5s w 20 cores
subs <- all_subsets_whhP('CH4_logn1', vars, Guild_CH4_dS_df, p_cut = 0.05, sort = "AIC", ncores = 20)

# dim(subs); 
subs[1:20,]  # show output
# Set vars for function: 
outcome  <- 'CH4_logn1'   # y var NAME
predictors  <- vars       # X vars NAME  -- taken from above output from LASSO
dataset <- Guild_CH4_dS_df  # data name# Compare GLM significances
testx <- glm(CH4_logn1 ~ CO2_mg_m2_h + Bulk_dens + H2O_FPS + N + SO4 + CH4_ac, data = Guild_CH4_dS_df)
#extractAIC(testx)
#BIC(testx)
# AICc(testx)
summary(testx)# compare SEM significances

# Warns that a solution can not be found w test, not test 2!

# "BEST MODEL" by GLM filter
SEM_las_CH4_mod.1 <- sem('CH4_logn1 ~ CO2_mg_m2_h + CH4_ac + MOB_IIa + AOB', data = Guild_CH4_dS_df,

#SEM_las_CH4_mod.1 <- sem('CH4_logn1 ~ CO2_mg_m2_h + Bulk_dens + H2O_FPS + N + SO4 + CH4_ac', data = Guild_CH4_dS_df,
                         fixed.x=FALSE, estimator = 'mlm', orthogonal = TRUE)

# Best MOB model
#SEM_las_CH4_mod.1 <- sem('CH4_logn1 ~ CO2_mg_m2_h + Bulk_dens + pH + CH4_ac + MOB_IIa', data = Guild_CH4_dS_df,
#                         fixed.x=FALSE, estimator = 'mlm', orthogonal = TRUE)

new_mod <- SEM_las_CH4_mod.1
get_SEM_R2s(new_mod)
sem_fit_sum(new_mod, sem_fit_params)
get_SEM_nonsig_vars(new_mod, p_cut =0.05)form_glm <- as.formula(form)Guild_CH4_dS_df <- data.frame(Guild_CH4_dS)# make linear model
lin_mod <- glm(test, data = Guild_CH4_dS_df)
summary(lin_mod)

# glmMulti model selection
# best_lin <- glmulti(lin_mod, level = 1, crit = 'bic')
#lin_mod$p_value
# make function for SEM selection 
sem_filter = function(subsets, data, top_n_models = 25, sort_by = 'bic', n_mod_ret = 10){
    
    # reduce candidate models by top n (subtests)
    mod_cands <- subsets[1:top_n_models,]                              # subset
    mods <- mod_cands$model                                            # get models
    
    # orig feats to keep 
    keep_orig <- c("n_pred", "AIC", "BIC", "AICc")
    orig_feats <- mod_cands[keep_orig]
    names(orig_feats) <- c("n_pred","aic_glm","bic_glm","aicc_glm")
    orig_feats$n_pred <- orig_feats$n_pred - 1                         # note init n_pred incl. y 
   
    # Get performance metrics for each
    mod_outs = {}

    for (i in seq(1:length(mods))){

        # Get SEM model
        form <- paste(mods[i])
        SEM_mod <- sem(form, data = data, fixed.x=FALSE, estimator = 'mlm', orthogonal = TRUE)
    
        # Get R2 data
        r2 <- get_SEM_R2s(SEM_mod)
        R2 <-r2[,2]
    
        # Get aic, bic data
        keep <- c('aic','bic','npar')
        fits <- sem_fit_sum(SEM_mod, sem_fit_params)
        fits <- fits[keep]
      
        # all feats significant?
        sig_test <- get_SEM_nonsig_vars(SEM_mod, p_cut =0.05)
        if (dim(sig_test)[1] == 0){
            all_sig <- 'yes'
            } else {
            all_sig <- 'no'
            }
        
        # get feats from original data
        keep_orig <- orig_feats[i,]        
                            
        # combine feats & assign to list
        out <- data.frame(fits, R2, form, all_sig, keep_orig)
        mod_outs[[i]] <- out #paste(mods[i])
    
        }

    # Assemble and sort list
    output <- data.frame(do.call(rbind, mod_outs))
    output <- output[order(output[sort_by]),]
    
    # clean and return output
    model <- seq(1:(length(mods)))
    output <- data.frame(model, output)
    
    
    output <- output[output$all_sig == 'yes',]
    output <- output[1:(n_mod_ret),]
    output <- output[complete.cases(output),]
    row.names(output) <- output$model
    return(output)
    
}semF <-sem_filter(subs, data = Guild_CH4_dS_df, top_n_models = 140, sort_by = 'aic', n_mod_ret = 140)
#semF

# 140 models in 20s
# Parallel SEM filter

# make function for SEM selection 
sem_filterP = function(subsets, data, top_n_models = 25, sort_by = 'bic', ncores = 20){       # n_mod_ret = 10, 
    
    # reduce candidate models by top n (subtests)
    mod_cands <- subsets[1:top_n_models,]                              # subset
    mods <- mod_cands$model                                            # get models
    
    # orig feats to keep 
    keep_orig <- c("model","n_pred", "AIC", "BIC", "AICc")
    orig_feats <- data.frame(mod_cands[keep_orig])
    names(orig_feats) <- c('model',"n_pred","aic_glm","bic_glm","aicc_glm")
    orig_feats$n_pred <- orig_feats$n_pred - 1                         # note init n_pred incl. y 
   
    # Get performance metrics for each
    mod_outs <- mclapply(mods, function(x){

        # Get SEM model
        form <- paste(x)
        SEM_mod <- sem(form, data = data, fixed.x=FALSE, estimator = 'mlm', orthogonal = TRUE)
    
        # Get R2 data
        r2 <- get_SEM_R2s(SEM_mod)
        R2 <-r2[,2]
    
        # Get aic, bic data
        keep <- c('aic','bic','npar')
        fits <- sem_fit_sum(SEM_mod, sem_fit_params)
        fits <- fits[keep]
      
        # all feats significant?
        sig_test <- get_SEM_nonsig_vars(SEM_mod, p_cut =0.05)
        if (dim(sig_test)[1] == 0){
            all_sig <- 'yes'
            } else {
            all_sig <- 'no'
            }
        
        # get feats from original data   
        keep_orig <- orig_feats[orig_feats$model == x,]  
        keep_orig <- keep_orig[,-1]
        
        # combine feats & assign to list
        data.frame(fits, R2, form, all_sig, keep_orig)

        }, mc.cores = ncores)

    # Assemble and sort list
    output <- data.frame(do.call(rbind, mod_outs))
    # output <- data.frame(rbindlist(mod_outs))         - rbindlist abandoned, no speedup & crash on NA
    output <- output[order(output[sort_by]),]
    
    # clean and return output
    model <- seq(1:(length(mods)))
    output <- data.frame(model, output)
    output <- output[output$all_sig == 'yes',]
    #output <- output[1:(n_mod_ret),]
    output <- output[complete.cases(output),]
    output <- data.frame(output)
    row.names(output) <- output$model
    
    # make cols numeric, fix forms       
    output2 <- data.frame(sapply(output, as.numeric))
    output2$form <-output$form
    
    return(output2)
    
}

semF <-sem_filterP(subs, data = Guild_CH4_dS_df, top_n_models = 140, sort_by = 'aic', ncores = 20) # depricate , n_mod_ret = 140,  
head(semF)
# 140 models in 22s (lapply, single core?  but looks like 5 in activity monitor?)
# 140 models in 3.5s (mcapply, 20 cores)
# test dev / debug
asemF <- sem_filterP(aCH4_lass0_subs, data = Guild_CH40, 
                                  top_n_models = 500, sort_by = 'bic', n_mod_ret = 500, ncores = 20)

#asemF <- aCH4_lass0_semFhead(asemF)
#asemF[asemF$n_pred == "...",]
# (is.na(asemF$n_pred) == FALSE)

#asemF[50:100,]


# function for R2 vs. AIC max by reg. resids
filter_byLM_resid = function(sem_tab, n_keep_split = 4, sort_by = "resid", split_by = "n_pred",
#filter_lSEM_tab = function(sem_tab, n_keep_split = 4, sort_by = "resid", split_by = "n_pred",
                             R2_filter = 0, resid_0filter = TRUE){

    # get linear model & add residuals to table
    lm <- lm(R2 ~ aic, data = sem_tab)
    resid <- round(resid(lm),3)
    sem_tab$resid <- unlist(resid)

    # sort, split, and filter list
    sem_tab <- sem_tab[order(-sem_tab[,sort_by]),]
    sem_tab_split <- split(sem_tab, sem_tab[,split_by])

    filtered_list <- lapply(sem_tab_split, function(x) {head(x, n_keep_split)}) 
    sem_tab_filt <- rbindlist(rev(filtered_list))

    if (resid_0filter == TRUE) {
        out_tab <- sem_tab_filt[sem_tab_filt$resid > 0,]
        } else {out_tab <- sem_tab_filt}
    
    out_tab <- data.frame(out_tab[out_tab$R2 > R2_filter])
    out_tab
    
    }


# demonstrate function
filt_models <- filter_byLM_resid(semF, n_keep_split = 3, R2_filter = 0.5, resid_0filter = TRUE)
#filt_models

# filter_lSEM_tab(semF, n_keep_split = 4)
# filter_lSEM_tab(sem_tab, n_keep_split = 4, sort_by = "resid", split_by = "n_pred", R2_filter = 0, resid_filter = TRUE)

# preview filter function on plot
#lasso_sel_sumplot <- plot_SEM_lasso_mods4_lin(aCH4_lass0, guild_colors, asemF, filt_models, exclude_mods = '')
#lasso_sel_sumplot <- plot_SEM_lasso_mods4_lin(CH4_lass0, guild_colors, semF, filt_models, exclude_mods = '')

options(repr.plot.width=12, repr.plot.height= 3)
lasso_sel_sumplot


########################################################################################
### 6) Plot LASSO-SEM model screening functions


# HELPER fxn -- get ggplot color wheel hue function -- Maybe use later?
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = (n + 1))
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# plot SEM LASSO models
plot_SEM_lasso_mods = function(sem_tab, label = 'FALSE', colors = ''){
    
    
    # Get color palette, default or custom
    options(warn=-1)                                            # need to supress warning in this function
    if (colors == ''){                                          # if (colors warns for vector, not needed
          cols <- gg_color_hue(max(sem_tab$n_pred))} else {     # default ggplot colors  OR
          cols <- colors}                                       # custom color input     generated in combined plots  
    
    # Get n_pred as factor 
    sem_tab$n_pred <- factor(sem_tab$n_pred)
      
    # plot w/o labels
    no_lab <- ggplot(sem_tab, aes(x=aic, y=R2, color = n_pred)) + geom_point() +
                    theme_minimal() + scale_color_manual(values = cols) #+ 
                    #geom_smooth(method=lm, aes(fill = n_pred))
    
    # plot w labels
    label_plot <- ggplot(sem_tab, aes(x=aic, y=R2, color = n_pred, label = model)) + geom_point() +
                     geom_text_repel() + theme_minimal() + scale_color_manual(values = cols)
    
    # is plot labeled?
    if (label == 'FALSE'){
          out_plot <- no_lab} else {
          out_plot <- label_plot}
    
    return(out_plot)    
}

#head(asemF)

# Demonstrate function
all_models <- plot_SEM_lasso_mods(semF)
#all_models <- plot_SEM_lasso_mods(asemF)

options(repr.plot.width=4.5, repr.plot.height=4)
all_models

# With restricted input, labels  -- Incorporated into g) below
# semF_top_n <- semF[order(semF$aic_glm),]
# semF_top_n <- semF_top_n[1:25,]

# plot selected models
# select_models <- plot_SEM_lasso_mods(semF_top_n, label = 'TRUE')
# select_models

# function for SEM model screening
plot_sem_lasso_selected = function(all_sem, sel_sem, x_pad = 0.9, y_pad = 0.9){
    
    # params from selected
    selected <- sel_sem$model                                                                   # selected mods
    x_min <- min(sel_sem$aic)*x_pad
    y_min <- min(sel_sem$R2)*y_pad
    
    # recode selected feats in all DF
    all_sem$selected <- ifelse(all_sem$model %in% selected, "true", "false")                    # code selected
    all_sem$predictors <- as.factor(ifelse(all_sem$selected == 'true', "true", all_sem$n_pred)) # recode n_pred
    all_sem$sel_lab <- ifelse(all_sem$selected == 'true', all_sem$model, "")                    # label selected 

    # Get colors for plot
    n_cols <- length(unique(all_sem$n_pred))     # n pred
    cols0 <- gg_color_hue(n_cols)                # gg pal
    cols <- c(cols0, '#000000')                  # add black for selected 

    # plot data
    ggplot(all_sem, aes(x=aic, y=R2, color = predictors, label = sel_lab)) +
        geom_point(aes(shape=selected))  + theme_minimal() + geom_text_repel() + 
        scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = cols) + 
        ylim(y_min, NA) + xlim(x_min, NA)
    
}


#g2

# Demonstrate function
get_models <- c(34, 35, 67, 69, 82, 107, 113, 116) 
sel_models <- semF[semF$model %in% get_models,]

options(repr.plot.width = 5, repr.plot.height = 4)
plot_sem_lasso_selected(semF, sel_models)
#save_plot(filename, plot, width, height)
# Function for splitting sem formulas & resorting form vars (by LASSO abs coef/ VIP)
coef_sort_sem_form_tab = function(lasso_model, sem_tab, fxn_operator = '~ ', sep_table = 'FALSE'){
    
    ### Prepare empty var inclusion table ###

    # sort lasso coefs by abs. vals
    coefs <- glmnet_coefs(lasso_model)                                                        # FXN get coefs  
    coefs$abs <- abs(coefs$coef)                                                              # get coef abs
    coefs <- coefs[order(-coefs$abs),]                                                        # sort

    # get vars & model matrix data
    ranked_vars <- coefs$var                                                                  # get ranked vars 
    n_coef <- length(ranked_vars)                                                             # n vars   (cols)
    n_mods <- dim(sem_tab)[1]                                                                 # n models (rows)             

    # make vars, model matrix DF                                                              # w/ empty model components
    coef_m <- data.frame(matrix(ncol=n_coef,nrow=n_mods, dimnames=list(NULL, ranked_vars)))   # make empty df
    sem_mods_vars <- cbind(sem_tab, coef_m)                                                   # join sem models

       
    ### Fill model matrix w. true / false ###
    for (i in ranked_vars)                                                                   # for each var
         sem_mods_vars  <- model_feat_ext(sem_mods_vars, i)                                  # use FXN for fill

 
    ### Rebuild sorted model formulas 
    # Get y var to redo formulas, (assume y same in all)
    form1 <- paste(sem_mods_vars$form[1])                                                    # first formula 
    #fxn_operator <- "~ "
    y_var <- unlist(strsplit(form1, fxn_operator))[1]                                        # get y var
    y_var <- paste0(y_var, fxn_operator)                                                     # add op back
    # y_var

    # Get x var formula componets from TF data
    sort_Xform <- apply(sem_mods_vars, 1, function(r) paste0(names(r)[r == TRUE], collapse = " + "))
    sem_mods_vars$form <- paste0(y_var, sort_Xform)                                         # add y, replace form                   
                                 
    ### Select output mode ###
    if (sep_table == 'FALSE'){
        out_tab <- sem_mods_vars[names(sem_tab)]} else {                                    # revert to input cols
        out_tab <- sem_mods_vars}                                                           # keep split cols (for plots)  
    
    out_tab
    
}# Demonstrate function:
coef_tab <-coef_sort_sem_form_tab(CH4_lass0, sel_models, sep_table = 'FALSE', fxn_operator = '~ ')  # Get sorted table
#coef_tab <- coef_sort_sem_form_tab(CH4_lass0, sel_models, sep_table = 'TRUE', fxn_operator = '~ ')   # Get full table (for plots)
coef_tab# HELPER function: single col TRUE / FALSE replace with value
TF_replace = function(coef_tab, var, value){
    var <- ifelse(coef_tab[,var] == TRUE, coef_tab[,value], 0)
    return(var)
}leaps_plot_sem_models = function(coef_tab, sort_var = "aic", color_var = "R2", exclude_mods = ""){
    
    ### replace T/F var columns with data (color_var) ###
    
    # Get variable sets for apply:
    coef_cols <- names(coef_tab)                             # all columns
    log_cols_test <- unlist(lapply(coef_tab, is.logical))    # which are T/F
    log_cols <- coef_cols[log_cols_test]                     # get T/F columns
    var_cols <- log_cols
    other_cols <- coef_cols[!coef_cols %in% log_cols]        # non-T/F columns 

    # make new aic2x var
    # coef_tab$aic2x <- rescale(coef_tab$aic)*rescale(coef_tab$aic_glm)
    # coef_tab$aic2x <- coef_tab$aic*coef_tab$aic_glm

    # Apply TF replace function to logical, passing others 
    var_col_list = {}                                        # empty list
    for (i in log_cols)                                      # for each T/F col 
        var_col_list[[i]] <- TF_replace(coef_tab, i, color_var)  # insert val or 0

    coef_tab[log_cols] <- data.frame(var_col_list)     # -- USE THIS ONE  # replace log cols in orig data
    
    ### get ggplot data ###
    
    # exclude models 
    if (exclude_mods == ""){
        coef_tab <- coef_tab} else {
        coef_tab <- coef_tab[!coef_tab$model %in% exclude_mods,]
    }
       

    
    # reorder model factor for plot
    coef_tab$model <- as.factor(as.character(coef_tab$model))
    coef_tab$model <- reorder(coef_tab$model, coef_tab[,sort_var])

    # Get vars & data
    gg_vars <- c('model', log_cols)
    gg_data <- coef_tab[gg_vars]
    gg_data_m <- melt(gg_data)
    
    # clean data
    vars <- gg_data_m$variable
    #gg_data_m$variable <- strtrim(vars, 9)          # trim var names, fails due to factor ordering issue, upstream issues before factor too
    gg_data_m[gg_data_m == 0] <-NA

    # plot data 
    plot <- ggplot(data = gg_data_m, aes(y=model, x=variable, fill=value)) + 
          geom_tile(color="white") +  scale_fill_gradient(na.value = "white", high = "#132B43", low = "#56B1F7",) + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + 
          labs(fill = color_var)
 
    plot
    
 }# coef_tab <-coef_tab[-1,]
options(repr.plot.width = 4, repr.plot.height = 4)

leaps_plot_sem_models(coef_tab, sort_var = "aic", color_var = "R2", exclude_mods = '')
# leaps_plot_sem_models(coef_tab, sort_var = "aic", color_var = "R2", exclude_mods = c('34', '113', '121')) # 



########################################################################################
### 7) SEM model / formula sorting by lasso variable importance


# extract model formula features into columns 
model_feat_ext <- function(mod_ext_df, var){

    exact <- paste0('\\<', var, '\\>')
    mod_ext_df[var] <- grepl(mod_ext_df$form, pattern = exact)
    return(mod_ext_df)

}

# HELPER function: single col TRUE / FALSE replace with value
TF_replace = function(coef_tab, var, value){
    var <- ifelse(coef_tab[,var] == TRUE, coef_tab[,value], 0)
    return(var)
}

# Function for splitting sem formulas & resorting form vars (by LASSO abs coef/ VIP)
coef_sort_sem_form_tab = function(lasso_model, sem_tab, fxn_operator = '~ ', sep_table = 'FALSE'){
    
    ### Prepare empty var inclusion table ###

    # sort lasso coefs by abs. vals
    coefs <- glmnet_coefs(lasso_model)                                                        # FXN get coefs  
    coefs$abs <- abs(coefs$coef)                                                              # get coef abs
    coefs <- coefs[order(-coefs$abs),]                                                        # sort

    # get vars & model matrix data
    ranked_vars <- coefs$var                                                                  # get ranked vars 
    n_coef <- length(ranked_vars)                                                             # n vars   (cols)
    n_mods <- dim(sem_tab)[1]                                                                 # n models (rows)             

    # make vars, model matrix DF                                                              # w/ empty model components
    coef_m <- data.frame(matrix(ncol=n_coef,nrow=n_mods, dimnames=list(NULL, ranked_vars)))   # make empty df
    sem_mods_vars <- cbind(sem_tab, coef_m)                                                   # join sem models

       
    ### Fill model matrix w. true / false ###
    for (i in ranked_vars)                                                                   # for each var
         sem_mods_vars  <- model_feat_ext(sem_mods_vars, i)                                  # use FXN for fill

 
    ### Rebuild sorted model formulas 
    # Get y var to redo formulas, (assume y same in all)
    form1 <- paste(sem_mods_vars$form[1])                                                    # first formula 
    y_var <- unlist(strsplit(form1, fxn_operator))[1]                                        # get y var
    y_var <- paste0(y_var, fxn_operator)                                                     # add op back

    # Get x var formula componets from TF data
    sort_Xform <- apply(sem_mods_vars, 1, function(r) paste0(names(r)[r == TRUE], collapse = " + "))
    sem_mods_vars$form <- paste0(y_var, sort_Xform)                                         # add y, replace form                   
                                 
    ### Select output mode ###
    if (sep_table == 'FALSE'){
        out_tab <- sem_mods_vars[names(sem_tab)]} else {                                    # revert to input cols
        out_tab <- sem_mods_vars}                                                           # keep split cols (for plots)  
    
    out_tab
    
}

# Demonstrate function:
coef_tab <-coef_sort_sem_form_tab(CH4_lass0, sel_models, sep_table = 'FALSE', fxn_operator = '~ ')  # Get sorted table
#coef_tab <- coef_sort_sem_form_tab(CH4_lass0, sel_models, sep_table = 'TRUE', fxn_operator = '~ ')   # Get full table (for plots)
coef_tab

# coef matrix function (trimmed from "coef_sort_sem_form_tab" function in Leaps like plot, then expanded)
# some legacy crud in here, no need for color var in TF replace, for ex...sort_var, etc...
# sem_model_coef_mtx = function(lasso_model, sel_models, sort_var = "aic", color_var = "R2", exclude_mods = ""){

sort_sem_model_coef_mtx = function(lasso_model, sel_models, VIP_sort = TRUE, rev_order = TRUE, color_var = "R2"){
   
    # sorts sem formula table rows by variables & their lasso VIP (for sem model trees, or heatmaps)
    
    ### replace T/F var columns with data (color_var) ###
    # get table of vars encoded as T/F, with other model scoring data
    coef_tab <- coef_sort_sem_form_tab(lasso_model, sel_models, sep_table = 'TRUE', fxn_operator = '~ ') 
   
    # Get T/F variable sets for apply:
    coef_cols <- names(coef_tab)                                 # all columns
    logic_cols_test <- unlist(lapply(coef_tab, is.logical))      # which are T/F
    log_cols <- coef_cols[logic_cols_test]                       # get vect ot T/F columns only 
 
    # Apply TF replace function to logical, passing others 
    var_col_list = {}                                            # empty list
    for (i in log_cols)                                          # for each T/F col 
        var_col_list[[i]] <- TF_replace(coef_tab, i, color_var)  # insert val or 0 -- HELPER TF_replace row-wise
    
    # Make DF with binary data
    coef_tab2 <- data.frame(var_col_list)                        # replace log cols in orig data
    row.names(coef_tab2) <- coef_tab$form                        # get row names as formulas
    
    # all vars 1, for VIP sort OR R2 instead of sort             # VIP = lasso var imp = scaled coefs 
    if (VIP_sort == TRUE) {
        coef_tab2[!coef_tab2== 0] <- 1}                          # all vars T = 1, else 0 
        else {}                                                  # else T = R2 for heatmap, breaks VIP sort
    
    ### Sort models by vars included, by lasso VIP rank ###
    n_col <- dim(coef_tab2)[2]                                   # get number of columns
    col_nos <- seq(1:n_col)                                      # get vector of col_nos
    coef_tab2 <- dfOrder(coef_tab2, col_nos)                     # sort rows by all columns, in order of columns 
    
    if (rev_order == TRUE){
        coef_tab2 <- coef_tab2[rev(seq_len(nrow(coef_tab2))),]}  # reverse row order
        else{}
    
    return(coef_tab2)

}

# get demo input data
#filt_models <- filter_lSEM_tab(semF, n_keep_split = 3, R2_filter = 0.5, resid_filter = TRUE)

# head(subs)

# demonstrate function 
matrix <-sort_sem_model_coef_mtx(CH4_lass0, filt_models)                 # , rev_order = F)
#matrix <-sort_sem_model_coef_mtx(CH4_lass0, subs2, color_var = "n_pred")
#matrix <-sort_sem_model_coef_mtx(CH4_lass0, asubs, color_var = "n_pred")
matrix <-sort_sem_model_coef_mtx(CH4_lass0, filt_models, VIP_sort = F)  # for heatmap, later
dim(matrix); head(matrix)
matU <- unique(matrix); dim(matU)
# show plot of data (DELETE later)
# library(Heatplus)
# suppressMessages(library(vegan))

# basic heatmap2 to show dendograms
# options(repr.plot.width=7, repr.plot.height=7)
# hm <- heatmap(as.matrix(matrix), Rowv = NA, Colv = NA, margins = c(8, 30))
# hm
vip_sort_sem_models = function(lasso_model, sem_model_matrix){
    
    # get coef weights from lasso
    coefs <- glmnet_coefs(lasso_model)                      # return only list of signif vars
    vars <- unlist(names(sem_model_matrix))                 # get var names 
    new_matrix <- sem_model_matrix                          # new matrix copy 

    # fill in lasso VIP coef weights for each time var included
    for (i in vars){                                            # for each var 
        value <- coefs[coefs$var == i, 1]                       # get coef
        new_matrix[i] <- ifelse(!sem_model_matrix[i]==0, value, 0)}        # if matrix = 1, coef, else 0
    
    return(new_matrix)
}

# demonstrate function
vip_matrix <- vip_sort_sem_models(CH4_lass0, matrix)
head(vip_matrix)

sem_model_vip_sorter = function(lasso_model, filt_models, color_var = 'n_pred'){
    
    
  ### 1) Get data needed from existing functions ###                                       # likely could more gracefully refactor
    # VIP sort formulas & rows into matrix
    matrix <- sort_sem_model_coef_mtx(lasso_model, filt_models, color_var = color_var)     # VIP formulas matrix   
    vip_matrix <- vip_sort_sem_models(lasso_model, matrix)                                 # VIP sort model rows  
    
    
  ### 2) Patch sem_table columns / values to align with VIP sorted ### 
    # Get all data with TF columns per var (unsorted)
    model_tab <-coef_sort_sem_form_tab(lasso_model, filt_models, sep_table = 'TRUE', fxn_operator = '~ ')

    # Apply TF replace function to logical columns, here
    log_cols <- names(vip_matrix)                                # get logical columns, here from VIP matrix ! 
    var_col_list = {}                                            # empty list
    for (i in log_cols)                                          # for each T/F col 
        var_col_list[[i]] <- TF_replace(model_tab, i, color_var) # insert val or 0 -- HELPER TF_replace row-wise
    
    # Get dataframe of vars, make binary                         # TF replace set up to fill from rows, modify  
    var_cols_mtx <- data.frame(var_col_list)                     # insert TF replace into DF 
    var_cols_mtx[!var_cols_mtx == 0] <-1                         # replace all !=0 w 1 
    
    # replace TF vars in model tab w binary
    model_tab[log_cols] <-var_cols_mtx                           # replace logical columns in model table  
 
    
  ### 3) Merge, sort & clean models sorted by VIP ###
    # prepare VIP matrix for merge
    vip_matrix[!vip_matrix == 0] <- 1                                            
    vip_matrix$vip_indx <- seq(1:nrow(vip_matrix))
       
    # merge model_tab & vip matrix
    combined_matrix <- merge(model_tab, vip_matrix, by = log_cols)
    combined_matrix <- combined_matrix[order(combined_matrix$vip_indx),]
    
    # drop var sorting cols, finalize table
    all_cols <- names(combined_matrix)
    keep_cols <- all_cols[!all_cols %in% log_cols]
    out_model_table <- combined_matrix[keep_cols]
    
    out_model_table
    
}

filt_models_VIP <- sem_model_vip_sorter(CH4_lass0, filt_models)
#filt_models_VIP#; filt_models


########################################################################################
### 8) Plots for model selection, heatmap, clusters -- descr needed below

model_cluster_dend = function(sem_mod_mtx, k_clusts = 4){
    
    # distance and hier clustering
    dd <- dist(sem_mod_mtx, method = "euclidean")                  # distance matrix
    hc <- hclust(dd, method = "ward.D2")                           # hier clustering

    # override cluster tree tip ordering, revert to input
    n_col <- dim(sem_mod_mtx)[1]                                   # get number of columns
    ord <- seq(1:n_col)                                            # get vector of col_nos       #ord <- seq(1:12)  
    ord_r <- rev(ord)                                              # reverse row ordering
    tree <- reorder(hc, ord_r)                                     # order cluster tree by input row order

    # make dend extend -> ggplot tree
    dend <- as.dendrogram(tree)
    
    # cols <- brewer.pal(n = k_clusts, name = "RdBu")               # could set colors manually, also list
    model_tree <- dend         # value = cols) 
    model_tree <- dend %>% set("branches_k_color", k = k_clusts)
    ggd1 <- as.ggdend(model_tree)
    
    plot <- ggplot(ggd1, horiz = TRUE, labels = FALSE)
    plot
}

# demonstrate function 
options(repr.plot.width=1, repr.plot.height=5)

model_tree <- model_cluster_dend(vip_matrix, k_clusts = 3)
model_tree
suppressMessages(library(ggtree))matrix2 <- vip_matrix 

#matrix2 <- data.matrix(matrix)

# Make clusters
dd <- dist(matrix2, method = "euclidean")                  # distance matrix
hc <- hclust(dd, method = "ward.D2")                       # hier clustering

# override tree tip ordering, revert to input
n_col <- dim(matrix2)[1]                                   # get number of columns
ord <- seq(1:n_col)                                        # get vector of col_nos       #ord <- seq(1:12)  
ord_r <- rev(ord)                                          # reverse row ordering

x <- reorder(hc, ord_r)                                    # order cluster tree by input row order

# plot data
ggtree <- ggtree(x, ladderize = F) + geom_tiplab() + coord_cartesian(clip = 'off') +
            theme_tree2(plot.margin=margin(6, 350, 6, 6)) #+ xlim(0, 0.08) #

options(repr.plot.width=10, repr.plot.height=5)
ggtree
#ggdendrogram(hc, rotate = F)# x# install.packages("ggdendro")
suppressMessages(library("ggdendro"))#install.packages('dendextend')
suppressMessages(library(dendextend))# simple ggdend...
options(repr.plot.width= 6, repr.plot.height=5)
ggdendrogram(x, rotate = TRUE)suppressMessages(library(ape))colors = c("red", "blue", "green", "black")
clus4 = cutree(hc, 4)

#options(repr.plot.width=20, repr.plot.height=20)
options(repr.plot.width=10, repr.plot.height=5)
plot(as.phylo(hc), cex = 1, label.offset = 0.5, tip.color = colors[clus4])
sem_model_feats_heatmap = function(lasso_model, filt_models, vip_matrix, color_var = "R2", high_col, low_col){

    # Get R2 matrix data
    R2_matrix <-sort_sem_model_coef_mtx(lasso_model, filt_models, VIP_sort = F) 
    R2_matrix <- data.frame(R2_matrix)
    R2_matrix$form <- row.names(R2_matrix)
    R2_vars <- c(unlist(names(R2_matrix)))                    # get list of vars in R2 data                       

    # replace orignal filt models data w R2 matrix
    in_mods_sorted <- sem_model_vip_sorter(lasso_model, filt_models)  # using new wrapper fxn
    mod_coef_tab <- merge(in_mods_sorted, R2_matrix, by="form")
    row.names(mod_coef_tab) <- mod_coef_tab$form              # make form row names (for sort by VIP)

    # Sort R2 by VIP model ordering 
    vip_rows <- c(unlist(row.names(vip_matrix)))              # Get VIP ordered models
    mod_coef_tab <- mod_coef_tab[vip_rows,]                   # sort coef data by VIP order

    # Get keep vars from mod coef table
    if(color_var == "R2"){                                    # if R2, get all R2 vars w. data from R2 matrix
        keep_cols <- R2_vars}
        else {keep_cols <- color_var}
    
    keep_vars <- c("model", keep_cols)                        # keep model NUMBER for ggplot, heat data 
    gg_data <- mod_coef_tab[,keep_vars]                       # get ggplot data

    # make model number factor for plotting
    gg_data$model <- as.factor(as.character(gg_data$model)) 

    # Get row ordering as model number for sorting melted 
    sort_mods <- unlist(as.character(paste(gg_data$model)))   # get sorted model numbers, still factor, pass to numeric
    sort_mods <- as.character(as.numeric(sort_mods))          # make char again (probably too many conv) 

    # melt data & reorder levels of model number factor in melted
    gg_data_m <- suppressMessages(melt(gg_data))              # Melt data
    
    gg_data_m$model <- as.factor(as.character(gg_data_m$model))
    gg_data_m$model <- reorder(gg_data_m$model, new.order = rev(sort_mods))    
    gg_data_m[gg_data_m == 0] <-NA                            # clean data for ggplot, 0 <- NA for blank
    #gg_data_m$variable <- strtrim(vars, 9)                   # trim var names, fails due to factor ordering 
    

    # plot data 
    plot <- ggplot(data = gg_data_m, aes(y=model, x=variable, fill=value)) + 
          geom_tile(color="white") +  scale_fill_gradient(na.value = "white", high = high_col, low = low_col) +             
          theme_minimal() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 0), axis.title.x = element_blank(), axis.title.y = element_blank()) + 
          scale_x_discrete(position = "top") +
          labs(fill = color_var) #+ theme(legend.position = "bottom")
 
    plot
    
    #mod_coef_tab
    #gg_data
}

# demonstrate function
options(repr.plot.width=5, repr.plot.height=4)

sem_model_feats_heatmap(CH4_lass0, filt_models, vip_matrix, color_var = "R2", high_col = "#132B43", low_col = "#56B1F7")

#sem_model_feats_heatmap(CH4_lass0, filt_models, vip_matrix, color_var = R2_vars, high_col = "#132B43", low_col = "#56B1F7")
# sem_model_feats_heatmap(CH4_lass0, filt_models, vip_matrix, color_var = "aic", high_col = "#ccff33", low_col = "#608000")
# sem_model_feats_heatmap(CH4_lass0, filt_models, vip_matrix, color_var = "resid", high_col = "#b3003b", low_col = "#ffff66")

combine_model_heatmaps = function(lasso_model, filt_models, vip_matrix, model_tree, alt_metric = "resid", 
                                  rel_widths = c(1,4,0.72, 0.25,1)){  # hclust = model_tree,
    
    # get heatmap plots, 3x
    R2_plot <- sem_model_feats_heatmap(lasso_model, filt_models, vip_matrix, color_var = "R2", high_col = "#132B43", low_col = "#56B1F7")
    aic_plot <- sem_model_feats_heatmap(lasso_model, filt_models, vip_matrix, color_var = "aic", high_col = "#ccff33", low_col = "#608000")
    resid_plot <- sem_model_feats_heatmap(lasso_model, filt_models, vip_matrix, color_var = "resid", high_col = "#b3003b", low_col = "#ffff66")
    
    # get legends only
    R2_legend <- as_ggplot(get_legend(R2_plot))
    aic_legend <- as_ggplot(get_legend(aic_plot))
    resid_legend <- as_ggplot(get_legend(resid_plot))

    # select alt metric in legends
    if(alt_metric == "resid"){
        alt_legend <- resid_legend
        } else {alt_legend <- aic_legend}

    legends <- plot_grid(R2_legend, alt_legend, nrow = 2, ncol = 1, align ="hv")

    # Strip legends, labels from main plots
    R2_plot <- R2_plot + theme(legend.position = "none")
    aic_plot <- aic_plot + theme(legend.position = "none", axis.text.y = element_blank())
    resid_plot <- resid_plot + theme(legend.position = "none", axis.text.y = element_blank())

    # select alt metric in main plots
    if(alt_metric == "resid"){
        alt_plot <- resid_plot
        } else {alt_plot <- aic_plot}

    options(warn = -1)
    if (model_tree == ''){
    plot <- plot_grid(R2_plot, alt_plot, NULL, legends, ncol = 4, nrow = 1, align = "h", rel_widths = c(6,1,0.5,1))
    } else {
    plot <- plot_grid(model_tree, R2_plot, alt_plot, NULL, legends, ncol = 5, nrow = 1, align = "h", rel_widths = rel_widths)
    }
   
    plot

}

# demonstrate function
options(repr.plot.width=5, repr.plot.height=4)

combine_model_heatmaps(CH4_lass0, filt_models, vip_matrix, model_tree, alt_metric = "resid")
# combine_model_heatmaps(CH4_lass0, filt_models, vip_matrix, model_tree, alt_metric = "aic") 
# combine_model_heatmaps(CH4_lass0, filt_models, vip_matrix, model_tree = '', alt_metric = "resid")

compare_sem_model_heat_dend = function(lasso_model, filt_models, k_clusts = 3, R2_filt = 0.5, alt_metric = resid,
                                      rel_widths = c(1,4,0.72, 0.25,1) ){
    
    # keep only models > R2_filt
    filt_models <- filt_models[filt_models$R2 > R2_filt, ]
    
    # break formulas to binary, sort formulas, cluster & dend
    matrix <-sort_sem_model_coef_mtx(lasso_model, filt_models)          # get formulas as binary by var
    vip_matrix <- vip_sort_sem_models(lasso_model, matrix)              # sort formulas by lasso VIP rank
    model_tree <- model_cluster_dend(vip_matrix, k_clusts = k_clusts)   # get dendogram
   
    # get output heat/dend plot
    plot <- combine_model_heatmaps(CH4_lass0, filt_models, vip_matrix, model_tree, 
                                   alt_metric = alt_metric)
    plot
    
    }

options(repr.plot.width=5, repr.plot.height=4)
compare_sem_model_heat_dend(CH4_lass0, filt_models, k_clusts = 4, R2_filt = 0.7, alt_metric = "resid")


########################################################################################
### 9) Combined lasSEM model selection plots

# This is actually breaking where n_pred aren't all in table
plot_SEM_lasso_mods2 = function(sem_tab, sel_sort = 'aic_glm', top_n = 25) {

    # Get color palettes
    cols <- gg_color_hue(max(sem_tab$n_pred))                      # get base colors
    n_pred <- sort(unique(sem_tab$n_pred))                         # don't make factor, crashes function
    col_df <- data.frame(cols, n_pred)                             # color x n_preds DF
    all_cols <- c(as.character(col_df$cols))                       # colors for all models
    
    ## Plot ALL models ##
    all_models <- plot_SEM_lasso_mods(sem_tab, label = 'FALSE', colors = all_cols) 
   
    ## Model subsets ##
    sem_top_n <- sem_tab[order(sem_tab[sel_sort]),]                # sort models by selected var    
    sem_top_n <- sem_top_n[1:top_n,]                               # keep only top mods
    
    # subset colors
    n_pred_sub <- data.frame(n_pred = unique(sem_top_n$n_pred))    # only cols in subset
    col_df_sub <- merge(col_df, n_pred_sub)                        # merge with color key (inner)
    sub_cols <- c(as.character(col_df_sub$cols))                   # get colors for plot
    
    ## Plot selected models ## 
    select_models <- plot_SEM_lasso_mods(sem_top_n, label = 'TRUE', colors = sub_cols) #+ 
                     #theme(legend.position = "none")

    ## Plot both figs ## 
    dual_plot <- plot_grid(all_models, select_models, ncol=2, nrow =1)

    return(dual_plot)

}options(repr.plot.width=9, repr.plot.height=4)
plot_SEM_lasso_mods2(semF, sel_sort = 'aic_glm', top_n = 25)
# head(asemF)
# unique(semF$n_pred)

# test edge case
# plot_SEM_lasso_mods2a(asemF, sel_sort = 'aic_glm', top_n = 200)
plot_SEM_lasso_mods3 = function(lass_mod, group_colors, sem_tab, sel_sort = 'aic_glm', top_n = 25){
    
    coef_plot <- plot_lasso_coefs(lass_mod, group_colors)
    models_plot <- plot_SEM_lasso_mods2(sem_tab, sel_sort = sel_sort, top_n = top_n)
        
    triple_plot <- plot_grid(coef_plot, models_plot, ncol=2, nrow =1, rel_widths = c(1, 2))
    
    #return(models_plot)
    return(triple_plot)
    
}# demonstrate function
lasso_sel_sumplot <- plot_SEM_lasso_mods3(CH4_lass0, guild_colors, semF, sel_sort = 'aic_glm', top_n = 25)

options(repr.plot.width=10, repr.plot.height=2.5)
lasso_sel_sumplot# edge case test (def well below)
# alasso_sel_sumplot <- plot_SEM_lasso_mods3(aCH4_lass0, guild_colors, asemF, sel_sort = 'aic_glm', top_n = 200)
# alasso_sel_sumplotplot_SEM_lasso_mods4_lin = function(lass_mod, group_colors, sem_tab, sel_models,
                               sort_var = "aic", color_var = "R2", exclude_mods = '', 
                               rel_widths =c(1, 1, 1, 1), ncol = 4, nrow =1){
    
    coef_plot   <-   plot_lasso_coefs(lass_mod, group_colors)
    models_plot <-   plot_SEM_lasso_mods(sem_tab) + geom_smooth(method = "lm", colour="darkgrey",
                         linetype = "dashed", size = 1, se = FALSE)    
    
    selected_plot <- plot_sem_lasso_selected(sem_tab, sel_models) + theme(legend.position = "none")
    
    coef_tab    <-   coef_sort_sem_form_tab(lass_mod, sel_models, sep_table = 'TRUE', fxn_operator = '~ ')  
    leaps_plot  <-   leaps_plot_sem_models(coef_tab, sort_var, color_var, exclude_mods) +
                        scale_x_discrete(position = "top") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 0), axis.title.x = element_blank())
    
    quad_plot <- plot_grid(coef_plot, models_plot, selected_plot, leaps_plot, 
                           ncol = ncol, nrow = nrow, labels = "AUTO", align = "h", rel_widths = rel_widths)
    
    return(quad_plot)
    
}# demonstrate function
lasso_sel_sumplot <- plot_SEM_lasso_mods4_lin(CH4_lass0, guild_colors, semF, sel_models, exclude_mods = '')
# lasso_sel_sumplot <- plot_SEM_lasso_mods4_lin(CH4_lass0, guild_colors, semF, sel_models, ncol =2, nrow = 2)

options(repr.plot.width=12, repr.plot.height= 3)
lasso_sel_sumplot
plot_lasso_semfilt_list = function(lasso_mod, group_colors, sem_tab, sel_models){
                               #sort_var = "aic", color_var = "R2", exclude_mods = '', 
                               #rel_widths =c(1, 1, 1, 1), ncol = 4, nrow =1){
    
    # get_plots for coefs, sem models & filt
    coef_plot   <-   plot_lasso_coefs(lasso_mod, group_colors)
                                   
    models_plot <-   plot_SEM_lasso_mods(sem_tab) + geom_smooth(method = "lm", colour="darkgrey",
                         linetype = "dashed", size = 1, se = FALSE)    
    
    selected_plot <- plot_sem_lasso_selected(sem_tab, sel_models) + theme(legend.position = "none") 
    
    # return list of plots
    plot_list  <- list(coef_plot, models_plot, selected_plot)
    plot_names <- c("coef_plot", "models_plot", "selected_plot")
    names(plot_list) <- plot_names

    return(plot_list)  
    
}

# demonstrate function

lasso_sel_sumplot <- plot_lasso_semfilt_list(CH4_lass0, guild_colors, semF, filt_models)
# lasso_sel_sumplot#[1]; #lasso_sel_sumplot$coef_plot
plot_grid(lasso_sel_sumplot[[1]], lasso_sel_sumplot[[2]], lasso_sel_sumplot[[3]], nrow =1, ncol = 3)
#                       ncol = ncol, nrow = nrow, labels = "AUTO", align = "h", rel_widths = rel_widths)

combine_model_heatmaps_list = function(lasso_model, filt_models, vip_matrix, alt_metric = "resid",                                       
                                  legend_row = 2, legend_col = 1){   
    
    # get heatmap plots, 3x
    R2_plot <- sem_model_feats_heatmap(lasso_model, filt_models, vip_matrix, color_var = "R2", high_col = "#132B43", low_col = "#56B1F7")
    aic_plot <- sem_model_feats_heatmap(lasso_model, filt_models, vip_matrix, color_var = "aic", high_col = "#ccff33", low_col = "#608000")
    resid_plot <- sem_model_feats_heatmap(lasso_model, filt_models, vip_matrix, color_var = "resid", high_col = "#b3003b", low_col = "#ffff66")
    
    # get legends only
    R2_legend <- as_ggplot(get_legend(R2_plot)) 
    aic_legend <- as_ggplot(get_legend(aic_plot))
    resid_legend <- as_ggplot(get_legend(resid_plot))

    # select alt metric in legends
    if(alt_metric == "resid"){
        alt_legend <- resid_legend
        } else {alt_legend <- aic_legend}

    legends <- plot_grid(R2_legend, alt_legend, nrow = 1, ncol = 2, align ="hv") 
    #legends <- plot_grid(R2_legend, alt_legend, nrow = 2, ncol = 1, align ="hv")

    # Strip legends, labels from main plots
    R2_plot <- R2_plot + theme(legend.position = "none")
    aic_plot <- aic_plot + theme(legend.position = "none", axis.text.y = element_blank())
    resid_plot <- resid_plot + theme(legend.position = "none", axis.text.y = element_blank())

    # select alt metric in main plots
    if(alt_metric == "resid"){
        alt_plot <- resid_plot
        } else {alt_plot <- aic_plot}

    
    # return list of plots
    plot_list  <- list(R2_plot, alt_plot, legends)
    plot_names <- c("R2_plot", "alt_plot", "legends")
    names(plot_list) <- plot_names

    return(plot_list)
}

# demonstrate function
multi_heat <- combine_model_heatmaps_list(CH4_lass0, filt_models, vip_matrix, alt_metric = "resid", legend_row = 2, legend_col = 1)
multi_heat
# plot_grid(multi_heat[[1]], multi_heat[[2]], multi_heat[[3]], nrow =1, ncol = 3, align = "h")

plot_lasso_sem_models = function(lasso_model, group_colors, sem_tab, filt_models, plot_title = '',
                                k_clusts = 3, R2_filt = 0, alt_metric = resid, 
                                rel_widths_main = c(4, 4.5, 4, .75, 4, 0.6)){
    
    # get lasso & sem model plots
    lasSem <- plot_lasso_semfilt_list(lasso_model, group_colors, sem_tab, filt_models)
    coef <- lasSem$coef_plot 
    sems <- lasSem$models_plot
    filt <- lasSem$selected_plot
    
    # filter & sort data for heatmaps #
    filt_models <- filt_models[filt_models$R2 > R2_filt, ]              # keep only models > R2_filt
    matrix <-sort_sem_model_coef_mtx(lasso_model, filt_models)          # get formulas as binary by var
    vip_matrix <- vip_sort_sem_models(lasso_model, matrix)              # sort formulas by lasso VIP rank
    
    
    # get heatmap components
    tree <- model_cluster_dend(vip_matrix, k_clusts = k_clusts)         # get dendogram
    heats <- combine_model_heatmaps_list(lasso_model, filt_models, vip_matrix,   # get heatmaps  
                alt_metric = "resid", legend_row = 2, legend_col = 1)  
    R2 <- heats$R2_plot
    res <- heats$alt_plot
    leg <- heats$legends #+ legend.text=element_text(size=6)

    # get all plots, no legends
    p1 <- plot_grid(coef, sems, filt, tree, R2, res, nrow = 1, ncol =6, 
                      labels = c('A', 'B', 'C', 'D', '',''),
                      align = 'h', rel_widths = rel_widths_main)
    
    # kludge to fit legends without smashing together (not same scale as other plots)
    p2 <- plot_grid(p1, leg, nrow =1, ncol = 2, rel_widths = c(10,1.5))
    
    # add plot title
    title = paste0('    ', plot_title)
    p2 + draw_figure_label(label = title, size = 14, fontface = "bold") 
    
}
    # plot <- plot_grid(coef, sems, filt, tree, R2, res, NULL, leg, nrow = 1, ncol =8, 
    #                  align = 'h', rel_widths = c(4, 4.5, 4, .75, 4, 0.6, 0.25, 1))   

semF[semF$model==107,]

#options(repr.plot.width=10, repr.plot.height= 4.25)
options(repr.plot.width=12, repr.plot.height= 3)
#options(repr.plot.width=16, repr.plot.height= 4.25)

lasso_sel_sumplot <- plot_lasso_sem_models(CH4_lass0, guild_colors, semF, filt_models, plot_title = "CH4_df_SRB",
                                           k_clusts = 3, R2_filt = 0.7, rel_widths_main = c(4, 4.5, 4, .75, 4, 0.7))
lasso_sel_sumplot



########################################################################################
### 10) Wrappers for combined lasso -> SEM models, plots, workflows


# likely these functions should be returning lists, including lasso_mod.
# NOTE changing lasso mod may hit downstream...

lasso_mod_plot = function(y_var, x_vars, data, var_colors, force_in = '', plot_title = ''){ 
    
    lasso_mod  <- glmnet_lassoCV(y_var, x_vars, data, force_in = force_in)              # LASSO models 
    lasso_vars <- glmnet_coefs(lasso_mod, output = 'vars', print_ns = F)                # Get LASSO sig vars
    coef_plot  <- plot_lasso_coefs(lasso_mod, var_colors) + ggtitle(plot_title)         # Plot coefs
    
    out        <- list(lasso_mod, lasso_vars, coef_plot)                                # list of feats
    names(out) <- c('lasso_mod', 'lasso_vars', 'coef_plot')
    out
}

#CH4_lass0 <- glmnet_lassoCV(y_var, x_vars, Guild_CH4_d0, force_in = c("Salinity.x",'SRB')) 

# Demonstrate function
options(repr.plot.width=4, repr.plot.height= 3)

CH4_lassO <- lasso_mod_plot(y_var, x_vars, Guild_CH4_dS_df, guild_colors, plot_title = 'CH4_delta: all vars')
CH4_lassO$coef_plot

# demonstrate variable forcing
CH4_lass.f <- lasso_mod_plot(y_var, x_vars, Guild_CH4_d0, guild_colors, 
                             force_in = c("Salinity.x",'SRB'),#, 'Bulk_dens', 'SO4', 'C'),
                             plot_title = 'CH4_delta: force Sal, SRB')
CH4_lass.f$coef_plot
#CH4_lass.f$lasso_mod

keep_only_models = function(lasso_model, sem_models, force_filt, all_any = "all" ){    

    # Get matrix data from models: fits/forms (coef_tab) & vars included in models (model_matrix)
    coef_tab <- coef_sort_sem_form_tab(lasso_model, sem_models, sep_table = 'TRUE', fxn_operator = '~ ')  # VIP sorted sem_models, TF enc
    model_matrix <-sort_sem_model_coef_mtx(lasso_model, sem_models, color_var = "n_pred")                                       # vars only VIP sorted, 01 enc 

    # Get forced vars
    force_n <-length(force_filt)                                          # get n vars forced
    if(all_any == "any"){force_n <-1}                                     # override n vars forced for "any"

    if (force_n ==1){                                                     # kludge for force_n = 1, for indexing & rowSums
        force_mtx <- model_matrix[force_filt]                             # get forced cols (binary T,F as 1,0)                      
        force_mtx$blank <- 0}   
        else {force_mtx <- model_matrix[,force_filt]}                     # get forced cols (binary T,F as 1,0) 
  
    # models to keep & their fit/formula data 
    keep_models <- row.names(force_mtx[rowSums(force_mtx) == force_n, ])  # get models passing (n TRUE = sum (0,1) = force n)
    filtered_mods <- coef_tab[coef_tab$form %in% keep_models,]            # get all data for these models

    # clean up/out variable matrix cols from data
    vars <- names(model_matrix)                                           # get vars in all models 
    scores <- names(coef_tab)                                             # get all model report cols
    keep_scores <- scores[!scores %in% vars]                              # get model scoring & formula cols (no var cols)
    filtered_out <- filtered_mods[keep_scores]                            # get model scores & formula data  (no var cols)

    return(filtered_out)
    
    #model_matrix; #force_mtx; #force_n; #keep_models
}

# extract model from mod_obj
#lasso_CH4.0  <- CH4_lass0$lasso_mod        # note need to extract object component if using this way !!
#lasso_CH4.f1 <- CH4_lass.f$lasso_mod       # note need to extract object component if using this way !!

# rename glm out vars to shoehorn into function made for sem out - later in wrapper
subs2 <- subs  # note formula is called model in glm output!, also need to rename n_pred & substitute into fxn... 

names(subs2) <- c('n_pred','AIC','BIC','AICc','form','all_sig')

# demonstrate function, note takes lasso model not object

#keep_only_models(CH4_lass0, subs2, force_filt = c("Salinity.x",'SRB'), all_any = "all" )
k <- keep_only_models(CH4_lass0, subs2, force_filt = c('SO4', 'N'), all_any = "any" )
k <- keep_only_models(CH4_lass0, subs2, force_filt = c('SO4'), all_any = "any" )

#o <-keep_only_models(lasso_model, sem_mods, force_filt = c("Salinity.x",'SRB'), all_any = "all")
#o <-keep_only_models(lasso_model, sem_mods, force_filt = c('SRB'), all_any = "any")
#o
#    coef_tab <- coef_sort_sem_form_tab(lasso_model, subs2, sep_table = 'TRUE', fxn_operator = '~ ')  # VIP sorted sem_models, TF enc
#    model_matrix <-sort_sem_model_coef_mtx(lasso_model, sem_models)                                       # vars only VIP sorted, 01 enc 


# might as well use above function?
# here, lazily passing over opportunity to only run SEM on top GLM models
# NEED TO PUT IN HARD FILTER for models only here?  Or separate function for filtering ?
# Note no VIP sorting here, messes up plot downstream (? not sure why)

lasso_to_SEM_mods = function(lasso_obj, y_var, data, n_cores = 20, 
            top_glm = '', p_cut = 0.05, sort_by = 'aic', GLM_print = FALSE,    # glm params - Default should suffice
            hard_filter = FALSE, force_filt = '', keep_any_all = "all") {      # var forcing params
                              
           
    # Get vars from lasso object
    lasso_vars <- lasso_obj$lasso_vars
    lasso_model <- lasso_obj$lasso_mod
    
    # GLM selection for ALL significant subsets
    glm_models <- all_subsets_whhP(y_var, lasso_vars, data, p_cut = p_cut, ncores = n_cores)
    
    
    # glm to SEM model pre-filters 
    if (force_filt == ''){hard_filter = FALSE}                                 # bypass hard filter if no filter vars
    
    if (hard_filter == TRUE){                                                  # keeps only models incl. listed vars
         names(glm_models)[5] <- "form"                                        # kludge- match hard filter to SEM fmt 
         glm_models <-keep_only_models(lasso_model, glm_models, force_filt = force_filt, all_any = keep_any_all)
         names(glm_models)[5] <- "model"}                                      # change back to orig 
     
    n_GLM <- dim(glm_models)[1]
    if (top_glm == ''){top_glm <- n_GLM}                                       # only top scoring models 

    
    # Get only significant SEM models                                  
    sem_models <- sem_filterP(glm_models, data = data, sort_by = sort_by, 
                              top_n_models = top_glm, ncores = n_cores)        # n_mod_ret = top_glm,
       
    # Print n sig vars?
    if (GLM_print == TRUE) {
        print(paste0(n_GLM, ' sig models'))}  
    
    return(sem_models)

}
 

# demonstrate function
sem_mods1 <- lasso_to_SEM_mods(CH4_lassO, 'CH4_logn1', Guild_CH4_dS_df)                       # no filter (default)
# sem_mods1 <- lasso_to_SEM_mods2(CH4_lassO, 'CH4_logn1', Guild_CH4_d0, hard_filter = TRUE)   # test no filt wrong flag
# sem_mods2 <- lasso_to_SEM_mods2(CH4_lass.f, 'CH4_logn1', Guild_CH4_dS_df, 
#              hard_filter = FALSE) #, force_filt = c(''), keep_any_all = "any")
# force in SRB, N, all the options 
sem_mods3 <- lasso_to_SEM_mods(CH4_lass.f, 'CH4_logn1', Guild_CH4_dS_df, n_cores = 20, 
              top_glm = '', p_cut = 0.05, sort_by = 'aic', GLM_print = TRUE, #)#,           # glm options - default fine
              hard_filter = TRUE, force_filt = c('SRB', 'N'), keep_any_all = "all")         # var filt options

# note keep_any_all = "all" requires all force filt vars in models; "any" allows mods with any forced vars to pass..

dim(sem_mods1); head(sem_mods1)
#dim(sem_mods2); head(sem_mods2)
#dim(sem_mods3); head(sem_mods3)

# Data prep is LASSO obj, then GLM -> SEM -> filter
# CH4_lassO <- lasso_mod_plot(y_var, x_vars, Guild_CH4_dS_df, guild_colors, plot_title = 'CH4_delta: all vars')

# options(repr.plot.width=4, repr.plot.height= 3)
# CH4_lassO$coef_plot

lasso_obj_to_sem_select_plots = function(lasso_obj, y_var, data, var_colors,
                                force_filt = '', hard_filter = TRUE, keep_any_all = "all",  # force_params 
                                nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,              # sem filts
                                heat_filt = 0.7, k_clusts = 3){                              # heatmap params
    
    # Get all sem models and LM resid filtered
    sem_models <- lasso_to_SEM_mods(lasso_obj, y_var, data, force_filt = force_filt, hard_filter = hard_filter)      
    filt_models <- filter_byLM_resid(sem_models, n_keep_split = nbest, R2_filter = keep_filt, 
                                     resid_0filter = resid_0filt)
    
    # Get composite model selection plot
    plot_title <- deparse(substitute(lasso_obj))
    lasso_mod <- lasso_obj$lasso_mod
    lasso_sel_sumplot <- plot_lasso_sem_models(lasso_mod, var_colors, sem_models, filt_models, 
                                           plot_title = plot_title, k_clusts = k_clusts, R2_filt = heat_filt) 
                                           #, rel_widths_main = c(4, 4.5, 4, .75, 4, 0.7))

    # VIP sort sem & filt models after plotting
    sem_models  <- sem_model_vip_sorter(lasso_mod, sem_models)     
    filt_models <- sem_model_vip_sorter(lasso_mod, filt_models)
    
    # Gather output components
    output = list()
    output$sem_models <- sem_models
    output$filt_models <- filt_models
    output$plot <- lasso_sel_sumplot
    
    #output <- lasso_sel_sumplot
    
    return(output)
    
}

# demonstrate function, makes lasSEM object unpacked below

CH4_lassO_sem <- lasso_obj_to_sem_select_plots(CH4_lassO, 'CH4_logn1', Guild_CH4_dS_df, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)

# all sem: 
sems = CH4_lassO_sem$sem_models       # sems

# filtered models
filts = CH4_lassO_sem$filt_models    #filts
# merge(sems, filts, by = "form")   # test subset works (filt in matching after VIP formula sorting)

# show combined plot
options(repr.plot.width=12, repr.plot.height= 3)

CH4_lassO_sem$plot

# Show with forced in SRB
CH4d_fSRB <- lasso_mod_plot(y_var, x_vars, Guild_CH4_dS_df, guild_colors, 
                             force_in = c('SRB'))

# forced SRB
CH4d_fSRB_sem <- lasso_obj_to_sem_select_plots(CH4d_fSRB, 'CH4_logn1', Guild_CH4_dS_df, guild_colors, 
                                         force_filt = 'SRB', nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.5, k_clusts = 4)

CH4d_fSRB_sem$plot


# CH4d_fSRB_sem$filt_models

# lasso <- CH4d_fSRB$lasso_mod
#lasso
# filt_mods <- CH4d_fSRB_sem$filt_models[,-13]

#compare_sem_model_heat_dend(lasso, filt_mods, k_clusts = 4, R2_filt = 0.7, alt_metric = "resid")

#filt_models

#compare_sem_model_heat_dend(CH4_lass0, filt_models, k_clusts = 4, R2_filt = 0.7, alt_metric = "resid")

########################################################################################
### 11) SEM auto-compositing functions

setwd("~/Desktop/Sal_SEM_models")

### corresponding uncomposited test input is:
composite_test_semEx7 <- read.table("SEM_models/fxn_example_tables/SEM_Exper7_res_CH4_lasso.txt", sep='\t') #in SEM models
semEx7_res <- composite_test_semEx7
semEx7_res

# new SUBSETS w completness scenarios:
# Doesn't work for auto-composite!! ?

keep_mods2 <- c('69','79','107','113','116','118','121')  # same as manual exp 7, all pass
keep_mods3 <- c(keep_mods2, '75')                         # only model not in manual with all passing
keep_mods5 <- c(keep_mods2, '10', '34')                   # add single comp fails
keep_mods6 <- c(keep_mods2, '8','16')                     # double comp fails

# Subsetting data, note just changing keep mods below
# semEx7_res2 <- keep_selected_mods_in_sem_obj(semEx7_res, keep_mods2, sem_object=FALSE)
semEx7_res2 <- semEx7_res[semEx7_res$model %in% keep_mods2,]
#keep_mods <- c('ch4L0C.69','ch4L0C.79','ch4L0C.107','ch4L0C.113','ch4L0C.116','ch4L0C.118','ch4L0C.121')


semEx7_res2



# Resulting output ALL lassSEM objects:
# c(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSO4_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem,  CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)
# sem_object_list <- list(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem, CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)
#sem_object_list#$plot
# head(filt_models)
# head(filt_mods)# filt_mod <- CH4_mod_all_sem$filt_models
# head(filt_mod)# filt_mod_ALL <- CH4_mod_ALL_sem$filt_models
# head(filt_mod_ALL)
# adding SEM response vars with string defs
CH4 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h'
CH4_H2 <- 'CH4_H2~ SRB_syn + CO2_mg_m2_h'
MOB_IIa <- 'MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac'
MOB_I <- 'MOB_I ~ CH4_H2 +CH4_ac'

# Get vector of response vars as strings
Other <- ""
response_vars <- c("CH4", "CH4_H2", "MOB_I","MOB_IIa", "Other")

# Test function
sem_models.2 <- update_sem_model_table(model_test, base_model = "ch4_mod0", new_mod_name = "ch4_mod0a.1", 
                      notes = "drop CH4~MOB_I",
                      formula_edits = c("CH4", "CH4_ug_m2_h ~ CH4_ac + CH4_H2 + MOB_IIa + CO2_mg_m2_h"))
                      # output_prefix = "sem."#non_feat_cols = c("model", "base", "notes", "formula"))

sem_models.2
##### get SEM models, needs to be parallel !!# Demonstrate function
#sem_models.2_out <- run_compare_sem_models(sem_models.2, data, estimator = "mlm")
# sem_models.2_out




remove_yvar_from_df_formula_column =function(df, col_name, fxn_operator = '~ '){
    
    # start function 
    row_data = paste(df[,col_name])                                              # get vect of rows in column

    form_list = list()                                                           # empty list

    for (i in (seq_along(row_data)))     {                                       # for each column
        form_list[i] <- paste(unlist(strsplit(row_data[i], '~ '))[2])            # split, x_vars only to list
        }

    col_replace <- unlist(form_list)                                             # unlist
    df[,col_name] <- col_replace                                                 # replace column data
    return(df)                                                                   # return df

}

# demonstrate function
remove_yvar_from_df_formula_column(sem_models.2, 'MOB_IIa')
#remove_yvar_from_df_formula_column(sem_models.2, 'formula') #-- just drop this one!

keep_selected_mods_in_sem_obj = function(SEM_mod_obj, keep_models_list, 
                                   model_col = "model", form_col = "form", sem_object =TRUE){
    
    # SEM object or model table
    if (sem_object == TRUE){
        sem_models <- SEM_mod_obj$sem_models                              # get SEM models from object
        } else {
        sem_models <- SEM_mod_obj                                         # input is SEM model table
    }
    
    # make model numbers char
    sem_models[model_col] <- as.character(sem_models[,model_col])         # model col for model number

    # Get keep_models subset of SEM_models
    sem_models_keep <- sem_models[sem_models[,model_col] %in% keep_models_list,]
    
    return(sem_models_keep)

    }

# Demonstrate fxn (with object, defaults) --- note object used appears later in workflow
# keep_mods <- c("25", "41", "40", "39")
# keep_selected_mods_in_sem_obj(MOBIIa_delt_ngt_sem, keep_mods)

# fxn with table as input, creating output for next fxns
#keep_mods = c('5','11','12','13','15','16')

# selM_MOBIIa <- keep_selected_mods_in_sem_obj(MOB_IIa_models, keep_mods, sem_object = FALSE,
#                                           model_col = "model", form_col = "form") # optional w defaults
# selM_MOBIIa

# fxn with table as input, creating output for next fxns

# keep_mods = c('ch4L0C.69','ch4L0C.107','ch4L0C.118','ch4L0C.121')

# selM_CH4_comp <- keep_selected_mods_in_sem_obj(ch4_lass0_comp_sem_tab, keep_mods, sem_object = FALSE,
#                                            model_col = "model", form_col = "form") # optional w defaults
# selM_CH4_comp



drop_vars_from_df = function(df, drop_vars){

    df_names <- unlist(names(df))
    keep_df_cols <- df_names[!df_names %in% drop_vars]
    df <- df[keep_df_cols]
    return(df)
}

# demonstrate function -- need more upstream 
# head(selM_CH4_comp)
# drop_vars <- c("formula", "model_out")
# drop_vars_from_df(selM_CH4_comp, drop_vars)
setwd("~/Desktop/Sal_SEM_models")# corresponding uncomposited test input is:
composite_test_semEx7 <- read.table("SEM_models/tables/SEM_Exper7_res_CH4_lasso.txt", sep='\t') #in SEM models
semEx7_res <- composite_test_semEx7
semEx7_res# new SUBSETS w completness scenarios:
# Doesn't work for auto-composite!! ?

keep_mods2 <- c('69','79','107','113','116','118','121')  # same as manual exp 7, all pass
keep_mods3 <- c(keep_mods2, '75')                         # only model not in manual with all passing
keep_mods5 <- c(keep_mods2, '10', '34')                   # add single comp fails
keep_mods6 <- c(keep_mods2, '8','16')                     # double comp fails

# Subsetting data, note just changing keep mods below
semEx7_res2 <- keep_selected_mods_in_sem_obj(semEx7_res, keep_mods2, sem_object=FALSE)

#keep_mods <- c('ch4L0C.69','ch4L0C.79','ch4L0C.107','ch4L0C.113','ch4L0C.116','ch4L0C.118','ch4L0C.121')

semEx7_res2

# Resulting output ALL lassSEM objects:
# c(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSO4_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem,  CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)
# sem_object_list <- list(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem, CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)
#sem_object_list#$plot
# head(filt_models)
# head(filt_mods)# filt_mod <- CH4_mod_all_sem$filt_models
# head(filt_mod)# filt_mod_ALL <- CH4_mod_ALL_sem$filt_models
# head(filt_mod_ALL)
CH4_gen <- c("CH4_ac", "CH4_H2", "CO2_mg_m2_h", "Salinity.x", "SO4", "SRB", "FeRB", "NP")

CH4_ox <- c("MOB_IIa", "MOB_I", "MOB_II", "Bulk_dens", "H2O_FPS", 
            "N", "AOB", "AOA", "NOB", "NO3_N", "NH4_N", "NO3_NH4")

composite_sets <- c("CH4_gen", "CH4_ox")
# Composite_base <- c('0 + CH4_gen + CH4_ox')

# Remove sample names from ALL data set here-- testing dysfunctional case
Guild_CH4_0 <- Guild_CH4
row.names(Guild_CH4_0) <- Guild_CH4_0$Sample
Guild_CH4_0 <- Guild_CH4_0[,-1]

# prepare data for function
# filt_mod <- CH4_mod_all_sem$filt_models
# head(filt_mod)

# head(filt_models)
# get_unique_Xvars_in_formula_dfcol(filt_mod[1,], form_col= "form")



get_unique_Xvars_in_formula_dfcol = function(model_table, form_col = "form"){
    
    # Remove y vars, get only formulae
    no_y_table <- remove_yvar_from_df_formula_column(model_table, form_col)  # drop y var 
    forms_no_y <- no_y_table[,form_col]
    forms_no_y <- str_replace_all(forms_no_y, " ", "")                       # KLUDGE, sometimes extra space introduced

    # Get unique X vars among ALL of these formulas
    all_forms <- paste(forms_no_y, collapse = "+")                          # collapse rows into single string
    all_vars <- strsplit(all_forms, "\\+")                                  # split forms into vars 
    unique_vars <- unique(all_vars[[1]])                                    # get unique vars
    
    return(unique_vars) 
}      

# Demonstrate function
get_unique_Xvars_in_formula_dfcol(semEx7_res, form_col = "form")

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
# get_unique_Xvars_in_formula_dfcol(filt_models, form_col = "form"
# get_unique_Xvars_in_formula_dfcol(filt_mod, form_col = "form")
# get_unique_Xvars_in_formula_dfcol(all_mods, form_col = "form")

# filt_mod2<- filt_mod
# filt_mod2$form <- as.factor(filt_mod2$form)


assign_to_composite_split_df = function(model_table, form_col = "form", composite_sets){
    
    # init accumulators
    var_list       <- list()                                                # X var list / formula
    var_group_list <- list()                                                # X vars in set / formula / set
    out_group_list <- list()                                                # X vars not in set ""
    set_list       <- list()                                                # X lists for each set
    out_list       <- list()                                                # X lists not in each set
    other_list     <- list()                                                # X vars in nether set 
    #other_form     <- list()                                                # X vars in nether set, " + " sep
    
    # Loop over sets (j), rows in model table
    for (j in composite_sets){                                              # get lists of vars from each set

        set_vars <- get(j)                                                  # get set vars data

        for (i in 1:nrow(model_table)){

            # remove y var and get list of X vars from formula
            var_list <- get_unique_Xvars_in_formula_dfcol(model_table[i,], form_col = form_col)
            var_group_list[[i]] <- set_vars[set_vars %in% var_list]               # vars in set   
            var_group_list[[i]] <- paste(var_group_list[[i]], collapse = " + ")   # make formula style (?) 
            
            out_group_list[[i]] <- var_list[!var_list %in% set_vars]              # vars not in set
        }
    
        # Gather results for each set    
        set_list[[j]] <- var_group_list
        out_list[[j]] <- out_group_list
    }    
  
    # Make DFs of vars in and !in sets
    var_sets_df <- data.frame(do.call(cbind, set_list))
    out_vars <- data.frame(do.call(cbind, out_list))                       # df of vars ! in sets
    
    # Get duplicate vars among not-in-set vars (ie not in any set)
    for (i in 1:nrow(out_vars)){

        other_vars <- unlist(c(out_vars[i,]))                             # get all vars in row of df
        names(other_vars) <- NULL                                         # drop their indicies
        other_list[[i]] <- other_vars[duplicated(other_vars)]             # get only duplicated 
        other_list[[i]] <- paste(other_list[[i]], collapse = " + ")
        }
    
    # Add other vars col to var_sets DF
    var_sets_df$other <- other_list
    #var_sets_df$other_fxn <- other_form
    
    return(var_sets_df)
    
}

# Demonstrate function
composite_base_df <- assign_to_composite_split_df(semEx7_res, form_col = "form", composite_sets)
composite_base_df

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
# x <- assign_to_composite(filt_models, form_col = "form", composite_sets)
# composite_base_df <- assign_to_composite_split_df(filt_models, form_col = "form", composite_sets)
#composite_base_df <- assign_to_composite_split_df(filt_mod_ALL, form_col = "form", composite_sets)
#head(filt_mod_ALL)
# composite_base_df# $CH4_gen

get_composite_split_form_table = function(model_table, composite_sets, assign_other = "y_var", other_name = "other",
                                      split_y = "  ~", new_y_op = "~ ", no_int = TRUE, set_ops = " <~ ",      
                                      form_col = "form", model_col = "model", model_out_pfx = "sem."){
    
    ### 1) Get main model data: IDs, y, x_vars, combine ###
    
    # Extract model numbers, out strings (to assign)
    model <- as.character(model_table[,model_col])
    model_out <- paste0(model_out_pfx, model)
    
    # Extract y var
    formulae <- model_table[,form_col]
    y_var <- str_split(formulae, split_y) [[1]][[1]]           # first row, column; assumes all are same
    
    # Build y formula base
    y_ops <- paste(y_var, new_y_op)                            # "y ~"  
    compositeX <- paste(composite_sets, collapse= " + ")       # "comp1 + comp2..."
    if (no_int == TRUE){
        compositeX <- paste0("0 + ", compositeX)               # "0 + comp1 + ..."
    }
    
    # Get y base formula column, name by y_var
    y_form_base <- data.frame(paste0(y_ops, compositeX), stringsAsFactors = FALSE)                   
    names(y_form_base) <- y_var
            
    # X vars : get df of composite var assignments -- prior function
    composite_base <- assign_to_composite_split_df(model_table, composite_sets, form_col = form_col)
    
    # Combine model, y and X 
    composite_df <- data.frame(model_out, model, y_form_base, composite_base, stringsAsFactors = FALSE)
    
    
    ### 2) Assigning "other" vars to formulae ###
    # check if assign_other to y_var or composite
    if (assign_other == "y_var"){                                                    # if assign_other = y_var                          
        assign_other = y_var}                                                        # substitute y_var above
   
    # combine formulae from assigned & other, where neither are blank
    mix_other <- ifelse(composite_df[,other_name] == '', composite_df[,assign_other], # if other '', '', else         
                 ifelse(composite_df[,assign_other] == '', composite_df[,other_name], # if assign blank, else 
                 paste(composite_df[,assign_other], composite_df[,other_name], sep = " + "))) # combine forms   
    
    mix_other_df <- data.frame(mix_other, stringsAsFactors=FALSE)
    
    # add to composite DF, kludge for assign_other sometimes being transposed (CH4_ox, not CH4_gen, y var ???)    
    if (nrow(mix_other_df) > 1) {                                                    # transp has only 1 row
        composite_df[assign_other] = mix_other_df} else {                            # replace w other DF, else
        composite_df[assign_other] = t(mix_other_df)                                 # replace w t(other) DF
    }
    
    
    ### 3) Clean formulas, validate composites, produce table ###
    # add y <~ 1* to composite x_var formulas, where not blank.  First var '1*' sets scale of composite var (+/-) 
    for (j in composite_sets){                                                       # for each composite
        composite_df[j] <- ifelse(composite_df[,j] == '', '',                        # where not blank
                                  paste0(j, set_ops, "1*",composite_df[,j]))         # add y, op, "1*" = scaling
    }
    
    return(composite_df)

}    

# Demonstrate function
get_composite_split_form_table(semEx7_res, form_col = "form", composite_sets, assign_other = "y_var")
get_composite_split_form_table(semEx7_res, form_col = "form", composite_sets, assign_other = "CH4_gen")

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
# get_composite_split_form_table(filt_mod, form_col = "form", composite_sets, assign_other = "y_var")
# get_composite_split_form_table(filt_mod_ALL, form_col = "form", composite_sets, assign_other = "y_var") #filt_mod_ALL
# get_composite_split_form_table(filt_mod_ALL, composite_sets, assign_other = "CH4_ox") #filt_mod_ALL


# second version, testing only bad composites w. paste var name
get_composites_table_failures = function(model_table, composite_sets, assign_other = "y_var", other_name = "other",
                                      split_y = "  ~", new_y_op = "~ ", no_int = TRUE, set_ops = " <~ ",      
                                      form_col = "form", model_col = "model", model_out_pfx = "sem."){
       
    # move this OUTSIDE FXN?  USING any of those options here?    
    # get composite forms model table from above fxn, pass all options
    composite_df <- get_composite_split_form_table(model_table, composite_sets, assign_other = assign_other, 
                                                  other_name = other_name, split_y = split_y, new_y_op = new_y_op, 
                                                  no_int = no_int, set_ops = set_ops, form_col = form_col, 
                                                  model_col = model_col, model_out_pfx = model_out_pfx)
  
    # Extract model numbers, out strings (to assign)
    model <- as.character(composite_df[,model_col])
      
    ### filter non-valid composite models (need > 1 var/composite) ###
    
    # method quirks: 1) str_split keeps y_var & x_var1 together, count as 1 var (not issue)
    #                2) length(vars_i...) counts empty strings as 1 var -- can work around:
    #                    - 2) workarounds: a) need > 2 vars for composite anyway; b) can drop blank w/o counter
    
    # get n vars in composites
    var_counts_plus1 <- list()                                      # row accumulator for length of split forms
    var_counts_compos <- list()                                     # compos. accum. for list of form. len/compos.
    
    for (j in composite_sets){    
       
        for (i in 1:nrow(composite_df)){
                       
            vars_i <- str_split(composite_df[,j], " \\+ ")          # split formulas into x's (keeps y~x1)
            var_counts_plus1[[i]] <-length(vars_i[i][[1]])          # get n elements (x's, BUT "" = 1, not 0)
          }
        
        var_counts_compos[[j]] <- var_counts_plus1
        }
       
    n_vars_compos <- data.frame(do.call(cbind, var_counts_compos))  # DF of variable counts/composite (w caveats above)
    
    
    # do n vars/model/composite exceed SEM min threshold? 
    var_thresh = 1                                                  # var threshold > 1, sem requires 2+)  
    
    for (j in composite_sets){ 
       
         n_vars_compos[j] <- ifelse(n_vars_compos[j] > var_thresh, "", j)
    }
    
    # something funny about data type (<chr[,1]>), simplify to char only
    n_vars_compos <- data.frame(lapply(n_vars_compos, as.character), stringsAsFactors = FALSE)
    n_vars_compos
    
    failed_list <- list()
    for (i in 1:nrow(n_vars_compos)){
        vars <- do.call(c, n_vars_compos[i,])
        names(vars) <- NULL
        vars <- vars[vars != ""]
        failed_list[[i]] <- vars
        }
    
    composite_df$comp_failed <- failed_list
    composite_df <- drop_vars_from_df(composite_df, other_name)
    composite_df
    
}  

# Demonstrate function
comp_f7 <- get_composites_table_failures(semEx7_res, form_col = "form", composite_sets, assign_other = "y_var")
comp_f7

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
#get_composites_table_failures(filt_mod, form_col = "form", composite_sets, assign_other = "y_var") #filt_mod_ALL
# y <- get_composites_table_failures(filt_mod_ALL, composite_sets, assign_other = "CH4_ox") #filt_mod_ALL

send_one_bad_composite_x_to_y_var = function(model_table, fail_comp, y_var, comp_string = " <~ 1\\*"){
    
    # strip formula
    fail_comp_data <- model_table[,fail_comp]                                    # vect of failed comp formulae
    fail_comp_form_y <- paste0(fail_comp, comp_string)                           # get y + op string to delete 
    fail_comp_keep_vars <- str_replace(fail_comp_data, fail_comp_form_y, "")     # replace it w "", leaving x
    model_table['pass_failed'] <- fail_comp_keep_vars                            # new col of failed comp x

    # replace failed composite name in y formula w/ failed composite x           # y form is base of composite
    model_table[y_var] <- ifelse(model_table[,'pass_failed'] == '',              # if failed x is blank, 
                    str_replace(model_table[,y_var], paste0(" \\+ ", fail_comp), model_table[,'pass_failed']),
                    str_replace(model_table[,y_var], fail_comp, model_table[,'pass_failed']))
    
    # drop all failed component formulae from table. Would have bad behavior in loop, but to be run by row 
    model_table[fail_comp] <- ""                            
    model_table <- drop_vars_from_df(model_table, 'pass_failed')                 # drop new pass_failed_col
    
    return(model_table)
   
}

# Demonstrate function
send_one_bad_composite_x_to_y_var(model_table = comp_f7, fail_comp = "CH4_ox", y_var = "CH4_logn1", 
                              comp_string = " <~ 1\\*")

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
# here y is output from "get_composites_table_failures" or equiv (has "comp_failed" col)
#z= y[1,]
#send_bad_composite_x_to_y_var(model_table = z, fail_comp = "CH4_ox", y_var = "CH4_logn1", 
#send_one_bad_composite_x_to_y_var(model_table = y, fail_comp = "CH4_ox", y_var = "CH4_logn1", 
#                              comp_string = " <~ 1\\*")


send_all_bad_composite_x_to_y_var = function(model_table, y_var, comp_string = " <~ 1\\*"){

    results = list()                             

    # Operate on each row, get component_lists and n 
    for (i in 1:nrow(model_table)){
    
        model_table_i <- model_table[i,]                          # Get model table row
        comp_failed_i <- model_table_i[,"comp_failed"][[1]]       # Get comp_failed / row
        n_failed_comp <- length(comp_failed_i)                    # Get length (comp_failed) = n_comp_failed /row
    
        # For row, send FIRST bad composite x to y_var, store as row output
        results_i <- send_one_bad_composite_x_to_y_var(model_table = model_table_i, fail_comp = comp_failed_i[1],
                                                      y_var = y_var, comp_string = comp_string) 
        
        # Accomodate ALL failed composites (if more than 1)
        if(length(comp_failed_i) < 2){                            # if 0 or 1 failed comps, 
            results[[i]] <- results_i                             # keep results_i
        } else {                                                  # if 2 or more failed comps 
            # iterate row function over other failed composites   # note model table = results_i
            for (k in comp_failed_i[-1]){                         # iterate over failed comps (not first):
            results_i <- send_one_bad_composite_x_to_y_var(model_table = results_i, fail_comp = k,
                              y_var = y_var, comp_string = comp_string)   
            }   
        }
        results[[i]] <- results_i                                 # accumulate row results
    }
    out <- do.call(rbind, results)                                       # make accumulated DF
    return(out)
}

# Demonstrate function  not working with Ex7 either, but works embedded in intended host function (g below)
# send_all_bad_composite_x_to_y_var(model_table = y, y_var = "CH4_logn1", comp_string = " <~ 1\\*")

# demonstrate function, BUT SOURCE DATA DOWNSTREAM !!
# here y is output from "get_composites_table_failures" or equiv (has "comp_failed" col)
# send_all_bad_composite_x_to_y_var(model_table = y, y_var = "CH4_logn1", comp_string = " <~ 1\\*")  
    

# send_all_bad_composite_x_to_y_var(model_table = y, y_var = "CH4_logn1", comp_string = " <~ 1\\*")

# BELOW, last "working", then edited to address problems here,
# doesn't work if no fails, part fails, needs if statements around those cases... 

get_composites_model_table = function(model_table, composite_sets, assign_other = "y_var", other_name = "other",
                                      compos_fail_collapse_all = TRUE, # new here, keep some composites?
                                      split_y = "  ~", new_y_op = "~ ", no_int = TRUE, set_ops = " <~ ",      
                                      form_col = "form", model_col = "model", model_out_pfx = "sem."){
    
   ### 1) Get composite fail list data & COUNT fails ###  
    # Get last composite DF with composite fails--    move this OUTSIDE FXN?  USING any of those options here?
    composite_df <- get_composites_table_failures(model_table, composite_sets, assign_other = assign_other, 
                                                  other_name = other_name, split_y = split_y, new_y_op = new_y_op, 
                                                  no_int = no_int, set_ops = set_ops, form_col = form_col, 
                                                  model_col = model_col, model_out_pfx = model_out_pfx)
    # Get model indices
    model <- composite_df[,model_col]
    
    # get data for composite failure classes -- note is list of failed composites/model
    failed_col = "comp_failed"          # col for failed composites, add to FXN options? hardcoded in last FXN
    fail_df <- composite_df[failed_col]
    fail_list <- composite_df[,failed_col]
        
    # count items in failed composites, make DF with model indexing
    fail_len_list <- list()   
    for (i in 1:length(fail_list)){                                       # for rows in composite failed list/col 
        fail_len_list[[i]] <- length((fail_list[i][[1]]))                 # get n items in failed list
    }
    n_fail <- do.call(rbind, fail_len_list)                               # list to rows 
    fails_df <- data.frame(model, n_fail, stringsAsFactors = FALSE)       # add model to fails DF
    
    
   ### 2) SUBSET models on composite pass, fail, OR "neither" (part failure) ###
    n_compos <- length(composite_sets)                                    # another potential FXN option  
    
    pass_models <- fails_df[fails_df["n_fail"] == 0,][,"model"]           # PASS all composites   
    failed_models <- fails_df[fails_df["n_fail"] == n_compos,][,"model"]  # FAIL all composites  
    pass_fail_mods <- c(pass_models, failed_models)                       # list both pass and fail  
    part_fail_models <- model[!model %in% pass_fail_mods]                 # SOME composite fails, others pass
    
    # Fate of failed composite vars (compos_fail_collapse_all == TRUE)
    if(compos_fail_collapse_all == TRUE){                                 # collapse all composites if one fails  
        failed_models <- c(failed_models, part_fail_models)               # add part_failed_models to failed_models
    }    
    
    # Get DATAFRAMES for passing, failed, some fail
    pass_df <- composite_df[model %in% pass_models, ]                     # Data for all passing, ** USE AS IS **
    #failed_df <- composite_df[model %in% failed_models, ]                 # all failing,          modify 
    #part_df <- composite_df[model %in% part_fail_models, ]                # part failed,          replace 

    
   ### 3) Clean failed composite models ###  
       
    if(length(failed_models) > 0){
        
        # Get failed DF
        failed_df <- composite_df[model %in% failed_models, ]                 # all failing,          modify 
    
        # Extract y var name 
        formulae <- model_table[,form_col]                                    # get formulae to split
        y_var <- str_split(formulae, split_y) [[1]][[1]]                      # first row, column; assumes all are same
        replace_orig <- c(model_col, y_var)                                   # c(model, y_var colnames) 
    
        # get orig input formulas for replace in failed out DF
        keep_orig_vars <- c(model_col, form_col)                              # model and formula cols
        keep_orig <- model_table[keep_orig_vars]                              # get mod & form from input table
        orig_fails <- keep_orig[keep_orig[, model_col] %in% failed_models,]   # keep only failed model rows
        orig_fails <- setNames(orig_fails, replace_orig)                      # rename model_col to y_var
 
        # clean failed df  
        failed_df[,composite_sets] <- ''                                      # clean out composite info 
        replace_orig <- c(model_col, y_var)                                   # cols to replace w orig data
        failed_df[replace_orig] <- orig_fails[replace_orig]                   # replace orig data
    
        #Add 0 intercept to composites failures ? 
        if(no_int == TRUE){
            failed_df[y_var] <- str_replace(failed_df[,y_var], new_y_op, paste0(new_y_op, "0 + "))
        }
    } else {failed_df <- data.frame()}
    
    
   ### 4) Clean part failed composite models ###  
    
    if(length(part_fail_models) > 0){        
        
        part_df <- composite_df[model %in% part_fail_models, ]                 # get part df from part_failed 
        
        # Fate of failed composite vars (compos_fail_collapse_all == TRUE)
        if(compos_fail_collapse_all == FALSE){                                 # collapse all composites if one fails  
            part_df <- send_all_bad_composite_x_to_y_var(model_table = part_df, y_var, comp_string = " <~ 1\\*") # LAST ARG TO OPTIONS!
        }
    }  else {part_df <- data.frame()}

    
   ### 5) Combine all results ###    
    if(compos_fail_collapse_all == TRUE){
        out <- rbind(pass_df, failed_df)
    } else {
        out <- rbind(pass_df, part_df, failed_df)
    }

    # RECONSTITUTE FORMULA COLUMN # ID which columns are formulas, get as DF 
    form_cols <- c(y_var, composite_sets)
    form_data <- out[form_cols]                                                       # df of formula cols
    
    # paste formula components for new combined formulae
    paste_args <- c(form_data, sep="\n")                                              # do.call paste args
    formula <- do.call(paste, paste_args)                                             # formulae to DF
    formula <- str_replace(formula, "\n\n", "")                                       # no n\n\ (from empty) 
    out$formula <- formula                                                            # formulae to output
    
    # final cleaning - all formula columns to character -- some glitches obs. otherwise
    for (f in form_cols){
        out[f] <- as.character(out[,f])
    }
    
    return(out) 
}       

# demonstrate function
get_composites_model_table(semEx7_res, composite_sets, assign_other = "CH4_ox", 
                           compos_fail_collapse_all = FALSE)

# Show main user behavior options
semEx7_autocomp <- get_composites_model_table(semEx7_res, composite_sets, assign_other = "CH4_gen", 
                           compos_fail_collapse_all = FALSE)

get_composites_model_table(semEx7_res, composite_sets, assign_other = "y_var", 
                           compos_fail_collapse_all = TRUE)

# Show main user behavior options
semEx7_autocomp <- get_composites_model_table(semEx7_res, composite_sets, assign_other = "CH4_gen", 
                           compos_fail_collapse_all = FALSE)
# # demonstrate function, BUT SOURCE DATA DOWNSTREAM !!
z <- get_composites_model_table(filt_mod, form_col = "form", composite_sets, assign_other = "y_var") #filt_mod_ALL
z <- get_composites_model_table(filt_mod, form_col = "form", composite_sets, assign_other = "CH4_ox") 
y <- get_composites_model_table(filt_mod_ALL, composite_sets, assign_other = "CH4_ox", 
                           compos_fail_collapse_all = FALSE) #filt_mod_ALL


setwd("~/Desktop/Sal_SEM_models")

# get external latent model table
ch4_lass0_comp_sems <- read.table("SEM_models/fxn_example_tables/CH4_flux_latent_mod_form_setup2.txt", header = T, na.strings=c(NA),)
ch4_lass0_comp_sems$other <- ""                   #head(ch4_lass0_comp_sems)

# remove base model formula
ncols_start <- dim(ch4_lass0_comp_sems)[2]
ch4_lass0_comp_sem_tab <- ch4_lass0_comp_sems[, -ncols_start]  # drop base formula

# add formula
#ch4_lass0_comp_sem_tab$formula <- with(ch4_lass0_comp_sem_tab, paste(CH4_flux, CH4_gen, CH4_ox, other, sep = '\n'))
#ch4_lass0_comp_sem_tab[,1:7]

ch4_lass0_comp_sem_tab <- ch4_lass0_comp_sem_tab[,1:7]  
ch4_lass0_comp_sem_tab

# ch4_L0C_mods_keep

MOB_IIa_models <- read.table("SEM_models/fxn_example_tables/SEM_Exper8_res_MOB_IIa_lasso.txt", header = T, sep = '\t')
head(MOB_IIa_models)

# demonstrate function
# keep_mods = c('5','11','12','13','15','16')

# selected_MOBIIa <- select_mods_from_sem_object(MOB_IIa_models, keep_mods,
#                                               model_col = "model", form_col = "form") # optional w defaults
# selected_MOBIIa

# MOBIIa_delt_ngt_sem

drop_vars_from_df = function(df, drop_vars){

    df_names <- unlist(names(df))
    keep_df_cols <- df_names[!df_names %in% drop_vars]
    df <- df[keep_df_cols]
    return(df)
}

# demonstrate function -- need more upstream 
# head(selM_CH4_comp)
# drop_vars <- c("formula", "model_out")
# drop_vars_from_df(selM_CH4_comp, drop_vars)

get_formula_col_list_from_df = function(df, fxn_operator = '~ '){
    
    # Get columns containing formulas
    cols = names(df)                                       # get column names

    col_list = list()                                      # empty list
 
    for (i in seq_along(cols))                             # for each column
        col_list[[i]] = unique(str_detect(df[,i], fxn_operator))    # does contain formula operator?

    keep_cols <- unlist(col_list)                          # unlist T/F
    form_cols <- cols[keep_cols]                                        # get columns with formulas
    form_cols
    
    
}

# Test data:
# Test function
sem_models.2 <- update_sem_model_table(model_test, base_model = "ch4_mod0", new_mod_name = "ch4_mod0a.1", 
                      notes = "drop CH4~MOB_I",
                      formula_edits = c("CH4", "CH4_ug_m2_h ~ CH4_ac + CH4_H2 + MOB_IIa + CO2_mg_m2_h"))
                      # output_prefix = "sem."#non_feat_cols = c("model", "base", "notes", "formula"))

# sem_models.2

# Demonstrate function
get_formula_col_list_from_df(sem_models.2, fxn_operator = '~ ')



get_model_formula_from_SEM_results = function(sem_results, set_prefix = '', 
                                          keep_cols = c('model','form'), rename_form = TRUE, fxn_op = " ~"){
    
    model_col <- keep_cols[1]
    form_col <- keep_cols[2]
    
    keep_data <- sem_results[keep_cols]
    keep_data[model_col] <- paste0(set_prefix, ".", keep_data[,model_col])
    
    if (rename_form ==TRUE) {
        
        forms <- as.character(keep_data[,form_col])
        y_vars <- unlist(strsplit(forms, fxn_op))[1]                                        # get y var
        y <- unique(y_vars)
        
        
        names(keep_data)[2] <-y
    }

    return(keep_data)
}

# Demonstrate function
CH4ex7_res2_forms <- get_model_formula_from_SEM_results(semEx7_res2, set_prefix = "ch4.7man")
CH4ex7_res2_forms

MOBIIa_forms <- get_model_formula_from_SEM_results(MOB_IIa_models, set_prefix = "mIIa.man")
MOBIIa_forms

# test data downstream
# MOBIIa_forms <- get_model_formula_from_SEM_results(selM_MOBIIa, set_prefix = "mIIa.man")
# MOBIIa_forms


combine_sem_model_layers = function(left_matrix, right_matrix, 
                                    left_mod = "model", right_mod = "model", out_prefix = "sem.", 
                                    Lform_drop = "formula"){

   # LEFT matrix prep #
    
    # drop any existing "formula" column
    left_matrix <- drop_vars_from_df(left_matrix, Lform_drop)
    
    #left_names <- names(left_matrix)
    #keep_L_cols <- left_names[!left_names %in% Lform_drop]
    #left_matrix <- left_matrix[keep_L_cols]
    

    
    # tag left model colname 
    new_lmod_col = paste0(left_mod, ".l")                                                     # L model name
    colnames(left_matrix)[colnames(left_matrix) == left_mod] <- new_lmod_col                  # mark L model name 

    # repeat left df, n = nrow right df
    left_expand <- left_matrix[rep(seq_len(nrow(left_matrix)), nrow(right_matrix)), ]         # expand left df
    
    
   # RIGHT matrix prep # 
    # tag right model colname
    new_rmod_col = paste0(right_mod, ".r")                                                    # R model name
    colnames(right_matrix)[colnames(right_matrix) == left_mod] <- new_rmod_col                # mark R model name
 
    # repeat right df, n = nrow left df; sort by index
    right_expand <- right_matrix[rep(seq_len(nrow(right_matrix)), nrow(left_matrix)), ]       # expand right df
    right_expand <- right_expand[order(as.numeric(row.names(right_expand))),]                 # sort by index

    
   # COMBINE L & R, make combined IDs and formulae #  
    LR_expanded <- cbind(left_expand, right_expand)                                           # combine L & R expanded 
    model <-  paste0(LR_expanded[,new_lmod_col], "_", LR_expanded[,new_rmod_col])             # combine model ids
    model_out <- paste0(out_prefix, model)
    
    # ID which columns are formulas, get as DF                                                # using newer fxn
    form_cols <- get_formula_col_list_from_df(LR_expanded, fxn_operator = '~')                # id formula cols
    form_data <- LR_expanded[form_cols]                                                       # df of formula cols

    # paste formula components for new combined formulae
    paste_args <- c(form_data, sep="\n")                                                      # do.call paste args
    formula <- do.call(paste, paste_args)                                                     # formulae to DF

    # add model IDs to formula data to finish, make all character
    combined_models <- cbind(model_out, model, form_data, formula) 
    combined_models <- data.frame(lapply(combined_models, as.character),stringsAsFactors = FALSE) 
}


# drop formula column from left matrix !!!!  -- no need, FIXED !!!
# selM_CH4_comp2 <- selM_CH4_comp[,1:7] 
# selM_CH4_comp

# Demonstrate function
CH4mob_comp_forms <- combine_sem_model_layers(ch4_lass0_comp_sem_tab, MOBIIa_forms) 
# CH4mob_comp_forms <- combine_sem_model_layers(selM_CH4_comp, MOBIIa_forms) 
                             #optional: # , left_mod = "model", right_mod = "model", 
                                        #  out_prefix = "sem.", Lform_drop = "formula")

#head(CH4mob_comp_forms)
#names(CH4mob_comp_forms)[6] <- 'MOB_IIa'
#head(CH4mob_comp_forms)

#CH4mob_comp_forms

# original test data (small scale, no fail)

# dim(model_test2)
# model_test2
#nrow(sem_model_table)
# sem_model_table$formula

### a) new TEST data from above

# DO NOT USED lasso SCALED DATA IN SEM !!!   log is good, scaled bad
# FURTHER run_compare needs non-default flags on keep_mod_descr  -- OR CHANGE DEFAULTS in orig. Fxn(s)
# Seems to run at 50% total CPU (part parallel??)

#CH4mob_comp_res <- suppressWarnings(run_compare_sem_models(CH4mob_comp_forms, Guild_CH4_d, estimator = "mlm",
#                                                            keep_mod_descr = c("model_out","model")))
#head(CH4mob_comp_res)

# test data from below

# combine with CH4_ac w base -- expect model 21_249 fails, but function runs/passes results
#CH4_ac_comp_f <- combine_sem_model_layers(CH4_mod_fch4.2_s, CH4_ac2)

# CH4_ac_comp_f2 <- CH4_ac_comp_f[18:19,]
# CH4_ac_comp_f2 <- CH4_ac_comp_f[15:25,]
# CH4_ac_comp_f2 <- CH4_ac_comp_f[18,]# one pass, one fail -- if not dropping fail above
#CH4_ac_comp_f2$model_out <- factor(CH4_ac_comp_f2$model_out)  -- needed?

# working function with new try catch, not parallel
sem_table_run_assign2 = function(sem_mod_table, data, estimator = 'mlm', prefix = "sem.",
                                fixed_x = FALSE, orthogonal = TRUE,
                                print_passing = FALSE, print_failed = TRUE){
    
        iters <- nrow(sem_mod_table)                                                  # use iters for loop
        passing <-list()                                                              # good SEM accumulator


        for(i in seq(1:iters)) {                                                      # for each row              
            skip <- FALSE                                                             # assume pass (for try)
            
            # run SEM model            
            curr_model <- sem_mod_table$formula[i]                                    # get ith sem formula
            sem <- tryCatch(sem(curr_model, data=data, fixed.x=fixed_x, estimator = estimator, 
                                orthogonal = orthogonal), #std.lv=TRUE), -- no diff   # TRY sem
                                error = function(e){skip <<- TRUE})                   # if SEM fails, skip = T  
            if(skip==TRUE){next}                                                      # if SEM failed, next iter
                                   
            # assign mod_out names to models (GLOBAL ENV)
            if(skip == FALSE){                                                        # if SEM pass
                sem_name <- paste(sem_mod_table$model_out[i])                         # get model_out
                assign(sem_name, sem, envir=.GlobalEnv)                               # assign out to global env 
                passing[[i]] <- sem_name                                              # accum model_out name
            }
        }
        passing <- unlist(passing)
    
        # Gather models run and (optional) print pass/fail lists to console
        models_run <- sem_mod_table$model_out
        if (print_passing == TRUE) {
            passing_pr <- paste(passing, collapse =', ')
            print(paste0("ran SEM to create models: ", passing_pr))}
        if (print_failed == TRUE){
            failed <- models_run[!models_run %in% passing]
            print(paste0("failed SEM models: ", failed))}
    
        return(passing)
    }   

# demonstrate function
sem_multi_test <- sem_table_run_assign2(model_test2, data=Guild_CH4_d, estimator ="mlm")
sem_multi_test

# Parallel mc apply version - fix assign to global problem in parallel, ugly not sure how functional try is...
sem_table_run_assign_P = function(sem_mod_table, dataset, estimator = 'mlm', prefix = "sem.", ncores = 20,
                                fixed_x = FALSE, orthogonal = TRUE,
                                print_passing = FALSE, print_failed = TRUE){
    
        iters <- 1:nrow(sem_mod_table)                                                # use iters for loop
        passing <- mclapply(iters, function(i){                                       # multicore apply 

            skip <- FALSE                                                             # assume pass (for try)
            
            # run SEM model, here try catch in case fails            
            curr_model <- sem_mod_table$formula[i]                                    # get ith sem formula
            sem <- tryCatch(sem(curr_model, data=dataset, fixed.x=fixed_x, estimator = estimator, 
                                orthogonal = orthogonal), #std.lv=TRUE), -- no diff   # TRY sem
                                error = function(e){skip <<- TRUE},#)                 # if SEM fails, skip = T  
            if(skip==TRUE){
                next})                                                                # if SEM failed, next iter
                                   
            # assign mod_out names to global env, failed here, just get SEM names
            if(skip == FALSE){                                                        # if SEM pass
                sem_name <- paste(sem_mod_table$model_out[i])                         # get model_out
                c(sem, sem_name)
            }
        }, mc.cores = ncores)                                                         # get n cores for parallel 
   
        # try catch produces nulls & won't assign. Extract list elements, then unlist
        sems <- list()
        sem_names <- list()
        
        for (i in 1:length(passing)){                         
            sems[i] <- (passing[i][[1]][[1]])
            sem_names[i] <- (passing[i][[1]][[2]])
        }
        
        # unlist new lists, no more nulls 
        sems <- unlist(sems)             
        sem_names <- unlist(sem_names)
    
        # assign SEMs to global env by name
        for (i in 1:length(sem_names)){                
                assign(sem_names[[i]], sems[[i]], envir=.GlobalEnv)
            }
            
        # Gather models run and (optional) print pass/fail lists to console
        models_run <- sem_mod_table$model_out
        if (print_passing == TRUE) {
            passing_pr <- paste(sem_names, collapse =', ')
            print(paste0("ran SEM to create models: ", passing_pr))}
    
        if (print_failed == TRUE){
            failed <- models_run[!models_run %in% sem_names]
            print(paste0("failed SEM models: ", failed))}
            else{}

        return(sem_names)
    }   

# test run and assign models from model table
sem_multi_test <- sem_table_run_assign_P(model_test2, data=Guild_CH4_d, estimator ="mlm", ncores = 20)
sem_multi_test; 
compare_sem_results(sem_multi_test, sem_fit_params)

# larger test, should be valid inline
sem_multi_test2 <- sem_table_run_assign_P(CH4mob_comp_forms, data=Guild_CH4_0, estimator ="mlm", ncores=20)
#                                        fixed_x = TRUE, orthogonal = TRUE)
# sem_multi_test2

# show test model summary, not helper output
table <- compare_sem_results(sem_multi_test2, sem_fit_params)
head(table[order(-table$pvalue),])

invisible(sem_table_run_assign_P(CH4mob_comp_forms, data=Guild_CH4_0, estimator ="mlm", ncores=20, print_failed == FALSE))
# test on larger dataset -- same as above?
sem_multi_test2 <- sem_table_run_assign_P(CH4mob_comp_forms, data=Guild_CH4_d, estimator ="mlm", 
                                        fixed_x = TRUE, orthogonal = TRUE)
# sem_multi_test2 -- should be list of models

# show results with model scoring (from assigned)
table <- compare_sem_results(sem_multi_test2, sem_fit_params)
head(table[order(-table$pvalue),])# test on larger dataset, is valid for examples or from below?  
sem_multi_test2 <- sem_table_run_assign_P(CH4_ac_comp_f, data=Guild_CH4_0, estimator ="mlm", 
                                        fixed_x = TRUE, orthogonal = TRUE)
# sem_multi_test2 -- should be list of models

# show results with model scoring (from assigned)
table <- compare_sem_results(sem_multi_test2, sem_fit_params)
head(table[order(-table$pvalue),])

run_compare_sem_models = function(mod_table, data, estimator = "mlm", ncores = 20,
                                  keep_mod_descr = c("model_out","model"), #"base", "notes"),
                                  sem_fit_params = c("pvalue", "chisq", "df", "npar","aic", "bic", 
                                                      "gfi", "cfi", "rni", "rmsea", "srmr"),
                                  print_failed = FALSE){
                                    
    
    # note mod_table requires columns: model, model_out, formula cols 
    
    # get model table features to include in summary
    mod_descr <- data.frame(mod_table[keep_mod_descr])
    mod_descr
    
    # run assign SEM models in table 
    sem_out_list <- sem_table_run_assign_P(mod_table, data, estimator, ncores = ncores, print_failed = print_failed)
    #sem_out_list

    # Get SEM results summary as previously
    fit_summary <- data.frame(compare_sem_results(sem_out_list, sem_fit_params))
    
    
    # Clean data summary
    names(fit_summary)[1] <- "model_out"
    fit_summary$chisq <- round(fit_summary$chisq, digits = 1)
    fit_summary$aic <- round(fit_summary$aic, digits = 1)
    fit_summary$bic <- round(fit_summary$bic, digits = 1)
    #fit_summary
    
    # combine model descriptions and fits
    out <- merge(mod_descr, fit_summary, by = "model_out")
    out
    #mod_descr
    #fit_summary
}

# Demonstrate function
out <- run_compare_sem_models(model_test2, data, estimator = "mlm")
out

# run models, from manual example -- why doesn't this work with 0 data (all???) it did at some pont...
# CH4mob_comp_res <- run_compare_sem_models(CH4mob_comp_forms, Guild_CH4_0, ncores = 20)  
                            #, estimator = "mlm", keep_mod_descr = c("model_out","model")) # defaults
#CH4mob_comp_res

# table <- CH4mob_comp_res
# table[order(-table$pvalue),]
# head(CH4mob_comp_res)
# get summary, full set of models -- MOVE DOWN to IV
#CH4mob_comp_t <- clean_SEM_model_table_for_output(CH4mob_comp_forms, CH4mob_comp_r, 
#                                drop_cols = c("model_out", "formula"), merge_col = c("model"))
#
#CH4mob_comp_t[order(-CH4mob_comp_t$pvalue),]# run models from later composite examples -- below
#CH4_ac_comp_r <- run_compare_sem_models(CH4_ac_comp_f, Guild_CH4_0, ncores = 20)
                            #, estimator = "mlm", keep_mod_descr = c("model_out","model")) # defaults
#CH4_ac_comp_r
# Parllel report: for this ex, run 39 models in CH4_ac_comp_f, data = Guild_CH4_0:
# non par: ~ 50-60% CPU load, 700% CPU), ca 70 sec
# parallel: ~ 90-100% CPU load, ca. 29s 
# get summary, full set of models -- MOVE DOWN to IV
#CH4_ac_comp_t <- clean_SEM_model_table_for_output(CH4_ac_comp_f, CH4_ac_comp_r, 
#                                drop_cols = c("model_out", "formula"), merge_col = c("model"))

#CH4_ac_comp_t[order(-CH4_ac_comp_t$pvalue),]# adding SEM response vars with string defs
CH4 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h'
CH4_H2 <- 'CH4_H2~ SRB_syn + CO2_mg_m2_h'
MOB_IIa <- 'MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac'
MOB_I <- 'MOB_I ~ CH4_H2 +CH4_ac'# Get vector of response vars as strings
Other <- ""
response_vars <- c("CH4", "CH4_H2", "MOB_I","MOB_IIa", "Other")head(sem# Test function
sem_models.2 <- update_sem_model_table(model_test, base_model = "ch4_mod0", new_mod_name = "ch4_mod0a.1", 
                      notes = "drop CH4~MOB_I",
                      formula_edits = c("CH4", "CH4_ug_m2_h ~ CH4_ac + CH4_H2 + MOB_IIa + CO2_mg_m2_h"))
                      # output_prefix = "sem."#non_feat_cols = c("model", "base", "notes", "formula"))

sem_models.2
sem_models.2

# Demonstrate function
sem_models.2_out <- run_compare_sem_models(sem_models.2, data, estimator = "mlm")
sem_models.2_out

remove_yvar_from_df_formula_column =function(df, col_name, fxn_operator = '~ '){
    
    # start function 
    row_data = paste(df[,col_name])                                              # get vect of rows in column

    form_list = list()                                                           # empty list

    for (i in (seq_along(row_data)))     {                                       # for each column
        form_list[i] <- paste(unlist(strsplit(row_data[i], '~ '))[2])            # split, x_vars only to list
        }

    col_replace <- unlist(form_list)                                             # unlist
    df[,col_name] <- col_replace                                                 # replace column data
    return(df)                                                                   # return df

}

# demonstrate function
remove_yvar_from_df_formula_column(sem_models.2, 'MOB_IIa')
#remove_yvar_from_df_formula_column(sem_models.2, 'formula') #-- just drop this one!

remove_yvars_from_model_table_formulae = function(model_table){
    
    form_cols <- get_formula_col_list_from_df(model_table)
    
    for (i in form_cols){
        model_table <- remove_yvar_from_df_formula_column(model_table, i)
    }

    model_table
}    

# demonstrate function
model_test2_out <- remove_yvars_from_model_table_formulae(sem_models.2)
head(model_test2_out)

clean_SEM_model_table_for_output = function(sem_models, sem_output, 
                                            drop_cols = c("model_out", "formula"), merge_col = "model"){
    
    # drop extra columns in sem models, outputs
    sem_models <- drop_vars_from_df(sem_models, drop_cols)
    sem_output <- drop_vars_from_df(sem_output, drop_cols)
    
    #model_cols <- names(sem_models)
    #keep_cols <- model_cols[!model_cols %in% drop_cols]
    #sem_models <- sem_models[keep_cols]
    
    # remove y vars from formulas
    clean_models <- remove_yvars_from_model_table_formulae(sem_models)  # need to drop formula, else
    
    # merge with results
    models_out <- merge(sem_output, clean_models, by = merge_col)       # merge
    models_out[models_out=='NA'] <-''                                   # drop mysterious NAs
    return(models_out)    
    
    # models_out$formula <- sem_models$form
    # repair column names -- probably need a for loop here TODO
    names(models_out) <- gsub(x = names(models_out), pattern = ".x", "R2")
    names(models_out) <- gsub(x = names(models_out), pattern = ".y", "f")
    #return(models_out)  
    
}

# demonstrate function
clean_SEM_model_table_for_output(sem_models.2, sem_models.2_out, 
                                drop_cols = c("model_out", "formula"), merge_col = c("model"))#, 'base', 'notes'))
# demonstrate function - ex post
CH4mob_comp_table <- clean_SEM_model_table_for_output(CH4mob_comp_forms, CH4mob_comp_res, 
                                drop_cols = c("model_out", "formula"), merge_col = c("model"))

table <- CH4mob_comp_table 
head(table[order(-table$pvalue),])

# head(CH4mob_comp_table)CH4mob_comp_table_sig <- CH4mob_comp_table[CH4mob_comp_table$pvalue > 0.05,]
CH4mob_comp_table_sig
# get summary, full set of models -- MOVE DOWN to IV
#CH4_ac_comp_t <- clean_SEM_model_table_for_output(CH4_ac_comp_f, CH4_ac_comp_r, 
#                                drop_cols = c("model_out", "formula"), merge_col = c("model"))

#CH4_ac_comp_t[order(-CH4_ac_comp_t$pvalue),]


########################################################################################
### 12) Workflow wrappers


run_sem_model_summary = function(sem_formula_table, data, sort_by = 'pvalue', sort_cutoff = 0, 
                                 estimator = 'mlm', output = 'summary_table', ncores = 20){
    
    # Run models for results
    sem_results <- suppressWarnings(run_compare_sem_models(sem_formula_table, data, ncores = ncores,
                                            estimator = "mlm", keep_mod_descr = c("model_out","model")))
    
    # Get summary table 
    sem_summary <- clean_SEM_model_table_for_output(sem_formula_table, sem_results, 
                                drop_cols = c("model_out", "formula"), merge_col = c("model"))

    # sort and filter results 
    sem_summary <- sem_summary[order(-sem_summary[sort_by]),]
    sem_summary <- sem_summary[sem_summary[sort_by] >= sort_cutoff,]
    
    #sem_summary <- sem_summary[order(-sem_summary$pvalue),]
    
    # Choose output type
    if (output == 'summary_table'){
        table_out <- sem_summary
    } else {table_out <- sem_results}
    
    return(table_out)
}       


# here add keep full formula, but doesn't seem to work
run_sem_model_summary = function(sem_formula_table, data, sort_by = 'pvalue', sort_cutoff = 0, 
                                 estimator = 'mlm', output = 'summary_table', full_formula = FALSE){
    
    # Run models for results
    sem_results <- suppressWarnings(run_compare_sem_models(sem_formula_table, data, 
                                            estimator = "mlm", keep_mod_descr = c("model_out","model")))
    
    # keep full formula in summary table?
    if(full_formula == FALSE){
        drop_cols = c("model_out", "formula")
        } else {
        drop_cols = c("model_out")
    }
    
    # Get summary table 
    sem_summary <- clean_SEM_model_table_for_output(sem_formula_table, sem_results, 
                                drop_cols = drop_cols, merge_col = c("model"))
                                #drop_cols = c("model_out", "formula"), merge_col = c("model"))

    # sort & filter results 
    sem_summary <- sem_summary[order(-sem_summary[sort_by]),]
    sem_summary <- sem_summary[sem_summary[sort_by] >= sort_cutoff,]
 
        
    # Choose output type
    if (output == 'summary_table'){
        table_out <- sem_summary
    } else {table_out <- sem_results}
    
    return(table_out)
}       
# Demonstrate function
# run_sem_model_summary(sem_models.2, Guild_CH4_0)    #  test fxn
# run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.001)   #  -- ex won't work here


# run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.001, full_formula = TRUE)
# run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.001, output = 'sem_results')
# run_sem_model_summary(CH4_ac_comp_f, Guild_CH4_0)   #  -- ex won't work here

# simple function produces data to run using run_sem_model_summary (above)
# more complex versions could be IO wrappers, run lists or input matrices,
# OR could run at each branch and filter results

sem_tree_formula_build = function(tree_base, branches, run = TRUE, warn_n_models = TRUE) {   # branches should be list of SEM_objects
    
    # add branch formulae to base model
    sem_formula_tree <- combine_sem_model_layers(tree_base, branches[[1]])   
    
    # iterate other branches
    if (length(branches) > 1) {
        for (i in 2:length(branches)){
           sem_formula_tree <- combine_sem_model_layers(sem_formula_tree, branches[[i]])
       }
     }
     #readline(prompt="Press [enter] to continue")
    n_models <- dim(sem_formula_tree)[1]
    print(paste(n_models, "sem models"))
    
    
    sem_formula_tree
     #print(dim(sem_formula_tree[1]))
    
    
}
# Examples below from later in dataset testing, need replacement with generics...last sequential version v.1.09# tree_base <- CH4_mod_fch4.2_s

# branch test data, from lasSEM object to filtered models
# branch_1 <- MOBIIa_mod_fAOANOB_sem$filt_models     # was MOBIIa_fng
# branch_2 <- CH4ac_mod.2_sem$filt_models            # CH4_ac2   
    names(branch_2)[6] <- "CH4_ac"                   #branch_2

# branchz <- list(branch_1)                          # branchz[[1]]
# branchz2 <- list(branch_1, branch_2)
# demonstrate function -- note prints n models in case too many         # 360 models, peak ram = 112 GB!
#sem_formula_tree <-sem_tree_formula_build(tree_base, branchz) # , run = TRUE, warn_n_models = TRUE)
# sem_formula_tree <-sem_tree_formula_build(tree_base, branchz2) # , run = TRUE, warn_n_models = TRUE)
# sem_formula_tree

# test primative fxn
#sem_formula_tree <- combine_sem_model_layers(tree_base, branch_1) 
# sem_formula_tree <- combine_sem_model_layers(tree_base, branchz[[1]])   

# run using SCALED data for some reason or crash, e.g. Guild_CH4_aS
# run it to summary
#run_sem_model_summary(sem_formula_tree, Guild_CH4_0)
# sem_run_test <- run_sem_model_summary(sem_formula_tree, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.00) # pass parallel!!
#run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0)


# branchz

# Get model formula and apply pre-fix from lassSEM
get_sem_models_add_prefix = function(lasso_object, models = "filt_models", prefix = "",
                                  filter = '', cutoff = 0){
    formulae <- lasso_object[[models]]                       # get data from object, could be "sem_models"
    formulae$model <- paste0(prefix, formulae$model)         # append model/ dataset prefix
    
    # filter data as desired (R2, p, etc)
    if (filter != ''){
        formulae <- formulae[formulae[filter] > cutoff,]         
    }
    return(formulae)
}


# Example use  - is from downstream, will need for upstream examples

#CH4_mod_fch4.2 <- get_sem_models_add_prefix(CH4_mod_fch4.2_sem, models = "filt_models", prefix = "fch4.2.")
# CH4_mod_fch4.2_A <- get_sem_models_add_prefix(CH4_mod_fch4.2_sem, models = "sem_models", prefix = "fch4.2.",
#                                          filter = 'R2', cutoff = 0.5)
# CH4_mod_fch4.2
# CH4_mod_fch4.2_sem

rename_formula_col_as_yvar = function(table, formula_col = 'form', fxn_operator = '~'){

    form1 <- paste(table[formula_col][[1]])                                                   # first formula 
    y_var <- unlist(strsplit(form1, fxn_operator))[1]                                         # get y var
    colnames(table)[colnames(table) == formula_col] <- y_var
    table
    #y_var #<- paste0(y_var, fxn_operator)                                                     # add op back
}

# Demonstrate function  -- will need new example, works but data below
# results_table <- CH4_mod_fch4.2_sem$filt_models
# rename_formula_col_as_yvar(results_table, formula_col = 'form', fxn_operator = '~')



########################################################################################
### II) Salinity gradient models (apply above functions
########################################################################################

# helper for dropping vars -- MOVE UP

drop_from_list = function(drop_list, parent_list){
    
    clean_list <- parent_list[!parent_list %in% drop_list]
    clean_list
}

# coef plot of lasso objects, input as str for obj name
plot_lasso_obj_list = function(lasso_list){     
    
    plot_list <- list()

    for (i in lasso_list){ 
        lasso <- get(i)
        plot_list[[i]] <- lasso$coef_plot

    }

    plot_grid(plotlist = plot_list, nrow=1)     
    
}

#mobIIa_lassos <- list("MOBIIa_delt_all", "MOBIIa_delt_ng", "MOBIIa_delt_fAOANOB", "MOBIIa_delt_ngt")

#options(repr.plot.width = 12, repr.plot.height = 3)

#plot_lasso_obj_list(mobIIa_lassos)

# Remove sample names
Guild_CH4_d0 <- Guild_CH4_d
row.names(Guild_CH4_d0) <- Guild_CH4_d0$Sample
Guild_CH4_d0 <- Guild_CH4_d0[,-1]

# SCALE ALL, needed for some functions (SEM screening)
Guild_CH4_dS <- scale(Guild_CH4_d0)              #  careful will scale y 
Guild_CH4_dS <- data.frame(Guild_CH4_dS)          #  make df
Guild_CH4_dS$CH4_logn1 <- Guild_CH4_d0$CH4_logn1  # fix CH4 scale to original

names(Guild_CH4_dS)


# prepare variable SUBSETS 

# get all vars
xvars0 <-names(Guild_CH4_dS)

# vars to exclude (maybe)
drop_gas <- c('CH4_ug_m2_h', 'CH4_logn1', 'CH4_CO2', 'CO2_soilC_mg_g_d') # keep this: ,'CO2_mg_m2_h',
drop_ratios <- c('CN','CP','NP','NO3_NH4','NP_ext')
drop_pw <- c('SO4_pw','Fe_pw','Fe','DOC_mg_L')
drop_guilds <- c('SOxB','MeOB','FeOB','AO_NOB','NOB_AO','mcr_pmo', 'pmo_mcr','Anamx')
drop_taxa <- c('Actino','Chlorf','Firmic')

# Get combined drop lists
drop_vars1 <- c(drop_gas, drop_pw, drop_guilds, drop_taxa, drop_ratios)

# Select guild & chem data:                # Remove drop vars from list
xvars1 <- xvars0[!xvars0 %in% drop_vars1]

# alternately, may want only guild or chem data lists, including ALL

# Def y var
y_var <- 'CH4_logn1'

# DEFINE various X data sets:

# Most VARS IN, no ratios etc.
x_vars <-xvars1                     # x_vars

# ALL VARS IN, w/ ratios, remove CH4 permutations
drop2 <- c(drop_gas, drop_pw)
drop_2cl <- c(drop2, "Cl")                  # also drop Cl!

max_vars <- xvars0[!xvars0 %in% drop2]
#max_vars <- xvars0[!xvars0 %in% drop_2cl]   # also drop Cl!

x_vars0 <- max_vars                                  


x_vars0

# x_CO2_C

options(repr.plot.width=4, repr.plot.height= 3)  # plot dims for coeff plots

# all vars with RATIOS, complete, nothing held back
CH4_delt_ALL <- lasso_mod_plot(y_var, x_vars0, Guild_CH4_dS, guild_colors, plot_title = 'CH4_delta: all vars')

# all vars
CH4_delt_all <- lasso_mod_plot(y_var, x_vars, Guild_CH4_dS, guild_colors, plot_title = 'CH4_delta: all vars')

# force in SO4
CH4_delt_fSO4 <- CH4_delt_all 

# force in SRB
CH4_delt_fSRB <- lasso_mod_plot(y_var, x_vars, Guild_CH4_dS, guild_colors, 
                             force_in = c('SRB'), plot_title = 'CH4_delta: force SRB')

# all vars
CH4_delt_fmobac <- lasso_mod_plot(y_var, x_vars, Guild_CH4_dS, guild_colors, 
                             force_in = c('MOB_IIa, CH4_ac'), plot_title = 'CH4_delta: all vars')


# all vars
CH4_delt_all_sem <- lasso_obj_to_sem_select_plots(CH4_delt_all, 'CH4_logn1', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)

# ALL vars w/ ratios, taxa (min call)
CH4_delt_ALL_sem <- lasso_obj_to_sem_select_plots(CH4_delt_ALL, 'CH4_logn1', Guild_CH4_dS, guild_colors)

# forced SO4
CH4_delt_fSO4_sem <- lasso_obj_to_sem_select_plots(CH4_delt_fSO4, 'CH4_logn1', Guild_CH4_dS, guild_colors,
                                         force_filt = 'SO4', nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)

# forced SRB
CH4_delt_fSRB_sem <- lasso_obj_to_sem_select_plots(CH4_delt_fSRB, 'CH4_logn1', Guild_CH4_dS, guild_colors,
                                         force_filt = 'SRB', nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)

                                         # see more models   
                                         #force_filt = 'SRB', nbest = 5, keep_filt = 0.5, resid_0filt = FALSE,  # sem filts
                                         #heat_filt = 0.5, k_clusts = 5)

# forced SRB
CH4_delt_f_mobac_sem <- lasso_obj_to_sem_select_plots(CH4_delt_fmobac, 'CH4_logn1', Guild_CH4_dS, guild_colors,
                                         force_filt = c('MOB_IIa','CH4_ac'), nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)

options(repr.plot.width=12, repr.plot.height= 3) # plot dims for combined plots

# CH4_delt_all_sem$plot

options(repr.plot.width=12, repr.plot.height= 3) # plot dims for combined plots

CH4_delt_ALL_sem$plot
CH4_delt_all_sem$plot
CH4_delt_fSO4_sem$plot
CH4_delt_fSRB_sem$plot 

CH4_delt_f_mobac_sem$plot

# Def y var
y_var <- 'CH4_logn1'

# DEFINE various X data sets:

# Most VARS IN, no ratios etc.
#x_vars <-xvars1                     # x_vars

# ALL VARS IN, w/ ratios, remove CH4 permutations
#max_vars <- xvars0[!xvars0 %in% drop_gas]
#x_vars0 <- max_vars                                  

# no CO2
CO2 <- 'CO2_mg_m2_h'
x_noCO2 <- x_vars[!x_vars %in% CO2]

# no CO2.2 - drop MOBII, wrong weight
x_noCO2.2 <- drop_from_list('MOB_II', x_noCO2)

# replace CO2 w CO2_C
x_CO2_C <- c(x_noCO2, 'CO2_soilC_mg_g_d')

x_noCO2.2
# x_CO2_C

options(repr.plot.width=4, repr.plot.height= 3)  # plot dims for coeff plots

# no CO2
CH4_delt_xCO2 <- lasso_mod_plot(y_var, x_noCO2, Guild_CH4_dS, guild_colors, plot_title = 'CH4_delta: noCO2')

# no CO2.2 -- drop MOBII
CH4_delt_xCO2.2 <- lasso_mod_plot(y_var, x_noCO2.2, Guild_CH4_dS, guild_colors, plot_title = 'CH4_delta: noCO2')

# CO2/C instead of CO2
CH4_delt_CO2C <- lasso_mod_plot(y_var, x_CO2_C, Guild_CH4_dS, guild_colors, plot_title = 'CH4_delta: CO2/C')

# CH4_delt_xCO2$coef_plot

 #CH4_delt_CO2C$coef_plot

# no CO2 - had to drop heat filt down to 0.65, NO resid filter!
CH4_delt_xCO2_sem <- lasso_obj_to_sem_select_plots(CH4_delt_xCO2, 'CH4_logn1', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = FALSE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

# no CO2.2 - had to drop heat filt down to 0.65
CH4_delt_xCO2.2_sem <- lasso_obj_to_sem_select_plots(CH4_delt_xCO2.2, 'CH4_logn1', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = FALSE,  # sem filts
                                         heat_filt = 0.5, k_clusts = 4)

# CO2/C
CH4_delt_CO2C_sem <- lasso_obj_to_sem_select_plots(CH4_delt_CO2C, 'CH4_logn1', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 4)

options(repr.plot.width=12, repr.plot.height= 3) # plot dims for combined plots

# CH4_delt_ALL_sem$plot
CH4_delt_xCO2_sem$plot
CH4_delt_xCO2.2_sem$plot 
CH4_delt_CO2C_sem$plot 

#x_vars <- x_vars[!x_vars %in% 'Cl']   # DROP Cl
x_vars <- x_vars[!x_vars %in% 'Salinity.x']   # DROP Cl
x_vars

# Initialize x vars for MOB_IIa
MOBIIa_drop.1 <- c('CO2_mg_m2_h','MOB_IIa','SRB_syn','SRB','FeRB')
x_MOBIIa.1 <- drop_from_list(MOBIIa_drop.1 , x_vars)                   # init list

MOBIIa_drop.2 <- c('MOB_I','MOB_II')                                   # corr not mech (wrong sign)
x_MOBIIa.2 <- drop_from_list(MOBIIa_drop.2, x_MOBIIa.1)

MOBIIa_drop.3 <- c('CH4_ac','CH4_H2')                                  # corr not mech (wrong sign)
x_MOBIIa.3 <- drop_from_list(MOBIIa_drop.3, x_MOBIIa.2)

MOBIIa_drop.4 <- c('N', 'C', 'P')                                      # overwhelm signal, sequential
x_MOBIIa.4 <- drop_from_list(MOBIIa_drop.4, x_MOBIIa.3)

# all vars
MOBIIa_delt_all <- lasso_mod_plot('MOB_IIa', x_MOBIIa.1, Guild_CH4_dS, guild_colors, plot_title = 'MOBIIa_delta1: all vars')

# drop other guilds
MOBIIa_delt_ng <- lasso_mod_plot('MOB_IIa', x_MOBIIa.3, Guild_CH4_dS, guild_colors, plot_title = 'MOBIIa_delta3: no_troph/gen')

# force in AOA, AOB to ng
MOBIIa_delt_fAOANOB <- lasso_mod_plot('MOB_IIa', x_MOBIIa.3, Guild_CH4_dS, guild_colors, 
                             force_in = c('AOA', 'NOB'), plot_title = 'MOBIIa_delta3: force AOA, NOB')

# drop soil totals
MOBIIa_delt_ngt <- lasso_mod_plot('MOB_IIa', x_MOBIIa.4, Guild_CH4_dS, guild_colors, plot_title = 'MOBIIa_delta4: ng/totals')

# coef plot of lasso objects, input as str for obj name
plot_lasso_obj_list = function(lasso_list){     
    
    plot_list <- list()

    for (i in lasso_list){ 
        lasso <- get(i)
        plot_list[[i]] <- lasso$coef_plot

    }

    plot_grid(plotlist = plot_list, nrow=1)     
    
}

mobIIa_lassos <- list("MOBIIa_delt_all", "MOBIIa_delt_ng", "MOBIIa_delt_fAOANOB", "MOBIIa_delt_ngt")

options(repr.plot.width = 12, repr.plot.height = 3)

plot_lasso_obj_list(mobIIa_lassos)

#"MOBIIa_delt_all", "MOBIIa_delt_ng", "MOBIIa_delt_fAOANOB", "MOBIIa_delt_ngt"

# MOB_IIa no CH4 guilds
MOBIIa_delt_ng_sem <- lasso_obj_to_sem_select_plots(MOBIIa_delt_ng, 'MOB_IIa', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

# MOB_IIa no CH4 guilds
MOBIIa_delt_fAOANOB_sem <- lasso_obj_to_sem_select_plots(MOBIIa_delt_fAOANOB, 'MOB_IIa', Guild_CH4_dS, guild_colors,
                                         force_filt = c('AOA','NOB'), nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 1)

# MOB_IIa no CH4 guilds
MOBIIa_delt_ngt_sem <- lasso_obj_to_sem_select_plots(MOBIIa_delt_ngt, 'MOB_IIa', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 1)

MOBIIa_delt_ng_sem$plot
MOBIIa_delt_fAOANOB_sem$plot
MOBIIa_delt_ngt_sem$plot

#x_vars
#x_vars <- x_vars[!x_vars %in% 'Cl']   # DROP Cl
x_vars <- x_vars[!x_vars %in% 'Salinity.x']   # DROP Cl
x_vars

# Initialize x vars for CH4_ac
CH4ac_drop.1 <- c('CO2_mg_m2_h','CH4_ac','CH4_H2','MOB_I','MOB_II','MOB_IIa')       # first drop CH4 related
x_CH4ac.1 <- drop_from_list(CH4ac_drop.1, x_vars)                   # init list

# next drop N cycling microbes
CH4ac_drop.2 <- c('AOA', 'AOB', 'NOB')                                 
x_CH4ac.2 <- drop_from_list(CH4ac_drop.2, x_CH4ac.1) 

# then add chem ratios
#x_CH4ac.3 <- c(x_CH4ac.2, drop_ratios)                        - didn't work, unlist??
x_CH4ac.3 <- c(x_CH4ac.2,'CN','NP','NO3_NH4','NP_ext')                             # NO CP in data!!
 
# last, add some PW components
x_CH4ac.4 <- c(x_CH4ac.3, "DOC_mg_L", "Fe_pw")
#x_CH4ac.3 <- c('Bulk_dens','H2O_FPS','pH','C','Salinity.x','SO4','SRB_syn','SRB','FeRB', 'DOC_mg_L','Fe_pw')

# drop hard to explain vars

CH4ac_drop.5 <- c('NP', 'N')   
x_CH4ac.5 <- drop_from_list(CH4ac_drop.5, x_CH4ac.4)


x_CH4ac.3

names(Guild_CH4_dS)

CH4ac_delta.1 <- lasso_mod_plot('CH4_ac', x_CH4ac.1, Guild_CH4_dS, guild_colors, plot_title = 'CH4ac_delta1: init vars')
CH4ac_delta.2 <- lasso_mod_plot('CH4_ac', x_CH4ac.2, Guild_CH4_dS, guild_colors, plot_title = 'CH4ac_delta2: no N guild')

CH4ac_delta.3 <- lasso_mod_plot('CH4_ac', x_CH4ac.3, Guild_CH4_dS, guild_colors, plot_title = 'CH4ac_delta3: add ratios')
CH4ac_delta.4 <- lasso_mod_plot('CH4_ac', x_CH4ac.4, Guild_CH4_dS, guild_colors, plot_title = 'CH4ac_delta4: add pw')

# Force SRB, FeRB into x4 (max vars)
CH4ac_d_fSrFeRB <- lasso_mod_plot('MOB_IIa', x_CH4ac.5, Guild_CH4_dS, guild_colors, 
                             force_in = c('FeRB','SRB_syn'), plot_title = 'CH4ac_delta4_fSrFeRB')

options(repr.plot.width = 15, repr.plot.height = 3)
ch4ac_lassos = c('CH4ac_delta.1','CH4ac_delta.2','CH4ac_delta.3','CH4ac_delta.4', 'CH4ac_d_fSrFeRB')

plot_lasso_obj_list(ch4ac_lassos)
options(repr.plot.width = 4, repr.plot.height = 3)
CH4ac_delta.1$coef_plot
CH4ac_delta.2$coef_plot
CH4ac_delta.3$coef_plot
CH4ac_delta.4$coef_plot 
# MOB_IIa no CH4 guilds
CH4ac_delta.1_sem <- lasso_obj_to_sem_select_plots(CH4ac_delta.1, 'CH4_ac', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

CH4ac_delta.2_sem <- lasso_obj_to_sem_select_plots(CH4ac_delta.2, 'CH4_ac', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

CH4ac_delta.3_sem <- lasso_obj_to_sem_select_plots(CH4ac_delta.3, 'CH4_ac', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

CH4ac_delta.4_sem <- lasso_obj_to_sem_select_plots(CH4ac_delta.4, 'CH4_ac', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

CH4ac_delta_f_sem <- lasso_obj_to_sem_select_plots(CH4ac_d_fSrFeRB, 'CH4_ac', Guild_CH4_dS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)


#CH4ac_delta_f_sem <- lasso_obj_to_sem_select_plots(CH4ac_d_fSrFeRB, 'CH4_ac', Guild_CH4_dS, guild_colors,
#                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
#                                         heat_filt = 0.65, k_clusts = 4)



CH4ac_delta.1_sem$plot
CH4ac_delta.2_sem$plot
CH4ac_delta.3_sem$plot
CH4ac_delta.4_sem$plot
CH4ac_delta_f_sem$plot





### i) select input data

### ii) create lasso objects, coef plots for diff scenarios

### iii) get SEM models, filtered, lists w plots

### iv) get model comparison plots





# Remove sample names
Guild_CH4_0 <- Guild_CH4
row.names(Guild_CH4_0) <- Guild_CH4_0$Sample
Guild_CH4_0 <- Guild_CH4_0[,-1]

# SCALE ALL, needed for some functions (SEM screening)
Guild_CH4_aS <- scale(Guild_CH4_0)              #  careful will scale y 
Guild_CH4_aS <- data.frame(Guild_CH4_aS)          #  make df
Guild_CH4_aS$CH4_logn1 <- Guild_CH4_0$CH4_logn1  # fix CH4 scale to original

names(Guild_CH4_aS)
# prepare variable SUBSETS -- SAME AS ABOVE! 

# get all vars
xvars0 <-names(Guild_CH4_dS)

# vars to exclude (maybe)
drop_gas <- c('CH4_ug_m2_h', 'CH4_logn1', 'CH4_CO2', 'CO2_soilC_mg_g_d') # keep this: ,'CO2_mg_m2_h',
drop_ratios <- c('CN','CP','NP','NO3_NH4','NP_ext')
drop_pw <- c('SO4_pw','Fe_pw','Fe','DOC_mg_L')
drop_guilds <- c('SOxB','MeOB','FeOB','AO_NOB','NOB_AO','mcr_pmo', 'pmo_mcr','Anamx')
drop_taxa <- c('Actino','Chlorf','Firmic')

# Get combined drop lists
drop_vars1 <- c(drop_gas, drop_pw, drop_guilds, drop_taxa, drop_ratios)# Select guild & chem data:                # Remove drop vars from list
xvars1 <- xvars0[!xvars0 %in% drop_vars1]

# alternately, may want only guild or chem data lists, including ALL
# Def y var
y_var <- 'CH4_logn1'

# DEFINE various X data sets:

# Most VARS IN, no ratios etc.
x_vars <-xvars1                     # x_vars

# ALL VARS IN, w/ ratios, remove CH4 permutations
drop2 <- c(drop_gas, drop_pw)

max_vars <- xvars0[!xvars0 %in% drop2]
x_vars0 <- max_vars                                  


# then add chem ratios
#x_CH4ac.3 <- c(x_CH4ac.2, drop_ratios)                        - didn't work, unlist??
x_CH4all.2 <- c(x_vars,'CN','NP','NO3_NH4','NP_ext') 

x_CH4all.2

#x_vars0

# x_CO2_C

options(repr.plot.width=4, repr.plot.height= 3)  # plot dims for coeff plots

# all vars with RATIOS, complete, nothing held back
CH4_mod_ALL <- lasso_mod_plot(y_var, x_vars0, Guild_CH4_aS, guild_colors, plot_title = 'CH4_delta: all vars')

# all vars (ratios, no taxa)
CH4_mod_ALL.2 <- lasso_mod_plot(y_var, x_CH4all.2, Guild_CH4_aS, guild_colors, plot_title = 'CH4_delta: all vars')

# all vars
CH4_mod_all <- lasso_mod_plot(y_var, x_vars, Guild_CH4_aS, guild_colors, plot_title = 'CH4_delta: all vars')

# force in SO4
CH4_mod_fSO4 <- CH4_mod_all 

# force in SRB
CH4_mod_fSRB <- lasso_mod_plot(y_var, x_vars, Guild_CH4_aS, guild_colors, 
                             force_in = c('SRB'), plot_title = 'CH4_delta: force SRB')


CH4_mod_fsal.1 <- lasso_mod_plot(y_var, x_vars0, Guild_CH4_aS, guild_colors, 
                             force_in = c('Salinity.x'), plot_title = 'CH4_delta: force Salin')

CH4_mod_fsal.2 <- lasso_mod_plot(y_var, x_CH4all.2, Guild_CH4_aS, guild_colors, 
                             force_in = c('Salinity.x', 'MOB_IIa', 'CH4_ac'), plot_title = 'CH4_delta: force Salin.2')

CH4_mod_fch4.2 <- lasso_mod_plot(y_var, x_CH4all.2, Guild_CH4_aS, guild_colors, 
                            plot_title = 'CH4_delta: filt ch4.2')



CH4all_lassos <- c('CH4_mod_ALL', 'CH4_mod_ALL.2', 'CH4_mod_all', 'CH4_mod_fSO4') 
CH4all_lassos2 <- c('CH4_mod_fSRB','CH4_mod_fsal.1','CH4_mod_fsal.2','CH4_mod_fch4.2')

options(repr.plot.width = 12, repr.plot.height = 3)

plot_lasso_obj_list(CH4all_lassos)
plot_lasso_obj_list(CH4all_lassos2)

# Resulting output lassSEM objects:
# c(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSO4_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem,  CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)

# all vars
CH4_mod_all_sem <- lasso_obj_to_sem_select_plots(CH4_mod_all, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)

# ALL vars w/ ratios, taxa (min call)
CH4_mod_ALL_sem <- lasso_obj_to_sem_select_plots(CH4_mod_ALL, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         nbest = 4, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)

# ALL vars w/ ratios, taxa (min call)
CH4_mod_ALL.2_sem <- lasso_obj_to_sem_select_plots(CH4_mod_ALL.2, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         nbest = 4, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)
# forced SO4
CH4_mod_fSO4_sem <- lasso_obj_to_sem_select_plots(CH4_mod_fSO4, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         force_filt = 'SO4', nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)
# forced SRB
CH4_mod_fSRB_sem <- lasso_obj_to_sem_select_plots(CH4_mod_fSRB, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         force_filt = 'SRB', nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 3)

                                         # see more models   
                                         #force_filt = 'SRB', nbest = 5, keep_filt = 0.5, resid_0filt = FALSE,  # sem filts
                                         #heat_filt = 0.5, k_clusts = 5)

# forced Salin
CH4_mod_fsal.1_sem <- lasso_obj_to_sem_select_plots(CH4_mod_fsal.1, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         force_filt = 'Salinity.x', nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.6, k_clusts = 3)

# forced Salin2
CH4_mod_fsal.2_sem <- lasso_obj_to_sem_select_plots(CH4_mod_fsal.2, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         force_filt = c('Salinity.x', 'CH4_ac'), nbest = 3, keep_filt = 0.4, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.4, k_clusts = 3)

# forced Salin2
CH4_mod_fch4.2_sem <- lasso_obj_to_sem_select_plots(CH4_mod_fch4.2, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         force_filt = c('CH4_ac','MOB_IIa','Bulk_dens'), 
                                         nbest = 6, keep_filt = 0.4, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.69, k_clusts = 3)

options(repr.plot.width=12, repr.plot.height= 3) # plot dims for combined plots

# CH4_delt_all_sem$plot

options(repr.plot.width=12, repr.plot.height= 3) # plot dims for combined plots

CH4_mod_ALL_sem$plot
CH4_mod_ALL.2_sem$plot
CH4_mod_all_sem$plot
#CH4_mod_fSO4_sem$plot
CH4_mod_fSRB_sem$plot 


CH4_mod_fsal.1_sem$plot

CH4_mod_fsal.2_sem$plot

CH4_mod_fch4.2_sem$plot



# Def y var
y_var <- 'CH4_logn1'

# DEFINE various X data sets:

# Most VARS IN, no ratios etc.
#x_vars <-xvars1                     # x_vars

# ALL VARS IN, w/ ratios, remove CH4 permutations
#max_vars <- xvars0[!xvars0 %in% drop_gas]
#x_vars0 <- max_vars                                  

# no CO2
CO2 <- 'CO2_mg_m2_h'
x_noCO2 <- x_vars[!x_vars %in% CO2]

# no CO2.2 - drop MOBII, wrong weight
x_noCO2.2 <- drop_from_list('MOB_II', x_noCO2)

## CO2_C 
# replace CO2 w CO2_C
x_CO2_C <- c(x_noCO2, 'CO2_soilC_mg_g_d')

# CO2_C with ratios
x_CO2_C.2 <- c(x_CH4all.2, 'CO2_soilC_mg_g_d')
x_CO2_C.2 <- x_CO2_C.2[!x_CO2_C.2 %in% CO2]


x_CO2_C.2

x_noCO2.2
# x_CO2_C

options(repr.plot.width=4, repr.plot.height= 3)  # plot dims for coeff plots

# no CO2
CH4_mod_xCO2 <- lasso_mod_plot(y_var, x_noCO2, Guild_CH4_aS, guild_colors, plot_title = 'CH4_all: noCO2')

# no CO2.2 -- drop MOBII
CH4_mod_xCO2.2 <- lasso_mod_plot(y_var, x_noCO2.2, Guild_CH4_aS, guild_colors, plot_title = 'CH4_all: noCO2')

# CO2/C instead of CO2
CH4_mod_CO2C <- lasso_mod_plot(y_var, x_CO2_C, Guild_CH4_aS, guild_colors, plot_title = 'CH4_all: CO2/C')

# # CO2/C w. soil ratios
CH4_mod_CO2C.2 <- lasso_mod_plot(y_var, x_CO2_C.2, Guild_CH4_aS, guild_colors, plot_title = 'CH4_all: CO2/C ratios')


# CH4_delt_xCO2$coef_plot

 #CH4_delt_CO2C$coef_plot

# no CO2 - had to drop heat filt down to 0.65, NO resid filter!
CH4_mod_xCO2_sem <- lasso_obj_to_sem_select_plots(CH4_mod_xCO2, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = FALSE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

# no CO2.2 - had to drop heat filt down to 0.65
CH4_mod_xCO2.2_sem <- lasso_obj_to_sem_select_plots(CH4_mod_xCO2.2, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = FALSE,  # sem filts
                                         heat_filt = 0.5, k_clusts = 4)

# CO2/C
CH4_mod_CO2C_sem <- lasso_obj_to_sem_select_plots(CH4_mod_CO2C, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 4)

CH4_mod_CO2C.2_sem <- lasso_obj_to_sem_select_plots(CH4_mod_CO2C.2, 'CH4_logn1', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.7, k_clusts = 4)

options(repr.plot.width=12, repr.plot.height= 3) # plot dims for combined plots

# CH4_delt_ALL_sem$plot
CH4_mod_xCO2_sem$plot
CH4_mod_xCO2.2_sem$plot 
CH4_mod_CO2C_sem$plot 
CH4_mod_CO2C.2_sem$plot 

# Initialize x vars for MOB_IIa
MOBIIa_drop.1 <- c('CO2_mg_m2_h','MOB_IIa','SRB_syn','SRB','FeRB')
x_MOBIIa.1 <- drop_from_list(MOBIIa_drop.1 , x_vars)                   # init list

MOBIIa_drop.2 <- c('MOB_I','MOB_II')                                   # corr not mech (wrong sign)
x_MOBIIa.2 <- drop_from_list(MOBIIa_drop.2, x_MOBIIa.1)

MOBIIa_drop.3 <- c('CH4_ac','CH4_H2')                                  # corr not mech (wrong sign)
x_MOBIIa.3 <- drop_from_list(MOBIIa_drop.3, x_MOBIIa.2)

MOBIIa_drop.4 <- c('N', 'C', 'P')                                      # overwhelm signal, sequential
x_MOBIIa.4 <- drop_from_list(MOBIIa_drop.4, x_MOBIIa.3)

# all vars
MOBIIa_mod_all <- lasso_mod_plot('MOB_IIa', x_MOBIIa.1, Guild_CH4_aS, guild_colors, plot_title = 'MOBIIa_mod.1: all vars')

# drop other guilds
MOBIIa_mod_ng <- lasso_mod_plot('MOB_IIa', x_MOBIIa.3, Guild_CH4_aS, guild_colors, plot_title = 'MOBIIa_mod.3: no_troph/gen')

# force in AOA, AOB to ng
MOBIIa_mod_fAOANOB <- lasso_mod_plot('MOB_IIa', x_MOBIIa.3, Guild_CH4_aS, guild_colors, 
                             force_in = c('AOA', 'NOB'), plot_title = 'MOBIIa_mod.3: force AOA, NOB')

# drop soil totals
MOBIIa_mod_ngt <- lasso_mod_plot('MOB_IIa', x_MOBIIa.4, Guild_CH4_aS, guild_colors, plot_title = 'MOBIIa_mod.4: ng/totals')

mobIIa_lassos <- list("MOBIIa_mod_all", "MOBIIa_mod_ng", "MOBIIa_mod_fAOANOB", "MOBIIa_mod_ngt")

options(repr.plot.width = 12, repr.plot.height = 3)

plot_lasso_obj_list(mobIIa_lassos)

#"MOBIIa_delt_all", "MOBIIa_delt_ng", "MOBIIa_delt_fAOANOB", "MOBIIa_delt_ngt"

# MOB_IIa no CH4 guilds
MOBIIa_mod_ng_sem <- lasso_obj_to_sem_select_plots(MOBIIa_mod_ng, 'MOB_IIa', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

# MOB_IIa no CH4 guilds
MOBIIa_mod_fAOANOB_sem <- lasso_obj_to_sem_select_plots(MOBIIa_mod_fAOANOB, 'MOB_IIa', Guild_CH4_aS, guild_colors,
                                         force_filt = c('AOA','NOB'), nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

# MOB_IIa no CH4 guilds
MOBIIa_mod_ngt_sem <- lasso_obj_to_sem_select_plots(MOBIIa_mod_ngt, 'MOB_IIa', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

MOBIIa_mod_ng_sem$plot
MOBIIa_mod_fAOANOB_sem$plot
MOBIIa_mod_ngt_sem$plot

x_vars

# Initialize x vars for CH4_ac
CH4ac_drop.1 <- c('CO2_mg_m2_h','CH4_ac','CH4_H2','MOB_I','MOB_II','MOB_IIa')       # first drop CH4 related
x_CH4ac.1 <- drop_from_list(CH4ac_drop.1, x_vars)                   # init list

# next drop N cycling microbes
CH4ac_drop.2 <- c('AOA', 'AOB', 'NOB', 'Olsen_P')                                 
x_CH4ac.2 <- drop_from_list(CH4ac_drop.2, x_CH4ac.1) 

# then add chem ratios
#x_CH4ac.3 <- c(x_CH4ac.2, drop_ratios)                        - didn't work, unlist??
x_CH4ac.3 <- c(x_CH4ac.2,'CN','NP','NO3_NH4','NP_ext')                             # NO CP in data!!
 
# last, add some PW components
x_CH4ac.4 <- c(x_CH4ac.3, "DOC_mg_L", "Fe_pw")
#x_CH4ac.3 <- c('Bulk_dens','H2O_FPS','pH','C','Salinity.x','SO4','SRB_syn','SRB','FeRB', 'DOC_mg_L','Fe_pw')

# drop hard to explain vars

CH4ac_drop.5 <- c('NP', 'N')   
x_CH4ac.5 <- drop_from_list(CH4ac_drop.5, x_CH4ac.4)


x_CH4ac.3

names(Guild_CH4_dS)

CH4ac_mod.1 <- lasso_mod_plot('CH4_ac', x_CH4ac.1, Guild_CH4_aS, guild_colors, plot_title = 'CH4ac_mod1: init vars')
CH4ac_mod.2 <- lasso_mod_plot('CH4_ac', x_CH4ac.2, Guild_CH4_aS, guild_colors, plot_title = 'CH4ac_mod2: no N guild')

CH4ac_mod.3 <- lasso_mod_plot('CH4_ac', x_CH4ac.3, Guild_CH4_aS, guild_colors, plot_title = 'CH4ac_mod3: add ratios')
CH4ac_mod.4 <- lasso_mod_plot('CH4_ac', x_CH4ac.4, Guild_CH4_aS, guild_colors, plot_title = 'CH4ac_mod4: add pw')

# Force SRB, FeRB into x4 (max vars)
CH4ac_a_fSrFeRB <- lasso_mod_plot('MOB_IIa', x_CH4ac.5, Guild_CH4_aS, guild_colors, 
                             force_in = c('FeRB','SRB_syn'), plot_title = 'CH4ac_mod4_fSrFeRB')

options(repr.plot.width = 15, repr.plot.height = 3)
ch4ac_lassos = c('CH4ac_mod.1','CH4ac_mod.2','CH4ac_mod.3','CH4ac_mod.4', 'CH4ac_a_fSrFeRB')

plot_lasso_obj_list(ch4ac_lassos)
options(repr.plot.width = 4, repr.plot.height = 3)
CH4ac_delta.1$coef_plot
CH4ac_delta.2$coef_plot
CH4ac_delta.3$coef_plot
CH4ac_delta.4$coef_plot 
# MOB_IIa no CH4 guilds
CH4ac_mod.1_sem <- lasso_obj_to_sem_select_plots(CH4ac_mod.1, 'CH4_ac', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

CH4ac_mod.2_sem <- lasso_obj_to_sem_select_plots(CH4ac_mod.2, 'CH4_ac', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

CH4ac_mod.3_sem <- lasso_obj_to_sem_select_plots(CH4ac_mod.3, 'CH4_ac', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

CH4ac_mod.4_sem <- lasso_obj_to_sem_select_plots(CH4ac_mod.4, 'CH4_ac', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.65, k_clusts = 4)

CH4ac_mod_f_sem <- lasso_obj_to_sem_select_plots(CH4ac_a_fSrFeRB, 'CH4_ac', Guild_CH4_aS, guild_colors,
                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
                                         heat_filt = 0.6, k_clusts = 3)


#CH4ac_delta_f_sem <- lasso_obj_to_sem_select_plots(CH4ac_d_fSrFeRB, 'CH4_ac', Guild_CH4_dS, guild_colors,
#                                         nbest = 3, keep_filt = 0.5, resid_0filt = TRUE,  # sem filts
#                                         heat_filt = 0.65, k_clusts = 4)



CH4ac_mod.1_sem$plot
CH4ac_mod.2_sem$plot
CH4ac_mod.3_sem$plot
CH4ac_mod.4_sem$plot
CH4ac_mod_f_sem$plot





CH4_gen <- c("CH4_ac", "CH4_H2", "CO2_mg_m2_h", "Salinity.x", "SO4", "SRB", "FeRB")#, "NP")

CH4_ox <- c("MOB_IIa", "MOB_I", "MOB_II", "Bulk_dens", "H2O_FPS", 
            "N", "AOB", "AOA", "NOB", "NO3_N", "NH4_N", "NO3_NH4")

composite_sets <- c("CH4_gen", "CH4_ox")

# test set - all data, forced CH4_ac & MOB_IIa

# rename models and "model" row prefix
CH4_mod_fch4.2_m <- CH4_mod_fch4.2_sem$filt_models                   # get filtered models
CH4_mod_fch4.2_m$model <- paste0("ch4f2." ,CH4_mod_fch4.2_m$model)
#CH4_mod_fch4.2_m


# all models, not just passing filter
#CH4_mod_fch4.2_a <- CH4_mod_fch4.2_sem$sem_models                   # get filtered models
#CH4_mod_fch4.2_a

# Get composite models
#CH4_mod_fch4.2_c <- get_composites_model_table(CH4_mod_fch4.2_m, composite_sets, assign_other = "y_var",#)"CH4_gen", 
#                           no_int = FALSE, compos_fail_collapse_all = TRUE)

CH4_mod_fch4.2_c <- get_composites_model_table(CH4_mod_fch4.2_m, composite_sets, assign_other = "CH4_gen", 
                           no_int = TRUE, compos_fail_collapse_all = TRUE)


# CH4_mod_fch4.2_c#$formula

# run composite models
CH4_mod_fch4.2_r <- suppressWarnings(run_compare_sem_models(CH4_mod_fch4.2_c, Guild_CH4_0, 
                                            estimator = "mlm", keep_mod_descr = c("model_out","model")))

# show composite models
# CH4_mod_fch4.2_r

# Get summary table and sort by p-value
CH4_mod_fch4.2_s <- clean_SEM_model_table_for_output(CH4_mod_fch4.2_c, CH4_mod_fch4.2_r, 
                                 drop_cols = c("model_out", "formula"), merge_col = c("model"))

head(CH4_mod_fch4.2_s[order(-CH4_mod_fch4.2_s$pvalue),])


# models are here:
CH4_mod_xCO2_sem$plot
CH4_mod_xCO2.2_sem$plot 
CH4_mod_CO2C_sem$plot 
CH4_mod_CO2C.2_sem$plot 
# Example: get filtered models from SEM object
CH4_mod_CO2C.2_fm <- CH4_mod_CO2C.2_sem$filt_models                   # get filtered models
CH4_mod_CO2C.2_a <- CH4_mod_CO2C.2_sem$sem_models                   # get filtered models
CH4_mod_CO2C.2_fm6 <- CH4_mod_CO2C.2_fm[CH4_mod_CO2C.2_fm$n_pred <7,] # get mods < 7 vars


# should make FXN HERE: 
# label model prefix
mod_prefix <- "ch4co2C2"
CH4_mod_CO2C.2_fm$model <- paste(mod_prefix, CH4_mod_CO2C.2_fm$model, sep = ".")
CH4_mod_CO2C.2_fm6$model <- paste(mod_prefix, CH4_mod_CO2C.2_fm6$model, sep = ".")

head(CH4_mod_CO2C.2_fm)
#CH4_mod_CO2C.2_a

CH4_mod_CO2C.2_c <- get_composites_model_table(CH4_mod_CO2C.2_fm, composite_sets, assign_other = "CH4_gen", 
                           compos_fail_collapse_all = TRUE, no_int = TRUE)

# cCH4_mod_CO2C.2_fm6 <- get_composites_model_table(CH4_mod_CO2C.2_fm, composite_sets, assign_other = "y_var", 
#                            compos_fail_collapse_all = TRUE, no_int = TRUE)


#cCH4_mod_CO2C.2_fm6 <- cCH4_mod_CO2C.2_fm6[1:7,] # test all comp pass subset... some mods not running


#cCH4_mod_CO2C.2_fm6
# x <-cCH4_mod_CO2C.2_fm6$comp_failed
# str(x)
# length(x)
# as.character(cCH4_mod_CO2C.2_fm6$comp_failed)
cCH4_mod_CO2C.2_fm6_p <- cCH4_mod_CO2C.2_fm6[as.character(cCH4_mod_CO2C.2_fm6$comp_failed) == 'character(0)',]
cCH4_mod_CO2C.2_fm6_p
# Run models for results -- not matching above!!
CH4_mod_CO2C.2_r <- suppressWarnings(run_compare_sem_models(CH4_mod_CO2C.2_c, Guild_CH4_0, 
                                            estimator = "mlm", keep_mod_descr = c("model_out","model")))

head(CH4_mod_CO2C.2_r)

CH4_mod_CO2C.2_s <- clean_SEM_model_table_for_output(CH4_mod_CO2C.2_c, CH4_mod_CO2C.2_r, 
                                 drop_cols = c("model_out", "formula"), merge_col = c("model"))

head(CH4_mod_CO2C.2_s[order(-CH4_mod_CO2C.2_s$pvalue),])
# order(cCH4_mod_CO2C.2_fm6_sum$pvalue)

# ggplot(x = CH4_logn1.x, y = aic, data = cCH4_mod_CO2C.2_fm6_sum) + geom_scatter()
# ggplot(data = cCH4_mod_CO2C.2_fm6_sum, aes(x = CH4_logn1.x, y = pvalue)) + geom_point()

# filter step likely should be part of workflow above, but doing it here...

# CH4 base with forced, comp models: CH4_mod_fch4.2_c
base_sig <- c("ch4f2.2", "ch4f2.21")
base_sim <- c('ch4f2.4','ch4f2.8','ch4f2.11','ch4f2.16')
base_pass <- c(base_sig, base_sim)
#base_pass

# passing composites
CH4_mod_fch4.2_s <- CH4_mod_fch4.2_c[CH4_mod_fch4.2_c$model %in% base_sig, ]
CH4_mod_fch4.2_p <- CH4_mod_fch4.2_c[CH4_mod_fch4.2_c$model %in% base_pass, ]
CH4_mod_fch4.2_p

# CH4 CO2/C with forced, comp models: CH4_mod_CO2C.2_c
co2C_sig <- c("ch4co2C2.623", "ch4co2C2.392")
co2C_sim <- c('ch4co2C2.628')
co2C_pass <- c(co2C_sig, co2C_sim)
co2C_pass

CH4_mod_CO2C.2_s <- CH4_mod_CO2C.2_c[CH4_mod_CO2C.2_c$model %in% co2C_sig, ]
CH4_mod_CO2C.2_p <- CH4_mod_CO2C.2_c[CH4_mod_CO2C.2_c$model %in% co2C_pass, ]
# CH4_mod_CO2C.2_p

# combine base composites (forced, CO2C)
CH4_mods_p <- rbind(CH4_mod_fch4.2_s, CH4_mod_CO2C.2_s)
CH4_mods_p

# MOB_IIa models are:
MOBIIa_ng <- MOBIIa_mod_ng_sem$filt_models         
MOBIIa_fng <- MOBIIa_mod_fAOANOB_sem$filt_models   # try first, forced AOA/NOB
MOBIIa_ngt <- MOBIIa_mod_ngt_sem$filt_models

# rename form
names(MOBIIa_fng)[6] <- "MOB_IIa"  # hmmm, should be part of workflow somewhere...
MOBIIa_fng 

# get combined formulae for tree  
# issues with not positive definite, see model 6, 26 (problem) -- need try catch here!
CH4mob_comp_f <- combine_sem_model_layers(CH4_mods_p, MOBIIa_fng)
#CH4mob_comp_f <- combine_sem_model_layers(CH4_mod_fch4.2_s, MOBIIa_fng)





# drop bad models  - DEPRECATED W TRY Catch
#bad_mods <- c("ch4f2.2_2", "ch4f2.2_11", "ch4f2.21_25")
# CH4mob_comp_f <- CH4mob_comp_f[!CH4mob_comp_f$model %in% bad_mods, ]
#CH4mob_comp_f[10,]
#CH4mob_comp_f[,-1] #<-  CH4mob_comp_f[1:10,]# note something off here, MOB_IIa col is "form", should be MOB_IIa

# Run models for results
CH4mob_comp_r <- suppressWarnings(run_compare_sem_models(CH4mob_comp_f, Guild_CH4_0, 
                                            estimator = "mlm", keep_mod_descr = c("model_out","model")))

#CH4mob_comp_r



# get summary

CH4mob_comp_t <- clean_SEM_model_table_for_output(CH4mob_comp_f, CH4mob_comp_r, 
                                drop_cols = c("model_out", "formula"), merge_col = c("model"))

head(CH4mob_comp_t[order(-CH4mob_comp_t$pvalue),])
# order(cCH4_mod_CO2C.2_fm6_sum$pvalue)

# get 2nd branch acetoclast models
CH4_ac2 <- CH4ac_mod.2_sem$filt_models
names(CH4_ac2)[6] <- "CH4_ac"

# combine with base model 
CH4_ac_comp_f <- combine_sem_model_layers(CH4mob_comp_f, CH4_ac2)
#CH4_ac_comp_f <- CH4_ac_comp_f[20:40,]

# drop bad models
#bad_mods <- c("ch4f2.21_249") # row 19 ! 
#CH4_ac_comp_f <- CH4_ac_comp_f[!CH4_ac_comp_f$model %in% bad_mods, ]
head(CH4_ac_comp_f) 



# below are "promising looking" models for CH4_ac composite branches, sent to semPlot in earlier version (dev 0.7?)
CH4_ac_comp_f[CH4_ac_comp_f$model =='ch4f2.2_24_133',]$formula;
CH4_ac_comp_f[CH4_ac_comp_f$model =='ch4f2.2_24_141',]$formula 
CH4_ac_comp_f[CH4_ac_comp_f$model =='ch4f2.2_24_39',]$formula 
CH4_ac_comp_f[CH4_ac_comp_f$model =='ch4f2.2_15_313',]$formula




CH4_ac2#CH4_mod_fch4.2_s
MOBIIa_fngCH4_ac_comp_f
# Run models for results
CH4_ac_comp_r <- suppressWarnings(run_compare_sem_models(CH4_ac_comp_f, Guild_CH4_0, 
                                            estimator = "mlm", keep_mod_descr = c("model_out","model")))

head(CH4_ac_comp_r)

# get summary
CH4_ac_comp_t <- clean_SEM_model_table_for_output(CH4_ac_comp_f, CH4_ac_comp_r, 
                                drop_cols = c("model_out", "formula"), merge_col = c("model"))

CH4_ac_comp_t[order(-CH4_ac_comp_t$pvalue),]  # sort on pvalue
# order(cCH4_mod_CO2C.2_fm6_sum$pvalue)

#summary(sem.ch4f2.2_24_236)

# Filter by pvalue, order by aic decending
CH4_ac_comp_tf <- CH4_ac_comp_t[CH4_ac_comp_t$pvalue > 0.15, ]
CH4_ac_comp_tf <- CH4_ac_comp_tf[order(CH4_ac_comp_tf$aic),] 
CH4_ac_comp_tf

#R.Version()
#library(semplot)
#install.packages("semPlot")

#library(devtools)

#install_github('SachaEpskamp/semPlot',  dependencies = T)

#library('semPlot')

# SEM tree base - Delta all
CH4_da_base_m <- get_sem_models_add_prefix(CH4_delt_all_sem, models = "filt_models", prefix = "dch4a.")

CH4_da_tree_c <- get_composites_model_table(CH4_da_base_m, composite_sets, assign_other = "y_var",#)"CH4_gen", 
                           no_int = FALSE, compos_fail_collapse_all = TRUE)


# SEM tree base - Delta all
#CH4_da_base_m <- CH4_delt_all_sem$filt_models         # source, LASS SEM FILTERED models (obj)
#CH4_da_base_m$model <- paste0("ch4da." ,CH4_da_base_m$model)

# Get composites for MODEL
#CH4_da_tree_c <- get_composites_model_table(CH4_da_base_m, composite_sets, assign_other = "y_var",#)"CH4_gen", 
#                           no_int = FALSE, compos_fail_collapse_all = TRUE)

#CH4_da_tree_c -- note all but below fail compositing

# pick base models
CH4_da_tree_s <- CH4_da_tree_c[1:5,]   # picking first 6 models, others look funny, don't run? 
tree_base <- CH4_da_tree_s
CH4_da_tree_s

# other base models possible? -- NO, next try CO2/C base
# sem_delta_base <- run_sem_model_summary(CH4_da_tree_c, Guild_CH4_dS, sort_by = 'pvalue', sort_cutoff = 0.00)

# CH4 base model with forced MOB_IIa/CH4_ac (not filtered)
CH4_fg_base_m <-CH4_delt_f_mobac_sem$sem_models                # all models
CH4_fg_base_m$model <- paste0("ch4fg." ,CH4_fg_base_m$model)

# Get composites for MODEL
CH4_fg_tree_c <- get_composites_model_table(CH4_fg_base_m, composite_sets, assign_other = "CH4_gen", 
                           no_int = FALSE, compos_fail_collapse_all = TRUE)

#CH4_fg_tree_c
CH4_fg_tree_s <- CH4_fg_tree_c[1:7,]
CH4_fg_tree_s

tree_base <- CH4_fg_tree_s

# CH4 base model with forced MOB_IIa/CH4_ac (not filtered)
CH4_co2c_base_m <-CH4_delt_CO2C_sem$sem_models                # all models
CH4_co2c_base_m$model <- paste0("ch4co2c." ,CH4_co2c_base_m$model)

# Get composites for MODEL
CH4_co2c_tree_c <- get_composites_model_table(CH4_co2c_base_m, composite_sets, assign_other = "CH4_gen", 
                           no_int = FALSE, compos_fail_collapse_all = FALSE)
#CH4_co2c_tree_c
CH4_co2c_tree_s <- CH4_co2c_tree_c[1:5,]
CH4_co2c_tree_c

ch4_co2_pick <- c(1, 6, 7, 8, 18)
CH4_co2c_tree_s2 <- CH4_co2c_tree_c[ch4_co2_pick,]
#tree_base <- CH4_co2c_tree_s2

#tree_base 

# newer workflow
d_mob2af <- get_sem_models_add_prefix(MOBIIa_delt_fAOANOB_sem, models = "filt_models", prefix = "mob2af.")
d_mob2afa <- get_sem_models_add_prefix(MOBIIa_delt_fAOANOB_sem, models = "sem_models", prefix = "mob2afa.",
                                     filter = 'R2', cutoff = 0.5)

d_mob2ngt <- get_sem_models_add_prefix(MOBIIa_delt_fAOANOB_sem, models = "filt_models", prefix = "mob2ngt.")
d_mob2ngta <- get_sem_models_add_prefix(MOBIIa_delt_fAOANOB_sem, models = "sem_models", prefix = "mob2ngta.",
                                     filter = 'R2', cutoff = 0.5)

# assign to branch and get y-var
branch_1 <- d_mob2af
branch_1 <- rename_formula_col_as_yvar(branch_1)

# branch_1
# prepare branch data

# NG (removes non-relevant guilds, not forcing AOA, NOB)
# branch_1 <-MOBIIa_delt_ng_sem$sem_models


# FORCE AOA NOB
branch_1 <- MOBIIa_delt_fAOANOB_sem$filt_models        # forced AOA / AOB
branch_1 <- MOBIIa_delt_fAOANOB_sem$sem_models        # forced AOA / AOB

# NGT
#branch_1 <- MOBIIa_delt_ngt_sem$filt_models           # tried w CH4 ac.2-.4
# branch_1 <- MOBIIa_delt_ngt_sem$sem_models           # ran comparison w CH4_ac.4

#branch_1 <- MOBIIa_delt_ngt_sem$sem_models 
# branch_1 <- MOBIIa_delt_ngt_sem$filt_models          # forced NO Guild
#    names(branch_1)[6] <- "MOB_IIa"

branch_1 <- rename_formula_col_as_yvar(branch_1, formula_col = 'form', fxn_operator = '~')

# Filter models by R2 (for all)
#branch_1 <- branch_1[branch_1$R2 > 0.5,]
#branch_1
d_ch4ac.2 <- get_sem_models_add_prefix(CH4ac_delta.2_sem, prefix = "dac2.")
d_ch4ac.2a <- get_sem_models_add_prefix(CH4ac_delta.2_sem, models = "sem_models", prefix = "dac2.",
                                     filter = 'R2', cutoff = 0.5)

d_ch4ac.3 <- get_sem_models_add_prefix(CH4ac_delta.3_sem, prefix = "dac3.")
d_ch4ac.3a <- get_sem_models_add_prefix(CH4ac_delta.3_sem, models = "sem_models", prefix = "dac3a.",
                                     filter = 'R2', cutoff = 0.5)
d_ch4ac.4 <- get_sem_models_add_prefix(CH4ac_delta.4_sem, prefix = "dac4.")
d_ch4ac.4a <- get_sem_models_add_prefix(CH4ac_delta.4_sem, models = "sem_models", prefix = "dac4a.",
                                     filter = 'R2', cutoff = 0.3)

d_ch4ac.4f <- get_sem_models_add_prefix(CH4ac_delta_f_sem, prefix = "dac4f.")
d_ch4ac.4fa <- get_sem_models_add_prefix(CH4ac_delta_f_sem, models = "sem_models", prefix = "dac4a.",
                                     filter = 'R2', cutoff = 0.3)

# assign to branch and get y-var
branch_2 <- d_ch4ac.2a
branch_2 <- rename_formula_col_as_yvar(branch_2)

#branch_2
branch_2 <- CH4ac_delta.2_sem$filt_models              # CH4_ac2 - works with SCALED data, crashes unscaled
#branch_2 <- CH4ac_delta.3_sem$filt_models              # CH4_ac2 - works with SCALED data
#branch_2 <- CH4ac_delta.4_sem$filt_models              # CH4_ac2 - works with SCALED data  
    names(branch_2)[6] <- "CH4_ac"

branch_2
# make list of branches
#head(branch_1)
branchz <- list(branch_1)
branchz2 <- list(branch_1, branch_2)

# combine branches, check n models-- note 720 models = 118 GB peak ram w delta DS
sem_formula_tree_d <-sem_tree_formula_build(tree_base, branchz2)
sem_formula_tree_d

# run simple delta models (4 min for 720? -- smaller data set) -- use SCALED dataset!! fails otherwise (& slower)!
#sem_run_test_d <- run_sem_model_summary(sem_formula_tree_d, Guild_CH4_dS, sort_by = 'pvalue', sort_cutoff = 0.001) # pass parallel!!

sem_run_test_d

# filter results before printing
# sem_run_test_d <- sem_run_test_d[sem_run_test_d$pvalue > 0.05,]

#write.table(sem_run_test, "auto_tables/delta_CO2C_mobaf_ac.2_noSal_675.txt", sep='\t', row.names=FALSE)
# write.table(sem_run_test_d, "auto_tables/delta_CH7f_mobngALL_ac2a_noSal_1953.txt", sep='\t', row.names=FALSE)

# CH4_mod_fch4.2 - could use for filter

# BASE OF TREE # ALL SITES, CO2 models
CH4_mod_fch4.2 <- get_sem_models_add_prefix(CH4_mod_fch4.2_sem, models = "filt_models", prefix = "fch4.2.")

# Get composites for MODEL
CH4_mod_fch4.2_c <- get_composites_model_table(CH4_mod_fch4.2, composite_sets, assign_other = "y_var",#)"CH4_gen", 
                           no_int = FALSE, compos_fail_collapse_all = FALSE)

# remove failed composites (not needed here)
# CH4_mod_fch4.2_c <- CH4_mod_fch4.2_c[as.character(CH4_mod_fch4.2_c$comp_failed) == 'character(0)',]

#head(CH4_mod_fch4.2_c)
# CH4_mod_fch4.2_c
# run composite models - NAH
CH4_mod_fch4.2_r <- suppressWarnings(run_compare_sem_models(CH4_mod_fch4.2_c, Guild_CH4_aS, 
                                            estimator = "mlm", keep_mod_descr = c("model_out","model")))

CH4_mod_fch4.2_r 
# CH4_mod_fch4.2_A <- get_sem_models_add_prefix(CH4_mod_fch4.2_sem, models = "sem_models", prefix = "fch4.2.",
#                                          filter = 'R2', cutoff = 0.5)

# BASE OF TREE # ALL SITES, CO2 models
# of tree prep, including composite sets define & assign.
#CH4_mod_fch4.2_m <- CH4_mod_fch4.2_sem$filt_models         # source, LASS SEM FILTERED models (obj)
#CH4_mod_fch4.2_m$model <- paste0("ch4f2." ,CH4_mod_fch4.2_m$model)

# Get composites for MODEL
#CH4_mod_fch4.2_c <- get_composites_model_table(CH4_mod_fch4.2_m, composite_sets, assign_other = "y_var",#)"CH4_gen", 
#                           no_int = FALSE, compos_fail_collapse_all = TRUE)

# cherry pick restults & get models (Using run compare, summarize NOT SHOWN HERE, picked only sig)
base_sig <- c("fch4.2.2", "fch4.2.21")
#base_sig <- c("ch4f2.2", "ch4f2.21")
CH4_mod_fch4.2_s <- CH4_mod_fch4.2_c[CH4_mod_fch4.2_c$model %in% base_sig, ]

tree_base <- CH4_mod_fch4.2_s
CH4_mod_fch4.2_s

# MOBIIa_fng <- MOBIIa_mod_fAOANOB_sem$filt_models  -- SEM lasso [filtered_models]?
# names(MOBIIa_fng)[6] <- "MOB_IIa"  # hmmm, should be part of workflow somewhere...  what is this?  Also need to label model # ?
# CH4_mod_fch4.2_c

# CO2 / C models, base
CH4_mod_CO2C.2_fm <- CH4_mod_CO2C.2_sem$filt_models  
CH4_mod_CO2C.2_fm$model <- paste0("ch4co2c." ,CH4_mod_CO2C.2_fm$model)

# Get composites
# CH4_mod_CO2C.2_c <- get_composites_model_table(CH4_mod_CO2C.2_fm, composite_sets, assign_other = "y_var",#)"CH4_gen", 
#                            no_int = FALSE, compos_fail_collapse_all = TRUE)

# Get composites -- NOTE SHIFT IN ASSIGN OTHER!
CH4_mod_CO2C.2_c <- get_composites_model_table(CH4_mod_CO2C.2_fm, composite_sets, assign_other = "CH4_gen", 
                           no_int = FALSE, compos_fail_collapse_all = TRUE)

# remove failed composites (not needed here)
CH4_mod_CO2C.2_c <- CH4_mod_CO2C.2_c[as.character(CH4_mod_CO2C.2_c$comp_failed) == 'character(0)',]
#CH4_mod_CO2C.2_c

#head(CH4_mod_CO2C.2_c)
# check results (incl. run compare?)

#keep_CO2C <- c(1, 10, 12)
#tree_base <- CH4_mod_CO2C.2_c[keep_CO2C, ]
#tree_base

# how to deal with failed composites?
# unlist(CH4_mod_CO2C.2_c$comp_failed)
# identical(CH4_mod_CO2C.2_c$comp_failed[[1]], character(0))
#grep('',form_tab)

# form_tab <- CH4_mod_CO2C.2_c

# lapply(colnames(form_tab), function(x) grep('', form_tab[,x]))

#CH4_mod_CO2C.2_c$comp_failed
#CH4_mod_CO2C.2_c$comp_failed2 <- if(identical(CH4_mod_CO2C.2_c$comp_failed, character(0)), "t",'f')

##ifelse(identical)

# newer workflow
mob2af <- get_sem_models_add_prefix(MOBIIa_mod_fAOANOB_sem, models = "filt_models", prefix = "mob2af.")
mob2afa <- get_sem_models_add_prefix(MOBIIa_mod_fAOANOB_sem, models = "filt_models", prefix = "mob2afa.",
                                     filter = 'R2', cutoff = 0.5)

# assign to branch and get y-var
branch_1 <- mob2af
branch_1 <- rename_formula_col_as_yvar(branch_1)
# branch test data, from lasSEM object to filtered models
branch_1 <- MOBIIa_mod_fAOANOB_sem$filt_models     # was MOBIIa_fng
branch_1$model <- paste0("mobF.",branch_1$model)

names(branch_1)[6] <- "MOB_IIa"


# branch_1

# newer workflow

# Prepare data for tree
ch4_ac.2 <- get_sem_models_add_prefix(CH4ac_mod.2_sem, models = "filt_models", prefix = "ac.2.")  # CH4_ac.2
ch4_ac.2a <- get_sem_models_add_prefix(CH4ac_mod.2_sem, models = "sem_models", prefix = "ac.2a.")  # CH4_ac.2a

ch4_ac.3 <- get_sem_models_add_prefix(CH4ac_mod.3_sem, models = "filt_models", prefix = "ac.3.")  # CH4_ac.2
ch4_ac.4 <- get_sem_models_add_prefix(CH4ac_mod.4_sem, models = "filt_models", prefix = "ac.4.")  # CH4_ac.2

# assign to branch and get y-var
branch_2 <- ch4_ac.2a
branch_2 <- rename_formula_col_as_yvar(branch_2)


#branch_2
# branch test data, from lasSEM object to filtered models

# CH4_ac.2
branch_2 <- CH4ac_mod.2_sem$filt_models            # CH4_ac2  
#branch_2 <- CH4ac_mod.2_sem$sem_models            # CH4_ac2  
branch_2$model <- paste0("ac.2",branch_2$model)

# CH4_ac.3
#branch_2 <- CH4ac_mod.3_sem$filt_models            # CH4_ac3 -- no fit measures if not converge -- run compare handling error!  
#branch_2$model <- paste0("ac.3",branch_2$model)


# CH4_ac.4
#branch_2 <- CH4ac_mod.4_sem$filt_models            # CH4_ac4 -- no converge error here too!  
#branch_2$model <- paste0("ac.4",branch_2$model)

branch_2 <- branch_2[branch_2$R2 > 0.5,]
    names(branch_2)[6] <- "CH4_ac."
#tree_base <- CH4_mod_fch4.2_s
tree_base <- CH4_mod_fch4.2_s

branchz <- list(branch_1)
branchz2 <- list(branch_1, branch_2)

sem_formula_tree <-sem_tree_formula_build(tree_base, branchz2) # , run = TRUE, warn_n_models = TRUE)
#sem_formula_tree

# 360 models, peak ram = 112 GB!

# below use scaled data? Guild_CH4_aS

# run it to summary
#run_sem_model_summary(sem_formula_tree, Guild_CH4_0)
sem_run_test <- run_sem_model_summary(sem_formula_tree, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.01) # pass parallel!!
#run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0)


# branchz

sem_run_test

#names(sem_run_test)

# sem_run_test[order(sem_run_test$npar),]

# filter results before printing
# sem_run_test <- sem_run_test_d[sem_run_test_d$pvalue > 0.05,]

write.table(sem_run_test, "auto_tables/ALL_sites_CH4_mobF_ac.2_360.txt", sep = '\t', row.names = FALSE)
#write.table(sem_run_test, "ALL_sites_CH4CO2C_mobF_ac.2.txt", sep = '\t', row.names = FALSE)


