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
# R_otu_modules/?
#source("Import_SalOTU_dat_Plot_test_v0.1.R")                             # load data and metadata
#source("OTU_subsetting_modules_v.0.2_strip.R")                           # gets guilds, contains colors 
#suppressMessages(source("Corr_ranks_module_v0.3.2_strip.R"))             # contains simple aggregation function

# For Cliff
source("~/Documents/GitHub/SF_microbe_methane/modules/Import_Silva_OTU_data4plots_v0.1.R")
source("~/Documents/GitHub/SF_microbe_methane/modules/3_OTU_subsetting_modules_v.0.4_strip.r")
source("~/Documents/GitHub/SF_microbe_methane/modules/6_Corr_ranks_module_v0.3.4_strip.R")