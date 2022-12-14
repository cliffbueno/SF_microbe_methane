# SEM - Lasso semi-auto for SF OTU guilds v1.12
#-  runs under r > 3.6 (glmnet for LASSO)

# ASSUME here that libraries loaded, incl. own R libs with data loading embedded

########################################################################################
### 1) Import data sets

## a) Import metadata -- Read data fom import module for testing convenience here
# source("Import_SalOTU_dat_Plot_test_v0.1.R")
# head(Meta_iTag)


## b) Imported OTU table
# head(otu_V)


## c) Separate Delta Sites (oligo and FW)                       # levels(Meta_iTag$SALTgroup)
Meta_iTag_FW <- Meta_iTag[Meta_iTag$SALTgroup =="FW",]
Meta_iTag_Oligo <- Meta_iTag[Meta_iTag$SALTgroup =="Oligo",]

Meta_iTag_Delta <-rbind(Meta_iTag_FW, Meta_iTag_Oligo)          # Meta_iTag_Delta
Delta_sites <-data.frame(Meta_iTag_Delta[,"Sample"])            # DF of Delta sites
colnames(Delta_sites) [1] <- "Sample"                           # Rename column "Sample"
#Delta_sites["Samp_Index"] <- seq(1:nrow(Delta_sites))          # Make sample index for reordering post-merges
# Delta_sites


########################################################################################
### 2) Process guild data and metadata

## 0) import guild counting from OTU table function from module 
# source("OTU_subsetting_modules_v.0.2_strip.R")                # IMPORTING R MODULE!! - MOVE 
# options(repr.plot.width=1.5, repr.plot.height=6) 


## a) get guild counts/sample data from OTU table
Guild_OTUs <- Get_16S_Guilds_alt(otu_V)                             # use Get_16S_Guilds to get guilds  
otu_V["OTU"] <- row.names(otu_V)                                # Make OTU number column
# dim(Guild_OTUs); head(Guild_OTUs)

# merge otu table and Guilds 
OTU_guilds <- merge(Guild_OTUs, otu_V, by="OTU")#, all.y=TRUE)
# head(OTU_guilds)


## b) Get C guilds from Phyla (based on shotgun data for C)
Firmic <- otu_V[subset(otu_V["Phylum"]=="Firmicutes"),]         # Firmicutes
Firmic['Guild'] <-"Firmic" 

Actino <- otu_V[subset(otu_V["Phylum"]=="Actinobacteriota"),] 
Actino['Guild'] <-"Actino" 

Chlorf <- otu_V[subset(otu_V["Phylum"]=="Chloroflexi"),] 
Chlorf['Guild'] <-"Chlorf" 

C_guilds <- rbind(Firmic, Actino, Chlorf)


## c) Process guild data and metadata 
# Get agg_by_cat aggregatopm function from corr module
# source("Corr_ranks_module_v0.3.2_strip.R")                    # IMPORTING R MODULE!! - MOVE  

# Aggregate OTU guilds, C guilds & combine  
Guild_agg <- agg_by_cat(OTU_guilds, "Guild")                    # aggregate by Guild
C_guild_agg <- agg_by_cat(C_guilds, "Guild")
Guild_agg2 <-rbind(Guild_agg, C_guild_agg)                      # Combine OTU and C guild data


## d) calculate new guild ratios, log transform guild data
# prepare data
Guild_agg2 <- data.frame(t(Guild_agg2[,-1]))                    # Transpose for new calcs
Guild_agg2 <- replace(Guild_agg2, Guild_agg2==0, 2)             # replace 0s with pseudo counts 

# Calculate guild ratios 
Guild_agg2$AO_NOB <- (Guild_agg2$AOA + Guild_agg2$AOB) / Guild_agg2$NOB
Guild_agg2$NOB_AO <- Guild_agg2$NOB / (Guild_agg2$AOA + Guild_agg2$AOB)
Guild_agg2$mcr_pmo <- ((Guild_agg2$CH4_ac + Guild_agg2$CH4_H2) / 
                            (Guild_agg2$MOB_I + Guild_agg2$MOB_II + Guild_agg2$MOB_IIa))
Guild_agg2$pmo_mcr <- ((Guild_agg2$MOB_I + Guild_agg2$MOB_II + Guild_agg2$MOB_IIa) / 
                       (Guild_agg2$CH4_ac + Guild_agg2$CH4_H2))

# Log 2 transform guild data
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

########################################################################################
## Metadata processing

## f) Process CH4 and other soil metadata
### Reimport metadata instead of using module loaded metadata
metaDB <- read.table("SF_sal_meta_FIX3b_gap_fill_MDS.txt", sep='\t', header=T)
# metaDB <- read.table("SF_sal_metaLOG_FIX2.txt", sep='\t', header=T)
# head(metaDB); names(metaDB)

# Get only numeric data, LOG 10 transform soil chemistry
metaCHEM <- metaDB[,15:ncol(metaDB)]                                             # should use lapply
metaCHEM_log <-log10(metaCHEM)                                                   # log10 chem
CH4_logn1 <- log10(metaDB[,"CH4_ug_m2_h"]- (min(metaDB[,"CH4_ug_m2_h"])*1.05))   # CH4 logn1
metaCHEM_log["CH4_logn1"] <- CH4_logn1                                           # metaCHEM_log      

# reattach non-numeric
metaDB = data.frame(metaDB[,1:14], metaCHEM_log)                                 # metaDB

# Merge site order and Samples
Meta_iTag <- merge(metaDB, OTU_samps, by='Sample')                               # colnames(metaDB)
rownames(Meta_iTag) <- Meta_iTag$Sample

# Reorder location factor
Meta_iTag$Location <-factor(Meta_iTag$Location, levels=c("Sandmound","WestPond","Mayberry","Browns","RushRanch","Joice","Goodyear","WhiteSlough","Tolay","ChinaCamp","Muzzi"))  #head(Meta_iTag)
Meta_iTag$Pl_Sp <-factor(Meta_iTag$Pl_Sp, levels=c("Cattail","Tule","ThreeSq","CattailNL","Phrag","PW","Cord"))

# Resort meta itag by index
indexer = 'EWsiteHyd_index'
Meta_iTag <- Meta_iTag[order(Meta_iTag[indexer]),]
# colnames(Meta_iTag) #head(Meta_iTag) #max(Meta_iTag$CH4_ug_m2_h) #plot(x=Meta_iTag$Sample, y=Meta_iTag$CH4_ug_m2_h)

# reduce metadata factors
CH4_samp <- c("Sample", "CH4_ug_m2_h","CH4_logn1",'CO2_mg_m2_h','CO2_soilC_mg_g_d','CH4_CO2', "Bulk_dens", 'H2O_FPS',
               'pH', 'C','N','P','CN','NP', 'NO3_N', 'NH4_N', 'Olsen_P', 'NO3_NH4','NP_ext',
               'Salinity.x', 'Cl', 'SO4', 'SO4_pw', 'Fe', 'Fe_pw','DOC_mg_L')
              
              # 'NO2_pw','NO3_pw', 'NH3_pw', 
              # "C_g_m2")  # not for log data, already been logn1 transf: "CH4_logn1", 
CH4 <- Meta_iTag[,CH4_samp]  
# CH4

########################################################################################
## Metadata processing, with Guilds data

## g) Merge metadata with guild data 
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
# dim(Guild_CH4_d); head(Guild_CH4_d); colnames(Guild_CH4_d)


## i) Write out convenience tables for later use 
# dir.create("SEM_data")

# write.table(Guild_CH4_d, "SEM_data/SEM_base_log2_10_delta_guild_soil_data.txt")
# write.table(Guild_CH4, "SEM_data/SEM_base_log2_10_all_guild_soil_data.txt")