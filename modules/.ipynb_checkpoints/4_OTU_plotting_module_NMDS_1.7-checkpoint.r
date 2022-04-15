
################################################################################
### OTU plotting module NMDS dev.v.1.7 

### Inputs 
#    - 1) otu table (ideally VST)
#    - 2) env metadata table 
#    - 3) color variable vector 

### Outputs NMDS plots (Function outline):
#  - F1): colors only -------------------------------- NMDS_plot_group(otu_V, Location)
#  - F2): colors, shapes ----------------------------- NMDS_group_shape(otu_V, Location, Restor)
#  - F3): colors only, env. fit arrows --------------- NMDS_group_fit(otu_V, group, Env_chem)
#  - F4): colors, shapes, env. fit arrows ------------ NMDS_group_shape_fit(otu_V, Location, Restor, Env_chem)

# Other: F5): gapfill NA data with minimum ----------- min_fill(metaDB, col) 

# Note: Fxns 1-4 were not easily combined with single call and options, could be later consolidated

################################################################################
#### Pre-formatting of input data may be desired, including:

# - Merge of metadata, otu table Samples (so number of samples matching)
# - Def of color and shape vars, incl. possible factor releveling
# - gapfilling of Env. data (sensitive to NA, min_fill fxn provided)
# - reducing Env. data, renaming cols. for envFit arrow clarity

# EXTENSIVE examples of preprocessing left commented out below in Load test data

################################################################################
# Import packagaes
library(ggplot2)
library(RColorBrewer)
library(reshape2)

suppressMessages(library(gtools))
suppressMessages(library(vegan))
suppressMessages(library(gplots))

################################################################################

# Load TEST data
# and preprocess examples

# Import OTU Table
#OTU_v <- read.table("SF_Sal_OTU_VSTcpm.txt", sep='\t', header=T, row.names=1)                   # dim(OTU_v); head(OTU_v)

# otu_V <- read.table("SF_Sal_OTU_VSTcpm.txt", sep='\t', header=T, row.names=1)                   # dim(OTU_v); head(OTU_v)

# Sort OTU table                                                                      
#otu_V <-OTU_v[order(OTU_v$Consensus.lineage),]                                  # sort by lineage  

#dim(otu_V)

# Drop Sandmound Cattail samps
#otu_V <- otu_V[, (colnames(otu_V) !='Sandmound_CattailA_D1')]
#otu_V <- otu_V[, (colnames(otu_V) !='Sandmound_CattailA_D2')]
                  
#dim(otu_V)
#c(names(otu_V))



# Import site colors
#site_colors <- read.table("Sal_siteColors_testR.txt", sep='\t', header=T, row.names=1)          # site_colors
#site_colours <- (site_colors$Salpal3_col)                                                       # only color

# Get list from vector, move inside functions?
#site_col <- levels(site_colours)[as.numeric(site_colours)]  # collect levels, as.numeric preserves ordering
#site_col

# Import Sample mapping
#metaDB <-read.table("SF_sal_meta_FIX3b_gap_fill_MDS.txt", sep="\t", header=TRUE)          # Import Metadata, keep all    
#metaDB <-read.table("SF_sal_meta_FIX3b.txt", sep="\t", header=TRUE)          # Import Metadata, keep all    
#row.names(metaDB) <- metaDB$Sample                                          # Row names are samples for phyloseq             
#metaDB = metaDB[,-1]      
# head(metaDB)

# Simple gapfilling function   NA <- 0.5 * min value
min_fill = function(metaDB, col){

    vect <- metaDB[col]
    vect[is.na(vect)] <- min(vect[!is.na(vect)])*0.5
    metaDB[col] <- vect
    return(metaDB)   
}

# Manually filled: NO3_N, Fe, Mn, Cu, Zn

#To be gapfilled (0.5*min):  NO2_pw, NO3_pw, NH3_pw, PO4_pw,#  Zn_pw, K_pw

# Apply gap_filling, recalc. stoich
# need an lapply here.
#metaDB <- min_fill("NO2_pw")
#metaDB <- min_fill("NO3_pw")
#metaDB <- min_fill("NH3_pw")
#metaDB <- min_fill("PO4_pw")

#fill_cols <- c("NO2_pw", "NO3_pw", "NH3_pw", "PO4_pw")
#metaDB[fill_cols]

#metaDB["NP_pw"] <- (metaDB["NO2_pw"] + metaDB["NO3_pw"] + metaDB["NH3_pw"]) /14/(metaDB["PO4_pw"]/31)
#metaDB["NP_pw"]
#metaDB["NC"] <- metaDB["N"]/14 /(metaDB["NO3_pw"]/12)

# Get matching OTU samples and metadata samples (now done by preprocess?)
# Get Sample names in OTU table               
#OTU_samps <- data.frame('Sample'=colnames(otu_V))                #OTU_samps

# Merge site order and Samples
#Meta_iTag <- merge(metaDB, OTU_samps, by='Sample')               #colnames(metaDB)
#rownames(Meta_iTag) <- Meta_iTag$Sample

# Resort meta itag by index
#indexer = 'EWsiteHyd_index'
#Meta_iTag <- Meta_iTag[order(Meta_iTag[indexer]),]
# colnames(Meta_iTag)

# Reorder factors as needed
# Reorder location factor
#Meta_iTag$Location <-factor(Meta_iTag$Location, levels=c("Sandmound","WestPond","Mayberry","Browns","RushRanch","Joice","Goodyear","WhiteSlough","Tolay","ChinaCamp","Muzzi"))  #head(Meta_iTag)
#Meta_iTag$Pl_Sp <-factor(Meta_iTag$Pl_Sp, levels=c("Cattail","Tule","ThreeSq","CattailNL","Phrag","PW","Cord"))

# Def NMDS vars
# Prepare data, cats for NMDS plot  -- move up to top or in function?
#Location <-Meta_iTag$Location                                                    # Get Location vector from meta_iTag
#Restor <- Meta_iTag$EWcoastGroup
#Plant <- Meta_iTag$Pl_Sp
# Plant

#chem_cols <- sapply(Meta_iTag, is.numeric)  
#Env_chem <- Meta_iTag[chem_cols]
#names(Env_chem)
#dim(Env_chem)
#names(Env_chem[5:40])
#names(Env_chem[5:30])

### Select and rename vars

#keep_vars <-c('Salinity.x', 'Bulk_dens', 'CO2_mg_m2_h', 'CH4_ug_m2_h', 'C', 'N', 'P', 'N_g_m2', #'C_g_m2', 
#              'NP', 'NO3_N', 'NH4_N', 'DOC_mg_L', 'NO2_pw','NO3_pw','NH3_pw','CN')#,
#Env_chem <-Env_chem[keep_vars]

### CAREFUL HERE, lazy renaming... comment out when testing newvars.
#colnames(Env_chem) <-c('Salinity', 'BD', 'CO2', 'CH4', 'C', 'N', 'P', 'N_g/m2',                  #'C_m2', 
#                       'N:P', 'NO3_s', 'NH4_s', 'DOC', 'NO2_pw','NO3_pw','NH3_pw','N:C')
#head(Env_chem)

#Env_chem <- Env_chem[,5:30] ; #chem_cols <- sapply(Meta_iTag, is.numeric); #Env_chem <- Meta_iTag[chem_cols[5:30]]



################################################################################
### 1) NMDS plots with groups only

NMDS_plot_group = function(otu_V, group){

    # Get OTU data for NMDS
    counts <- sapply(otu_V, is.numeric)                                               # Select numeric columns
    otu_d <- data.matrix(otu_V[,counts])                                              # Filter from OTU tab          #length(Restor) #length(Plant) #dim(otu_d)
    biomT <-t(otu_d)                                                                  # Transpose so cols are vars

    # Get NMDS scores
    bT_mds<-metaMDS(biomT, distance="bray",k=3, trymax=10);                          # run NMDS
    bT_mds_DF = data.frame(NMDS1=bT_mds$points[,1], NMDS2=bT_mds$points[,2], group)  # Make data frame of scores
    stress=round(bT_mds$stress, digits = 3)                                                                      
    # return(bT_mds_DF)  
    
    # Test ADONIS models, extract params
    bT_adonis <- adonis(biomT  ~ group, permutations=99, method="bray"); # data=env,
    R2 <- round(bT_adonis$aov.tab$R2[1], digits=3)
    P <- bT_adonis$aov.tab$"Pr(>F)"[1]
    
    # PLOT NMDS by group color, shape
    pA <- ggplot(bT_mds_DF) +
        geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=group, shape=shape_var), alpha = 0.8) + #+ geom_point(alpha = 0.8) +
        stat_ellipse(mapping = aes(x=NMDS1, y=NMDS2, color=group), alpha = 0.8) #  level=0.95, + ggtitle(paste("ADONIS R2 =", R2, "Stress:", stress));
    pB <- pA + scale_color_manual(values = c(site_col)) + scale_shape_manual(values=c(20, 17), labels = c("Reference", "Restored"))
    pC <- pB + theme(legend.title=element_blank())       + 
       theme_bw() + theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 6),
                          legend.key.height = unit(0.5, "cm"),
                          legend.key.width = unit(0.5, "cm"),
                          legend.spacing.y = unit(0, "cm"),
                          legend.margin = margin(5,0,0,0),
                          legend.box.margin = margin(0,0,0,-5),
                          axis.title.y = element_text(size = 10, margin = margin(t = 0, r = -2.5, b = 0, l = 0)),
                          axis.title.x = element_text(size = 10))
    
    return(pC)
}

# test function
# options(repr.plot.width=6, repr.plot.height=4)
# NMDS_plot_group(otu_V, Location)

################################################################################
### 2) NMDS plots with groups, shapes

NMDS_group_shape = function(otu_V, group, shape_var){

    # Get OTU data for NMDS
    counts <- sapply(otu_V, is.numeric)                                               # Select numeric columns
    otu_d <- data.matrix(otu_V[,counts])                                              # Filter from OTU tab          #length(Restor) #length(Plant) #dim(otu_d)
    biomT <-t(otu_d)                                                                  # Transpose so cols are vars

    # Get NMDS scores
    bT_mds<-metaMDS(biomT, distance="bray",k=3, trymax=10);                          # run NMDS
    bT_mds_DF = data.frame(NMDS1=bT_mds$points[,1], NMDS2=bT_mds$points[,2], group)  # Make data frame of scores
    stress=round(bT_mds$stress, digits = 3)                                                                      
    # return(bT_mds_DF)  
    
    # Test ADONIS models, extract params
    bT_adonis <- adonis(biomT  ~ group + shape_var, permutations=99, method="bray"); # data=env,
    R2 <- round(bT_adonis$aov.tab$R2[1], digits=3)
    P <- bT_adonis$aov.tab$"Pr(>F)"[1]
    
    # PLOT NMDS by group color, shape
    pA <- ggplot(bT_mds_DF) +
        geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=group, shape=shape_var), alpha = 0.8) + #+ geom_point(alpha = 0.8) +
        stat_ellipse(mapping = aes(x=NMDS1, y=NMDS2, color=group), alpha = 0.8) #  level=0.95, + ggtitle(paste("ADONIS R2 =", R2, "Stress:", stress));
    pB <- pA + scale_color_manual(values = c(site_col)) + scale_shape_manual(values=c(20, 17), labels = c("Reference", "Restored"))
    pC <- pB + theme(legend.title=element_blank())       + 
       theme_bw() + theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 6),
                          legend.key.height = unit(0.5, "cm"),
                          legend.key.width = unit(0.5, "cm"),
                          legend.spacing.y = unit(0, "cm"),
                          legend.margin = margin(5,0,0,0),
                          legend.box.margin = margin(0,0,0,-5),
                          axis.title.y = element_text(size = 10, margin = margin(t = 0, r = -2.5, b = 0, l = 0)),
                          axis.title.x = element_text(size = 10))
    
    return(pC)
}

# test function
# options(repr.plot.width=6, repr.plot.height=4)
# NMDS_group_shape(otu_V, Location, Restor)

################################################################################
### 3) NMDS plots with groups only, env fit arrows

NMDS_group_fit = function(otu_V, group, Env_chem){

    # Get OTU data for NMDS
    counts <- sapply(otu_V, is.numeric)                                               # Select numeric columns
    otu_d <- data.matrix(otu_V[,counts])                                              # Filter from OTU tab          #length(Restor) #length(Plant) #dim(otu_d)
    biomT <-t(otu_d)                                                                  # Transpose so cols are vars

    # Get NMDS scores
    bT_mds<-metaMDS(biomT, distance="bray",k=3, trymax=10);                          # run NMDS
    bT_mds_DF = data.frame(NMDS1=bT_mds$points[,1], NMDS2=bT_mds$points[,2], group)  # Make data frame of scores
    stress=round(bT_mds$stress, digits = 3)                                                                      
    # return(bT_mds_DF)  
    
    # Test ADONIS models, extract params
    bT_adonis <- adonis(biomT  ~ group, permutations=99, method="bray"); # data=env,
    R2 <- round(bT_adonis$aov.tab$R2[1], digits=3)
    P <- bT_adonis$aov.tab$"Pr(>F)"[1]
    
    # PLOT NMDS by group color, shape
    pA <- ggplot(bT_mds_DF) +
        geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=group, shape=shape_var), alpha = 0.8) + #+ geom_point(alpha = 0.8) +
        stat_ellipse(mapping = aes(x=NMDS1, y=NMDS2, color=group), alpha = 0.8) #  level=0.95, + ggtitle(paste("ADONIS R2 =", R2, "Stress:", stress));
    pB <- pA + scale_color_manual(values = c(site_col)) + scale_shape_manual(values=c(20, 17), labels = c("Reference", "Restored"))
    pC <- pB + theme(legend.title=element_blank())       + 
       theme_bw() + theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 6),
                          legend.key.height = unit(0.5, "cm"),
                          legend.key.width = unit(0.5, "cm"),
                          legend.spacing.y = unit(0, "cm"),
                          legend.margin = margin(5,0,0,0),
                          legend.box.margin = margin(0,0,0,-5),
                          axis.title.y = element_text(size = 10, margin = margin(t = 0, r = -2.5, b = 0, l = 0)),
                          axis.title.x = element_text(size = 10))
    
    # Get Env. fits 
    fitE <- envfit(bT_mds, Env_chem, perm=0, na.rm=TRUE)
    scrsE <-as.data.frame(scores(fitE, display="vectors")) # Get scores
    scrsE <-cbind(scrsE, var=row.names(scrsE))            # add Env var name   
    scrsE <-data.frame(scrsE)                             # make DF

    ##### Rescale env. scores
    range_x <- max(scrsO[,1]) - min(scrsO[,1])
    range_y <- max(scrsO[,2]) - min(scrsO[,2])
    scrsE$sNMDS1<-scrsE[,1]*(range_x/2)
    scrsE$sNMDS2<-scrsE[,2]*(range_y/2)
    
    pD <- pC + coord_fixed()+
      geom_segment(data = scrsE, aes(x = 0, xend = sNMDS1, y = 0, yend = sNMDS2), alpha=0.6,
            arrow = arrow(length = unit(0.2, "cm")), colour = "orange")                          + # +
      geom_text(data = scrsE, aes(x = sNMDS1, y = sNMDS2, label = var), size = 2.9) 
    
    
    return(pD)
}

# test function
# options(repr.plot.width=6, repr.plot.height=4)
# NMDS_group_fit(otu_V, group, Env_chem)

################################################################################
### 4) NMDS plots with groups only, env fit arrows

NMDS_group_shape_fit = function(otu_V, group, shape_var, envChem){

    # Get OTU data for NMDS
    counts <- sapply(otu_V, is.numeric)                                               # Select numeric columns
    otu_d <- data.matrix(otu_V[,counts])                                              # Filter from OTU tab          #length(Restor) #length(Plant) #dim(otu_d)
    biomT <-t(otu_d)                                                                  # Transpose so cols are vars

    # Get NMDS scores
    bT_mds<-metaMDS(biomT, distance="bray",k=3, trymax=10);                          # run NMDS
    bT_mds_DF = data.frame(NMDS1=bT_mds$points[,1], NMDS2=bT_mds$points[,2], group)  # Make data frame of scores
    stress=round(bT_mds$stress, digits = 3)                                                                      
    # return(bT_mds_DF)  
    
    # Test ADONIS models, extract params
    bT_adonis <- adonis(biomT  ~ group + shape_var, permutations=99, method="bray"); # data=env,
    R2 <- round(bT_adonis$aov.tab$R2[1], digits=3)
    P <- bT_adonis$aov.tab$"Pr(>F)"[1]
    
    # PLOT NMDS by group color, shape
    pA <- ggplot(bT_mds_DF) +
        geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=group, shape=shape_var), alpha = 0.8) + #+ geom_point(alpha = 0.8) +
        stat_ellipse(mapping = aes(x=NMDS1, y=NMDS2, color=group), alpha = 0.8) #  level=0.95, + ggtitle(paste("ADONIS R2 =", R2, "Stress:", stress));
    pB <- pA + scale_color_manual(values = c(site_col)) + scale_shape_manual(values=c(20, 17), labels = c("Reference", "Restored"))
    pC <- pB + theme(legend.title=element_blank())       + 
       theme_bw() + theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 6),
                          legend.key.height = unit(0.5, "cm"),
                          legend.key.width = unit(0.5, "cm"),
                          legend.spacing.y = unit(0, "cm"),
                          legend.margin = margin(5,0,0,0),
                          legend.box.margin = margin(0,0,0,-5),
                          axis.title.y = element_text(size = 10, margin = margin(t = 0, r = -2.5, b = 0, l = 0)),
                          axis.title.x = element_text(size = 10))
    
    # Get Env. fits 
    fitE <- envfit(bT_mds, Env_chem, perm=0, na.rm=TRUE)
    scrsE <-as.data.frame(scores(fitE, display="vectors")) # Get scores
    scrsE <-cbind(scrsE, var=row.names(scrsE))            # add Env var name   
    scrsE <-data.frame(scrsE)                             # make DF

    ##### Rescale env. scores
    range_x <- max(bT_mds_DF[,1]) - min(bT_mds_DF[,1])
    range_y <- max(bT_mds_DF[,2]) - min(bT_mds_DF[,2])
    scrsE$sNMDS1<-scrsE[,1]*(range_x/2)
    scrsE$sNMDS2<-scrsE[,2]*(range_y/2)
    scrsE$var <- gsub('_s', '', scrsE$var) # Clean up NH4 and NO3 text
    scrsE <- subset(scrsE, var != "C:N") # Remove C.N from plot
    
    pD <- pC + 
        # coord_fixed()+
      geom_segment(data = scrsE, aes(x = 0, xend = sNMDS1, y = 0, yend = sNMDS2), alpha=0.6,
            arrow = arrow(length = unit(0.2, "cm")), colour = "orange")                          + # +
      geom_text(data = scrsE, aes(x = sNMDS1, y = sNMDS2, label = var), size = 2.9) 
    
    
    return(pD)
    return(scrsE)
}

# test function
# options(repr.plot.width=6, repr.plot.height=4)
# NMDS_group_shape_fit(otu_V, Location, Restor, Env_chem)
