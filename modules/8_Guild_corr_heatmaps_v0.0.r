
# Plot heatmap of GUILDS corr by Tax Rank & corr filtered

# This module could use some WORK & cleanup
# could use better var naming & encapsulation of complete user functions
# currently these are shown in a Jupyter notebook as sequence of steps
# also depends on other module imports, user selects correlation cutoff.

# esp helper functions are poorly named & neeed to be combined: b, c, d 
# unfinished in several ways, user pastes plots in f,g.
# Plots have goofy OTU labels, flexible variables are hardcoded, etc.


# Guild heatmap module OUTLINE
#     a) Merge guild OTUs, taxonomy, colors     -- Guild_color_data(Guild_OTUs, Guild_colors, otu_V, sort_tax ="Consensus.lineage", agg_tax="Genus")
#     b) Merge guild color bar w Taxcorr        -- Guild_corr_d(Guild_bar_dU, Tax_corrU)
#     c) Make guild barplot                     -- Guild_color_bar(Corr_filt_Guild_d)
#     d) guild data only in guild corr data     -- filter_corr_guild_data(Corr_filt_Guild_d)
#     e) Extract Heatmap log2 abundance data    -- log2_HM_Fabund_data(otu_V, Corr_filt_Guild_Gd, "Genus", "G_order") 
#     f) Plot log2 abundance heatmap            -- log2_abundHM(HM_d, ylab=F) 
#  *  g) Plot log2 with colorbar    -- plot_guild_log2_abund_site_colorbar(otu_V, Meta_iTag, site_colors, Corr_filt_Guild_Gd)
#  *  h) Plot guild corr tax heat abundance    plot_guild_corr_tax_heat_abund(Taxons3d, Corr_filt_Guild_Gd, Tax_corrU_r, heat_cols)


#  * h) should be user function, but also need g.  Currently user also does a-d.


# library(ggplot2, hmisc, reshape,  etc...)
# these are user imported with other modules: subsetting & barplots



##############################
# test data
# source("../modules/Import_Silva_OTU_data4plots_v0.1.R")
# SILVA OTU table format needs repair here! But probably not iTagger
# row.names(otu_V) <- otu_V[,"OTU"]                                # Make OTU row names -- required!!!
# otu_V <- otu_V[,-1]                                              # Dropping OTU column (1st, some )


# Get guild OTUs
# source("../modules/3_OTU_subsetting_modules_v.0.4_strip.r")
# Guild_OTUs <- Get_16S_Guilds(otu_V)
# dim(Guild_OTUs); head(Guild_OTUs)
# head(Guild_colors)

# Reimport color legend from barplots
#Taxons3d <-read.table("data/colors/Silva_taxonomy_color_pallete.txt", sep="\t", header = T)

# Get correlation ranks
# CH4corrRanksOTU <- read.table("taxon_analysis/correlations/CH4corrRanks_Silva.txt", sep = '\t')

# Genus level only corr -- Drop OTU vars or heatmap will never render
# CH4corrRanks <-CH4corrRanksOTU[,1:(ncol(CH4corrRanksOTU)-2)]           # if need be, drop OTU data for unique at genus / Consensus.lineage
# Tax_corrU <- unique(CH4corrRanks)                                      # get only unique rows       

# correlation filtered -- needs to be Genus
# Corr_filt_r <- unique(corrRanks_filt(Tax_corrU, "Genus", 0.5))


##############################


#######################################################################################
#     a) Merge guild OTUs, taxonomy, colors     -- Guild_color_data(Guild_OTUs, Guild_colors, otu_V, sort_tax ="Consensus.lineage", agg_tax="Genus")
#######################################################################################

Guild_color_data = function(Guild_OTUs, Guild_colors, otu_t, sort_tax ="Consensus.lineage", agg_tax="Genus") {

    # Get OTU, consensus lineage from input OTU table
    # otu_t <- otu_V                                                              # abstraction for function
    # keep_cols <-c("Consensus.lineage", "Genus")                                 # keep Consens. Lineage and Genus   
    keep_cols <-c(sort_tax, agg_tax)                                 # keep Consens. Lineage and Genus   
    All_OTUs <- otu_t[keep_cols]                                                # keep " in input OTU table
    All_OTUs["OTU"]<-row.names(All_OTUs)   #All_OTUs                            # get OTU Ids from rownames

    # Merge OTU data with Guilds and Guild colors
    Allotu_Guilds <- merge(All_OTUs, Guild_OTUs, all.x =TRUE)                   # Merge Guilds with all OTUs
    Allotu_Guilds_colors <- merge(Allotu_Guilds, Guild_colors, all.x=TRUE)      # Merge with Guild color palette

    # Replace missing colors (NA from merge)
    Allotu_Guilds_colors["G_color"] <- ifelse(is.na(Allotu_Guilds_colors$color)==TRUE,    
        "#FFFFFF", as.character(Allotu_Guilds_colors$color)) 

    # drop OTUs in Guild data, Get genus level data
    Guild_d <-c("Guild", "G_color","Consensus.lineage", "Genus")                          # get columns to keep
    Guild_bar_d <- Allotu_Guilds_colors[Guild_d]                                          # seletct columns  
    Guild_bar_dU <- unique(Guild_bar_d)                                                   # get unique at Genus

    # Prepare data for plot
    Guild_bar_dU["counts"] <- 1
    return(Guild_bar_dU)
}    

# Test Guild color data function, 
# note sort_tax and agg_tax are defaults to sort on, and level of aggregation used
# Guild_bar_dU <- Guild_color_data(Guild_OTUs, Guild_colors, otu_V, sort_tax ="Consensus.lineage", agg_tax="Genus")
# head(Guild_bar_dU)



#######################################################################################
#     b) Merge guild color bar w Taxcorr     -- Guild_corr_data(Guild_bar_dU, Tax_corrU)
#######################################################################################

# Merge guild bar data with TaxCorr object (filtered or unfiltered)
# Plotting below looks to work using Tax_corrU OR Corr_filt_r
# REALLY NEEDS BETTER FUNCTION NAME


Guild_corr_d = function(Guild_bar_d, Tax_corrU) {

    # Merge Guild data with Tax Corr data (here corr filtered)  
    Corr_Guild_d <- merge(Guild_bar_d, Tax_corrU)                                       # merge by Genus

    # Reorder data so Corr_filt_r data is first (for using plotting functions)
    col_order <-c(names(Tax_corrU), names(Guild_bar_d))                                 # Get column ordering 
    Corr_Guild_d <-Corr_Guild_d[col_order]                                              # Reorder columns
    #length(unique(Corr_filt_Guild_d$Genus))

    # make Genus names unique (to match plotting in corrRanks), sort and order Genus factor
    Corr_Guild_d$Genus <- make.unique(as.character(Corr_Guild_d$Genus))                 # Make unique
    Corr_Guild_d <- Corr_Guild_d[rev(order(Corr_Guild_d$Consensus.lineage)),]           # reverse sort by lineage
    Corr_Guild_d["order"] <-seq(1:nrow(Corr_Guild_d))                                   # get Index for Genus factor
    Corr_Guild_d$Genus <- reorder(Corr_Guild_d$Genus, Corr_Guild_d$order)               # Reorder Genus factor by order
    return(Corr_Guild_d)
}


# Test function Guild_corr_d
# Corr_Guild_d <- Guild_corr_d(Guild_bar_dU, Tax_corrU)
# Corr_filt_Guild_d <- Guild_corr_d(Guild_bar_dU, Corr_filt_r)

# dim(Corr_Guild_d); head(Corr_Guild_d)
# # dim(Corr_filt_Guild_d); head(Corr_filt_Guild_d)

# Does updated dataset work with older functions  -- YES
# options(repr.plot.width=4, repr.plot.height=6) 
# Tax_heat_corrRanks(Taxons3d, Corr_Guild_d, "PiYG")
# Tax_heat_corrRanks(Taxons3d, Corr_filt_Guild_d, heat_cols=grey_red)
#corr_heat_colorIn(Taxons3d, Corr_filt_Guild_d)


#######################################################################################
#     c) Make guild barplot     -- Guild_color_bar(Corr_filt_Guild_d)
#######################################################################################

Guild_color_bar = function(Corr_Guild_d){

    # get only guilds and color data from merged data                        # Here not full compliment of guilds 
    keep<- c("Guild", "G_color")                                             # Keep columns
    Guild_colr_d <-Corr_Guild_d[keep]                                   # Filter merged data on keep

    # Merge data with Guild ordering
    Guild_color_O <- unique(merge(Guild_colr_d, Guild_colors, by="Guild"))   # Merge Guild color_d (fewer), Guild_colors
    Guild_color_O <- Guild_color_O[order(Guild_color_O$G_index),]            # Sort by ordering index
    Guild_cols <-Guild_color_O$G_color                                       # Get final colors from ordered

    # Make barplot using geom_tile method
    g <-ggplot(Corr_Guild_d, aes(x=counts, y=Genus, fill=Guild))                            # fill loc. 
    g2 <- g + geom_tile(position="identity", alpha = 0.8) + scale_fill_manual(values=c(Guild_cols))  +   
         # Axis params
         theme(axis.title.y = element_blank(), axis.text.y=element_blank(),                # hide y - axis labels  
                axis.ticks.y=element_blank())                                +
         theme(axis.title.x = element_blank())                                +            # hide x axis title                      
         scale_x_discrete(position = "top") +
         theme(axis.text.x = element_text(angle = 90, hjust = 1))             +

         # Panel params
         theme(panel.background = element_blank(),                                         # Remove panel borders and grid lines  #
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())                                #+                                                # Remove legend entirely
     
         # Legend params
         #theme(legend.position = "bottom")                                    + 
         #theme(legend.title = element_blank())
g2
}

# options(repr.plot.width=2, repr.plot.height=4) 
# Guild_color_bar(Corr_Guild_d)
# Guild_color_bar(Corr_filt_Guild_d)




#######################################################################################
#     d) guild data only in guild corr data     -- filter_corr_guild_data(Corr_filt_Guild_d)
#######################################################################################

# makes sure only Guild matching data is in corr-filtered data.
# NEEDED, not wrapped, seems like a patch for oversight in c) ? 


filter_corr_guild_data = function(Corr_filt_Guild_d){
    
    Corr_filt_Guild_Gd <-Corr_filt_Guild_d[!is.na(Corr_filt_Guild_d$Guild), ]   # Filter by guilds non-NA only
    #head(Corr_filt_Guild_d)

    Corr_filt_Guild_Gd <-Corr_filt_Guild_Gd[order(Corr_filt_Guild_Gd$Guild),]     # Sort df by Guilds 
    Corr_filt_Guild_Gd["G_order"]<-rep(1:nrow(Corr_filt_Guild_Gd))               # Make guild ordering index 
    #head(Corr_filt_Guild_d)

    # reorder genus levels by guild ordering
    Corr_filt_Guild_Gd$Genus <- reorder(Corr_filt_Guild_Gd$Genus, rev(Corr_filt_Guild_Gd$G_order)) 
    return(Corr_filt_Guild_Gd)
}

# test the function
# Corr_filt_Guild_Gd <- filter_corr_guild_data(Corr_filt_Guild_d)




#######################################################################################
#     e) Extract Heatmap log2 abundance data     -- log2_HM_Fabund_data(otu_V, Corr_filt_Guild_Gd, "Genus", "G_order") 
#######################################################################################

# Function to extract Heatmap log2 abundance data
log2_HM_Fabund_data = function(otu_t, filt_otu, agg_var, order_var){ 

    # aggregate OTU table, get filter list by agg_var
    otu_t[agg_var]<-make.unique(as.character(otu_t[,agg_var]))
    
    otu_agg <- agg_by_cat(otu_t, agg_var)                                          # Aggregate to agg_var           (e.g. Genus) 
    otu_f <- data.frame(filt_otu[agg_var])                                         # agg_var levels from OTU filter (e.g. Genus)
    colnames(otu_f) <- agg_var                                                    # rename DF col for agg_var

    # Filter OTU table down 
    otu_agg_f <- merge(otu_f, otu_agg, all.x=TRUE)                                            # slice OTU table using filter df (otu_f)                     
   
    # Separate numeric data and Taxonomy in OTU_table
    count_cols <- sapply(otu_agg_f, is.numeric)                                   # Select numeric columns
    otu_d <- data.frame(otu_agg_f[,count_cols])                                   # Keep numeric cols                      
    Tax_data <- otu_agg_f[(count_cols == FALSE)]                                  # Save taxonomy data                     

    # Fill zero counts with two prior to log2 transform
    otu_d[otu_d==0] <- 2                                                           # replace 0 with 2-- now log2(2) = 0
    row.names(otu_d) <- Tax_data[,1]                                              # add back rownames

    # Log 2 and scale (z-score) data frame
    otu_dL2 <- log2(otu_d)                                                         # log2 all data 
    otu_dL2z <- t(scale(t(otu_dL2)))                                              # Center and scale rows (transpose twice)
    otu_dL2zd <- data.frame(otu_dL2z)                                             # make DF
    otu_dL2zd[agg_var] <- row.names(otu_dL2zd)                                    # add back Genus before merge

    # Get Genus ordering cols, merge
    keepC <- c(agg_var, order_var)#, "G_color")                                    # Get ordering vars to keep
    G_sort <- filt_otu[keepC]                                                      # Get ordering DF
    #G_sort["rev_ord"] <- rev(G_sort[,order_var])
    otu_dL2z_S <- merge(G_sort, otu_dL2zd)                                        # merge ordering, Z-score data

    # Reorder Genus factor after re-paste                                         # drops non-included genera from factor)
    otu_dL2z_S[agg_var] <- factor(otu_dL2z_S[,agg_var])                           # re-paste factor
    #otu_dL2z_S[agg_var] <- reorder(otu_dL2z_S[,agg_var], otu_dL2z_S[,"rev_ord"])  # Reorder Genus factor by G_order
    otu_dL2z_S <- otu_dL2z_S[,-2]                                                 # drop G_order                 
    otu_dL2z_S <- otu_dL2z_S[,-2]                                                 # drop rev_order                 

    return(otu_dL2z_S)
}

# test function to prepare heatmap data
# HM_d <-log2_HM_Fabund_data(otu_V, Corr_filt_Guild_Gd, "Genus", "G_order")          # log2_HM_Fabund_data(otu_t, filt_otu, agg_var, order_var)

#head(Corr_filt_Guild_d)
#HM_d <-log2_HM_Fabund_data(otu_V, Corr_filt_Guild_d, "Genus", "order")          # log2_HM_Fabund_data(otu_t, filt_otu, agg_var, order_var)

# dim(HM_d); dim(Corr_filt_Guild_Gd)



#######################################################################################
#     f) Plot log2 abundance heatmap     -- log2_abundHM(HM_d, ylab=F) 
#######################################################################################


log2_abundHM = function(HM_d, heat_cols=def_cols, ylab=T, xlab=F){

    # prepare data (melt, default colors)
    HM_dM <-melt(HM_d)                                                  # melt data before plotting                                                                  
    def_cols <- c("steelblue", "white", "darkred")                      # Get default heatmap colors 
    #def_cols <- c("blue", "white", "red")                      # Get default heatmap colors 

    
    # make ggplot
    p <- ggplot(HM_dM, aes(variable, Genus)) + geom_tile(aes(fill=value))     + 
          scale_fill_gradientn(colours=heat_cols, lim=c(-3,3), na.value = "white") +            

         # Axis params
         theme(axis.title.y = element_blank())                                +            # hide y axis title  
         scale_y_discrete(position="right")  + #,                                                # y labs to right
          #                limits = rev(levels(HM_d[,"Genus"])))               +            # reverse y ordering
         theme(axis.title.x = element_blank())                                +            # hide x axis title                      
         scale_x_discrete(position = "bottom")                                +            # x labs on top
         theme(axis.text.x = element_text(angle = 90, hjust = 1))             +            # rotate x labs 
       
         # Panel params
         theme(panel.background = element_blank(),                                         # Remove panel borders 
           panel.grid.major = element_blank(),                                             # and grid lines  
           panel.grid.minor = element_blank())                                +            # Remove legend entirely
     
         # Legend params
         theme(legend.position = "bottom")                                    +            # legend on bottom
         theme(legend.title = element_blank())                                             # no legend title
        
         ifelse(ylab==T, p <- p,                                                           # hide y labels ? 
                p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()))
         ifelse(xlab==T, p <- p,                                                           # hide y labels ? 
                p <- p + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()))
    p
}

# Test plot

# options(repr.plot.width=6, repr.plot.height=4) 
# log2_abundHM(HM_d)         #, ylab=F) drops y labels 
#log2_abundHM(HM_d, ylab=F) # drops y labels


#######################################################################################
#     g) Plot log2 with colorbar    -- plot_guild_log2_abund_site_colorbar(otu_V, Meta_iTag, site_colors, Corr_filt_Guild_Gd)
#######################################################################################

# Note this is hardcoded sorting, but ability to make more flexible...
# E.g. could change site ordering, etc.  Color bar as graph, etc.


plot_guild_log2_abund_site_colorbar = function(otu_V, Meta_iTag, site_colors, Corr_filt_Guild_Gd){
    # Get site color bar
    sb <-site_colbar(Meta_iTag,"Sample", "Location", site_colors, "EWsiteHyd_index", plot="")
    #sb <-site_colbar(Meta_iTag,"Sample", "Location", site_colors, "EWsiteHyd_index", plot="graph")

    # Get abundance hm
    HM_d <-log2_HM_Fabund_data(otu_V, Corr_filt_Guild_Gd, "Genus", "G_order")          # log2_HM_Fabund_data(otu_t, filt_otu, agg_var, order_var)
    #HM_dS <- sort_otu_agg_by_meta(HM_d, Meta_iTag, "Sample", "CH4_logn1")             # Sort agg by (Fxn. n) on agg

    heat_abund <-log2_abundHM(HM_d, ylab=T)

    # Composite plotting of components
    pg_hm <- plot_grid(sb, heat_abund, align="v", nrow=2, rel_heights = c(1,5), axis = "rlbt") 
    return(pg_hm)
}

# options(repr.plot.width=6, repr.plot.height=5)
# plot_guild_log2_abund_site_colorbar(otu_V, Meta_iTag, site_colors, Corr_filt_Guild_Gd)


#######################################################################################
#     h) Plot guild corr tax heat abundance    plot_guild_corr_tax_heat_abund(Taxons3d, Corr_filt_Guild_Gd, Tax_corrU_r, heat_cols)
#######################################################################################


plot_guild_corr_tax_heat_abund = function(Taxons3d, Corr_filt_Guild_Gd, Tax_corrU_r, heat_cols){

    # Get heatmap of corrRanks
    Tax_corrU_r <- ggHeat_corrRanks_data(Taxons3d, Corr_filt_Guild_Gd, tax_sort="G_order")           # use function to get corrRanks data
    heatRanks <- heatmap_corrRanks(Tax_corrU_r, heat_cols, ylab=F)       # plotting function for heatmap of corrRanks

    # Get Taxonomy colorBar
    Color_bar_d <- TaxColorBar_dat(Taxons3d, Corr_filt_Guild_Gd, tax_sort="G_order")                 # use function to get Tax colorBar data 
    TaxColors <- plot_TaxColorbar(Color_bar_d)                          # plotting function for Tax colorBar

    # Split Taxononmy colorBar into no-legend; legend
    TaxColors_NL <- TaxColors + theme(legend.position="none")
    TaxColors_Legend <-get_legend(TaxColors)

    # Get guilds bar
    Guild_bar <- Guild_color_bar(Corr_filt_Guild_Gd)
    Guild_bar_NL <- Guild_bar + theme(legend.position="none")
    Guild_bar_Legend <-get_legend(Guild_bar)

    # Get abundance hm
    HM_d <-log2_HM_Fabund_data(otu_V, Corr_filt_Guild_Gd, "Genus", "G_order")          # log2_HM_Fabund_data(otu_t, filt_otu, agg_var, order_var)
    heat_abund <-log2_abundHM(HM_d, ylab=T)

    # Composite plotting of components
    pg_hm <- plot_grid(Guild_bar_NL, heatRanks, TaxColors_NL, heat_abund,  align="h", ncol=4, rel_widths = c(1.5,3,1,12), axis = "rlbt") 

    # pg_hm <- plot_grid(Guild_bar_NL, heatRanks, heat_abund, align="h", ncol=3, rel_widths = c(1,2,8), axis = "rlbt") 
    # pg_hm <- plot_grid(Guild_bar_NL, heatRanks, TaxColors_NL, heat_abund, align="h", ncol=4, rel_widths = c(1,2,.7,6), axis = "rlbt") 

    # pg_hm <- plot_grid(Guild_bar_NL, heatRanks, heat_abund, TaxColors_NL, align="h", ncol=4, rel_widths = c(1.5,3,8,1), axis = "rlbt") 
    # pg_hm <- plot_grid(Guild_bar_NL, heatRanks, TaxColors_NL, heat_abund,  align="h", ncol=4, rel_widths = c(1.5,3,1,12), axis = "rlbt") 
    
    return(pg_hm)
}

# test function
# heat_cols <-grey_red

# options(repr.plot.width=6, repr.plot.height=5) 
# plot_guild_corr_tax_heat_abund(Taxons3d, Corr_filt_Guild_Gd, Tax_corrU_r, heat_cols)