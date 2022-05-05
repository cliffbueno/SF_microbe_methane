# Plot heatmap of 16s corr by Tax Rank 

# Heatmap module OUTLINE

# 1) Show rank corrs on heatmap, w/ site colors 
#     a) Get data for heatmap of corrRanks                 Tax_corrU_r <- ggHeat_corrRanks_data(Taxons3d, Tax_corrU)
#     b) Plot corrRanks heatmap function                   heatmap_corrRanks(Tax_corrU_r, "PiYG")
#     c) Get taxonomy color bar data                       Color_bar_d <- TaxColorBar_dat(Taxons3d, Tax_corrU) 
#     d) Plot taxonomy color bar                           plot_TaxColorbar(Color_bar_d)
#  *  e) Combine Heatmap and Taxon color bar plot          Tax_heat_corrRanks(Taxons3d, Tax_corrU, "PiYG")   
#  * User function wrapper

# 2) Filter corrRanks by cuttoff for a tax rank    



### Test data ### 
# source("Import_Silva_OTU_data4plots_v0.1.R")
# SILVA OTU table format needs repair here! But probably not iTagger
# row.names(otu_V) <- otu_V[,"OTU"]                                # Make OTU row names -- required!!!
# otu_V <- otu_V[,-1]                                              # Dropping OTU column (1st, some )

# Reimport color legend from barplots
# Taxons3d <-read.table("data/colors/Silva_taxonomy_color_pallete.txt", sep="\t", header = T)

# Import CH4corrRanks at OTU level -- some glitch otherwise 
# CH4corrRanksOTU <- read.table("taxon_analysis/correlations/CH4corrRanks_Silva.txt", sep = '\t')
# head(CH4corrRanksOTU)

# Drop OTU vars or heatmap will never render - this appears to be finnicky part re last rank shown
# CH4corrRanks <-CH4corrRanksOTU[,1:(ncol(CH4corrRanksOTU)-2)]           # if need be, drop OTU data for unique at genus / Consensus.lineage
# Tax_corrU <- unique(CH4corrRanks)                                      # get only unique rows                  
# dim(Tax_corrU); 
# head(Tax_corrU) 

# NOTE Tax_corrU is GENUS LEVEL tax corr ranks (OTU data removed)



# Import packagaes
suppressMessages(library(Hmisc))
library(ggplot2)
library(RColorBrewer)
library(reshape2)
suppressMessages(library(cowplot))
theme_set(theme_grey())             



# import needed modules
# source("modules/6_Corr_ranks_module_v0.3.4_strip.R")
source("../modules/3_OTU_subsetting_modules_v.0.4_strip.r")
source("../modules/5_OTU_barplots_module_v0.4.R")  # uses ColorMap_other() to map colors

# Define heatmap colors -- defaults to slateblue-> red
grey_red <- c("slategray","white","red4")


#######################################################################################
#     a) Get data for heatmap of corrRanks                 Tax_corrU_r <- ggHeat_corrRanks_data(Taxons3d, Tax_corrU)
#######################################################################################
# uses ColorMap_other() imported from 5_OTU_barplots_module


ggHeat_corrRanks_data = function(Taxons3d, Tax_corrU, tax_sort="Consensus.lineage") {                       # uses colorMap_other function 
                                                                              #  add row_var, Tax_sort?
    ## Input data using colorMap function
    TaxCol_otu <-ColorMap_other(Taxons3d, Tax_corrU)                          # get Colormap merged corrRanks 
    TaxCol_otu <-TaxCol_otu[rev(order(TaxCol_otu[,tax_sort])),]    # sort by Cons.lin, rev for stacked bar
    row.names(TaxCol_otu) <- make.unique(as.character(TaxCol_otu[,"Genus"]))  # Make unique rownames for genera


    # Get only count data from OTU color table
    Count_cols <- sapply(TaxCol_otu, is.numeric)                               # are columns numeric?
    TaxCol_otu_d <-TaxCol_otu[Count_cols]                                      # get df with only numeric columns
    TaxCol_otu_d <-TaxCol_otu_d[,-1]                                           # drops Index

    # get heatmap plotting data
    Tax_corrU_r <-data.matrix(TaxCol_otu_d[,2:6])                              # select Tax Ranks to display
    #return(TaxCol_otu)
    return(Tax_corrU_r)
}

# Apply function test
# Tax_corrU_r <-ggHeat_corrRanks_data(Taxons3d, Tax_corrU)
# head(Tax_corrU_r)





#######################################################################################
#     b) Plot corrRanks heatmap function                   heatmap_corrRanks(Tax_corrU_r, "PiYG")
#######################################################################################
# only heatmap, no taxa decoration


heatmap_corrRanks = function(Tax_corrU_r, heat_cols = def_cols, ylab=F, xlab=T) {

    # Melt data for ggplot                                                                 # here if compositing subsets
    Tax_corrU_r_Melt<- melt(Tax_corrU_r)                                                   # Melt data
    colnames(Tax_corrU_r_Melt)[3] <- "Correlation"                                         # name corr. col "Correlation"

    # Get colors from heatmap
    def_cols<-c("steelblue", "white", "darkred")    
     
    # Make ggplot heatmap
    p <- ggplot(Tax_corrU_r_Melt, aes(Var2, Var1)) + geom_tile(aes(fill=Correlation)) + 
          scale_fill_gradientn(colours=heat_cols, na.value = "white") +  
          #scale_fill_gradientn(colours=heat_cols)                              +            # Fill with heat_cols  #breaks = c(-0.6,-0.3, 0, 0.3, 0.6),  # guide="legend"   

         # Axis params
         theme(axis.title.y = element_blank())                                +            # hide y axis title  
         scale_y_discrete(position="right")                                   +            # y labs to right
         theme(axis.title.x = element_blank())                                +            # hide x axis title                      
         scale_x_discrete(position = "bottom", labels = c("Phylum", "Class", "Order", "Family", "Genus"), expand = c(0, 0))                               +            # x labs on top
         theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5))             +            # rotate x labs 
         
         # Panel params
         theme(panel.background = element_blank(),                                         # Remove panel borders 
           panel.grid.major = element_blank(),                                             # and grid lines  
           panel.grid.minor = element_blank())                                +            # Remove legend entirely
     
         # Legend params
        labs(fill = "r") +
         theme(legend.position = "top", legend.justification = "right")                                    +            # legend on bottom
         guides(fill = guide_colourbar(barwidth = 3.5, barheight = 0.5, title.position = "right")) + # legend size and position
         theme(legend.title = element_text(size = 8),
               legend.text = element_text(size = 6),
               legend.margin = margin(0,5,0,5),
               legend.box.margin = margin(-5,0,-5,0),
               plot.margin = unit(c(0.1, -0.2, 0.1, 0.1), "cm"))                                             # no legend title
        
         ifelse(ylab==T, p <- p,                                                           # hide y labels ? 
                p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()))
         ifelse(xlab==T, p <- p,                                                           # hide y labels ? 
                p <- p + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()))
 
    p
}

# Test function
# options(repr.plot.width=1.5, repr.plot.height=6) 
#dim(Tax_corrU_r)

# Define heatmap colors -- defaults to slateblue-> red
# grey_red <- c("slategray","white","red4")

# heatmap_corrRanks(Tax_corrU_r, heat_cols=grey_red)       # Don't need to define heat_cols, but can
                                                         # Can also use y_lab=T to get Genus names   

#### Guild version
# X axis on top
# Legend centered and on bottom, title on top
# Plot margins adjusted
heatmap_corrRanks_guild = function(Tax_corrU_r, heat_cols = def_cols, ylab=F, xlab=T) {

    # Melt data for ggplot                                                                 # here if compositing subsets
    Tax_corrU_r_Melt<- melt(Tax_corrU_r)                                                   # Melt data
    colnames(Tax_corrU_r_Melt)[3] <- "Correlation"                                         # name corr. col "Correlation"

    # Get colors from heatmap
    def_cols<-c("steelblue", "white", "darkred")    
     
    # Make ggplot heatmap
    p <- ggplot(Tax_corrU_r_Melt, aes(Var2, Var1)) + geom_tile(aes(fill=Correlation)) + 
          scale_fill_gradientn(colours=heat_cols, na.value = "white") +  
          #scale_fill_gradientn(colours=heat_cols)                              +            # Fill with heat_cols  #breaks = c(-0.6,-0.3, 0, 0.3, 0.6),  # guide="legend"   

         # Axis params
         theme(axis.title.y = element_blank())                                +            # hide y axis title  
         scale_y_discrete(position="right")                                   +            # y labs to right
         theme(axis.title.x = element_blank())                                +            # hide x axis title                      
         scale_x_discrete(position = "top", labels = c("Phylum", "Class", "Order", "Family", "Genus"), expand = c(0, 0))                               +            # x labs on top
         theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0),
               axis.text.x.top = element_text(vjust = 0.5))             +            # rotate x labs 
         
         # Panel params
         theme(panel.background = element_blank(),                                         # Remove panel borders 
           panel.grid.major = element_blank(),                                             # and grid lines  
           panel.grid.minor = element_blank())                                +            # Remove legend entirely
     
         # Legend params
        labs(fill = "r") +
         theme(legend.position = "bottom")                                    +            # legend on top
         guides(fill = guide_colourbar(barwidth = 3.5, barheight = 0.5, title.position = "top", title.hjust = 0.5)) + # legend size and position
         theme(legend.title = element_text(size = 8),
               legend.text = element_text(size = 6),
               legend.margin = margin(0,5,0,5),
               legend.box.margin = margin(-5,0,0,0),
               plot.margin = unit(c(0.1, -0.2, 0.1, -0.2), "cm"))                                             # no legend title
        
         ifelse(ylab==T, p <- p,                                                           # hide y labels ? 
                p <- p + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()))
         ifelse(xlab==T, p <- p,                                                           # hide y labels ? 
                p <- p + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()))
 
    p
}


#######################################################################################
#     c) Get taxonomy color bar data                       Color_bar_d <- TaxColorBar_dat(Taxons3d, Tax_corrU) 
#######################################################################################
# uses ColorMap_other() imported from 5_OTU_barplots_module


TaxColorBar_dat = function(Taxons3d, Tax_corrU, tax_sort="Consensus.lineage"){
                                                 
    # Get genus scale Taxonomy color data, from using colorMap_other merge 
    TaxCol <-ColorMap_other(Taxons3d, Tax_corrU)                                    # Colormap filtered corrRanks 
    TaxCol <-TaxCol[rev(order(TaxCol[,tax_sort])),]                      # sort by Consensus.lineage, reverse for bar stacking

    # Slice Taxonomy, color data 
    keep_cols <-c("newTax", "newColor")                                    # Get columns to keep from TaxColors, for plotting

    #keep_cols <-c("newTax", "Index", "newColor")                                    # Get columns to keep from TaxColors, for plotting
    Color_bar_d <- TaxCol[,keep_cols]                                               # data, only Keep cols   
    colnames(Color_bar_d)[1] <- "Taxonomy" 

    # Plotting data, y, x, colors 
    Color_bar_d["Genus"] <- as.factor(make.unique(rownames(TaxCol)))                # "y": Get genus names from rownames 
    Color_bar_d["Gen_ind"] <- seq(1:nrow(Color_bar_d))                              # Make sequential genus names to reorder factor
    Color_bar_d["Genus"] <- reorder(Color_bar_d[,"Genus"], Color_bar_d[,"Gen_ind"]) # reorder "y" factor by index
    Color_bar_d["counts"] <- 1                                                      # "x": Fake 1 counts for barplot 
    return(Color_bar_d)
}

# Use TaxColorBar_dat function 
# Color_bar_d <- TaxColorBar_dat(Taxons3d, Tax_corrU)

#Color_bar_d <- TaxColorBar_dat(Taxons3d, Corr_filt_r)
# head(Color_bar_d)




#######################################################################################
#     d) Plot taxonomy color bar                           plot_TaxColorbar(Color_bar_d)
#######################################################################################

plot_TaxColorbar = function(Color_bar_d, ylab=F, alpha = 1){
    
    # Get colors from levels
    Tax_colU <-levels(Color_bar_d[,"newColor"]) 

    # Make barplot using geom_tile method
    t <-ggplot(Color_bar_d, aes(x=counts, y=Genus, fill=Taxonomy))                            # fill loc. 
    t2 <- t + geom_tile(position="identity", alpha = alpha) + scale_fill_manual(values=c(Tax_colU))           + # manual fill col2
          
          # Axis params      
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # Remove x labs
          theme(axis.title.y=element_blank()) +#, axis.text.y=element_blank(), axis.ticks.y=element_blank()) + # Remove y labs
          
          # Panel params
          theme(panel.background = element_blank(),                                           # Remove panel borders and grid lines  #
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())                                              +
          
          # Legend params   
          theme(legend.position="left")          +
          theme(legend.key.size = unit(0.5, "cm"))                                         +   # Shrink legend
          guides(fill = guide_legend(ncol=1))                                                  # 1 col legend
    
          ifelse(ylab==T, t2 <- t2,                                                           # hide y labels ? 
                t2 <- t2 + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()))
  
    t2 
}

# Note artifacts from using alpha = 0.8 above, or any alpha-- can erase but colors not matching orig...

# use plot_TaxColorbar function
# options(repr.plot.width=3, repr.plot.height=6) 
# plot_TaxColorbar(Color_bar_d)





#######################################################################################
#     e) Combine Heatmap and Taxon color bar plot          Tax_heat_corrRanks(Taxons3d, Tax_corrU, "PiYG")   
#######################################################################################
# User wrapper function for all of the above

Tax_heat_corrRanks = function(Taxons3d, Tax_corrU, heat_cols=default){           # uses cowplot::plot_grid 

    # Get heatmap of corrRanks
    Tax_corrU_r <- ggHeat_corrRanks_data(Taxons3d, Tax_corrU)           # use function to get corrRanks data
    heatRanks <- heatmap_corrRanks(Tax_corrU_r, heat_cols)                 # plotting function for heatmap of corrRanks

    # Get Taxonomy colorBar
    Color_bar_d <- TaxColorBar_dat(Taxons3d, Tax_corrU)                 # use function to get Tax colorBar data 
    TaxColors <- plot_TaxColorbar(Color_bar_d)                          # plotting function for Tax colorBar

    # Split Taxononmy colorBar into no-legend; legend
    TaxColors_NL <- TaxColors + theme(legend.position="none")
    TaxColors_Legend <-get_legend(TaxColors)

    # Composite plotting of components
    pg_hm <- plot_grid(heatRanks, TaxColors_NL, align="h", ncol=2, rel_widths = c(3.5, 1), axis = "rlbt")
    # pg_hmL <- plot_grid(pg_hm, TaxColors_Legend, ncol=2, rel_widths = c(1.3, 2))#, axis = "rlbt")  align="h"
    # return(pg_hmL)
    return(pg_hm)
}

# Test corrRanks heatmap w taxon color bar
# options(repr.plot.width=5, repr.plot.height=6) 
# Tax_heat_corrRanks(Taxons3d, Tax_corrU, heat_cols=grey_red)               # grey_red <- c("slategray","white","red4")
# ggsave('CH4_otu_corrheatmapUMX_0.1.pdf')
                        





#######################################################################################
# 2) Filter corrRanks by cuttoff for a tax rank                 Tax_corrU_r <- ggHeat_corrRanks_data(Taxons3d, Tax_corrU)
#######################################################################################
# 

corrRanks_filt = function(corrRanks, rank_var, r_cut){          # Simplified from v.3.1

    #corrRanks <- CH4corrRanksOTU_Delta 
    #chunk_var <- "Genus"
    rank_r <- paste0(rank_var,"_r")
    abs_rank_r <-paste0("abs_", rank_r)
    #r_cut <- 0.5

    # Get abs val of ranks
    corrRanks[abs_rank_r] <-abs(corrRanks[,rank_r])

    # Get corr filtered ranks 
    corrRanks_F <- corrRanks[corrRanks[,abs_rank_r] > r_cut,]         
    corrRanks_F <- corrRanks_F[,1:(ncol(corrRanks_F)-1)]
    row.names(corrRanks_F) <-  make.unique(as.character(corrRanks_F[,rank_var]))  # Rename rows
    return(corrRanks_F)
    
    
    # Select column of interest, make DF
    #r_Filt_rank <-data.frame(corrRanks_F[, rank_var])
    #colnames(r_Filt_rank)[1] <- rank_var
    #r_Filt_rank[rank_r]<-corrRanks_F[,rank_r]
    #return(r_Filt_rank)
    return(corrRanks_F)
}

# Use corrRanks_Filt function - GENUS
# Corr_filt_r <- unique(corrRanks_filt(Tax_corrU, "Genus", 0.5))
# Tax_heat_corrRanks(Taxons3d, Corr_filt_r, heat_cols=grey_red)               # grey_red <- c("slategray","white","red4")

# No option to print taxa names?  Seems like there should be.  maybe in guilds.
# Corr_filt_r <- unique(corrRanks_filt(CH4corrRanksOTU, "OTU", 0.6))