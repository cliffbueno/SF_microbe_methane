 Lasso semi-auto for SF OTU guilds v1.12
#-  runs under r > 3.6 (glmnet for LASSO)

# Assume have loaded modules 0-7
# likely pollulted with examples, possibly cruft


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
