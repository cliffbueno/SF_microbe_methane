#  Lasso semi-auto for SF OTU guilds v1.12
#  runs under r > 3.6 (glmnet for LASSO)

# Assume have loaded modules 0-8
# likely pollulted with examples, possibly cruft


########################################################################################
### 9) Combined lasSEM model selection plots

########################################################################################
### a) initial screening functions -- DEPRECATE ALL!!

### XXX DELETE XXX ###

# This is actually breaking where n_pred aren't all in table

# i) 
#plot_SEM_lasso_mods2 = function(sem_tab, sel_sort = 'aic_glm', top_n = 25) {
#
#    # Get color palettes
#    cols <- gg_color_hue(max(sem_tab$n_pred))                      # get base colors
#    n_pred <- sort(unique(sem_tab$n_pred))                         # don't make factor, crashes function
#    col_df <- data.frame(cols, n_pred)                             # color x n_preds DF
#    all_cols <- c(as.character(col_df$cols))                       # colors for all models
#     
#    ## Plot ALL models ##
#    all_models <- plot_SEM_lasso_mods(sem_tab, label = 'FALSE', colors = all_cols) 
#   
#    ## Model subsets ##
#    sem_top_n <- sem_tab[order(sem_tab[sel_sort]),]                # sort models by selected var    
#    sem_top_n <- sem_top_n[1:top_n,]                               # keep only top mods
    
#    # subset colors
#    n_pred_sub <- data.frame(n_pred = unique(sem_top_n$n_pred))    # only cols in subset
#    col_df_sub <- merge(col_df, n_pred_sub)                        # merge with color key (inner)
#    sub_cols <- c(as.character(col_df_sub$cols))                   # get colors for plot
    
#    ## Plot selected models ## 
#    select_models <- plot_SEM_lasso_mods(sem_top_n, label = 'TRUE', colors = sub_cols) #+ 
                     #theme(legend.position = "none")

#    ## Plot both figs ## 
#    dual_plot <- plot_grid(all_models, select_models, ncol=2, nrow =1)
#    return(dual_plot)
#}

# options(repr.plot.width=9, repr.plot.height=4)
# plot_SEM_lasso_mods2(semF, sel_sort = 'aic_glm', top_n = 25)
# head(asemF)
# unique(semF$n_pred)
# test edge case
# plot_SEM_lasso_mods2a(asemF, sel_sort = 'aic_glm', top_n = 200)


## ii) 
# plot_SEM_lasso_mods3 = function(lass_mod, group_colors, sem_tab, sel_sort = 'aic_glm', top_n = 25){
#     
#     coef_plot <- plot_lasso_coefs(lass_mod, group_colors)
#     models_plot <- plot_SEM_lasso_mods2(sem_tab, sel_sort = sel_sort, top_n = top_n)       
#     triple_plot <- plot_grid(coef_plot, models_plot, ncol=2, nrow =1, rel_widths = c(1, 2))  
#     #return(models_plot)
#     return(triple_plot)  
# }

# demonstrate function
# lasso_sel_sumplot <- plot_SEM_lasso_mods3(CH4_lass0, guild_colors, semF, sel_sort = 'aic_glm', top_n = 25)

# options(repr.plot.width=10, repr.plot.height=2.5)
# lasso_sel_sumplot# edge case test (def well below)
# alasso_sel_sumplot <- plot_SEM_lasso_mods3(aCH4_lass0, guild_colors, asemF, sel_sort = 'aic_glm', top_n = 200)


## ii) 
# alasso_sel_sumplotplot_SEM_lasso_mods4_lin = function(lass_mod, group_colors, sem_tab, sel_models,
#                                sort_var = "aic", color_var = "R2", exclude_mods = '', 
#                               rel_widths =c(1, 1, 1, 1), ncol = 4, nrow =1){
#   
#    coef_plot   <-   plot_lasso_coefs(lass_mod, group_colors)
#    models_plot <-   plot_SEM_lasso_mods(sem_tab) + geom_smooth(method = "lm", colour="darkgrey",
#                         linetype = "dashed", size = 1, se = FALSE)    
#    
#    selected_plot <- plot_sem_lasso_selected(sem_tab, sel_models) + theme(legend.position = "none")
#    
#    coef_tab    <-   coef_sort_sem_form_tab(lass_mod, sel_models, sep_table = 'TRUE', fxn_operator = '~ ')  
#    leaps_plot  <-   leaps_plot_sem_models(coef_tab, sort_var, color_var, exclude_mods) +
#                        scale_x_discrete(position = "top") +
#                        theme(axis.text.x = element_text(angle = 45, hjust = 0), axis.title.x = element_blank())
#    
#    quad_plot <- plot_grid(coef_plot, models_plot, selected_plot, leaps_plot, 
#                           ncol = ncol, nrow = nrow, labels = "AUTO", align = "h", rel_widths = rel_widths)
#   
#    return(quad_plot)    
#}

# demonstrate function
# lasso_sel_sumplot <- plot_SEM_lasso_mods4_lin(CH4_lass0, guild_colors, semF, sel_models, exclude_mods = '')
# lasso_sel_sumplot <- plot_SEM_lasso_mods4_lin(CH4_lass0, guild_colors, semF, sel_models, ncol =2, nrow = 2)

# options(repr.plot.width=12, repr.plot.height= 3)
# lasso_sel_sumplot


########################################################################################
### b) list of lasso & SEM model plots

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
# lasso_sel_sumplot <- plot_lasso_semfilt_list(CH4_lass0, guild_colors, semF, filt_models)
# lasso_sel_sumplot#[1]; #lasso_sel_sumplot$coef_plot
# plot_grid(lasso_sel_sumplot[[1]], lasso_sel_sumplot[[2]], lasso_sel_sumplot[[3]], nrow =1, ncol = 3)
#                       ncol = ncol, nrow = nrow, labels = "AUTO", align = "h", rel_widths = rel_widths)

########################################################################################
### c) list of heatmap plots

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
# multi_heat <- combine_model_heatmaps_list(CH4_lass0, filt_models, vip_matrix, alt_metric = "resid", legend_row = 2, legend_col = 1)
# multi_heat
# plot_grid(multi_heat[[1]], multi_heat[[2]], multi_heat[[3]], nrow =1, ncol = 3, align = "h")


########################################################################################
### d) plot lasso, sem, model inspection function (wrapper, user fxn)

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

