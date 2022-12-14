# SEM - Lasso semi-auto for SF OTU guilds v1.12
#-  runs under r > 3.6 (glmnet for LASSO)

# Assume have loaded modules 0-5
# likely pollulted with examples, possibly cruft

########################################################################################
### 6) Plot LASSO-SEM model screening functions


# HELPER fxn -- get ggplot color wheel hue function -- Maybe use later?
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = (n + 1))
  hcl(h = hues, l = 65, c = 100)[1:n]
}

########################################################################################
### a) plot SEM LASSO models
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


########################################################################################
### b) Plot selected SEM models vs. all

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

# Demonstrate function
get_models <- c(34, 35, 67, 69, 82, 107, 113, 116) 
sel_models <- semF[semF$model %in% get_models,]

options(repr.plot.width = 5, repr.plot.height = 4)
plot_sem_lasso_selected(semF, sel_models)
#save_plot(filename, plot, width, height)

########################################################################################
### c) sort variable in SEM model table -- MOVED BELOW -- DEPRECATED??

###  XXX DELETE THIS XXX ###

# Function for splitting sem formulas & resorting form vars (by LASSO abs coef/ VIP)
# coef_sort_sem_form_tab = function(lasso_model, sem_tab, fxn_operator = '~ ', sep_table = 'FALSE'){
#    
#     ### Prepare empty var inclusion table ###
#
#    # sort lasso coefs by abs. vals
#    coefs <- glmnet_coefs(lasso_model)                                                        # FXN get coefs  
#    coefs$abs <- abs(coefs$coef)                                                              # get coef abs
#    coefs <- coefs[order(-coefs$abs),]                                                        # sort

#    # get vars & model matrix data
#    ranked_vars <- coefs$var                                                                  # get ranked vars 
#    n_coef <- length(ranked_vars)                                                             # n vars   (cols)
#    n_mods <- dim(sem_tab)[1]                                                                 # n models (rows)             

#    # make vars, model matrix DF                                                              # w/ empty model components
#    coef_m <- data.frame(matrix(ncol=n_coef,nrow=n_mods, dimnames=list(NULL, ranked_vars)))   # make empty df
#    sem_mods_vars <- cbind(sem_tab, coef_m)                                                   # join sem models
     
#    ### Fill model matrix w. true / false ###
#    for (i in ranked_vars)                                                                   # for each var
#          sem_mods_vars  <- model_feat_ext(sem_mods_vars, i)                                  # use FXN for fill
 
#    ### Rebuild sorted model formulas 
#    # Get y var to redo formulas, (assume y same in all)
#    form1 <- paste(sem_mods_vars$form[1])                                                    # first formula 
#    #fxn_operator <- "~ "
#    y_var <- unlist(strsplit(form1, fxn_operator))[1]                                        # get y var
#    y_var <- paste0(y_var, fxn_operator)                                                     # add op back
#    # y_var

#    # Get x var formula componets from TF data
#    sort_Xform <- apply(sem_mods_vars, 1, function(r) paste0(names(r)[r == TRUE], collapse = " + "))
#    sem_mods_vars$form <- paste0(y_var, sort_Xform)                                         # add y, replace form                   
                                 
#    ### Select output mode ###
#    if (sep_table == 'FALSE'){
#        out_tab <- sem_mods_vars[names(sem_tab)]} else {                                    # revert to input cols
#        out_tab <- sem_mods_vars}                                                           # keep split cols (for plots)    
#    out_tab
    
#}
# Demonstrate function:
# coef_tab <-coef_sort_sem_form_tab(CH4_lass0, sel_models, sep_table = 'FALSE', fxn_operator = '~ ')  # Get sorted table
#coef_tab <- coef_sort_sem_form_tab(CH4_lass0, sel_models, sep_table = 'TRUE', fxn_operator = '~ ')   # Get full table (for plots)
# coef_tab# 


########################################################################################
### KEEP HELPER FUNCTION !!!
# HELPER function: single col TRUE / FALSE replace with value

TF_replace = function(coef_tab, var, value){
    var <- ifelse(coef_tab[,var] == TRUE, coef_tab[,value], 0)
    return(var)
}

########################################################################################
### d) LEAPS style SEM model plot by variable used -- DEPRECATED  

###  XXX DELETE THIS XXX ###

#leaps_plot_sem_models = function(coef_tab, sort_var = "aic", color_var = "R2", exclude_mods = ""){
    
#    ### replace T/F var columns with data (color_var) ###
    
#    # Get variable sets for apply:
#    coef_cols <- names(coef_tab)                             # all columns
#    log_cols_test <- unlist(lapply(coef_tab, is.logical))    # which are T/F
#   log_cols <- coef_cols[log_cols_test]                     # get T/F columns
#    var_cols <- log_cols
#    other_cols <- coef_cols[!coef_cols %in% log_cols]        # non-T/F columns 

#    # make new aic2x var
#    # coef_tab$aic2x <- rescale(coef_tab$aic)*rescale(coef_tab$aic_glm)
#    # coef_tab$aic2x <- coef_tab$aic*coef_tab$aic_glm

#    # Apply TF replace function to logical, passing others 
#    var_col_list = {}                                        # empty list
#    for (i in log_cols)                                      # for each T/F col 
#        var_col_list[[i]] <- TF_replace(coef_tab, i, color_var)  # insert val or 0

#    coef_tab[log_cols] <- data.frame(var_col_list)     # -- USE THIS ONE  # replace log cols in orig data
    
#    ### get ggplot data ###
    
#    # exclude models 
#    if (exclude_mods == ""){
#        coef_tab <- coef_tab} else {
#        coef_tab <- coef_tab[!coef_tab$model %in% exclude_mods,]
#    }
           
#    # reorder model factor for plot
#    coef_tab$model <- as.factor(as.character(coef_tab$model))
#    coef_tab$model <- reorder(coef_tab$model, coef_tab[,sort_var])

#    # Get vars & data
#    gg_vars <- c('model', log_cols)
#    gg_data <- coef_tab[gg_vars]
#    gg_data_m <- melt(gg_data)
    
#    # clean data
#    vars <- gg_data_m$variable
#    #gg_data_m$variable <- strtrim(vars, 9)          # trim var names, fails due to factor ordering issue, upstream issues before factor too
#    gg_data_m[gg_data_m == 0] <-NA

#    # plot data 
#    plot <- ggplot(data = gg_data_m, aes(y=model, x=variable, fill=value)) + 
#          geom_tile(color="white") +  scale_fill_gradient(na.value = "white", high = "#132B43", low = "#56B1F7",) + 
#          theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + 
#          labs(fill = color_var)
#    plot
# }
 
# coef_tab <-coef_tab[-1,]
# options(repr.plot.width = 4, repr.plot.height = 4)

# leaps_plot_sem_models(coef_tab, sort_var = "aic", color_var = "R2", exclude_mods = '')
# leaps_plot_sem_models(coef_tab, sort_var = "aic", color_var = "R2", exclude_mods = c('34', '113', '121')) # 


