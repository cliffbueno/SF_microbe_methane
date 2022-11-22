# Lasso semi-auto for SF OTU guilds v1.12
# runs under r > 3.6 (glmnet for LASSO)

# Assume have loaded modules 0-9
# likely pollulted with examples, possibly cruft

########################################################################################
### 10) Wrappers for combined lasso -> SEM models, plots, workflows


# likely these functions should be returning lists, including lasso_mod.
# NOTE changing lasso mod may hit downstream...


########################################################################################
### a) initial Lasso and coefficient plot

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

########################################################################################
### b) hard filter to force only selected variables into GLM / SEM

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

########################################################################################
### c) lasso to SEM wrapper w/hard filter

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

########################################################################################
### d) lasSEM object to SEM + filter + model selection plots workflow

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
