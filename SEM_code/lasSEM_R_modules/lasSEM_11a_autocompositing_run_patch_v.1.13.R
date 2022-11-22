####################################################################################################

# lasSEM_11a_autocompositing_run_patch_v.1.13

# to load inline after lasSEM_11 for testing
# then replace relevant functions

# patch given bugs during runs in lasSEM branching model runs (all sites)
# problems were wrapper script failures.
# try catch was built into some functions in 11), [re-written from 3)], not complete error coverage

# here issue presenting is some SEM models do not converge, but don't flag primary error
# thus, they are still assigned to output lists & global enviroment.
# then model summary extractors look for fit parameters in models with no fit data.


# Chain of functions causing error: (top down)
#   run_compare_sem_models -->
#		compare_sem_results -->
#			compare_sem_fits -->
#				sem_fit_sum
# 


####################################################################
## Patched functions 
#   note all functions renamed (number) EXCEPT top level wrapper same
#   run_compare_sem_models #  so that this can be called as previously

#  TODO revise function source code in 11) to fully resolve.


####################################################################
# a) compare_sem_fits -- added try catch here, left mess

compare_sem_fits3 = function(models, sem_fit_params){
    
    out_list <- {}
    
    
    for (i in seq_along(models)) {
        skip <- FALSE     
        curr_mod <- get(models[[i]])
        fits <- tryCatch(sem_fit_sum(curr_mod, sem_fit_params),
            error = function(e){skip <<- TRUE})                   # if SEM fails, skip = T  
            if(skip==TRUE){next} 
        
        fits$model <- models[i]
        out_list[[i]] <- fits
        #names(out_list[i]) <- models[i]
        }

    if (is.null(out_list)==TRUE) {
        print('SEM models failed!')
        
        } else {
        out <- do.call("rbind", out_list)
        #out$model <- models
        return(out)
    }
        
    #out <- do.call("rbind", out_list)
    #out$model <- models
    #return(out_list)
    
}

# demonstrate function (on output from ...)
# fits_test <- compare_sem_fits3(sem_out_list, sem_fit_params)
# head(fits_test[order(-fits_test$pvalue),])


####################################################################
# b) combine fit and R2 summaries -- NO REAL CHANGE, just embedding above change

compare_sem_results2 = function(models, sem_fit_params){
    
        fits <- compare_sem_fits3(models, sem_fit_params)        # v3 has try catch
    
        models_passed <- fits$model                              # include only models passing catch 
        R2s <- compare_sem_R2s(models_passed)                    # compare R2
        out <- data.frame(merge(fits, R2s, by="model"))
        return(out)

}

# fit_summary <- data.frame(compare_sem_results2(sem_out_list, sem_fit_params))
# head(fit_summary[order(-fit_summary$pvalue),])


####################################################################
# c)  ### 11) II c) run and compare SEM models in table (Updated)
# -- NO REAL CHANGE, just embedding above changes, DON'T need after remove renumbered fxn

#      - implements parallel sem_table_run_assign_P, passes n_cores to options

run_compare_sem_models2 = function(mod_table, data, estimator = "mlm", ncores = 20,
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
    fit_summary <- data.frame(compare_sem_results2(sem_out_list, sem_fit_params))
    
    
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

# out <- run_compare_sem_models2(sem_formula_tree, Guild_CH4_aS, estimator = "mlm")
# head(out[order(-out$pvalue),])

####################################################################
# d) ### 12 c) run sem model summary
# -- NO REAL CHANGE, just embedding above changes, DON'T need after remove renumbered fxn
#   - takes an SEM formula table and runs it, creates summary table with filtering & output options (& parallel cores

run_sem_model_summary = function(sem_formula_table, data, sort_by = 'pvalue', sort_cutoff = 0, 
                                 estimator = 'mlm', output = 'summary_table', ncores = 20){
    
    # Run models for results
    sem_results <- suppressWarnings(run_compare_sem_models2(sem_formula_table, data, ncores = ncores,
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

# Demonstrate function
# run_sem_model_summary(sem_models.2, Guild_CH4_0)    #  test fxn
# run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.001)   #  -- ex won't work here