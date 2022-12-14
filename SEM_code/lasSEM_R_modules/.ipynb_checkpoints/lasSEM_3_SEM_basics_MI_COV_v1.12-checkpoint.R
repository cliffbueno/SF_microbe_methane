########################################################################################
### 3) SEM functions for extracting fits, MI covariates

# SEM - Lasso semi-auto for SF OTU guilds v1.12
# runs under r > 3.6 (glmnet for LASSO)

# ASSUME have loaded modules 0-2
# TODO-- too much example data cut to functions only (but retain examples too, elsewhere)

# good information about acccessing fit parameters can be found here: http://www.understandingdata.net/2017/03/22/cfa-in-lavaan/ (http://www.understandingdata.net/2017/03/22/cfa-in-lavaan/)
# How to get features from Lavaan output, used in functions -- Deactivated learning example
# parameterEstimates(meth_mod.sem, standardized=TRUE)  # gets P value for each Regression in summary
# inspect(meth_mod.sem)                                # gets various SEM model params 
# inspect(meth_mod.sem, 'r2')                          # gets R2 for each fitted var
# fitMeasures(meth_mod.sem)                            # gets all model fit measures
# fitMeasures(meth_mod.sem, c("aic", "bic"))           # gets selected model fit measures

########################################################################################
## 0) test data - uses Guild_CH4_d, loaded in 0-2, formula defs below

# "Simple model" - example for 
meth_mod <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h 
             
             CH4_H2 ~ SRB_syn + CO2_mg_m2_h
             MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac
             MOB_I ~ CH4_H2 +CH4_ac'

# Calculate SEM, basic method
meth_mod.sem <- sem(meth_mod, data=Guild_CH4_d, fixed.x=FALSE, estimator ="mlm")    # Generates observed var > 1000 x warning, divide otu vas try again
sem_fit_params <- c( "pvalue", "chisq", "df", "npar","aic", "bic", "gfi", "cfi", "rni", "rmsea", "srmr")


########################################################################################
### A) check model results
# Purpose: new functions to extract flat table features from a Lavaan SEM model summary 

## i)   get_SEM_nonsig_vars(meth_mod.sem, p_cut =0.05) 
## ii)  get_SEM_R2s(meth_mod.sem) 
## iii) sem_fit_sum(meth_mod.sem, sem_fit_params)


########################################
## i) check for SEM variable significance

get_SEM_nonsig_vars = function(sem_mod, p_cut) {

    # get significant only 
    sem_params <- parameterEstimates(sem_mod)                          # Get estimates
    sem_params_bad <- sem_params[sem_params$pvalue > p_cut,]           # cut out significant

    # clean up
    keep <- c('lhs', 'op','rhs','pvalue')                              # params to keep
    sem_params_bad <- sem_params_bad[keep]                             # keep
    sem_params_bad$pvalue <- round(sem_params_bad$pvalue, 2)           #  round p value
    sem_params_bad <- sem_params_bad[complete.cases(sem_params_bad),]  # drop all NA values (comes up in latent mods)
    sem_params_bad <- sem_params_bad[!sem_params_bad$op == "~~",]       # drop covariate fails 
    
    return(sem_params_bad)
    
}

# Demonstrate function
get_SEM_nonsig_vars(meth_mod.sem, p_cut =0.05)

########################################
## ii) get R2 for each variable in model

get_SEM_R2s = function(sem_mod){
    
    sem_name <- paste(as.character(deparse(substitute(sem_mod))))
    model <- data.frame(model = sem_name)
    
    R2s <- inspect(sem_mod, 'r2')
    R2s <- data.frame(round(R2s, 3))
    out <- t(R2s)
    
    out <- data.frame(model, out)
    row.names(out) <-NULL
    return(out)
}

# demonstrate function
get_SEM_R2s(meth_mod.sem)

########################################
## iii) function for custom SEM fit parameters summary

sem_fit_sum = function(sem_mod, sem_fit_params) {
    
    sem_name <- paste(as.character(deparse(substitute(sem_mod))))
    model <- data.frame(model = sem_name)
    
    fits <- data.frame(fitMeasures(sem_mod, sem_fit_params))
    fits <- round(fits, digits = 3)
    out <- t(fits)
    
    out <- data.frame(model, out)
    row.names(out) <-NULL
   # names(out)[1] <- sem_name
    return(out)      
}

# Test function
sem_fit_sum(meth_mod.sem, sem_fit_params)

################################################################################
### B) compare models

## i) compare_sem_fits(models2, sem_fit_params) 
## ii) compare_sem_R2s(models2) 
## iii) compare_sem_results(models2, sem_fit_params)

########################################
## i) compare fits of more than one model
# models = c(meth_mod.sem, meth_mod.sem)
meth_mod.sem2 <- meth_mod.sem
models2 = c("meth_mod.sem","meth_mod.sem2")

compare_sem_fits = function(models, sem_fit_params){
    
    out_list <- {}
    
    for (i in seq_along(models)) {
        curr_mod <- get(models[[i]])
        out_list[[i]] <- sem_fit_sum(curr_mod, sem_fit_params)
        }

    out <- do.call("rbind", out_list)
    out$model <- models
    return(out)
    
}

# Demonstrate function
compare_sem_fits(models2, sem_fit_params)

########################################
## ii) compare R2 fits for more than one model

compare_sem_R2s = function(models){
    
    out_list <- {}
    
    for (i in seq_along(models)) {

        curr_mod <- get(models[[i]])
        R2c <- get_SEM_R2s(curr_mod)
        R2c$model <- models[i]
        out_list[[i]] <- R2c
        #out_list[[i]] <- get_SEM_R2s(curr_mod)
        #names(out_list[i]) <- models[i]
        }
    
    out <- Reduce(function(...) merge(..., all=T), out_list)
    
    #out <- Reduce(merge, out_list)
    #out <- do.call("rbind", out_list)
    #out$model <- models
    #out_list
    out
}

# Demonstrate function
#models2 <- c('ch4_mod0_f','ch4_mod0a_f', 'ch4_mod0b_f')
compare_sem_R2s(models2)  # note this doesn't look right with two of the same models names

########################################
# iii) combine fit and R2 summaries

compare_sem_results = function(models, sem_fit_params){
    
        fits <- compare_sem_fits(models, sem_fit_params)
        R2s <- compare_sem_R2s(models)
        out <- data.frame(merge(fits, R2s, by="model"))
        return(out)

}

# models2
# meth_mod.sem

# Demonstrate function
compare_sem_results(models2, sem_fit_params)

# DEACTIVATE edge testing, FIXED prior issues way downstream -- non-matching outputs with running auto composites
# man_mod_7 <- c('sem.ch4L0C.107','sem.ch4L0C.113','sem.ch4L0C.116','sem.ch4L0C.118',
#               'sem.ch4L0C.121','sem.ch4L0C.69','sem.ch4L0C.79')
# auto_mod_7m <- c('sem.107','sem.113','sem.116','sem.118','sem.121','sem.69','sem.79')
# auto_mod_7 <- c('sem.8','sem.10','sem.16','sem.34','sem.107','sem.113','sem.116','sem.118','sem.121','sem.69','sem.75','sem.79')

# compare_sem_results(man_mod_7, sem_fit_params)
# compare_sem_results(auto_mod_7m, sem_fit_params)
# compare_sem_results(auto_mod_7, sem_fit_params)

################################################################################
### C) extract modification index output

sem_mi_table <- function(sem_mod, mi_lower = 4, mi_sort = "T", head = "T"){
    
    mi_tab <- data.frame(modindices(sem_mod))  # get mod indices table
    mi_tab <- mi_tab[mi_tab$mi > mi_lower,]         # impose cutoff on lower MI 
    
    # sort by mi?  or defaults to covariances first
    if (mi_sort == "T") {out_tab <- mi_tab[order(-mi_tab$mi),]}
        else {out_tab <- mi_tab}
    
    # crop to top 5 most important effects?
    if (head == "T") {out <- head(out_tab)}
        else {out <- out_tab}
    
    return(out)
        
}

# demonstrate function
sem_mi_table(meth_mod.sem)
#sem_mi_table(meth_mod.sem, mi_lower = 4, mi_sort = "F", head = "T")


################################################################################
### D) access model residual covariances 

## i) function to get /trim residual covariate data frame

get_sem_covariances = function(sem_mod, n_res =5){  # output_factors = F - consider output as list not table?
    
    # Get covariance matrix, melt
    resid <- residuals(sem_mod)$cov
    resid[lower.tri(resid, diag = T)] <-0         # get upper triangle only, lower as 0
    resid_M <- melt(resid)                        # melt DF

    # Clean, sort table
    resid_M <- resid_M[!resid_M$value == 0,]      # Drop 0/lower
    resid_M$abs <- abs(resid_M$value)             # Get abs for sorting effect size 
    resid_M <- resid_M[order(-resid_M$abs),]      # sort by abs
    resid_M <- resid_M[,1:3]                      # drop abs col 
    resid_M <- resid_M[1:n_res,]                  # get first n results
    resid_M
              
} 
    

# Demonstrate function
get_sem_covariances(meth_mod.sem, n_res =5)


## ii) function to plot residual covariances in model
plot_matrix <- function(matrix_toplot){
corrplot(matrix_toplot, is.corr = FALSE, 
               type = 'lower', 
               order = "original", 
               tl.col='black', tl.cex=.75)
}


options(repr.plot.width=4, repr.plot.height=4) 
plot_matrix(residuals(meth_mod.sem)$cov)

# residuals(meth_mod.sem)$cov
# residuals(meth_mod.sem, type = "cor")$cor   # "keep an eye out for residual corr > 0.1" from http://www.understandingdata.net/2017/03/22/cfa-in-lavaan/

################################################################################
## E) Suggest model updates with external covariates (based on MI and model covariates) 

# prepare data matrix from input
data <- Guild_CH4_d[,-1]
data <- data[!is.na(data$CH4_ug_m2_h),]
data <-data.matrix(data)


## i) function to get external covariates (from internal pair)

get_mi_covars_pair = function(data, r_cut = 0.5, keep_pair){
    
    corr_m <- cor(data)                                             # get corr
    corr_m[abs(corr_m) < r_cut] <- NA                               # get abs value over threshold only, or NA

    select_cols <- data.frame(corr_m[,keep_pair])                   # Get only cols from MI inputs
    # select_cols <- data.frame(corr_m[,keep_vars])                 -- too many vars at once not helping
    common_data <- t(select_cols[complete.cases(select_cols),])        # get only complete cases -- no NA
    common_data <- round(common_data,2)
    #select_cols
    common_data
}

# Note behavior change with many variables, only common feats kept
# keep_pair = c("CH4_H2", "SRB_syn")
keep_pair = c("CH4_H2", "AOA")
# keep_pair = cov_uniq[1:4]   

# Demonstrate function
get_mi_covars_pair(data, r_cut = 0.6, keep_pair)


### ii) function to get mod index covariates for many

get_mi_covars_many = function(data, r_cut = 0.5, keep_vars){
    
    corr_m <- cor(data)                                             # get corr
    corr_m[abs(corr_m) < r_cut] <- NA                               # get abs value over threshold only, or NA

    # select_cols <- data.frame(corr_m[,keep_pair])                 # Get only cols from MI inputs
    select_cols <- data.frame(corr_m[,keep_vars])                   #  too many vars at once not helping
    # common_data <- select_cols[complete.cases(select_cols),]      # get only complete cases -- no NA
    common_data <- t(select_cols[rowSums(is.na(select_cols)) != ncol(select_cols), ])
    common_data <- round(common_data,2)
    common_data
}

keep_vars <- c("CH4_H2", "MOB_I", "AOA", "NOB", "AO_NOB")
#keep_vars <-cov_uniq

# Demonstrate function
get_mi_covars_many(data, r_cut = 0.6, keep_vars)

################################################################################
### F) SEM formula builder

########################################
## 0) Initial variable and list definition (example)
# Define simple model" - 
ch4_mod0 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h 
             
             CH4_H2 ~ SRB_syn + CO2_mg_m2_h
             MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac
             MOB_I ~ CH4_H2 +CH4_ac'
# try adding SEM response vars with string defs
CH4 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h'
CH4_H2 <- 'CH4_H2~ SRB_syn + CO2_mg_m2_h'
MOB_IIa <- 'MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac'
MOB_I <- 'MOB_I ~ CH4_H2 +CH4_ac'

# Get vector of response vars as strings
Other <- ""
response_vars <- c("CH4", "CH4_H2", "MOB_I","MOB_IIa", "Other")

########################################
### i) make SEM model building data frame 

#define_new_sem_table <- function(mod_name, response_vars, output_prefix = "sem."){

    # Build data frame of sub-formulas from names
#    feat_list <- {}
    
#    for (i in seq_along(response_vars)) {
#
 #       feat_list[[i]] <- get(response_vars[i])

  #  }
    
   # feats <- data.frame(feat_list)
    #row.names(feats) <- response_vars
    #feats <-data.frame(t(feats))
              
    # Get model formula from feats
    #formula <- do.call(paste, c(feats, sep = '\n'))
    #feats$formula <- formula
     
    # add model info
    #model_out <- paste0(output_prefix, mod_name)
    #mod_info <- data.frame(model_out = model_out, model = mod_name, base = mod_name, notes = "base")
    
    #out <- data.frame(mod_info, feats)
    #row.names(out) <- NULL
    #return(out)
        
#}
                  
define_new_sem_table <- function(mod_name, response_vars, output_prefix = "sem."){

    # Build data frame of sub-formulas from names (not working so manually do). Note need version for 5, 7 and 9 vars
    #feat_list <- {}
    #for (i in length(response_vars)) {
#
 #       feat_list[[i]] <- get(response_vars[i])
#
 #   }
    feat_list <- list()
    feat_list[[1]] <- get(response_vars[1])
    feat_list[[2]] <- get(response_vars[2])
    feat_list[[3]] <- get(response_vars[3])
    feat_list[[4]] <- get(response_vars[4])
    feat_list[[5]] <- get(response_vars[5])
    if (length(response_vars) == 7) {
        feat_list[[6]] <- get(response_vars[6])
        feat_list[[7]] <- get(response_vars[7])
}
    if (length(response_vars) == 9) {
        feat_list[[6]] <- get(response_vars[6])
        feat_list[[7]] <- get(response_vars[7])
        feat_list[[8]] <- get(response_vars[8])
        feat_list[[9]] <- get(response_vars[9])
}
    
    feats <- data.frame(feat_list)
    names(feats) <- response_vars
#    row.names(feats) <- response_vars
#    feats <-data.frame(t(feats))
              
    # Get model formula from feats
    formula <- do.call(paste, c(feats, sep = '\n'))
    feats$formula <- formula
     
    # add model info
    model_out <- paste0(output_prefix, mod_name)
    mod_info <- data.frame(model_out = model_out, model = mod_name, base = mod_name, notes = "base")
    
    out <- data.frame(mod_info, feats)
    row.names(out) <- NULL
    return(out)
        
}

# Demonstrate function
model_test <- define_new_sem_table(mod_name = "ch4_mod0", response_vars)
model_test

# Inspect formula
# cat(model_test$formula)    #model_test$formula

# check formula extraction
formula_test <- model_test$formula
# cat(formula_test)

# Confirm this works in SEM model
ch4_mod0_build_f <- sem(formula_test, data=Guild_CH4_d, fixed.x=FALSE, estimator ="mlm")    # Generates observed var > 1000 x warning, divide otu vas try again
sem_fit_sum(ch4_mod0_build_f, sem_fit_params)

########################################
### ii) update / edit models in table

update_sem_model_table = function(mod_table, base_model, new_mod_name, notes, formula_edits, 
                                  output_prefix = "sem.", 
                                  non_feat_cols = c("model_out","model", "base", "notes", "formula")){
    
     # copy base model & update descriptions
     new_model <- mod_table[mod_table$model == base_model,]
     new_model$base <- base_model
     new_model$model <- new_mod_name
     new_model$model_out <- paste0(output_prefix, new_mod_name) 
     new_model$notes <- notes
    
     # edit formula components
     form_feat <- formula_edits[1]
     form_edit <- formula_edits[2]
    
     new_model[form_feat] <- form_edit
    
     # rebuild full formula
     cols <- names(new_model)
     feat_cols <- cols[!cols %in% non_feat_cols]
     feats <- data.frame(new_model[feat_cols])
                         
     formula <- do.call(paste, c(feats, sep = '\n'))
     new_model$formula <- formula
    
     # append to original table 
     mod_table_updated <- rbind(mod_table, new_model)
                 
     mod_table_updated
     #new_model
}

# Test function
model_test2 <- update_sem_model_table(model_test, base_model = "ch4_mod0", new_mod_name = "ch4_mod0a.1", 
                      notes = "drop CH4~MOB_I",
                      formula_edits = c("CH4", "CH4_ug_m2_h ~ CH4_ac + CH4_H2 + MOB_IIa + CO2_mg_m2_h"))
                      # output_prefix = "sem."#non_feat_cols = c("model", "base", "notes", "formula"))

model_test2
# Test while dropping entire internal feature
model_test3 <- update_sem_model_table(model_test2, base_model = "ch4_mod0a.1", new_mod_name = "ch4_mod0a.2", 
                      notes = "drop MOB_I ~ .",
                      formula_edits = c("MOB_I", ""))

model_test3
dim(model_test2)
#nrow(sem_model_table)
# sem_model_table$formula

########################################
## iii) run SEM models & assign from model table

sem_table_run_assign = function(sem_mod_table, data, estimator = 'mlm', prefix = "sem."){

        iters <- nrow(sem_mod_table)                                                  # seq_along misbehaving, use iters for loop
    
        for(i in seq(1:iters)) {
         
            # run SEM model            
            curr_model <- sem_mod_table$formula[i]                                    # get ith sem formula
            sem <- sem(curr_model, data=data, fixed.x=FALSE, estimator = estimator, orthogonal =TRUE)   # run SEM; fixed.x not working in fxn
            
            # assign mod_out names to models (GLOBAL ENV)
            sem_name <- paste(sem_mod_table$model_out[i])                             # get model name from model_out
            assign(sem_name, sem, envir=.GlobalEnv)                                   # assign to global env 
        }
        
        # Gather models run and print list to console
        models_run <- paste(sem_mod_table$model_out, collapse =', ')
        print(paste0("ran SEM to create models: ", models_run))
        
        models_run_out_list <- c(paste(sem_mod_table$model_out))
        models_run_out_list
    
    }
    
# Note, could have run through output as list, but not necessary as only assigning vars here  # mods <- {}

# test run and assign models from model table
sem_multi_test <- sem_table_run_assign(model_test2, data=Guild_CH4_d, estimator ="mlm")
sem_multi_test
#cat(x)
#data.frame(x)

# test model summary
compare_sem_results(sem_multi_test, sem_fit_params)

########################################
## iv) new run and summary function from model table

run_compare_sem_models = function(mod_table, data, estimator = "mlm", 
                                  keep_mod_descr = c("model_out","model", "base", "notes"),
                                  sem_fit_params = c("pvalue", "chisq", "df", "npar","aic", "bic", 
                                                      "gfi", "cfi", "rni", "rmsea", "srmr")){
    
    # note mod_table requires columns: model, model_out, formula cols 
    
    # get model table features to include in summary
    mod_descr <- data.frame(mod_table[keep_mod_descr])
    mod_descr
    
    # run assign SEM models in table 
    sem_out_list <- sem_table_run_assign(mod_table, data, estimator)
    #sem_out_list

    # Get SEM results summary as previously
    fit_summary <- data.frame(compare_sem_results(sem_out_list, sem_fit_params))
    
    
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
out <- run_compare_sem_models(model_test2, data, estimator = "mlm")
out


########################################
### v) Report SEM models & output summary table (specific to this data, non-latent)

clean_SEM_models_for_output = function(models, models_out){
    
    ### Clean up models ####  
    models <- models[, -ncol(models)]                           # drop formula column  
    models <- models[, -1]                                      # drop model out columns

    # Hide predictors in formulas
    models$CH4 <- gsub("  ", " ", models$CH4)                   # Note there are double spaces in some, likely result of replace ~.
    models$CH4 <- gsub("CH4_ug_m2_h ~ ", "", models$CH4)        # hide predictor

    models$CH4_H2 <- gsub("  ", " ", models$CH4_H2)             # remove double spaces
    models$CH4_H2 <- gsub("CH4_H2~ ", "", models$CH4_H2)        # hide predictor     

    models$MOB_I <- gsub("MOB_I ~ ", "", models$MOB_I)        

    models$MOB_IIa <- gsub("  ", " ", models$MOB_IIa)           # remove double spaces
    models$MOB_IIa <- gsub("MOB_IIa ~ ", " ", models$MOB_IIa)   # hide predictor          
    # return(models)
    
    # Merge with results
    results_tab <- merge(models_out, models, by = c('model', 'base', 'notes'))
    return(results_tab)
}

# Demonstrate function
# Exper1_results <- clean_SEM_models_for_output(Ex1_sem_model_table, Ex1_sem_out)
