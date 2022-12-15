########################################################################################
### 5) ML feature selection functions for extracting fits, MI covariates

# SEM - Lasso semi-auto for SF OTU guilds v1.12
# runs under r > 3.6 (glmnet for LASSO)

# ASSUME have loaded modules 0-4
# TODO-- too much example data cut to functions only (but retain examples too, elsewhere)

#### Rationale
# select with LASSO CV (glmmnet); BUT not all preds individually significant
# screen significant with all subsets GLM (modified, custom for all vars sig); BUT diff penalty than SEM
# get only SEM significant models with SEM screening

# Alternately, had tried regSEM but failed on larger models (incl. stability selection)
# on smaller models passing stability selection to adaptive lasso worked well in regSEM
# note stability selection for linear models available in pacakges 'stabs' and 'stabsel' 


##### Functions (needs update!)
# - glmnet_lassoCV
# - glmnet_coefs
# - all_subsets_whh
# - sem_filter
# - plot_sem_filter

########################################################################################
## 0) prepare test data

# SCALE DATA, careful will scale y
# Guild_CH4_dS <- scale(Guild_CH4_d[,-1])

# Remove sample names
Guild_CH4_d0 <- Guild_CH4_d
row.names(Guild_CH4_d0) <- Guild_CH4_d0$Sample
Guild_CH4_d0 <- Guild_CH4_d0[,-1]

# names(Guild_CH4_d)

# prepare data subsets

# get vars
xvars0 <- names(Guild_CH4_d0)

# vars to exclude (maybe)
drop_gas <- c('CH4_ug_m2_h', 'CH4_logn1','CH4_CO2', 'CO2_soilC_mg_g_d') # keep this: ,'CO2_mg_m2_h',
drop_ratios <- c('CN','CP','NP','NP_ext')
drop_pw <- c('SO4_pw', 'Fe_pw','Fe','DOC_mg_L')
drop_guilds <- c('SOxB','MeOB','FeOB','AO_NOB','NOB_AO','mcr_pmo', 'pmo_mcr','Anamx')
drop_taxa <- c('Actino','Chlorf','Firmic')

# Get combined drop lists
drop_vars1 <- c(drop_gas, drop_pw, drop_guilds, drop_taxa, drop_ratios)

# Select guild & chem data:                # Remove drop vars from list
xvars1 <- xvars0[!xvars0 %in% drop_vars1]

# alternately, may want only guild or chem data lists, including ALL

# make glmnet data
data <- Guild_CH4_d0
y_var <- 'CH4_logn1'
x_vars <- xvars1

# SCALE ALL, needed for some functions (SEM screening)
Guild_CH4_dS <- scale(Guild_CH4_d[,-1])              #  careful will scale y 
Guild_CH4_dS_df <- data.frame(Guild_CH4_dS)          #  make df
Guild_CH4_dS_df$CH4_logn1 <- Guild_CH4_d0$CH4_logn1  # fix CH4 scale to original

# scale data, get separate x & y -- not needed, moved into FUNC
# Get SCALED X vars
# CH4_pred1 <- Guild_CH4_d0[xvars1]
# CH4_predS1 <- scale(CH4_pred1)
# names(CH4_pred1)# ; head(CH4_predS1)

# get Y var (CH4_flux)
# CH4_flux <- Guild_CH4_d0$CH4_logn1


########################################################################################
### a) glmnet_lassoCV function      # Revised function with force in
# note data split & scale inside function, might want to move it out

glmnet_lassoCV = function(y_var, x_vars, data, force_in = '', max_pred = 13){  #, max_pred = 
    
    # Get x & y data
    y <- data[,y_var]
    x <- data.matrix(data[x_vars])
    xS <- scale(x)                                       # lasso scales X by default?  results differ,
    
    # Get force in vars
    force_in_test <- x_vars %in% force_in
    force_in_vec <- ifelse(force_in_test == FALSE, 1, 0)

    # cross val for best lambda
    set.seed(123) 
    cv <- cv.glmnet(xS, y, alpha = 1, penalty.factor = force_in_vec, dfmax = max_pred)
    print(paste0("best_lambda: ", round(cv$lambda.min, 4)))
    
    # lasso w best lambda
    pfit <- glmnet(xS, y, alpha = 1, lambda = cv$lambda.min, dfmax = max_pred,
                  penalty.factor = force_in_vec) # , exclude = force_out_vec)
    
    # get predictions / fit scores (print)
    predictions <- pfit %>% predict(xS)
    fit <- data.frame(R2(predictions, y), RMSE(predictions, y))
    names(fit) <- c('R2','RMSE')
    fit <-round(fit, 3)
    print(fit)
    
    return(pfit)
    }

# demonstrate function
CH4_lass0 <- glmnet_lassoCV(y_var, x_vars, Guild_CH4_d0)
#CH4_lass0.2 <- glmnet_lassoCV3(y_var, x_vars, Guild_CH4_d0, force_in = c("C","Salinity.x","SRB"))

# function testing, look ahead
# dev screen advanced
# CH4_lass0 <- glmnet_lassoCV2(y_var, x_vars, Guild_CH40)
# glmnet_coefs(CH4_lass0)            
### dev advanced check -- coeffs
# glmnet_coefs(CH4_lass00)
# glmnet_coefs(CH4_lass0.2)
# dev advanced check -- plots
# CH4_lass0_cplot <- plot_lasso_coefs(CH4_lass0, guild_colors)

 
 ########################################################################################
 ### b) extract LASSO coeffs function
 
# extract lasso coeffs
glmnet_coefs = function(model, output = 'coefs', print_ns = FALSE){
    
    # Extract model coeffs
    coef_m <- data.frame(data.matrix(coef(model)))   # extract coeffs
    names(coef_m) <- 'coef'                          # label coef  
    coef_m$var <- row.names(coef_m)                  # get feats
    
    # get only model selected features
    sig_ests <- coef_m[!coef_m$coef==0,]             # drop NS feats  
    sig_ests <- sig_ests[-1,]                        # drop intercept
    
    sig_feats <- sig_ests$var                        # features only
    
    # non-sig feats
    all_feats <- coef_m$var
    ns_int <- all_feats[!all_feats %in% sig_feats]
    ns <- unlist(paste(ns_int[-1], collapse = ', ' ))
    ns_out <- paste('non-sig: ', ns)
    
    if (print_ns == TRUE){
    print(ns_out)}
    
    if (output == 'vars') {
        data_out <- sig_feats
        } else {data_out <- sig_ests}
      
    return(data_out)
}

# demonstrate function
# glmnet_coefs(CH4_lass0)                             # return table of coefs (e.g. for plotting)
vars <- glmnet_coefs(CH4_lass0, output = 'vars')      # return only list of signif vars
vars     #'CO2_mg_m2_h''Bulk_dens''H2O_FPS''pH''N''SO4''CH4_ac''MOB_II''MOB_IIa''AOB'

########################################################################################
### c) plot lasso coefficients
# rename col to var
names(guild_colors)[1] <- "var"
#colors <- guild_colors
#model <- CH4_lass0

plot_lasso_coefs = function(model_out, colors, nafill = '#708090'){

    # dev lasso effects plot
    coefs <- glmnet_coefs(model_out)              # Get significant effects

    # merge w colors
    coef_color <- data.frame(merge(coefs, colors, all.x=T, by = 'var'))     # coef_color
    colors2 <- as.character(coef_color$color)                               # only levels of var used  
    colors2[is.na(colors2)] <- '#708090'                                    # replace NA color 
    coef_color$color <- as.factor(colors2)                                  # replace colors (levels)
                         
    # reorder data / factors by coef
    coef_color <- coef_color[order(coef_color$coef),]                       # sort by coef
    coef_color$var <- reorder(coef_color$var, coef_color$coef)              # reorder var by coef
    coef_color$color <- reorder(coef_color$color, coef_color$coef)          # reorder color by coef

    # Make ggplot 
    plot_colors <- c(paste(coef_color$color))                               # Get colors for ggplot

    plot <- ggplot(data=coef_color, 
                   aes(x=var, y=coef, fill = var)) +       # base plot
                geom_bar(stat="identity") + 
                scale_fill_manual(values= plot_colors) +  # fill by 
                #geom_hline(yintercept=0) +                                            # add y = 0 line
                coord_flip() + 
                theme_minimal() + 
                theme(legend.position = "none")       # flip; clean plot

    return(plot)
}

# Demonstrate function
names(guild_colors)[1] <- 'var'

CH4_lass0_cplot <- plot_lasso_coefs(CH4_lass0, guild_colors)

options(repr.plot.width=6, repr.plot.height=4)
CH4_lass0_cplot

########################################################################################
### d) all subsets regression GLM screening function -- DEPRECATE FOR PARALLEL BELOW (e)?

# to determine models with _all significant features_, not done by leaps, meifly, bestglm packages
# modified from: https://rpubs.com/kaz_yos/exhaustive with coeffs & significance tests added
# could further get R2 if run LM in parallel, not sure its needed.

all_subsets_whh = function(outcome, predictors, dataset, p_cut, sort = "AIC"){
    
    ## Create list of models
    list.of.models <- lapply(seq_along((predictors)), function(n) {

        left.hand.side  <- outcome
        right.hand.side <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")

        paste(left.hand.side, right.hand.side, sep = "  ~  ")
    })

    ## Convert to a vector
    vector.of.models <- unlist(list.of.models)

    ## Fit glm to all models
    list.of.fits <- lapply(vector.of.models, function(x) {

        formula    <- as.formula(x)
        fit        <- glm(formula, data = dataset)
        result.AIC <- extractAIC(fit)

        # all sig filter place
        coeffsM <- data.frame(coef(summary(fit))[,1:4])            # Get coeffs & significance
        coeffsM <- coeffsM[-1,]                                    # drops intercept from NS test

        # test for significance of all feats
        bad_model <- ifelse(any(coeffsM$Pr...t.. > p_cut) == TRUE, "no", "yes")
               
        # get BIC
        BIC <- BIC(fit)
        AICc <- AICc(fit)
        
        ## combine data
        data.frame(num.predictors = result.AIC[1],
               AIC            = result.AIC[2],
               BIC            = BIC,
               AICc           = AICc,
               model          = x,
               all_sig            = bad_model)
    })

    ## Collapse to a data frame
    result <- data.frame(do.call(rbind, list.of.fits))
    result <- result[order(result[,sort]),]
    result <- result[result$all_sig == 'yes',]
    return(result)
}

# use function, 10s? / 1023 permutes
# subs <- all_subsets_whh('CH4_logn1', vars, Guild_CH4_dS_df, p_cut = 0.05, sort = "AIC")
# dim(subs); subs[1:20,]  # show output

# MOB_mods <- subs[grep("MOB", subs$model),]
# head(MOB_mods)

########################################################################################
### e) Parallel all subsets regression GLM screening function
# library(parallel)

all_subsets_whhP = function(outcome, predictors, dataset, p_cut, sort = "AIC", ncores = 20 ){
    
    ## Create list of models
    list.of.models <- lapply(seq_along((predictors)), function(n) {

        left.hand.side  <- outcome
        right.hand.side <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")

        paste(left.hand.side, right.hand.side, sep = "  ~  ")
    })

    ## Convert to a vector
    vector.of.models <- unlist(list.of.models)

    ## Fit glm to all models
    list.of.fits <- mclapply(vector.of.models, function(x) {

        formula    <- as.formula(x)
        fit        <- glm(formula, data = dataset)
        result.AIC <- extractAIC(fit)

        # all sig filter place
        coeffsM <- data.frame(coef(summary(fit))[,1:4])            # Get coeffs & significance
        coeffsM <- coeffsM[-1,]                                    # drops intercept from NS test

        # test for significance of all feats
        bad_model <- ifelse(any(coeffsM$Pr...t.. > p_cut) == TRUE, "no", "yes")
               
        # get BIC
        BIC <- BIC(fit)
        AICc <- AICc(fit)
        
        ## combine data
        data.frame(n_pred = result.AIC[1],
               AIC            = result.AIC[2],
               BIC            = BIC,
               AICc           = AICc,
               model           = x,
               all_sig            = bad_model)
    }, mc.cores = ncores)

    ## Collapse to a data frame
    #result <- data.frame(do.call(rbind, list.of.fits))
    result <- data.frame(rbindlist(list.of.fits))                 # Much faster than do.call!
    result <- result[order(result[,sort]),]
    result <- result[result$all_sig == 'yes',]
    return(result)
}

# use function, 10s? / 1023 permutes single core; 1.5s w 20 cores
subs <- all_subsets_whhP('CH4_logn1', vars, Guild_CH4_dS_df, p_cut = 0.05, sort = "AIC", ncores = 20)

# dim(subs); 
subs[1:20,]  # show output

########################################################################################
### f) filter significant SEM models (after GLM screening) -- -- DEPRECATE FOR PARALLEL BELOW (e)?
# make function for SEM selection 

sem_filter = function(subsets, data, top_n_models = 25, sort_by = 'bic', n_mod_ret = 10){
    
    # reduce candidate models by top n (subtests)
    mod_cands <- subsets[1:top_n_models,]                              # subset
    mods <- mod_cands$model                                            # get models
    
    # orig feats to keep 
    keep_orig <- c("n_pred", "AIC", "BIC", "AICc")
    orig_feats <- mod_cands[keep_orig]
    names(orig_feats) <- c("n_pred","aic_glm","bic_glm","aicc_glm")
    orig_feats$n_pred <- orig_feats$n_pred - 1                         # note init n_pred incl. y 
   
    # Get performance metrics for each
    mod_outs = {}

    for (i in seq(1:length(mods))){

        # Get SEM model
        form <- paste(mods[i])
        SEM_mod <- sem(form, data = data, fixed.x=FALSE, estimator = 'mlm', orthogonal = TRUE)
    
        # Get R2 data
        r2 <- get_SEM_R2s(SEM_mod)
        R2 <-r2[,2]
    
        # Get aic, bic data
        keep <- c('aic','bic','npar')
        fits <- sem_fit_sum(SEM_mod, sem_fit_params)
        fits <- fits[keep]
      
        # all feats significant?
        sig_test <- get_SEM_nonsig_vars(SEM_mod, p_cut =0.05)
        if (dim(sig_test)[1] == 0){
            all_sig <- 'yes'
            } else {
            all_sig <- 'no'
            }
        
        # get feats from original data
        keep_orig <- orig_feats[i,]        
                            
        # combine feats & assign to list
        out <- data.frame(fits, R2, form, all_sig, keep_orig)
        mod_outs[[i]] <- out #paste(mods[i])
    
        }

    # Assemble and sort list
    output <- data.frame(do.call(rbind, mod_outs))
    output <- output[order(output[sort_by]),]
    
    # clean and return output
    model <- seq(1:(length(mods)))
    output <- data.frame(model, output)
    
    
    output <- output[output$all_sig == 'yes',]
    output <- output[1:(n_mod_ret),]
    output <- output[complete.cases(output),]
    row.names(output) <- output$model
    return(output)
    
}

# Demonstrate function
semF <-sem_filter(subs, data = Guild_CH4_dS_df, top_n_models = 140, sort_by = 'aic', n_mod_ret = 140)
#semF

# 140 models in 20s
# Parallel SEM filter

########################################################################################
## g) Parallel filter significant SEM models (after GLM screening) 

# make function for SEM selection 
sem_filterP = function(subsets, data, top_n_models = 25, sort_by = 'bic', ncores = 20){       # n_mod_ret = 10, 
    
    # reduce candidate models by top n (subtests)
    mod_cands <- subsets[1:top_n_models,]                              # subset
    mods <- mod_cands$model                                            # get models
    
    # orig feats to keep 
    keep_orig <- c("model","n_pred", "AIC", "BIC", "AICc")
    orig_feats <- data.frame(mod_cands[keep_orig])
    names(orig_feats) <- c('model',"n_pred","aic_glm","bic_glm","aicc_glm")
    orig_feats$n_pred <- orig_feats$n_pred - 1                         # note init n_pred incl. y 
   
    # Get performance metrics for each
    mod_outs <- mclapply(mods, function(x){

        # Get SEM model
        form <- paste(x)
        SEM_mod <- sem(form, data = data, fixed.x=FALSE, estimator = 'mlm', orthogonal = TRUE)
    
        # Get R2 data
        r2 <- get_SEM_R2s(SEM_mod)
        R2 <-r2[,2]
    
        # Get aic, bic data
        keep <- c('aic','bic','npar')
        fits <- sem_fit_sum(SEM_mod, sem_fit_params)
        fits <- fits[keep]
      
        # all feats significant?
        sig_test <- get_SEM_nonsig_vars(SEM_mod, p_cut =0.05)
        if (dim(sig_test)[1] == 0){
            all_sig <- 'yes'
            } else {
            all_sig <- 'no'
            }
        
        # get feats from original data   
        keep_orig <- orig_feats[orig_feats$model == x,]  
        keep_orig <- keep_orig[,-1]
        
        # combine feats & assign to list
        data.frame(fits, R2, form, all_sig, keep_orig)

        }, mc.cores = ncores)

    # Assemble and sort list
    output <- data.frame(do.call(rbind, mod_outs))
    # output <- data.frame(rbindlist(mod_outs))         - rbindlist abandoned, no speedup & crash on NA
    output <- output[order(output[sort_by]),]
    
    # clean and return output
    model <- seq(1:(length(mods)))
    output <- data.frame(model, output)
    output <- output[output$all_sig == 'yes',]
    #output <- output[1:(n_mod_ret),]
    output <- output[complete.cases(output),]
    output <- data.frame(output)
    row.names(output) <- output$model
    
    # make cols numeric, fix forms       
    output2 <- data.frame(sapply(output, as.numeric))
    output2$form <-output$form
    
    return(output2)
    
}

semF <-sem_filterP(subs, data = Guild_CH4_dS_df, top_n_models = 140, sort_by = 'aic', ncores = 20) # depricate , n_mod_ret = 140,  
head(semF)
# 140 models in 22s (lapply, single core?  but looks like 5 in activity monitor?)
# 140 models in 3.5s (mcapply, 20 cores)


# test dev / debug
# asemF <- sem_filterP(aCH4_lass0_subs, data = Guild_CH40, 
#                                   top_n_models = 500, sort_by = 'bic', n_mod_ret = 500, ncores = 20)
#asemF <- aCH4_lass0_semFhead(asemF)
#asemF[asemF$n_pred == "...",]
# (is.na(asemF$n_pred) == FALSE)
#asemF[50:100,]


########################################################################################
### h) filter best SEM models by LM residuals
# function for R2 vs. AIC max by reg. resids


filter_byLM_resid = function(sem_tab, n_keep_split = 4, sort_by = "resid", split_by = "n_pred",
#filter_lSEM_tab = function(sem_tab, n_keep_split = 4, sort_by = "resid", split_by = "n_pred",
                             R2_filter = 0, resid_0filter = TRUE){

    # get linear model & add residuals to table
    lm <- lm(R2 ~ aic, data = sem_tab)
    resid <- round(resid(lm),3)
    sem_tab$resid <- unlist(resid)

    # sort, split, and filter list
    sem_tab <- sem_tab[order(-sem_tab[,sort_by]),]
    sem_tab_split <- split(sem_tab, sem_tab[,split_by])

    filtered_list <- lapply(sem_tab_split, function(x) {head(x, n_keep_split)}) 
    sem_tab_filt <- rbindlist(rev(filtered_list))

    if (resid_0filter == TRUE) {
        out_tab <- sem_tab_filt[sem_tab_filt$resid > 0,]
        } else {out_tab <- sem_tab_filt}
    
    out_tab <- data.frame(out_tab[out_tab$R2 > R2_filter])
    out_tab
    
    }


# demonstrate function
filt_models <- filter_byLM_resid(semF, n_keep_split = 3, R2_filter = 0.5, resid_0filter = TRUE)
#filt_models

# filter_lSEM_tab(semF, n_keep_split = 4)
# filter_lSEM_tab(sem_tab, n_keep_split = 4, sort_by = "resid", split_by = "n_pred", R2_filter = 0, resid_filter = TRUE)

# preview filter function on plot
#lasso_sel_sumplot <- plot_SEM_lasso_mods4_lin(aCH4_lass0, guild_colors, asemF, filt_models, exclude_mods = '')
#lasso_sel_sumplot <- plot_SEM_lasso_mods4_lin(CH4_lass0, guild_colors, semF, filt_models, exclude_mods = '')
# options(repr.plot.width=12, repr.plot.height= 3)
# lasso_sel_sumplot

