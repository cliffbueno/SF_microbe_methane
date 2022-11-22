# SEM - Lasso semi-auto for SF OTU guilds v1.12
#-  runs under r > 3.6 (glmnet for LASSO)

# Assume have loaded modules 0-6
# likely pollulted with examples, possibly cruft


########################################################################################
### 7) SEM model / formula sorting by lasso variable importance


# extract model formula features into columns 
model_feat_ext <- function(mod_ext_df, var){

    exact <- paste0('\\<', var, '\\>')
    mod_ext_df[var] <- grepl(mod_ext_df$form, pattern = exact)
    return(mod_ext_df)

}

# HELPER function: single col TRUE / FALSE replace with value
TF_replace = function(coef_tab, var, value){
    var <- ifelse(coef_tab[,var] == TRUE, coef_tab[,value], 0)
    return(var)
}

# Function for splitting sem formulas & resorting form vars (by LASSO abs coef/ VIP)
coef_sort_sem_form_tab = function(lasso_model, sem_tab, fxn_operator = '~ ', sep_table = 'FALSE'){
    
    ### Prepare empty var inclusion table ###

    # sort lasso coefs by abs. vals
    coefs <- glmnet_coefs(lasso_model)                                                        # FXN get coefs  
    coefs$abs <- abs(coefs$coef)                                                              # get coef abs
    coefs <- coefs[order(-coefs$abs),]                                                        # sort

    # get vars & model matrix data
    ranked_vars <- coefs$var                                                                  # get ranked vars 
    n_coef <- length(ranked_vars)                                                             # n vars   (cols)
    n_mods <- dim(sem_tab)[1]                                                                 # n models (rows)             

    # make vars, model matrix DF                                                              # w/ empty model components
    coef_m <- data.frame(matrix(ncol=n_coef,nrow=n_mods, dimnames=list(NULL, ranked_vars)))   # make empty df
    sem_mods_vars <- cbind(sem_tab, coef_m)                                                   # join sem models

       
    ### Fill model matrix w. true / false ###
    for (i in ranked_vars)                                                                   # for each var
         sem_mods_vars  <- model_feat_ext(sem_mods_vars, i)                                  # use FXN for fill

 
    ### Rebuild sorted model formulas 
    # Get y var to redo formulas, (assume y same in all)
    form1 <- paste(sem_mods_vars$form[1])                                                    # first formula 
    y_var <- unlist(strsplit(form1, fxn_operator))[1]                                        # get y var
    y_var <- paste0(y_var, fxn_operator)                                                     # add op back

    # Get x var formula componets from TF data
    sort_Xform <- apply(sem_mods_vars, 1, function(r) paste0(names(r)[r == TRUE], collapse = " + "))
    sem_mods_vars$form <- paste0(y_var, sort_Xform)                                         # add y, replace form                   
                                 
    ### Select output mode ###
    if (sep_table == 'FALSE'){
        out_tab <- sem_mods_vars[names(sem_tab)]} else {                                    # revert to input cols
        out_tab <- sem_mods_vars}                                                           # keep split cols (for plots)  
    
    out_tab
    
}

# Demonstrate function:
coef_tab <-coef_sort_sem_form_tab(CH4_lass0, sel_models, sep_table = 'FALSE', fxn_operator = '~ ')  # Get sorted table
#coef_tab <- coef_sort_sem_form_tab(CH4_lass0, sel_models, sep_table = 'TRUE', fxn_operator = '~ ')   # Get full table (for plots)
coef_tab

# coef matrix function (trimmed from "coef_sort_sem_form_tab" function in Leaps like plot, then expanded)
# some legacy crud in here, no need for color var in TF replace, for ex...sort_var, etc...
# sem_model_coef_mtx = function(lasso_model, sel_models, sort_var = "aic", color_var = "R2", exclude_mods = ""){

sort_sem_model_coef_mtx = function(lasso_model, sel_models, VIP_sort = TRUE, rev_order = TRUE, color_var = "R2"){
   
    # sorts sem formula table rows by variables & their lasso VIP (for sem model trees, or heatmaps)
    
    ### replace T/F var columns with data (color_var) ###
    # get table of vars encoded as T/F, with other model scoring data
    coef_tab <- coef_sort_sem_form_tab(lasso_model, sel_models, sep_table = 'TRUE', fxn_operator = '~ ') 
   
    # Get T/F variable sets for apply:
    coef_cols <- names(coef_tab)                                 # all columns
    logic_cols_test <- unlist(lapply(coef_tab, is.logical))      # which are T/F
    log_cols <- coef_cols[logic_cols_test]                       # get vect ot T/F columns only 
 
    # Apply TF replace function to logical, passing others 
    var_col_list = {}                                            # empty list
    for (i in log_cols)                                          # for each T/F col 
        var_col_list[[i]] <- TF_replace(coef_tab, i, color_var)  # insert val or 0 -- HELPER TF_replace row-wise
    
    # Make DF with binary data
    coef_tab2 <- data.frame(var_col_list)                        # replace log cols in orig data
    row.names(coef_tab2) <- coef_tab$form                        # get row names as formulas
    
    # all vars 1, for VIP sort OR R2 instead of sort             # VIP = lasso var imp = scaled coefs 
    if (VIP_sort == TRUE) {
        coef_tab2[!coef_tab2== 0] <- 1}                          # all vars T = 1, else 0 
        else {}                                                  # else T = R2 for heatmap, breaks VIP sort
    
    ### Sort models by vars included, by lasso VIP rank ###
    n_col <- dim(coef_tab2)[2]                                   # get number of columns
    col_nos <- seq(1:n_col)                                      # get vector of col_nos
    coef_tab2 <- dfOrder(coef_tab2, col_nos)                     # sort rows by all columns, in order of columns 
    
    if (rev_order == TRUE){
        coef_tab2 <- coef_tab2[rev(seq_len(nrow(coef_tab2))),]}  # reverse row order
        else{}
    
    return(coef_tab2)

}

# get demo input data
#filt_models <- filter_lSEM_tab(semF, n_keep_split = 3, R2_filter = 0.5, resid_filter = TRUE)

# head(subs)

# demonstrate function 
matrix <-sort_sem_model_coef_mtx(CH4_lass0, filt_models)                 # , rev_order = F)
#matrix <-sort_sem_model_coef_mtx(CH4_lass0, subs2, color_var = "n_pred")
#matrix <-sort_sem_model_coef_mtx(CH4_lass0, asubs, color_var = "n_pred")
matrix <-sort_sem_model_coef_mtx(CH4_lass0, filt_models, VIP_sort = F)  # for heatmap, later
dim(matrix); head(matrix)
matU <- unique(matrix); dim(matU)
# show plot of data (DELETE later)
# library(Heatplus)
# suppressMessages(library(vegan))

# basic heatmap2 to show dendograms
# options(repr.plot.width=7, repr.plot.height=7)
# hm <- heatmap(as.matrix(matrix), Rowv = NA, Colv = NA, margins = c(8, 30))
# hm
vip_sort_sem_models = function(lasso_model, sem_model_matrix){
    
    # get coef weights from lasso
    coefs <- glmnet_coefs(lasso_model)                      # return only list of signif vars
    vars <- unlist(names(sem_model_matrix))                 # get var names 
    new_matrix <- sem_model_matrix                          # new matrix copy 

    # fill in lasso VIP coef weights for each time var included
    for (i in vars){                                            # for each var 
        value <- coefs[coefs$var == i, 1]                       # get coef
        new_matrix[i] <- ifelse(!sem_model_matrix[i]==0, value, 0)}        # if matrix = 1, coef, else 0
    
    return(new_matrix)
}

# demonstrate function
vip_matrix <- vip_sort_sem_models(CH4_lass0, matrix)
head(vip_matrix)

sem_model_vip_sorter = function(lasso_model, filt_models, color_var = 'n_pred'){
    
    
  ### 1) Get data needed from existing functions ###                                       # likely could more gracefully refactor
    # VIP sort formulas & rows into matrix
    matrix <- sort_sem_model_coef_mtx(lasso_model, filt_models, color_var = color_var)     # VIP formulas matrix   
    vip_matrix <- vip_sort_sem_models(lasso_model, matrix)                                 # VIP sort model rows  
    
    
  ### 2) Patch sem_table columns / values to align with VIP sorted ### 
    # Get all data with TF columns per var (unsorted)
    model_tab <-coef_sort_sem_form_tab(lasso_model, filt_models, sep_table = 'TRUE', fxn_operator = '~ ')

    # Apply TF replace function to logical columns, here
    log_cols <- names(vip_matrix)                                # get logical columns, here from VIP matrix ! 
    var_col_list = {}                                            # empty list
    for (i in log_cols)                                          # for each T/F col 
        var_col_list[[i]] <- TF_replace(model_tab, i, color_var) # insert val or 0 -- HELPER TF_replace row-wise
    
    # Get dataframe of vars, make binary                         # TF replace set up to fill from rows, modify  
    var_cols_mtx <- data.frame(var_col_list)                     # insert TF replace into DF 
    var_cols_mtx[!var_cols_mtx == 0] <-1                         # replace all !=0 w 1 
    
    # replace TF vars in model tab w binary
    model_tab[log_cols] <-var_cols_mtx                           # replace logical columns in model table  
 
    
  ### 3) Merge, sort & clean models sorted by VIP ###
    # prepare VIP matrix for merge
    vip_matrix[!vip_matrix == 0] <- 1                                            
    vip_matrix$vip_indx <- seq(1:nrow(vip_matrix))
       
    # merge model_tab & vip matrix
    combined_matrix <- merge(model_tab, vip_matrix, by = log_cols)
    combined_matrix <- combined_matrix[order(combined_matrix$vip_indx),]
    
    # drop var sorting cols, finalize table
    all_cols <- names(combined_matrix)
    keep_cols <- all_cols[!all_cols %in% log_cols]
    out_model_table <- combined_matrix[keep_cols]
    
    out_model_table
    
}

filt_models_VIP <- sem_model_vip_sorter(CH4_lass0, filt_models)
#filt_models_VIP#; filt_models
