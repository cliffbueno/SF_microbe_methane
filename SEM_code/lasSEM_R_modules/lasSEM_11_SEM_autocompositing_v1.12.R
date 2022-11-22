# Lasso semi-auto for SF OTU guilds v1.12
# runs under r > 3.6 (glmnet for LASSO)

# Assume have loaded modules 0-10
# likely pollulted with examples, possibly cruft

########################################################################################
### 11) SEM auto-compositing functions


########################################################################################
### 00) TEST DATA

# setwd("~/Desktop/Sal_SEM_models")

## i) uncomposited model table example
### corresponding uncomposited test input is:
composite_test_semEx7 <- read.table("SEM_models/fxn_example_tables/SEM_Exper7_res_CH4_lasso.txt", sep='\t') #in SEM models
semEx7_res <- composite_test_semEx7
semEx7_res

# new SUBSETS w completness scenarios:
# Doesn't work for auto-composite!! ?

keep_mods2 <- c('69','79','107','113','116','118','121')  # same as manual exp 7, all pass
keep_mods3 <- c(keep_mods2, '75')                         # only model not in manual with all passing
keep_mods5 <- c(keep_mods2, '10', '34')                   # add single comp fails
keep_mods6 <- c(keep_mods2, '8','16')                     # double comp fails

# Subsetting data, note just changing keep mods below
# semEx7_res2 <- keep_selected_mods_in_sem_obj(semEx7_res, keep_mods2, sem_object=FALSE)
semEx7_res2 <- semEx7_res[semEx7_res$model %in% keep_mods2,]
#keep_mods <- c('ch4L0C.69','ch4L0C.79','ch4L0C.107','ch4L0C.113','ch4L0C.116','ch4L0C.118','ch4L0C.121')
semEx7_res2


# ii) original dev test data from CH4_all lasSEM models/ objects
# Resulting output ALL lassSEM objects:
# c(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSO4_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem,  CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)
# sem_object_list <- list(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem, CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)
#sem_object_list#$plot
# head(filt_models)
# head(filt_mods)# filt_mod <- CH4_mod_all_sem$filt_models
# head(filt_mod)# filt_mod_ALL <- CH4_mod_ALL_sem$filt_models
# head(filt_mod_ALL)


# iii) composite model test data -- REPEATED below, find / delete 
# adding SEM response vars with string defs
CH4 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h'
CH4_H2 <- 'CH4_H2~ SRB_syn + CO2_mg_m2_h'
MOB_IIa <- 'MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac'
MOB_I <- 'MOB_I ~ CH4_H2 +CH4_ac'

# Get vector of response vars as strings
Other <- ""
response_vars <- c("CH4", "CH4_H2", "MOB_I","MOB_IIa", "Other")


# Test function -- from module 3) 
# model update function, copied from elsewhere (
sem_models.2 <- update_sem_model_table(model_test, base_model = "ch4_mod0", new_mod_name = "ch4_mod0a.1", 
                      notes = "drop CH4~MOB_I",
                      formula_edits = c("CH4", "CH4_ug_m2_h ~ CH4_ac + CH4_H2 + MOB_IIa + CO2_mg_m2_h"))
                      # output_prefix = "sem."#non_feat_cols = c("model", "base", "notes", "formula"))

sem_models.2
##### get SEM models, needs to be parallel !!# Demonstrate function
#sem_models.2_out <- run_compare_sem_models(sem_models.2, data, estimator = "mlm")
# sem_models.2_out


########################################################################################
### 0) Helper functions

### a) 
remove_yvar_from_df_formula_column =function(df, col_name, fxn_operator = '~ '){
    
    # start function 
    row_data = paste(df[,col_name])                                              # get vect of rows in column

    form_list = list()                                                           # empty list

    for (i in (seq_along(row_data)))     {                                       # for each column
        form_list[i] <- paste(unlist(strsplit(row_data[i], '~ '))[2])            # split, x_vars only to list
        }

    col_replace <- unlist(form_list)                                             # unlist
    df[,col_name] <- col_replace                                                 # replace column data
    return(df)                                                                   # return df

}

# demonstrate function
remove_yvar_from_df_formula_column(sem_models.2, 'MOB_IIa')
#remove_yvar_from_df_formula_column(sem_models.2, 'formula') #-- just drop this one!


### b) 
keep_selected_mods_in_sem_obj = function(SEM_mod_obj, keep_models_list, 
                                   model_col = "model", form_col = "form", sem_object =TRUE){
    
    # SEM object or model table
    if (sem_object == TRUE){
        sem_models <- SEM_mod_obj$sem_models                              # get SEM models from object
        } else {
        sem_models <- SEM_mod_obj                                         # input is SEM model table
    }
    
    # make model numbers char
    sem_models[model_col] <- as.character(sem_models[,model_col])         # model col for model number

    # Get keep_models subset of SEM_models
    sem_models_keep <- sem_models[sem_models[,model_col] %in% keep_models_list,]
    
    return(sem_models_keep)

    }

# Demonstrate fxn (with object, defaults) --- note object used appears later in workflow
# keep_mods <- c("25", "41", "40", "39")
# keep_selected_mods_in_sem_obj(MOBIIa_delt_ngt_sem, keep_mods)

# fxn with table as input, creating output for next fxns
#keep_mods = c('5','11','12','13','15','16')
# selM_MOBIIa <- keep_selected_mods_in_sem_obj(MOB_IIa_models, keep_mods, sem_object = FALSE,
#                                           model_col = "model", form_col = "form") # optional w defaults
# selM_MOBIIa
# fxn with table as input, creating output for next fxns
# keep_mods = c('ch4L0C.69','ch4L0C.107','ch4L0C.118','ch4L0C.121')
# selM_CH4_comp <- keep_selected_mods_in_sem_obj(ch4_lass0_comp_sem_tab, keep_mods, sem_object = FALSE,
#                                            model_col = "model", form_col = "form") # optional w defaults
# selM_CH4_comp


### c) 
drop_vars_from_df = function(df, drop_vars){

    df_names <- unlist(names(df))
    keep_df_cols <- df_names[!df_names %in% drop_vars]
    df <- df[keep_df_cols]
    return(df)
}

# demonstrate function -- need more upstream 
# head(selM_CH4_comp)
# drop_vars <- c("formula", "model_out")
# drop_vars_from_df(selM_CH4_comp, drop_vars)

########################################################################################
## I) assign vars to composite formulae 

#  user defines composite_sets, e.g. = c("CH4_gen", "CH4_ox")... & each is defined as var list:
#  CH4_gen = c('CH4_ac', 'CH4_H2', 'CO2'...); CH4_ox = c('MOB_I','MOB_II', 'MOB_IIa',..)
#  functions allow unassinged vars to be assigned to main y var (root response) or any named composite (below)


##### main user function at end of chain is (g)
# get_composites_model_table(semEx7_res, composite_sets, assign_other = "CH4_gen", compos_fail_collapse_all = FALSE) 
# opt: assign_other = 'y_var' is default, to y (root), may wish to use 'CH4_gen' for example
# opt: composite_fail_collapse_all = TRUE: any failure of composite var flattens root model; FALSE leaves passing (>1 var) intact...

##################################
## 0) test data -- MOVED TO TOP, DEPRECATE 
### XXX DELETE XXXX

# setwd("~/Desktop/Sal_SEM_models")
# corresponding uncomposited test input is:
#composite_test_semEx7 <- read.table("SEM_models/tables/SEM_Exper7_res_CH4_lasso.txt", sep='\t') #in SEM models
#semEx7_res <- composite_test_semEx7
#semEx7_res# new SUBSETS w completness scenarios:
# Doesn't work for auto-composite!! ?

#keep_mods2 <- c('69','79','107','113','116','118','121')  # same as manual exp 7, all pass
#keep_mods3 <- c(keep_mods2, '75')                         # only model not in manual with all passing
#keep_mods5 <- c(keep_mods2, '10', '34')                   # add single comp fails
#keep_mods6 <- c(keep_mods2, '8','16')                     # double comp fails

# Subsetting data, note just changing keep mods below
#semEx7_res2 <- keep_selected_mods_in_sem_obj(semEx7_res, keep_mods2, sem_object=FALSE)
#keep_mods <- c('ch4L0C.69','ch4L0C.79','ch4L0C.107','ch4L0C.113','ch4L0C.116','ch4L0C.118','ch4L0C.121')
#semEx7_res2

# Resulting output ALL lassSEM objects:
# c(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSO4_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem,  CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)
# sem_object_list <- list(CH4_mod_ALL_sem, CH4_mod_all_sem, CH4_mod_fSRB_sem, CH4_mod_fsal.1_sem, CH4_mod_fsal.2_sem, CH4_mod_fch4.2_sem)
#sem_object_list#$plot
# head(filt_models)
# head(filt_mods)# filt_mod <- CH4_mod_all_sem$filt_models
# head(filt_mod)# filt_mod_ALL <- CH4_mod_ALL_sem$filt_models
# head(filt_mod_ALL)


##################################
## 0) USER DEFINE variable assignments for composite functions -- MOVE TO EXAMPLES AT TOP, out to usr notebooks

CH4_gen <- c("CH4_ac", "CH4_H2", "CO2_mg_m2_h", "Salinity.x", "SO4", "SRB", "FeRB", "NP")

CH4_ox <- c("MOB_IIa", "MOB_I", "MOB_II", "Bulk_dens", "H2O_FPS", 
            "N", "AOB", "AOA", "NOB", "NO3_N", "NH4_N", "NO3_NH4")

composite_sets <- c("CH4_gen", "CH4_ox")
# Composite_base <- c('0 + CH4_gen + CH4_ox')

# Remove sample names from ALL data set here-- testing dysfunctional case
Guild_CH4_0 <- Guild_CH4
row.names(Guild_CH4_0) <- Guild_CH4_0$Sample
Guild_CH4_0 <- Guild_CH4_0[,-1]

# prepare data for function
# filt_mod <- CH4_mod_all_sem$filt_models
# head(filt_mod)

# head(filt_models)
# get_unique_Xvars_in_formula_dfcol(filt_mod[1,], form_col= "form")

##################################
### a)  get unique X vars from formulae in DF col
#   note, depends on helper FXN "remove y_var_from_formula_column"

get_unique_Xvars_in_formula_dfcol = function(model_table, form_col = "form"){
    
    # Remove y vars, get only formulae
    no_y_table <- remove_yvar_from_df_formula_column(model_table, form_col)  # drop y var 
    forms_no_y <- no_y_table[,form_col]
    forms_no_y <- str_replace_all(forms_no_y, " ", "")                       # KLUDGE, sometimes extra space introduced

    # Get unique X vars among ALL of these formulas
    all_forms <- paste(forms_no_y, collapse = "+")                          # collapse rows into single string
    all_vars <- strsplit(all_forms, "\\+")                                  # split forms into vars 
    unique_vars <- unique(all_vars[[1]])                                    # get unique vars
    
    return(unique_vars) 
}      

# Demonstrate function
get_unique_Xvars_in_formula_dfcol(semEx7_res, form_col = "form")

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
# get_unique_Xvars_in_formula_dfcol(filt_models, form_col = "form"
# get_unique_Xvars_in_formula_dfcol(filt_mod, form_col = "form")
# get_unique_Xvars_in_formula_dfcol(all_mods, form_col = "form")

# filt_mod2<- filt_mod
# filt_mod2$form <- as.factor(filt_mod2$form)


##################################
### b) assign models to composite split df

assign_to_composite_split_df = function(model_table, form_col = "form", composite_sets){
    
    # init accumulators
    var_list       <- list()                                                # X var list / formula
    var_group_list <- list()                                                # X vars in set / formula / set
    out_group_list <- list()                                                # X vars not in set ""
    set_list       <- list()                                                # X lists for each set
    out_list       <- list()                                                # X lists not in each set
    other_list     <- list()                                                # X vars in nether set 
    #other_form     <- list()                                                # X vars in nether set, " + " sep
    
    # Loop over sets (j), rows in model table
    for (j in composite_sets){                                              # get lists of vars from each set

        set_vars <- get(j)                                                  # get set vars data

        for (i in 1:nrow(model_table)){

            # remove y var and get list of X vars from formula
            var_list <- get_unique_Xvars_in_formula_dfcol(model_table[i,], form_col = form_col)
            var_group_list[[i]] <- set_vars[set_vars %in% var_list]               # vars in set   
            var_group_list[[i]] <- paste(var_group_list[[i]], collapse = " + ")   # make formula style (?) 
            
            out_group_list[[i]] <- var_list[!var_list %in% set_vars]              # vars not in set
        }
    
        # Gather results for each set    
        set_list[[j]] <- var_group_list
        out_list[[j]] <- out_group_list
    }    
  
    # Make DFs of vars in and !in sets
    var_sets_df <- data.frame(do.call(cbind, set_list))
    out_vars <- data.frame(do.call(cbind, out_list))                       # df of vars ! in sets
    
    # Get duplicate vars among not-in-set vars (ie not in any set)
    for (i in 1:nrow(out_vars)){

        other_vars <- unlist(c(out_vars[i,]))                             # get all vars in row of df
        names(other_vars) <- NULL                                         # drop their indicies
        other_list[[i]] <- other_vars[duplicated(other_vars)]             # get only duplicated 
        other_list[[i]] <- paste(other_list[[i]], collapse = " + ")
        }
    
    # Add other vars col to var_sets DF
    var_sets_df$other <- other_list
    #var_sets_df$other_fxn <- other_form
    
    return(var_sets_df)
    
}

# Demonstrate function
composite_base_df <- assign_to_composite_split_df(semEx7_res, form_col = "form", composite_sets)
composite_base_df

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
# x <- assign_to_composite(filt_models, form_col = "form", composite_sets)
# composite_base_df <- assign_to_composite_split_df(filt_models, form_col = "form", composite_sets)
#composite_base_df <- assign_to_composite_split_df(filt_mod_ALL, form_col = "form", composite_sets)
#head(filt_mod_ALL)
# composite_base_df# $CH4_gen


##################################
### c) get composite split formulas table

get_composite_split_form_table = function(model_table, composite_sets, assign_other = "y_var", other_name = "other",
                                      split_y = "  ~", new_y_op = "~ ", no_int = TRUE, set_ops = " <~ ",      
                                      form_col = "form", model_col = "model", model_out_pfx = "sem."){
    
    ### 1) Get main model data: IDs, y, x_vars, combine ###
    
    # Extract model numbers, out strings (to assign)
    model <- as.character(model_table[,model_col])
    model_out <- paste0(model_out_pfx, model)
    
    # Extract y var
    formulae <- model_table[,form_col]
    y_var <- str_split(formulae, split_y) [[1]][[1]]           # first row, column; assumes all are same
    
    # Build y formula base
    y_ops <- paste(y_var, new_y_op)                            # "y ~"  
    compositeX <- paste(composite_sets, collapse= " + ")       # "comp1 + comp2..."
    if (no_int == TRUE){
        compositeX <- paste0("0 + ", compositeX)               # "0 + comp1 + ..."
    }
    
    # Get y base formula column, name by y_var
    y_form_base <- data.frame(paste0(y_ops, compositeX), stringsAsFactors = FALSE)                   
    names(y_form_base) <- y_var
            
    # X vars : get df of composite var assignments -- prior function
    composite_base <- assign_to_composite_split_df(model_table, composite_sets, form_col = form_col)
    
    # Combine model, y and X 
    composite_df <- data.frame(model_out, model, y_form_base, composite_base, stringsAsFactors = FALSE)
    
    
    ### 2) Assigning "other" vars to formulae ###
    # check if assign_other to y_var or composite
    if (assign_other == "y_var"){                                                    # if assign_other = y_var                          
        assign_other = y_var}                                                        # substitute y_var above
   
    # combine formulae from assigned & other, where neither are blank
    mix_other <- ifelse(composite_df[,other_name] == '', composite_df[,assign_other], # if other '', '', else         
                 ifelse(composite_df[,assign_other] == '', composite_df[,other_name], # if assign blank, else 
                 paste(composite_df[,assign_other], composite_df[,other_name], sep = " + "))) # combine forms   
    
    mix_other_df <- data.frame(mix_other, stringsAsFactors=FALSE)
    
    # add to composite DF, kludge for assign_other sometimes being transposed (CH4_ox, not CH4_gen, y var ???)    
    if (nrow(mix_other_df) > 1) {                                                    # transp has only 1 row
        composite_df[assign_other] = mix_other_df} else {                            # replace w other DF, else
        composite_df[assign_other] = t(mix_other_df)                                 # replace w t(other) DF
    }
    
    
    ### 3) Clean formulas, validate composites, produce table ###
    # add y <~ 1* to composite x_var formulas, where not blank.  First var '1*' sets scale of composite var (+/-) 
    for (j in composite_sets){                                                       # for each composite
        composite_df[j] <- ifelse(composite_df[,j] == '', '',                        # where not blank
                                  paste0(j, set_ops, "1*",composite_df[,j]))         # add y, op, "1*" = scaling
    }
    
    return(composite_df)

}    

# Demonstrate function
get_composite_split_form_table(semEx7_res, form_col = "form", composite_sets, assign_other = "y_var")
get_composite_split_form_table(semEx7_res, form_col = "form", composite_sets, assign_other = "CH4_gen")

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
# get_composite_split_form_table(filt_mod, form_col = "form", composite_sets, assign_other = "y_var")
# get_composite_split_form_table(filt_mod_ALL, form_col = "form", composite_sets, assign_other = "y_var") #filt_mod_ALL
# get_composite_split_form_table(filt_mod_ALL, composite_sets, assign_other = "CH4_ox") #filt_mod_ALL


##################################
### d) get composite table failures

# second version, testing only bad composites w. paste var name
get_composites_table_failures = function(model_table, composite_sets, assign_other = "y_var", other_name = "other",
                                      split_y = "  ~", new_y_op = "~ ", no_int = TRUE, set_ops = " <~ ",      
                                      form_col = "form", model_col = "model", model_out_pfx = "sem."){
       
    # move this OUTSIDE FXN?  USING any of those options here?    
    # get composite forms model table from above fxn, pass all options
    composite_df <- get_composite_split_form_table(model_table, composite_sets, assign_other = assign_other, 
                                                  other_name = other_name, split_y = split_y, new_y_op = new_y_op, 
                                                  no_int = no_int, set_ops = set_ops, form_col = form_col, 
                                                  model_col = model_col, model_out_pfx = model_out_pfx)
  
    # Extract model numbers, out strings (to assign)
    model <- as.character(composite_df[,model_col])
      
    ### filter non-valid composite models (need > 1 var/composite) ###
    
    # method quirks: 1) str_split keeps y_var & x_var1 together, count as 1 var (not issue)
    #                2) length(vars_i...) counts empty strings as 1 var -- can work around:
    #                    - 2) workarounds: a) need > 2 vars for composite anyway; b) can drop blank w/o counter
    
    # get n vars in composites
    var_counts_plus1 <- list()                                      # row accumulator for length of split forms
    var_counts_compos <- list()                                     # compos. accum. for list of form. len/compos.
    
    for (j in composite_sets){    
       
        for (i in 1:nrow(composite_df)){
                       
            vars_i <- str_split(composite_df[,j], " \\+ ")          # split formulas into x's (keeps y~x1)
            var_counts_plus1[[i]] <-length(vars_i[i][[1]])          # get n elements (x's, BUT "" = 1, not 0)
          }
        
        var_counts_compos[[j]] <- var_counts_plus1
        }
       
    n_vars_compos <- data.frame(do.call(cbind, var_counts_compos))  # DF of variable counts/composite (w caveats above)
    
    
    # do n vars/model/composite exceed SEM min threshold? 
    var_thresh = 1                                                  # var threshold > 1, sem requires 2+)  
    
    for (j in composite_sets){ 
       
         n_vars_compos[j] <- ifelse(n_vars_compos[j] > var_thresh, "", j)
    }
    
    # something funny about data type (<chr[,1]>), simplify to char only
    n_vars_compos <- data.frame(lapply(n_vars_compos, as.character), stringsAsFactors = FALSE)
    n_vars_compos
    
    failed_list <- list()
    for (i in 1:nrow(n_vars_compos)){
        vars <- do.call(c, n_vars_compos[i,])
        names(vars) <- NULL
        vars <- vars[vars != ""]
        failed_list[[i]] <- vars
        }
    
    composite_df$comp_failed <- failed_list
    composite_df <- drop_vars_from_df(composite_df, other_name)
    composite_df
    
}  

# Demonstrate function
comp_f7 <- get_composites_table_failures(semEx7_res, form_col = "form", composite_sets, assign_other = "y_var")
comp_f7

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
#get_composites_table_failures(filt_mod, form_col = "form", composite_sets, assign_other = "y_var") #filt_mod_ALL
# y <- get_composites_table_failures(filt_mod_ALL, composite_sets, assign_other = "CH4_ox") #filt_mod_ALL


##################################
### e) helper: send one bad composite x to y var
#      - note helper requires input like comp_f7 above (get_composite_table_failures)
#      - is run embedded in other function so no major concern re this format

send_one_bad_composite_x_to_y_var = function(model_table, fail_comp, y_var, comp_string = " <~ 1\\*"){
    
    # strip formula
    fail_comp_data <- model_table[,fail_comp]                                    # vect of failed comp formulae
    fail_comp_form_y <- paste0(fail_comp, comp_string)                           # get y + op string to delete 
    fail_comp_keep_vars <- str_replace(fail_comp_data, fail_comp_form_y, "")     # replace it w "", leaving x
    model_table['pass_failed'] <- fail_comp_keep_vars                            # new col of failed comp x

    # replace failed composite name in y formula w/ failed composite x           # y form is base of composite
    model_table[y_var] <- ifelse(model_table[,'pass_failed'] == '',              # if failed x is blank, 
                    str_replace(model_table[,y_var], paste0(" \\+ ", fail_comp), model_table[,'pass_failed']),
                    str_replace(model_table[,y_var], fail_comp, model_table[,'pass_failed']))
    
    # drop all failed component formulae from table. Would have bad behavior in loop, but to be run by row 
    model_table[fail_comp] <- ""                            
    model_table <- drop_vars_from_df(model_table, 'pass_failed')                 # drop new pass_failed_col
    
    return(model_table)
   
}

# Demonstrate function
send_one_bad_composite_x_to_y_var(model_table = comp_f7, fail_comp = "CH4_ox", y_var = "CH4_logn1", 
                              comp_string = " <~ 1\\*")

# Demonstrate fxn, BUT SOURCE DATA DOWNSTREAM !!
# here y is output from "get_composites_table_failures" or equiv (has "comp_failed" col)
#z= y[1,]
#send_bad_composite_x_to_y_var(model_table = z, fail_comp = "CH4_ox", y_var = "CH4_logn1", 
#send_one_bad_composite_x_to_y_var(model_table = y, fail_comp = "CH4_ox", y_var = "CH4_logn1", 
#                              comp_string = " <~ 1\\*")

##################################
### f) helper: send ALL bad composite x to y var (iterates previous fxn)
#      - iterate over multiple failed composites -- not likely w 2 compos, but > 2)... test
#      - note helper requires input like comp_f7 above (get_composite_table_failures), though this FAILED here
#      - is run embedded in other function so no major concern re this format


send_all_bad_composite_x_to_y_var = function(model_table, y_var, comp_string = " <~ 1\\*"){

    results = list()                             

    # Operate on each row, get component_lists and n 
    for (i in 1:nrow(model_table)){
    
        model_table_i <- model_table[i,]                          # Get model table row
        comp_failed_i <- model_table_i[,"comp_failed"][[1]]       # Get comp_failed / row
        n_failed_comp <- length(comp_failed_i)                    # Get length (comp_failed) = n_comp_failed /row
    
        # For row, send FIRST bad composite x to y_var, store as row output
        results_i <- send_one_bad_composite_x_to_y_var(model_table = model_table_i, fail_comp = comp_failed_i[1],
                                                      y_var = y_var, comp_string = comp_string) 
        
        # Accomodate ALL failed composites (if more than 1)
        if(length(comp_failed_i) < 2){                            # if 0 or 1 failed comps, 
            results[[i]] <- results_i                             # keep results_i
        } else {                                                  # if 2 or more failed comps 
            # iterate row function over other failed composites   # note model table = results_i
            for (k in comp_failed_i[-1]){                         # iterate over failed comps (not first):
            results_i <- send_one_bad_composite_x_to_y_var(model_table = results_i, fail_comp = k,
                              y_var = y_var, comp_string = comp_string)   
            }   
        }
        results[[i]] <- results_i                                 # accumulate row results
    }
    out <- do.call(rbind, results)                                       # make accumulated DF
    return(out)
}

# Demonstrate function  not working with Ex7 either, but works embedded in intended host function (g below)
# send_all_bad_composite_x_to_y_var(model_table = y, y_var = "CH4_logn1", comp_string = " <~ 1\\*")

# demonstrate function, BUT SOURCE DATA DOWNSTREAM !!
# here y is output from "get_composites_table_failures" or equiv (has "comp_failed" col)
# send_all_bad_composite_x_to_y_var(model_table = y, y_var = "CH4_logn1", comp_string = " <~ 1\\*")  
    

# send_all_bad_composite_x_to_y_var(model_table = y, y_var = "CH4_logn1", comp_string = " <~ 1\\*")

# BELOW, last "working", then edited to address problems here,
# doesn't work if no fails, part fails, needs if statements around those cases... 

##################################
### g) make model composite table
#    - final wrapper, now reverts non-composite functions

get_composites_model_table = function(model_table, composite_sets, assign_other = "y_var", other_name = "other",
                                      compos_fail_collapse_all = TRUE, # new here, keep some composites?
                                      split_y = "  ~", new_y_op = "~ ", no_int = TRUE, set_ops = " <~ ",      
                                      form_col = "form", model_col = "model", model_out_pfx = "sem."){
    
   ### 1) Get composite fail list data & COUNT fails ###  
    # Get last composite DF with composite fails--    move this OUTSIDE FXN?  USING any of those options here?
    composite_df <- get_composites_table_failures(model_table, composite_sets, assign_other = assign_other, 
                                                  other_name = other_name, split_y = split_y, new_y_op = new_y_op, 
                                                  no_int = no_int, set_ops = set_ops, form_col = form_col, 
                                                  model_col = model_col, model_out_pfx = model_out_pfx)
    # Get model indices
    model <- composite_df[,model_col]
    
    # get data for composite failure classes -- note is list of failed composites/model
    failed_col = "comp_failed"          # col for failed composites, add to FXN options? hardcoded in last FXN
    fail_df <- composite_df[failed_col]
    fail_list <- composite_df[,failed_col]
        
    # count items in failed composites, make DF with model indexing
    fail_len_list <- list()   
    for (i in 1:length(fail_list)){                                       # for rows in composite failed list/col 
        fail_len_list[[i]] <- length((fail_list[i][[1]]))                 # get n items in failed list
    }
    n_fail <- do.call(rbind, fail_len_list)                               # list to rows 
    fails_df <- data.frame(model, n_fail, stringsAsFactors = FALSE)       # add model to fails DF
    
    
   ### 2) SUBSET models on composite pass, fail, OR "neither" (part failure) ###
    n_compos <- length(composite_sets)                                    # another potential FXN option  
    
    pass_models <- fails_df[fails_df["n_fail"] == 0,][,"model"]           # PASS all composites   
    failed_models <- fails_df[fails_df["n_fail"] == n_compos,][,"model"]  # FAIL all composites  
    pass_fail_mods <- c(pass_models, failed_models)                       # list both pass and fail  
    part_fail_models <- model[!model %in% pass_fail_mods]                 # SOME composite fails, others pass
    
    # Fate of failed composite vars (compos_fail_collapse_all == TRUE)
    if(compos_fail_collapse_all == TRUE){                                 # collapse all composites if one fails  
        failed_models <- c(failed_models, part_fail_models)               # add part_failed_models to failed_models
    }    
    
    # Get DATAFRAMES for passing, failed, some fail
    pass_df <- composite_df[model %in% pass_models, ]                     # Data for all passing, ** USE AS IS **
    #failed_df <- composite_df[model %in% failed_models, ]                 # all failing,          modify 
    #part_df <- composite_df[model %in% part_fail_models, ]                # part failed,          replace 

    
   ### 3) Clean failed composite models ###  
       
    if(length(failed_models) > 0){
        
        # Get failed DF
        failed_df <- composite_df[model %in% failed_models, ]                 # all failing,          modify 
    
        # Extract y var name 
        formulae <- model_table[,form_col]                                    # get formulae to split
        y_var <- str_split(formulae, split_y) [[1]][[1]]                      # first row, column; assumes all are same
        replace_orig <- c(model_col, y_var)                                   # c(model, y_var colnames) 
    
        # get orig input formulas for replace in failed out DF
        keep_orig_vars <- c(model_col, form_col)                              # model and formula cols
        keep_orig <- model_table[keep_orig_vars]                              # get mod & form from input table
        orig_fails <- keep_orig[keep_orig[, model_col] %in% failed_models,]   # keep only failed model rows
        orig_fails <- setNames(orig_fails, replace_orig)                      # rename model_col to y_var
 
        # clean failed df  
        failed_df[,composite_sets] <- ''                                      # clean out composite info 
        replace_orig <- c(model_col, y_var)                                   # cols to replace w orig data
        failed_df[replace_orig] <- orig_fails[replace_orig]                   # replace orig data
    
        #Add 0 intercept to composites failures ? 
        if(no_int == TRUE){
            failed_df[y_var] <- str_replace(failed_df[,y_var], new_y_op, paste0(new_y_op, "0 + "))
        }
    } else {failed_df <- data.frame()}
    
    
   ### 4) Clean part failed composite models ###  
    
    if(length(part_fail_models) > 0){        
        
        part_df <- composite_df[model %in% part_fail_models, ]                 # get part df from part_failed 
        
        # Fate of failed composite vars (compos_fail_collapse_all == TRUE)
        if(compos_fail_collapse_all == FALSE){                                 # collapse all composites if one fails  
            part_df <- send_all_bad_composite_x_to_y_var(model_table = part_df, y_var, comp_string = " <~ 1\\*") # LAST ARG TO OPTIONS!
        }
    }  else {part_df <- data.frame()}

    
   ### 5) Combine all results ###    
    if(compos_fail_collapse_all == TRUE){
        out <- rbind(pass_df, failed_df)
    } else {
        out <- rbind(pass_df, part_df, failed_df)
    }

    # RECONSTITUTE FORMULA COLUMN # ID which columns are formulas, get as DF 
    form_cols <- c(y_var, composite_sets)
    form_data <- out[form_cols]                                                       # df of formula cols
    
    # paste formula components for new combined formulae
    paste_args <- c(form_data, sep="\n")                                              # do.call paste args
    formula <- do.call(paste, paste_args)                                             # formulae to DF
    formula <- str_replace(formula, "\n\n", "")                                       # no n\n\ (from empty) 
    out$formula <- formula                                                            # formulae to output
    
    # final cleaning - all formula columns to character -- some glitches obs. otherwise
    for (f in form_cols){
        out[f] <- as.character(out[,f])
    }
    
    return(out) 
}       

# demonstrate function
get_composites_model_table(semEx7_res, composite_sets, assign_other = "CH4_ox", 
                           compos_fail_collapse_all = FALSE)

# Show main user behavior options
semEx7_autocomp <- get_composites_model_table(semEx7_res, composite_sets, assign_other = "CH4_gen", 
                           compos_fail_collapse_all = FALSE)

get_composites_model_table(semEx7_res, composite_sets, assign_other = "y_var", 
                           compos_fail_collapse_all = TRUE)

# Show main user behavior options
semEx7_autocomp <- get_composites_model_table(semEx7_res, composite_sets, assign_other = "CH4_gen", 
                           compos_fail_collapse_all = FALSE)
# # demonstrate function, BUT SOURCE DATA DOWNSTREAM !!
# z <- get_composites_model_table(filt_mod, form_col = "form", composite_sets, assign_other = "y_var") #filt_mod_ALL
# z <- get_composites_model_table(filt_mod, form_col = "form", composite_sets, assign_other = "CH4_ox") 
# y <- get_composites_model_table(filt_mod_ALL, composite_sets, assign_other = "CH4_ox", 
#                            compos_fail_collapse_all = FALSE) #filt_mod_ALL

######################################################################################################
### II) combine SEM tables in hierarchy functions

##################################
### 0) get test data

# setwd("~/Desktop/Sal_SEM_models")


## i) CH4 manual composite model table (Experiment 7)
# get external latent model table
ch4_lass0_comp_sems <- read.table("SEM_models/fxn_example_tables/CH4_flux_latent_mod_form_setup2.txt", header = T, na.strings=c(NA),)
ch4_lass0_comp_sems$other <- ""                   #head(ch4_lass0_comp_sems)

# remove base model formula
ncols_start <- dim(ch4_lass0_comp_sems)[2]
ch4_lass0_comp_sem_tab <- ch4_lass0_comp_sems[, -ncols_start]  # drop base formula

# add formula
#ch4_lass0_comp_sem_tab$formula <- with(ch4_lass0_comp_sem_tab, paste(CH4_flux, CH4_gen, CH4_ox, other, sep = '\n'))
#ch4_lass0_comp_sem_tab[,1:7]

ch4_lass0_comp_sem_tab <- ch4_lass0_comp_sem_tab[,1:7]  
ch4_lass0_comp_sem_tab
# ch4_L0C_mods_keep

# ii) MOBI_IIa models (example branch
MOB_IIa_models <- read.table("SEM_models/fxn_example_tables/SEM_Exper8_res_MOB_IIa_lasso.txt", header = T, sep = '\t')
head(MOB_IIa_models)

# demonstrate function
# keep_mods = c('5','11','12','13','15','16')
# selected_MOBIIa <- select_mods_from_sem_object(MOB_IIa_models, keep_mods,
#                                               model_col = "model", form_col = "form") # optional w defaults
# selected_MOBIIa

# iii) MOB_IIa (SEM OBJECT)
# MOBIIa_delt_ngt_sem

##################################
### a) helper functions


## i) 
drop_vars_from_df = function(df, drop_vars){

    df_names <- unlist(names(df))
    keep_df_cols <- df_names[!df_names %in% drop_vars]
    df <- df[keep_df_cols]
    return(df)
}

# demonstrate function -- need more upstream 
# head(selM_CH4_comp)
# drop_vars <- c("formula", "model_out")
# drop_vars_from_df(selM_CH4_comp, drop_vars)

# ii) 

get_formula_col_list_from_df = function(df, fxn_operator = '~ '){
    
    # Get columns containing formulas
    cols = names(df)                                       # get column names

    col_list = list()                                      # empty list
 
    for (i in seq_along(cols))                             # for each column
        col_list[[i]] = unique(str_detect(df[,i], fxn_operator))    # does contain formula operator?

    keep_cols <- unlist(col_list)                          # unlist T/F
    form_cols <- cols[keep_cols]                                        # get columns with formulas
    form_cols
    
    
}

# Test data:
# Test function
sem_models.2 <- update_sem_model_table(model_test, base_model = "ch4_mod0", new_mod_name = "ch4_mod0a.1", 
                      notes = "drop CH4~MOB_I",
                      formula_edits = c("CH4", "CH4_ug_m2_h ~ CH4_ac + CH4_H2 + MOB_IIa + CO2_mg_m2_h"))
                      # output_prefix = "sem."#non_feat_cols = c("model", "base", "notes", "formula"))

# sem_models.2

# Demonstrate function
get_formula_col_list_from_df(sem_models.2, fxn_operator = '~ ')

##################################
### b) get model formula from SEM results

get_model_formula_from_SEM_results = function(sem_results, set_prefix = '', 
                                          keep_cols = c('model','form'), rename_form = TRUE, fxn_op = " ~"){
    
    model_col <- keep_cols[1]
    form_col <- keep_cols[2]
    
    keep_data <- sem_results[keep_cols]
    keep_data[model_col] <- paste0(set_prefix, ".", keep_data[,model_col])
    
    if (rename_form ==TRUE) {
        
        forms <- as.character(keep_data[,form_col])
        y_vars <- unlist(strsplit(forms, fxn_op))[1]                                        # get y var
        y <- unique(y_vars)
        
        
        names(keep_data)[2] <-y
    }

    return(keep_data)
}

# Demonstrate function
CH4ex7_res2_forms <- get_model_formula_from_SEM_results(semEx7_res2, set_prefix = "ch4.7man")
CH4ex7_res2_forms

MOBIIa_forms <- get_model_formula_from_SEM_results(MOB_IIa_models, set_prefix = "mIIa.man")
MOBIIa_forms

# test data downstream
# MOBIIa_forms <- get_model_formula_from_SEM_results(selM_MOBIIa, set_prefix = "mIIa.man")
# MOBIIa_forms

##################################
### c) combine model layers (branches) in SEM hierarchy

combine_sem_model_layers = function(left_matrix, right_matrix, 
                                    left_mod = "model", right_mod = "model", out_prefix = "sem.", 
                                    Lform_drop = "formula"){

   # LEFT matrix prep #
    
    # drop any existing "formula" column
    left_matrix <- drop_vars_from_df(left_matrix, Lform_drop)
    
    #left_names <- names(left_matrix)
    #keep_L_cols <- left_names[!left_names %in% Lform_drop]
    #left_matrix <- left_matrix[keep_L_cols]
    

    
    # tag left model colname 
    new_lmod_col = paste0(left_mod, ".l")                                                     # L model name
    colnames(left_matrix)[colnames(left_matrix) == left_mod] <- new_lmod_col                  # mark L model name 

    # repeat left df, n = nrow right df
    left_expand <- left_matrix[rep(seq_len(nrow(left_matrix)), nrow(right_matrix)), ]         # expand left df
    
    
   # RIGHT matrix prep # 
    # tag right model colname
    new_rmod_col = paste0(right_mod, ".r")                                                    # R model name
    colnames(right_matrix)[colnames(right_matrix) == left_mod] <- new_rmod_col                # mark R model name
 
    # repeat right df, n = nrow left df; sort by index
    right_expand <- right_matrix[rep(seq_len(nrow(right_matrix)), nrow(left_matrix)), ]       # expand right df
    right_expand <- right_expand[order(as.numeric(row.names(right_expand))),]                 # sort by index

    
   # COMBINE L & R, make combined IDs and formulae #  
    LR_expanded <- cbind(left_expand, right_expand)                                           # combine L & R expanded 
    model <-  paste0(LR_expanded[,new_lmod_col], "_", LR_expanded[,new_rmod_col])             # combine model ids
    model_out <- paste0(out_prefix, model)
    
    # ID which columns are formulas, get as DF                                                # using newer fxn
    form_cols <- get_formula_col_list_from_df(LR_expanded, fxn_operator = '~')                # id formula cols
    form_data <- LR_expanded[form_cols]                                                       # df of formula cols

    # paste formula components for new combined formulae
    paste_args <- c(form_data, sep="\n")                                                      # do.call paste args
    formula <- do.call(paste, paste_args)                                                     # formulae to DF

    # add model IDs to formula data to finish, make all character
    combined_models <- cbind(model_out, model, form_data, formula) 
    combined_models <- data.frame(lapply(combined_models, as.character),stringsAsFactors = FALSE) 
}


# drop formula column from left matrix !!!!  -- no need, FIXED !!!
# selM_CH4_comp2 <- selM_CH4_comp[,1:7] 
# selM_CH4_comp

# Demonstrate function
CH4mob_comp_forms <- combine_sem_model_layers(ch4_lass0_comp_sem_tab, MOBIIa_forms) 
# CH4mob_comp_forms <- combine_sem_model_layers(selM_CH4_comp, MOBIIa_forms) 
                             #optional: # , left_mod = "model", right_mod = "model", 
                                        #  out_prefix = "sem.", Lform_drop = "formula")

#head(CH4mob_comp_forms)
#names(CH4mob_comp_forms)[6] <- 'MOB_IIa'
#head(CH4mob_comp_forms)

#CH4mob_comp_forms


######################################################################################################
### III) Run and compare functions for composite formulas
#        - original functions from 3F above
#        - add try catch for misspecified hierarchical models 
#        - add parallel execution
#        - here investigate behavior of SEM options on composite runs -- added options to fxns



##################################
### 0) test data -- DEPRECATED, DELETE XXX ###


# original test data (small scale, no fail)
# dim(model_test2)
# model_test2
#nrow(sem_model_table)
# sem_model_table$formula

### a) new TEST data from above
# DO NOT USED lasso SCALED DATA IN SEM !!!   log is good, scaled bad
# FURTHER run_compare needs non-default flags on keep_mod_descr  -- OR CHANGE DEFAULTS in orig. Fxn(s)
# Seems to run at 50% total CPU (part parallel??)

#CH4mob_comp_res <- suppressWarnings(run_compare_sem_models(CH4mob_comp_forms, Guild_CH4_d, estimator = "mlm",
#                                                            keep_mod_descr = c("model_out","model")))
#head(CH4mob_comp_res)

# test data from below -- don't expect this to run!!
# combine with CH4_ac w base -- expect model 21_249 fails, but function runs/passes results
# CH4_ac_comp_f <- combine_sem_model_layers(CH4_mod_fch4.2_s, CH4_ac2)
# CH4_ac_comp_f2 <- CH4_ac_comp_f[18:19,]
# CH4_ac_comp_f2 <- CH4_ac_comp_f[15:25,]
# CH4_ac_comp_f2 <- CH4_ac_comp_f[18,]# one pass, one fail -- if not dropping fail above
#CH4_ac_comp_f2$model_out <- factor(CH4_ac_comp_f2$model_out)  -- needed?


##################################
### a) run & assign SEM hierarchical model table 
# 		-- working function with new try catch (error handling), NOT parallel

sem_table_run_assign2 = function(sem_mod_table, data, estimator = 'mlm', prefix = "sem.",
                                fixed_x = FALSE, orthogonal = TRUE,
                                print_passing = FALSE, print_failed = TRUE){
    
        iters <- nrow(sem_mod_table)                                                  # use iters for loop
        passing <-list()                                                              # good SEM accumulator


        for(i in seq(1:iters)) {                                                      # for each row              
            skip <- FALSE                                                             # assume pass (for try)
            
            # run SEM model            
            curr_model <- sem_mod_table$formula[i]                                    # get ith sem formula
            sem <- tryCatch(sem(curr_model, data=data, fixed.x=fixed_x, estimator = estimator, 
                                orthogonal = orthogonal), #std.lv=TRUE), -- no diff   # TRY sem
                                error = function(e){skip <<- TRUE})                   # if SEM fails, skip = T  
            if(skip==TRUE){next}                                                      # if SEM failed, next iter
                                   
            # assign mod_out names to models (GLOBAL ENV)
            if(skip == FALSE){                                                        # if SEM pass
                sem_name <- paste(sem_mod_table$model_out[i])                         # get model_out
                assign(sem_name, sem, envir=.GlobalEnv)                               # assign out to global env 
                passing[[i]] <- sem_name                                              # accum model_out name
            }
        }
        passing <- unlist(passing)
    
        # Gather models run and (optional) print pass/fail lists to console
        models_run <- sem_mod_table$model_out
        if (print_passing == TRUE) {
            passing_pr <- paste(passing, collapse =', ')
            print(paste0("ran SEM to create models: ", passing_pr))}
        if (print_failed == TRUE){
            failed <- models_run[!models_run %in% passing]
            print(paste0("failed SEM models: ", failed))}
    
        return(passing)
    }   

# demonstrate function
sem_multi_test <- sem_table_run_assign2(model_test2, data=Guild_CH4_d, estimator ="mlm")
sem_multi_test

##################################
### b) run & assign SEM hierarchical model table -- Parallel
#      - fixed assign fault in mcapply, but ugly

# Parallel mc apply version - fix assign to global problem in parallel, ugly not sure how functional try is...
sem_table_run_assign_P = function(sem_mod_table, dataset, estimator = 'mlm', prefix = "sem.", ncores = 20,
                                fixed_x = FALSE, orthogonal = TRUE,
                                print_passing = FALSE, print_failed = TRUE){
    
        iters <- 1:nrow(sem_mod_table)                                                # use iters for loop
        passing <- mclapply(iters, function(i){                                       # multicore apply 

            skip <- FALSE                                                             # assume pass (for try)
            
            # run SEM model, here try catch in case fails            
            curr_model <- sem_mod_table$formula[i]                                    # get ith sem formula
            sem <- tryCatch(sem(curr_model, data=dataset, fixed.x=fixed_x, estimator = estimator, 
                                orthogonal = orthogonal), #std.lv=TRUE), -- no diff   # TRY sem
                                error = function(e){skip <<- TRUE},#)                 # if SEM fails, skip = T  
            if(skip==TRUE){
                next})                                                                # if SEM failed, next iter
                                   
            # assign mod_out names to global env, failed here, just get SEM names
            if(skip == FALSE){                                                        # if SEM pass
                sem_name <- paste(sem_mod_table$model_out[i])                         # get model_out
                c(sem, sem_name)
            }
        }, mc.cores = ncores)                                                         # get n cores for parallel 
   
        # try catch produces nulls & won't assign. Extract list elements, then unlist
        sems <- list()
        sem_names <- list()
        
        for (i in 1:length(passing)){                         
            sems[i] <- (passing[i][[1]][[1]])
            sem_names[i] <- (passing[i][[1]][[2]])
        }
        
        # unlist new lists, no more nulls 
        sems <- unlist(sems)             
        sem_names <- unlist(sem_names)
    
        # assign SEMs to global env by name
        for (i in 1:length(sem_names)){                
                assign(sem_names[[i]], sems[[i]], envir=.GlobalEnv)
            }
            
        # Gather models run and (optional) print pass/fail lists to console
        models_run <- sem_mod_table$model_out
        if (print_passing == TRUE) {
            passing_pr <- paste(sem_names, collapse =', ')
            print(paste0("ran SEM to create models: ", passing_pr))}
    
        if (print_failed == TRUE){
            failed <- models_run[!models_run %in% sem_names]
            print(paste0("failed SEM models: ", failed))}
            else{}

        return(sem_names)
    }   

# test run and assign models from model table
sem_multi_test <- sem_table_run_assign_P(model_test2, data=Guild_CH4_d, estimator ="mlm", ncores = 20)
sem_multi_test; 
compare_sem_results(sem_multi_test, sem_fit_params)

# larger test, should be valid inline
sem_multi_test2 <- sem_table_run_assign_P(CH4mob_comp_forms, data=Guild_CH4_0, estimator ="mlm", ncores=20)
#                                        fixed_x = TRUE, orthogonal = TRUE)
# sem_multi_test2

# show test model summary, not helper output
table <- compare_sem_results(sem_multi_test2, sem_fit_params)
head(table[order(-table$pvalue),])

#invisible(sem_table_run_assign_P(CH4mob_comp_forms, data=Guild_CH4_0, estimator ="mlm", ncores=20, print_failed == FALSE))


# Further testing, DEPRECATE
# test on larger dataset -- same as above?
#sem_multi_test2 <- sem_table_run_assign_P(CH4mob_comp_forms, data=Guild_CH4_d, estimator ="mlm", 
#                                         fixed_x = TRUE, orthogonal = TRUE)
# sem_multi_test2 -- should be list of models

# show results with model scoring (from assigned)
# table <- compare_sem_results(sem_multi_test2, sem_fit_params)
# head(table[order(-table$pvalue),])# test on larger dataset, is valid for examples or from below?  
# sem_multi_test2 <- sem_table_run_assign_P(CH4_ac_comp_f, data=Guild_CH4_0, estimator ="mlm", 
#                                         fixed_x = TRUE, orthogonal = TRUE)
# sem_multi_test2 -- should be list of models

# show results with model scoring (from assigned)
# table <- compare_sem_results(sem_multi_test2, sem_fit_params)
# head(table[order(-table$pvalue),])

##################################
### c) run and compare SEM models in table (Updated)
#      - implements parallel sem_table_run_assign_P, passes n_cores to options

run_compare_sem_models = function(mod_table, data, estimator = "mlm", ncores = 20,
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

# run models, from manual example -- why doesn't this work with 0 data (all???) it did at some pont...
# CH4mob_comp_res <- run_compare_sem_models(CH4mob_comp_forms, Guild_CH4_0, ncores = 20)  
                            #, estimator = "mlm", keep_mod_descr = c("model_out","model")) # defaults
#CH4mob_comp_res

# table <- CH4mob_comp_res
# table[order(-table$pvalue),]
# head(CH4mob_comp_res)
# get summary, full set of models -- MOVE DOWN to IV
#CH4mob_comp_t <- clean_SEM_model_table_for_output(CH4mob_comp_forms, CH4mob_comp_r, 
#                                drop_cols = c("model_out", "formula"), merge_col = c("model"))
#
#CH4mob_comp_t[order(-CH4mob_comp_t$pvalue),]# run models from later composite examples -- below
#CH4_ac_comp_r <- run_compare_sem_models(CH4_ac_comp_f, Guild_CH4_0, ncores = 20)
                            #, estimator = "mlm", keep_mod_descr = c("model_out","model")) # defaults
#CH4_ac_comp_r

# Parllel report: for this ex, run 39 models in CH4_ac_comp_f, data = Guild_CH4_0:
# non par: ~ 50-60% CPU load, 700% CPU), ca 70 sec
# parallel: ~ 90-100% CPU load, ca. 29s 
# get summary, full set of models -- MOVE DOWN to IV
#CH4_ac_comp_t <- clean_SEM_model_table_for_output(CH4_ac_comp_f, CH4_ac_comp_r, 
#                                drop_cols = c("model_out", "formula"), merge_col = c("model"))
#CH4_ac_comp_t[order(-CH4_ac_comp_t$pvalue),]# adding SEM response vars with string defs


######################################################################################################
### IV) Create clean SEM model summary tables
#      - cleans sem model summaries (results and formula) for compact output
#      - intended for complex models (multiple y, R2, diff. workflow than LassSEM models
#      - may wish to create object, wrapper, export list functions...

##### user function is d): 
# clean_SEM_model_table_for_output(CH4mob_comp_forms, CH4mob_comp_res, 
#                                 drop_cols = c("model_out", "formula"), merge_col = c("model"))
                                
# - other helper functions may have utility elsewhere (get formula columns, remove y_vars from formulae)

##################################
## 0) Simple test data

### XXX DELETE XXX -- Deprecated, redundant with above

# CH4 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h'
# CH4_H2 <- 'CH4_H2~ SRB_syn + CO2_mg_m2_h'
# MOB_IIa <- 'MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac'
# MOB_I <- 'MOB_I ~ CH4_H2 +CH4_ac'# Get vector of response vars as strings
# Other <- ""
# response_vars <- c("CH4", "CH4_H2", "MOB_I","MOB_IIa", "Other")head(sem# Test function
# sem_models.2 <- update_sem_model_table(model_test, base_model = "ch4_mod0", new_mod_name = "ch4_mod0a.1", 
#                       notes = "drop CH4~MOB_I",
#                       formula_edits = c("CH4", "CH4_ug_m2_h ~ CH4_ac + CH4_H2 + MOB_IIa + CO2_mg_m2_h"))
#                       # output_prefix = "sem."#non_feat_cols = c("model", "base", "notes", "formula"))

# sem_models.2
# Demonstrate function
sem_models.2_out <- run_compare_sem_models(sem_models.2, data, estimator = "mlm")
# sem_models.2_out


##################################
## a) HELPER-- remove y vars from model table formula (single)

remove_yvar_from_df_formula_column =function(df, col_name, fxn_operator = '~ '){
    
    # start function 
    row_data = paste(df[,col_name])                                              # get vect of rows in column

    form_list = list()                                                           # empty list

    for (i in (seq_along(row_data)))     {                                       # for each column
        form_list[i] <- paste(unlist(strsplit(row_data[i], '~ '))[2])            # split, x_vars only to list
        }

    col_replace <- unlist(form_list)                                             # unlist
    df[,col_name] <- col_replace                                                 # replace column data
    return(df)                                                                   # return df

}

# demonstrate function
remove_yvar_from_df_formula_column(sem_models.2, 'MOB_IIa')
#remove_yvar_from_df_formula_column(sem_models.2, 'formula') #-- just drop this one!

##################################
## b) remove y vars from model table formula (multiple, any)

remove_yvars_from_model_table_formulae = function(model_table){
    
    form_cols <- get_formula_col_list_from_df(model_table)
    
    for (i in form_cols){
        model_table <- remove_yvar_from_df_formula_column(model_table, i)
    }

    model_table
}    

# demonstrate function
model_test2_out <- remove_yvars_from_model_table_formulae(sem_models.2)
head(model_test2_out)

##################################
## b) get model table summary data

clean_SEM_model_table_for_output = function(sem_models, sem_output, 
                                            drop_cols = c("model_out", "formula"), merge_col = "model"){
    
    # drop extra columns in sem models, outputs
    sem_models <- drop_vars_from_df(sem_models, drop_cols)
    sem_output <- drop_vars_from_df(sem_output, drop_cols)
    
    #model_cols <- names(sem_models)
    #keep_cols <- model_cols[!model_cols %in% drop_cols]
    #sem_models <- sem_models[keep_cols]
    
    # remove y vars from formulas
    clean_models <- remove_yvars_from_model_table_formulae(sem_models)  # need to drop formula, else
    
    # merge with results
    models_out <- merge(sem_output, clean_models, by = merge_col)       # merge
    models_out[models_out=='NA'] <-''                                   # drop mysterious NAs
    return(models_out)    
    
    # models_out$formula <- sem_models$form
    # repair column names -- probably need a for loop here TODO
    names(models_out) <- gsub(x = names(models_out), pattern = ".x", "R2")
    names(models_out) <- gsub(x = names(models_out), pattern = ".y", "f")
    #return(models_out)  
    
}

# demonstrate function
clean_SEM_model_table_for_output(sem_models.2, sem_models.2_out, 
                                drop_cols = c("model_out", "formula"), merge_col = c("model"))#, 'base', 'notes'))
# demonstrate function - ex post
# CH4mob_comp_table <- clean_SEM_model_table_for_output(CH4mob_comp_forms, CH4mob_comp_res, 
#                                drop_cols = c("model_out", "formula"), merge_col = c("model"))

# table <- CH4mob_comp_table 
# head(table[order(-table$pvalue),])

# head(CH4mob_comp_table)CH4mob_comp_table_sig <- CH4mob_comp_table[CH4mob_comp_table$pvalue > 0.05,]
# CH4mob_comp_table_sig
# get summary, full set of models -- MOVE DOWN to IV
# CH4_ac_comp_t <- clean_SEM_model_table_for_output(CH4_ac_comp_f, CH4_ac_comp_r, 
#                                drop_cols = c("model_out", "formula"), merge_col = c("model"))

#CH4_ac_comp_t[order(-CH4_ac_comp_t$pvalue),]

