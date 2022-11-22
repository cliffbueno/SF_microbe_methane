#  Lasso semi-auto for SF OTU guilds v1.12
#  runs under r > 3.6 (glmnet for LASSO)

# Assume have loaded modules 0-6
# likely pollulted with examples, possibly cruft


########################################################################################
### 12) lasSEM branching model permutation workflow wrappers


##################################
## a) stray helper functions 

# i) drop from list (for variable input list filtering) ## 

drop_from_list = function(drop_list, parent_list){
    
    clean_list <- parent_list[!parent_list %in% drop_list]
    clean_list
}


## ii) plot lasso object list (should be MOVED to lasso plot wrapper module 10)) ##
#     - coef plot of lasso objects, input as str for obj name

plot_lasso_obj_list = function(lasso_list){     
    
    plot_list <- list()

    for (i in lasso_list){ 
        lasso <- get(i)
        plot_list[[i]] <- lasso$coef_plot

    }

    plot_grid(plotlist = plot_list, nrow=1)     
    
}


##################################
# b) branching workflow helper functions

## i) rename y var to column name (helper) ##

rename_formula_col_as_yvar = function(table, formula_col = 'form', fxn_operator = '~'){

    form1 <- paste(table[formula_col][[1]])                                                   # first formula 
    y_var <- unlist(strsplit(form1, fxn_operator))[1]                                         # get y var
    colnames(table)[colnames(table) == formula_col] <- y_var
    table
    #y_var #<- paste0(y_var, fxn_operator)                                                     # add op back
}

# Demonstrate function  -- will need new example, works but data below
# results_table <- CH4_mod_fch4.2_sem$filt_models
# rename_formula_col_as_yvar(results_table, formula_col = 'form', fxn_operator = '~')


## ii) prepare model input data ## 
#   -- Get model formula and apply pre-fix from lassSEM

get_sem_models_add_prefix = function(lasso_object, models = "filt_models", prefix = "",
                                  filter = '', cutoff = 0){
    formulae <- lasso_object[[models]]                       # get data from object, could be "sem_models"
    formulae$model <- paste0(prefix, formulae$model)         # append model/ dataset prefix
    
    # filter data as desired (R2, p, etc)
    if (filter != ''){
        formulae <- formulae[formulae[filter] > cutoff,]         
    }
    return(formulae)
}


# Example use  - is from downstream, will need for upstream examples
# CH4_mod_fch4.2 <- get_sem_models_add_prefix(CH4_mod_fch4.2_sem, models = "filt_models", prefix = "fch4.2.")
# CH4_mod_fch4.2_A <- get_sem_models_add_prefix(CH4_mod_fch4.2_sem, models = "sem_models", prefix = "fch4.2.",
#                                          filter = 'R2', cutoff = 0.5)
# CH4_mod_fch4.2
# CH4_mod_fch4.2_sem





# run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.001, full_formula = TRUE)
# run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.001, output = 'sem_results')
# run_sem_model_summary(CH4_ac_comp_f, Guild_CH4_0)   #  -- ex won't work here

# simple function produces data to run using run_sem_model_summary (above)
# more complex versions could be IO wrappers, run lists or input matrices,
# OR could run at each branch and filter results


##################################
### d) sem_tree_formula_build
#   - makes formulae for sem trees, adding n branches to base (from lasSEM object results, composites)

sem_tree_formula_build = function(tree_base, branches, run = TRUE, warn_n_models = TRUE) {   # branches should be list of SEM_objects
    
    # add branch formulae to base model
    sem_formula_tree <- combine_sem_model_layers(tree_base, branches[[1]])   
    
    # iterate other branches
    if (length(branches) > 1) {
        for (i in 2:length(branches)){
           sem_formula_tree <- combine_sem_model_layers(sem_formula_tree, branches[[i]])
       }
     }
     #readline(prompt="Press [enter] to continue")
    n_models <- dim(sem_formula_tree)[1]
    print(paste(n_models, "sem models"))
    
    
    sem_formula_tree
     #print(dim(sem_formula_tree[1]))
        
}


# Examples below from later in dataset testing, need replacement with generics...last sequential version v.1.09# tree_base <- CH4_mod_fch4.2_s
# branch test data, from lasSEM object to filtered models
# branch_1 <- MOBIIa_mod_fAOANOB_sem$filt_models     # was MOBIIa_fng
# branch_2 <- CH4ac_mod.2_sem$filt_models            # CH4_ac2   
#    names(branch_2)[6] <- "CH4_ac"                   #branch_2

# branchz <- list(branch_1)                          # branchz[[1]]
# branchz2 <- list(branch_1, branch_2)
# demonstrate function -- note prints n models in case too many         # 360 models, peak ram = 112 GB!
#sem_formula_tree <-sem_tree_formula_build(tree_base, branchz) # , run = TRUE, warn_n_models = TRUE)
# sem_formula_tree <-sem_tree_formula_build(tree_base, branchz2) # , run = TRUE, warn_n_models = TRUE)
# sem_formula_tree

# test primative fxn
#sem_formula_tree <- combine_sem_model_layers(tree_base, branch_1) 
# sem_formula_tree <- combine_sem_model_layers(tree_base, branchz[[1]])   

# run using SCALED data for some reason or crash, e.g. Guild_CH4_aS
# run it to summary
#run_sem_model_summary(sem_formula_tree, Guild_CH4_0)
# sem_run_test <- run_sem_model_summary(sem_formula_tree, Guild_CH4_0, sort_by = 'pvalue', sort_cutoff = 0.00) # pass parallel!!
#run_sem_model_summary(CH4mob_comp_f, Guild_CH4_0)





