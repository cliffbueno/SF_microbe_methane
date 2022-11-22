# SEM - Lasso semi-auto for SF OTU guilds v1.12
#-  runs under r > 3.6 (glmnet for LASSO)

# Assume have loaded modules 0-3

########################################################################################
### 4) Simple SEM model testing -- example workflow


## a) Define simple model
ch4_mod0 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h 
             
             CH4_H2 ~ SRB_syn + CO2_mg_m2_h
             MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac
             MOB_I ~ CH4_H2 +CH4_ac'
# Instead of full function, add SEM response vars with string defs
CH4 <- 'CH4_ug_m2_h ~  CH4_ac + CH4_H2 + MOB_IIa + MOB_I + CO2_mg_m2_h'
CH4_H2 <- 'CH4_H2~ SRB_syn + CO2_mg_m2_h'
MOB_IIa <- 'MOB_IIa ~ NOB + AOA + AOB +CH4_H2 +CH4_ac'
MOB_I <- 'MOB_I ~ CH4_H2 +CH4_ac'

# Get vector of response vars as strings
Other <- ""
response_vars <- c("CH4", "CH4_H2", "MOB_I","MOB_IIa", "Other")

# Build model
sem_ch4_0_base_table <- define_new_sem_table(mod_name = "ch4_mod0", response_vars)
sem_ch4_0_base_table


## b) evaluate base model
ch4_mod0 <- sem_ch4_0_base_table$formula

# run simple model
ch4_mod0_f <- sem(ch4_mod0, data=Guild_CH4_d, fixed.x=FALSE, estimator ="mlm")

#compare_sem_results(ch4_mod0_f, sem_fit_params)

# Get summary of fit
sem_fit_sum(ch4_mod0_f, sem_fit_params)

# Get R2: 
get_SEM_R2s(ch4_mod0_f)

# Get non-significant vars
get_SEM_nonsig_vars(ch4_mod0_f, p_cut =0.05)