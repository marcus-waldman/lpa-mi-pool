rm(list = ls())

df2results<-function(df, pop_params, dtype = "", p = NA, rep = NA, z = NA, pm = NA, pva = NA){
  
  J_Y_z = pop_params$J_Y_z
  K_z = pop_params$K_z
  
  intercept = paste("Y", 1:J_Y_z, "~1", sep = "")
  vcovar = mat.or.vec(nr = J_Y_z, nc = J_Y_z)
  for (i in 1:J_Y_z){
    for (j in 1:i){
      vcovar[i,j] = paste("Y",j,"~~","Y",i, sep = "")
    }
  }
  covariances = vcovar[lower.tri(vcovar, diag = TRUE)]
  
  params = as.character(c(intercept,covariances))
  rep_df = data.frame(class = sort(rep(1:K_z, length(params))))
  rep_df = transform(rep_df, parameter = rep(params, K_z))
  rep_df$parameter = as.character(rep_df$parameter)
  
  rep_df = transform(rep_df, estimate = NA, se = NA, value = NA)
  for(k in unique(rep_df$class)){
    
    inds_k = which(df$class ==k)
    rep_df$estimate[which(endsWith(rep_df$parameter,"~1") & rep_df$class==k)] = apply(df[inds_k,1:J_Y_z],2,mean, na.rm = TRUE)
    
    vc_k = cov(df[inds_k, 1:J_Y_z], use = "pairwise.complete.obs")
    rep_df$estimate[which(!endsWith(rep_df$parameter,"~1") & rep_df$class==k)] = vc_k[lower.tri(vc_k, diag = TRUE)]

    
    mu_z = pop_params$mu_z
    rep_df$value[which(endsWith(rep_df$parameter,"~1") & rep_df$class==k)] = pop_params$mu_z[1:J_Y_z,k]
    
    S_z = pop_params$S_z[1:J_Y_z,1:J_Y_z,k]
    rep_df$value[which(!endsWith(rep_df$parameter,"~1") & rep_df$class==k)] = S_z[lower.tri(S_z, diag = TRUE)]
  }
  
  rep_df = transform(rep_df, data_type = dtype, p = p, rep = rep, z = z, pm = pm, pva = pva)
  return(rep_df)  
  
}

lavaan2results<-function(obj_lavaan, pop_params, dtype = "", 
                         p = NA, rep = NA, z = NA,  pm = NA, pva = NA){
  J_Y_z = pop_params$J_Y_z
  summ_lavaan = summary(obj_lavaan);
  est_lavaan = summ_lavaan$est; se_lavaan = summ_lavaan$se; 
  lhs_lavaan = summ_lavaan$lhs; op_lavaan = summ_lavaan$op; 
  rhs_lavaan = summ_lavaan$rhs; gr_lavaan = summ_lavaan$group;
  rep_df = data.frame(
    class = summ_lavaan$group,
    parameter = paste(summ_lavaan$lhs, summ_lavaan$op, summ_lavaan$rhs, sep = ""), 
    estimate = est_lavaan,
    se = se_lavaan 
  )
  rep_df = transform(rep_df, value = NA)
  for(k in unique(rep_df$class)){
    
    mu_z = pop_params$mu_z
    rep_df$value[which(op_lavaan=="~1" & rep_df$class==k)] = pop_params$mu_z[1:J_Y_z,k]
    
    S_z = pop_params$S_z[1:J_Y_z,1:J_Y_z,k]
    var_k = diag(S_z)
    cov_k = as.vector(S_z[lower.tri(S_z)])
    rep_df$value[which(op_lavaan=="~~" & rep_df$class==k)] = c(var_k,cov_k)
  }
  rep_df = transform(rep_df, data_type = dtype, p = p, rep = rep, z = z, pm = pm, pva = pva)
  return(rep_df)
}

model_lavaan = '
  # Means
  Y1~1 
  Y2~1
  Y3~1
  
  # Variances
  Y1~~Y1
  Y2~~Y2
  Y3~~Y3

  # Covariance
  Y1~~Y2
  Y1~~Y3
  Y2~~Y3
'

library(plyr)
library(tidyverse)
library(lpa.mi.src)
library(semTools)
library(doParallel)
library(foreach)


Processors = 5
Replications = 400
save_it = FALSE

#dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
#temp_dir = "C:/Users/marcu/Documents/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)
dropbox_wd =  "E:/Dropbox/Dropbox"  #  Dropbox folder on local machine
#temp_dir = "X:/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)        
#dropbox_wd = "C:/Users/mwaldman1/Dropbox"
#temp_dir = "C:/Users/mwaldman1/Documents/lpa-mi-pool"
#dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine
#temp_dir = "D:/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)        

# Fully-crossed simulation conditions
pctmiss_vec = c(0.5)     #  Percent of observations with *at least* one missing value (indicator or missing data correlate) -> if other specification needed one must modify list_inds.miss generation in get_obs_data.R
N_vec = c(5E2,2E3)                #  Sample sizes within each replication
J_Y_vec = c(3)                #  Number of latent class indicators
J_Xcom_vec = c(1)             #  Number of complete data missing data correlates
J_Xinc_vec = c(0)             #  Number of missing data correlates which themselves contain missing data
MD_vec = c(3)               #  WARNING! MD = 1 RESULTS IN ERROR! Mahalanobis distances between classes
pi_list = list(               #  Number of classes, K (maximum of 3)
  A = rep(0.5,2) , 
  B = rep(1/3,3) 
  #C = rep(1/4,4)
)
rho_YX_vec = c(0.4)          #  Strength of missing data correlates
C_modifies_YX_vec = c(TRUE)  # Whether or not the Y-X relationships are different between latent classes 

# Define the imputation alogrithms to compare
methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
methods_list$procedure[[1]] = "stratamelia";  methods_list$name[[1]] = "mg-Amelia";     methods_list$args[[1]] = list(m = 100, p2s = 1, tolerance = 1E-4)
#methods_list$procedure[[2]] = "mice";         methods_list$name[[2]] = "mg-mice";       methods_list$args[[2]] = list(method = "bygroup", imputationFunction = "norm", m = 20, maxit = 50)
#methods_list$procedure[[3]] = "mice";         methods_list$name[[3]] = "mg-mvn";        methods_list$args[[3]] = list(method = "bygroup", imputationFunction = "norm", m = 20, maxit = 50, blocks = NULL)


#### Pre-processing ####

#Sort and check pctmiss is between 0 and 1
pctmiss_vec = sort(pctmiss_vec)
if( (sum(pctmiss_vec<=0) + sum(pctmiss_vec>=1))>0){stop("elements in pctmiss_vec must be between 0 and 1 (exclusive).")}

# Define working directory
main_wd = paste(dropbox_wd, "/Dissertation/lpa-mi-pool", sep = "")
temp_fname = paste("sim-lpa-mi-pool", format(Sys.time(), "%b %d %Y %H %M"), sep = " ")

# Create temporary folders if saved
if (save_it){
  
  # Check that the dropbox and temp_dir are found on machine
  if (dir.exists(dropbox_wd)*dir.exists(temp_dir)==0){
    stop("Provided dropbox_wd or temp_dir paths not found." )
  }
  
  # Create temporary folder
  temp_wd =  paste(temp_dir, temp_fname, sep = "/")
  dir.create(temp_wd)

  # Create processor-specific temporary folder
  setwd(temp_wd)
  for (p in 1:Processors){
    dir.create(paste("p",p, sep = ""))
  }
  temp_wd_p_vec = dir(path = temp_wd, full.names = TRUE)
  
  # Within the processor-specific temporary folder, create subfolders
  for (p in 1:Processors){
    setwd(temp_wd_p_vec[p])
    dir.create("Observed data") 
    dir.create("Complete data")
    dir.create("Imputed data")
  }
  
}

# Create folder to store results in main directory
setwd(paste(dropbox_wd, "/Dissertation/results", sep = ""))
dir.create(temp_fname)
setwd(temp_fname)
results_wd = getwd()

# Identify all conditions to be evaluated across reps = 1:Replications for each method in method_list 
data_conditions = data.frame(expand.grid(N = N_vec, J_Y = J_Y_vec, J_Xcom = J_Xcom_vec, 
                                         J_Xinc = J_Xinc_vec, MD = MD_vec, class_size = names(pi_list), rho_YX = rho_YX_vec, 
                                         C_modifies_YX = C_modifies_YX_vec))

Kmax = max(simplify2array(lapply(pi_list, FUN = length)))

# Add in the number of classes corresponding to each condition and the marginal probabilities (pi_k)
pi_df = data.frame(mat.or.vec(nr = length(pi_list), nc = Kmax + 2))+NA
names(pi_df) = c("class_size","K", paste("pi_", 1:Kmax, sep = ""))
for(p in 1:nrow(pi_df)){
  pi_df$class_size[p] = names(pi_list)[p]
  pi_df$K[p] = length(pi_list[[p]])
  pi_df[p,seq(3,pi_df$K[p]+2)] = pi_list[[p]]
}
data_conditions = merge(x = data_conditions, y = pi_df)
rm(pi_df)


# Ensure that there is at least one data condition
Ztot = nrow(data_conditions)
if (Ztot==0){stop("Total number supported of data conditions is zero.")}


# Initiate a list of folder files
dffolderfiles = NULL

# # writeout the data
# setwd(main_wd); save.image(file = "image lpa-mi-pool.RData")

# create a list of packages loaded (needed for the "foreach" function)
temp = strsplit(subset(search(), startsWith(search(),"package")),"package:", fixed = TRUE)
packages_vec = rep(NA,length(temp))
for (t in 1:length(temp)){
  packages_vec[t] = temp[[t]][2]
}
rm(temp)

#p = 1
#### Conduct parallel processing
# Register processors
cl<-makePSOCKcluster(Processors)
registerDoParallel(cl)



foreach (p = 1:Processors, .packages = packages_vec) %dopar%{
    
    results_df = NULL
    
    for(z in 1:Ztot){
        
        pop_params = get_FMM_params(z = z, data_conditions = data_conditions)
    
        for (rep in 1:Replications){
          
            # Generate data
            
           
            list_complete = lpa.mi.src::get_complete_data(z = z, data_conditions = data_conditions,
                                                          save_it = save_it) 
    
            
           
            list_observed = lpa.mi.src::get_obs_data(z=z, df =  list_complete$dfcom,
                                                     pctmiss_vec = pctmiss_vec,
                                                     data_conditions = data_conditions,
                                                     save_it = save_it)
   
          

            list_imputed = lpa.mi.src::get_imputed_data(z = z, list_get_obs = list_observed,
                                                       list_get_complete = list_complete,
                                                       methods_list = methods_list,
                                                       data_conditions = data_conditions,
                                                       pctmiss_vec = pctmiss_vec, save_it = save_it, 
                                                       p = p, rep = rep)

            
            # Evaluate complete data model
            temp1_df = df2results(df = list_complete$dfcom, dtype = "complete",p = p, rep = rep, z = z)
            
            temp2_df <- temp3_df<-NULL
            
            for(pm in seq(1,length(pctmiss_vec))){
              
              # Evaluate observed data model
              temp2_pm = df2results(df = list_observed$list_obsdf[[pm]], dtype = "observed",p = p, rep = rep, z = z, pm = pm)
              temp2_df = rbind(temp2_df, temp2_pm)
              
              for(pva in seq(1, length(methods_list$procedure))){
                # Evaluate imputation model
                mids_obj = list_imputed$obj_call[[pm]][[pva]]
                mice.imp2<-lapply(seq(mids_obj$m),function(im) complete(mids_obj,im))
                obj_lavaan = lavaan.mi(model = model_lavaan, data = mice.imp2, group = "strata", 
                                       information = "observed")
                temp3_pmpva = lavaan2results(obj_lavaan = obj_lavaan, pop_params = pop_params, dtype = "imputed", p = p, rep = rep, z = z, pm = pm, pva = pva)
                temp3_df = rbind(temp3_df, temp3_pmpva)
              } #end for pva = 
            }# end for pm = 
            temp_df = rbind(temp1_df, temp2_df, temp3_df)
            results_df = rbind(results_df, temp_df)
        } #end for rep = 1:Replications
    } #end for z = 1:Ztot
    setwd(results_wd)
    save(results_df, file = paste("p", p, ".RData", sep = "") )

} #for each p = 1,...,P

stopImplicitCluster()
save.image("environment.RData")
