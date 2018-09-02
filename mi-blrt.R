rm(list = ls())

library(plyr)
library(tidyverse)
library(lpa.mi.src)
library(doParallel)
library(foreach)
library(doRNG)

Processors = 20
Replications = 50
save_it = TRUE
KKmax = 4;
starts_txt = "10 4;";
M = 10

dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
temp_dir = "C:/Users/marcu/Documents/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)
#dropbox_wd =  "E:/Dropbox/Dropbox"  #  Dropbox folder on local machine
#temp_dir = "X:/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)        
#dropbox_wd = "C:/Users/mwaldman1/Dropbox"
#temp_dir = "C:/Users/mwaldman1/Documents/lpa-mi-pool"
#dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine
#temp_dir = "D:/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)        

# Fully-crossed simulation conditions
pctmiss_vec = c(0.35)     #  Percent of observations with *at least* one missing value (indicator or missing data correlate) -> if other specification needed one must modify list_inds.miss generation in get_obs_data.R
N_vec = c(2E3)                #  Sample sizes within each replication
J_Y_vec = c(3)                #  Number of latent class indicators
J_Xcom_vec = c(1)             #  Number of complete data missing data correlates
J_Xinc_vec = c(0)             #  Number of missing data correlates which themselves contain missing data
MD_vec = c(2)               #  WARNING! MD = 1 RESULTS IN ERROR! Mahalanobis distances between classes
pi_list = list(               #  Number of classes, K (maximum of 3)
  #A = c(0.425,0.425,0.15), 
  B = rep(1/3,3) 
  #C = rep(1/4,4)
)
rho_YX_vec = c(0.4)          #  Strength of missing data correlates
C_modifies_YX_vec = c(TRUE)  # Whether or not the Y-X relationships are different between latent classes 

# Define the imputation alogrithms to compare
methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
methods_list$procedure[[1]] = "stratamelia";  methods_list$name[[1]] = "mg-Amelia";     methods_list$args[[1]] = list(m = M, p2s = 0, tolerance = 1E-4)
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
#setwd(main_wd); save.image(file = "image lpa-mi-pool.RData")

# create a list of packages loaded (needed for the "foreach" function)
temp = strsplit(subset(search(), startsWith(search(),"package")),"package:", fixed = TRUE)
packages_vec = rep(NA,length(temp))
for (t in 1:length(temp)){
  packages_vec[t] = temp[[t]][2]
}
rm(temp)


#### Conduct parallel processing
# Register processors
#cl<-makePSOCKcluster(Processors)
#registerDoParallel(cl)

p = 1; rep = 1; z = 1;
#out_foreach<- foreach(p = 1:Processors,
#                      .packages = packages_vec, .options.RNG=42) %dorng%{
    dff_complete <- dff_imputed <- NULL;
#    for(rep in 1:Replications){
      results_df = NULL
#        for (z in 1:Ztot){
          
#          sink(file = paste0(temp_wd_p_vec[p],"/sink-p", p, ".txt"))
          
          file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), full.names = TRUE))
          file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), full.names = TRUE))
          file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), full.names = TRUE))
          
          
          pop_params = get_FMM_params(z = z, data_conditions = data_conditions)

          #### Construct Data Files ####
            list_complete = lpa.mi.src::get_complete_data(z = z, 
                                                          data_conditions = data_conditions,
                                                          rep = rep, 
                                                          p = p, 
                                                          save_it = TRUE, 
                                                          temp_wd_p_vec = temp_wd_p_vec)

            list_observed = lpa.mi.src::get_obs_data(z=z, 
                                                     df =  list_complete$dfcom,
                                                     pctmiss_vec = pctmiss_vec,
                                                     data_conditions = data_conditions,
                                                     save_it = FALSE)

            list_imputed = lpa.mi.src::get_imputed_data(z = z, 
                                                       list_get_obs = list_observed,
                                                       list_get_complete = list_complete,
                                                       methods_list = methods_list,
                                                       data_conditions = data_conditions,
                                                       pctmiss_vec = pctmiss_vec, 
                                                       save_it = TRUE,
                                                       p = p, 
                                                       rep = rep, 
                                                       temp_wd_p_vec = temp_wd_p_vec)
            
            dff_complete = list_complete$dffolderfiles
            dff_imputed = subset(list_imputed$dffolderfiles, is.na(m))
            
            
            LL_df = data.frame(k = 1:KKmax, params = rep(NA, KKmax), free = rep(NA,KKmax), fixed = rep(NA,KKmax))
            for(kk in 1:KKmax){
            
                # Delete all inp files
                file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".inp", full.names = TRUE))
                file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".out", full.names = TRUE))
                file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".inp", full.names = TRUE))
                file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".out", full.names = TRUE))
                file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".inp", full.names = TRUE))
                file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".out", full.names = TRUE))
                
                # STarting values for Mplus
                if(kk == pop_params$K_z){
                  pop_params_kk = pop_params;
                } else {
                  pop_params_kk = pop_params;
                  pop_params_kk$K_z = kk
                  pop_params_kk$mu_z = NULL;
                  pop_params_kk$S_z = NULL
                  pop_params_kk$pi_z = NULL;
                }
                
                inp_complete = lpa.mi.src::create_naiveMplus_inpfile(z = z, 
                                                                     out_get_FMM = pop_params_kk, 
                                                                     dffolderfiles = dff_complete, 
                                                                     temp_wd_p = temp_wd_p_vec[p], 
                                                                     starts_txt = starts_txt)
                
                inp_imputed = lpa.mi.src::create_naiveMplus_inpfile(z = z, 
                                                                    out_get_FMM = pop_params_kk, 
                                                                    dffolderfiles = dff_imputed, 
                                                                    temp_wd_p = temp_wd_p_vec[p], 
                                                                    starts_txt = starts_txt, 
                                                                    type_imputation = TRUE)
                
                runModels(target = paste0(temp_wd_p_vec[p], "/Imputed data"))
                out_imputed <- readModels(target = paste0(temp_wd_p_vec[p], "/Imputed data"), 
                                          what = c("warn_err", "summaries","parameters"), 
                                          quiet = TRUE)
                # Extract the loglikelihoods for each imputed data set
                print(kk)
                print(out_imputed$summaries$Parameters)
                LL_df$params[kk] = out_imputed$summaries$Parameters
                LL_df$free[kk] = out_imputed$summaries$LL_Mean
                
                
                svals_txt<-extract_svals(file = list.files(path = paste0(temp_wd_p_vec[p], "/Imputed data"), pattern = ".out"), 
                              path =  paste0(temp_wd_p_vec[p], "/Imputed data"))
                model_txt = c("MODEL:\n",gsub("*","@",svals_txt, fixed = TRUE))
                
                # Difference between below and above is that this fixes the value at the pooled estimate
                inp_imputed = lpa.mi.src::create_naiveMplus_inpfile(z = z, 
                                                                    out_get_FMM = pop_params_kk, 
                                                                    dffolderfiles = dff_imputed, 
                                                                    temp_wd_p = temp_wd_p_vec[p], 
                                                                    type_imputation = TRUE, 
                                                                    Model_txt = model_txt)
                runModels(target = paste0(temp_wd_p_vec[p], "/Imputed data"))
                out_imputed <- readModels(target = paste0(temp_wd_p_vec[p], "/Imputed data"), 
                                          what = c("warn_err", "summaries","parameters"), 
                                          quiet = TRUE)
                LL_df$fixed[kk]  = out_imputed$summaries$LL_Mean
              
            }
            LRT_df = data.frame(diff(as.matrix(LL_df[,-1]))); 
            names(LRT_df) = c("df","LL_free","LL_fixed")
            LRT_df = transform(LRT_df, ARIV = (2*(LL_free-LL_fixed)*(M+1))/(df*(M-1)) )
            LRT_df = transform(LRT_df, D = (2*LL_free)/(df*(1+ARIV)))
            LRT_df = transform(LRT_df, Test = paste0(as.character(seq(1,KKmax-1)), " class vs. ", as.character(seq(2,KKmax)), " class"))
            
          
              
              
#          sink()
#        } #end for z = 1:Ztot
      
      #Save the progress

      # Save the result
#    } #end for rep = 1:Replications
#} #for each p = 1,...,P

#stopImplicitCluster()
#save.image("environment.RData")
