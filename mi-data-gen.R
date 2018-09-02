rm(list = ls())

library(plyr)
library(tidyverse)
library(lpa.mi.src)
library(doParallel)
library(foreach)
library(doRNG)
library(semTools)


Processors = detectCores()
Replications = 64
KKmax = 4;



# Acer-laptop
  #dropbox_wd =  "C:/Users/marcu/Dropbox"             #  Dropbox folder on local machine
  #temp_dir = "C:/Users/marcu/Documents/lpa-mi-pool"  #  Temporary directory for local maching (will automaticall add new temporary folder)

# Queen Mary's Revenge
  dropbox_wd =  "D:/Dropbox"            #  Dropbox folder on local machine
  temp_dir = paste0(LETTERS[seq(26-Processors+1,26)],":") #  Temporary directory for local maching (will automaticall add new temporary folder)        

# Xi-GSU
  #dropbox_wd = "C:/Users/mwaldman1/Dropbox"
  #temp_dir = "C:/Users/mwaldman1/Documents/lpa-mi-pool"

# Big Bertha
  #dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine
  #temp_dir = "D:/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)        


# Nuke directory
  lpa.mi.src::nuke_dirs(temp_dir)
    
# Fully-crossed simulation conditions
pctmiss_vec = c(0.35)     #  Percent of observations with *at least* one missing value (indicator or missing data correlate) -> if other specification needed one must modify list_inds.miss generation in get_obs_data.R
N_vec = c(500,2E3)                #  Sample sizes within each replication
J_Y_vec = c(3)                #  Number of latent class indicators
J_Xcom_vec = c(1)             #  Number of complete data missing data correlates
J_Xinc_vec = c(0)             #  Number of missing data correlates which themselves contain missing data
MD_vec = c(2,2.5,3,3.5)       #  WARNING! MD = 1 RESULTS IN ERROR! Mahalanobis distances between classes
pi_list = list(               #  Number of classes, K (maximum of 3)
  A = c(0.5,0.5),
  B = c(0.8, 0.2) , 
  C = c(0.4,0.4,0.2), 
  D = rep(1/3,3)
)
rho_YX_vec = c(0.4)          #  Strength of missing data correlates
C_modifies_YX_vec = c(TRUE)  # Whether or not the Y-X relationships are different between latent classes 

# Define the imputation alogrithms to compare
methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
methods_list$procedure[[1]] = "stratamelia";  methods_list$name[[1]] = "mg-Amelia";     methods_list$args[[1]] = list(m = 1, p2s = 0, tolerance = 1E-4)
#methods_list$procedure[[2]] = "mice";         methods_list$name[[2]] = "mg-mice";       methods_list$args[[2]] = list(method = "bygroup", imputationFunction = "norm", m = 5, maxit = 20)
#methods_list$procedure[[3]] = "mice";         methods_list$name[[3]] = "mg-mvn";        methods_list$args[[3]] = list(method = "bygroup", imputationFunction = "norm", m = 20, maxit = 50, blocks = NULL)


#### Pre-processing ####

# Identify all conditions to be evaluated across reps = 1:Replications for each method in method_list 
data_conditions = data.frame(expand.grid(N = N_vec, J_Y = J_Y_vec, J_Xcom = J_Xcom_vec, 
                                         J_Xinc = J_Xinc_vec, MD = MD_vec, class_size = names(pi_list), rho_YX = rho_YX_vec, 
                                         C_modifies_YX = C_modifies_YX_vec))

K_list = lapply(pi_list, FUN = length)
Kmax = max(simplify2array(K_list))

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


#Sort and check pctmiss is between 0 and 1
pctmiss_vec = sort(pctmiss_vec)
if( (sum(pctmiss_vec<=0) + sum(pctmiss_vec>=1))>0){stop("elements in pctmiss_vec must be between 0 and 1 (exclusive).")}

# Define working directory
main_wd = paste(dropbox_wd, "/Dissertation/lpa-mi-pool", sep = "")
temp_fname = paste("sim-lpa-mi-pool", format(Sys.time(), "%b %d %Y %H %M"), sep = " ")

# #  Create folder to store results in main directory
 setwd(paste(dropbox_wd, "/Dissertation/results", sep = ""))
 dir.create(temp_fname)
 setwd(temp_fname)
 results_wd = getwd()

# Check that the dropbox and temp_dir are found on machine
if(length(temp_dir)==1){
  if (dir.exists(dropbox_wd)*dir.exists(temp_dir)==0){
    stop("Provided dropbox_wd or temp_dir paths not found." )
  }
} else {
  if (!dir.exists(dropbox_wd)){
    stop("Provided dropbox_wd paths not found." )
  }
  wd_tmp = getwd()
  message("Testing access to each temp_dir:")
  for (p in 1:length(temp_dir)){
    setwd(temp_dir[p])
    message(paste0("  ",temp_dir[p], " successfully accessed"))
    setwd(wd_tmp)
  }
}

# Create temporary folder
temp_wd =  paste(temp_dir, temp_fname, sep = "/")
for (p in 1:length(temp_wd)){
  dir.create(temp_wd[p])
}

# Create processor-specific temporary folder

if (Processors%%length(temp_wd)!=0){
  stop("Number of working directories is not a multiple of the number of processors. Make sure length(temp_wd)%%Processors == 0.")
}

if (length(temp_wd)==1){
  setwd(temp_wd)
  for (p in 1:Processors){
    dir.create(paste("p",p, sep = ""))
  }
  temp_wd_p_vec = dir(path = temp_wd, full.names = TRUE)
} else {
  wt_tmp = getwd()
  temp_wd_p_vec = rep("",Processors)
  n_twd = 1
  for (p in 1:Processors){
    setwd(temp_wd[n_twd])
    dir.create(paste("p",p, sep = ""))
    setwd(paste("p",p, sep = ""))
    
    temp_wd_p_vec[p] = getwd()
    n_twd = ifelse(n_twd==length(temp_wd), 1, n_twd+1)
    setwd(wd_tmp)
  }
  
} 


# Within the processor-specific temporary folder, create subfolders
for (p in 1:Processors){
  setwd(temp_wd_p_vec[p])
  dir.create("Complete data")
  dir.create("Observed data")
  dir.create("Imputed data")
  dir.create("Stacked data")
  
  for (pm in 1:length(pctmiss_vec)){
    dir.create(paste0(temp_wd_p_vec[p],"/Observed data/pm", pm))
    dir.create(paste0(temp_wd_p_vec[p],"/Imputed data/pm", pm))
    dir.create(paste0(temp_wd_p_vec[p],"/Stacked data/pm", pm))
    
    for (pva in 1:length(methods_list$procedure)){
      dir.create(paste0(temp_wd_p_vec[p],"/Imputed data/pm", pm,"/pva",pva))
      dir.create(paste0(temp_wd_p_vec[p],"/Stacked data/pm", pm,"/pva",pva))
    }
  }
  
}




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


#rep = 1
#z = 1
#p = 1

# # #### Conduct parallel processing
# # Register processors
 cl<-makePSOCKcluster(Processors)
 registerDoParallel(cl)

 out_foreach<- foreach(p = 1:Processors,
                       .packages = packages_vec, .options.RNG=42) %dorng%{
    for(rep in 1:Replications){
      summary_df <- parameters_df <- NULL
        for (z in 1:Ztot){
          
#          sink(file = paste0(temp_wd_p_vec[p],"/sink-p", p, ".txt"))
          
          file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), full.names = TRUE, recursive = TRUE))
          file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), full.names = TRUE, recursive = TRUE))
          file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), full.names = TRUE, recursive = TRUE))
          file.remove(list.files(paste0(temp_wd_p_vec[p],"/Stacked data"), full.names = TRUE, recursive = TRUE))
          
          
          
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
                                                       save_it = TRUE, 
                                                       temp_wd_p_vec = temp_wd_p_vec,
                                                       p =p , rep = rep)

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
            dff_observed= list_observed$dffolderfiles
            dff_imputed = subset(list_imputed$dffolderfiles, !is.na(m))

            
            
            ### Enumerate and save results ###
            KKmax_z = K_list[[which(data_conditions$class_size[z] == names(K_list))]]+1
#            for (kk in 3){
              
                # pop_params_kk = pop_params;
                # pop_params_kk$K_z = kk;
                # pop_params_kk$mu_z = NULL;
                # pop_params_kk$S_z = NULL;
                # pop_params_kk$pi_z = NULL;
                
                kk = pop_params$K_z
                
                print("")
                print("")
                print("")
                print("Complete Data")
                  # Complete data
                    source_tmp = "Complete data"
                    target_wd = paste0(temp_wd_p_vec[p], "/Complete data")
                    # Fit the model
                    out_ftc = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                  dff_target = dff_complete, target_wd = target_wd, pop_params_kk = pop_params, 
                                                  starts0 = 20, mX_max = 6)
                    # Extract the parameters
                    if (out_ftc$problem == FALSE){
                        out_Mplus = readModels(target = target_wd, 
                                             what = c("summaries", "parameters"))
                        summary_complete = out_Mplus$summaries %>% transform(p = p, rep = rep, z = z, kk = kk, pm = NA, pva = NA, Source = source_tmp) %>% 
                                                               transform(Converged = out_ftc$Converged, Rcond = out_ftc$Rcond, Starts = out_ftc$Starts)
                        parameters_complete = out_Mplus$parameters$unstandardized %>% transform(p = p, rep = rep, z = z, kk = kk, pm = NA, pva = NA, Source = source_tmp)
                        
                        # Resolve Label Switch
                        mu_df = parameters_complete[parameters_complete$paramHeader == "Means" & startsWith(parameters_complete$param,"Y"), c("est","LatentClass")]
                        mu_est_mat = matrix(mu_df$est, nrow = pop_params$J_Y_z, ncol = kk)
                        S_est_array = array(0, dim = c(3,3,kk))
                        for (k in 1:kk){
                          for (i in 1:3){
                            S_est_array[i,i,k] = 1
                          }
                        }
                        out_label_switch = resolve_label_switch(mu_est_mat, S_est_array, z, data_conditions, name = "glpk", parameters_df = parameters_complete)
                        parameters_complete = out_label_switch$parameters_df
                        summary_complete = transform(summary_complete, Labels_Switched = out_label_switch$class_switched)
                    }
                    
                print("")
                print("")
                print("")
                print("Observed data")
                  # Observed data
                    source_tmp = "Observed data"
                    summary_observed <- parameters_observed <- NULL
                    for (pm in 1:length(pctmiss_vec)){
                      target_wd = paste0(temp_wd_p_vec[p], "/Observed data/pm", pm)
                      # Fit the model
                        out_ftc = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, dff_target = dff_observed, 
                                                      target_wd = target_wd, pop_params_kk = pop_params, 
                                                      starts0 = 20, mX_max = 6)
                        if (out_ftc$problem == FALSE){
                          # Extract the parameters
                          out_Mplus = readModels(target = target_wd, 
                                                 what = c("summaries", "parameters"))
                          summary_pm = out_Mplus$summaries %>% transform(p = p, rep = rep, z = z, kk = kk, pm = pm, pva = NA, Source = source_tmp) %>% 
                            transform(Converged = out_ftc$Converged, Rcond = out_ftc$Rcond, Starts = out_ftc$Starts)
                          parameters_pm = out_Mplus$parameters$unstandardized %>% transform(p = p, rep = rep, z = z, kk = kk, pm = pm, pva = NA, Source = source_tmp)
                          # Augment
                          summary_observed = rbind(summary_observed, summary_pm)
                          parameters_observed = rbind(parameters_observed, parameters_pm)
                          
                          # Resolve Label Switch
                          mu_df = parameters_observed[parameters_observed$paramHeader == "Means" & startsWith(parameters_observed$param,"Y"), c("est","LatentClass")]
                          mu_est_mat = matrix(mu_df$est, nrow = pop_params$J_Y_z, ncol = kk)
                          S_est_array = array(0, dim = c(3,3,kk))
                          for (k in 1:kk){
                            for (i in 1:3){
                              S_est_array[i,i,k] = 1
                            }
                          }
                          out_label_switch = resolve_label_switch(mu_est_mat, S_est_array, z, data_conditions, name = "glpk", parameters_df = parameters_observed)
                          parameters_observed = out_label_switch$parameters_df
                          summary_observed = transform(summary_observed, Labels_Switched = out_label_switch$class_switched)
                        }
                    }

                    
                print("")
                print("")
                print("")
                print("Imputed data")
                  # Imputed data
                    source_tmp = "Imputed data"
                    summary_imputed <- parameters_imputed <- NULL
                    for (pm in 1:length(pctmiss_vec)){
                      for(pva in 1:length(methods_list$procedure)){ 
                        target_wd = paste0(temp_wd_p_vec[p], "/Imputed data/pm", pm,"/pva",pva)
                        # Fit the model
                        out_ftc = fit_til_convergence(p = p, z = z, temp_wd_p_vec = temp_wd_p_vec, dff_target = dff_imputed, 
                                                      target_wd = target_wd, pop_params_kk = pop_params, 
                                                      starts0 = 20, mX_max = 6)
                        if (out_ftc$problem == FALSE){
                          # Extract the parameters
                          out_Mplus = readModels(target = target_wd, 
                                                 what = c("summaries", "parameters"))
                          summary_pm_pva = out_Mplus$summaries %>% transform(p = p, rep = rep, z = z, kk = kk, pm = pm, pva = pva, Source = source_tmp)  %>% 
                                                                   transform(Converged = out_ftc$Converged, Rcond = out_ftc$Rcond, Starts = out_ftc$Starts)
                          parameters_pm_pva = out_Mplus$parameters$unstandardized %>% transform(p = p, rep = rep, z = z, kk = kk, pm = pm, pva = pva, Source = source_tmp)
                          # Augment
                          summary_imputed = rbind(summary_imputed, summary_pm_pva)
                          parameters_imputed = rbind(parameters_imputed, parameters_pm_pva)
                          # Resolve Label Switch
                          mu_df = parameters_imputed[parameters_imputed$paramHeader == "Means" & startsWith(parameters_imputed$param,"Y"), c("est","LatentClass")]
                          mu_est_mat = matrix(mu_df$est, nrow = pop_params$J_Y_z, ncol = kk)
                          S_est_array = array(0, dim = c(3,3,kk))
                          for (k in 1:kk){
                            for (i in 1:3){
                              S_est_array[i,i,k] = 1
                            }
                          }
                          out_label_switch = resolve_label_switch(mu_est_mat, S_est_array, z, data_conditions, name = "glpk", parameters_df = parameters_imputed)
                          parameters_imputed = out_label_switch$parameters_df
                          summary_imputed = transform(summary_imputed, Labels_Switched = out_label_switch$class_switched)
                        }
                      }
                    }
                  
                    
                  # Combine all data sources
                  summary_df = rbind(summary_df, summary_imputed, summary_observed, summary_complete)
                  parameters_df = rbind(parameters_df, parameters_imputed, parameters_observed, parameters_complete)
                    
                    
#            } # End for kk
            
          sink()
        } #end for z = 1:Ztot
      
      #Save the progress
      progress_df = data.frame(progress = round(100*(rep/Replications),2))
      write.csv(x = progress_df, file = paste0(results_wd,"/P", p, "-progress.csv"), row.names = FALSE)
      
      # Save the result
      save(summary_df, file = paste0(results_wd,"/summary-p", p, "-rep",rep,".RData"))
      save(parameters_df, file = paste0(results_wd,"/parameters-p", p, "-rep",rep,".RData"))
    } #end for rep = 1:Replications
} #for each p = 1,...,P

 stopImplicitCluster()
 setwd(results_wd)
 save.image("environment.RData")
