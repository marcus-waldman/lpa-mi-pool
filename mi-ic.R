rm(list = ls())

library(plyr)
library(tidyverse)
library(lpa.mi.src)
library(doParallel)
library(foreach)
library(doRNG)

Processors = detectCores()
Replications = 1
KKmax = 4;




# Acer-laptop
  #dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
  #temp_dir = "C:/Users/marcu/Documents/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)
# Queen Mary's Revenge
  dropbox_wd =  "D:/Dropbox"  #  Dropbox folder on local machine
  temp_dir = paste0(LETTERS[seq(26-Processors+1,26)],":") #  Temporary directory for local maching (will automaticall add new temporary folder)        
# Xi-GSU
  #dropbox_wd = "C:/Users/mwaldman1/Dropbox"
  #temp_dir = "C:/Users/mwaldman1/Documents/lpa-mi-pool"
# Big Bertha
  #dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine
  #temp_dir = "D:/lpa-mi-pool" #  Temporary directory for local maching (will automaticall add new temporary folder)        


# Nuke directory
  lpa.mi.src::nuke_dirs(temp_dir)
    
#  Create folder to store results in main directory
   setwd(paste(dropbox_wd, "/Dissertation/results", sep = ""))
   temp_fname = paste("sim-lpa-mi-pool", format(Sys.time(), "%b %d %Y %H %M"), sep = " ")
   dir.create(temp_fname)
   setwd(temp_fname)
   results_wd = getwd()
  

# Fully-crossed simulation conditions
pctmiss_vec = c(0.5)     #  Percent of observations with *at least* one missing value (indicator or missing data correlate) -> if other specification needed one must modify list_inds.miss generation in get_obs_data.R
N_vec = c(500, 2E3)                #  Sample sizes within each replication
J_Y_vec = c(3)                #  Number of latent class indicators
J_Xcom_vec = c(1)             #  Number of complete data missing data correlates
J_Xinc_vec = c(0)             #  Number of missing data correlates which themselves contain missing data
MD_vec = c(2.75, 3.25, 4)               #  WARNING! MD = 1 RESULTS IN ERROR! Mahalanobis distances between classes
pi_list = list(               #  Number of classes, K (maximum of 3)
  A = rep(1/3,3) #,
  #B = c(0.4,0.4,0.2), 
  #C = rep(1/4,4)
)
rho_YX_vec = c(0.4)          #  Strength of missing data correlates
C_modifies_YX_vec = c(TRUE)  # Whether or not the Y-X relationships are different between latent classes 

# Define the imputation alogrithms to compare
methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
methods_list$procedure[[1]] = "stratamelia";  methods_list$name[[1]] = "mg-Amelia";     methods_list$args[[1]] = list(m = 25, p2s = 0, tolerance = 1E-4)
#methods_list$procedure[[2]] = "mice";         methods_list$name[[2]] = "mg-mice";       methods_list$args[[2]] = list(method = "bygroup", imputationFunction = "norm", m = 5, maxit = 20)
#methods_list$procedure[[3]] = "mice";         methods_list$name[[3]] = "mg-mvn";        methods_list$args[[3]] = list(method = "bygroup", imputationFunction = "norm", m = 20, maxit = 50, blocks = NULL)


#### Pre-processing ####

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


#Sort and check pctmiss is between 0 and 1
pctmiss_vec = sort(pctmiss_vec)
if( (sum(pctmiss_vec<=0) + sum(pctmiss_vec>=1))>0){stop("elements in pctmiss_vec must be between 0 and 1 (exclusive).")}

# Define working directory
main_wd = paste(dropbox_wd, "/Dissertation/lpa-mi-pool", sep = "")


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


# #### Conduct parallel processing
# Register processors
 cl<-makePSOCKcluster(Processors)
 registerDoParallel(cl)



out_foreach<- foreach(p = 1:Processors,
                      .packages = packages_vec, .options.RNG=42) %dorng%{
    dff_complete<-dff_imputed<-dff_majority<-dff_averaged<-dff_stacked<-NULL;
    for(rep in 1:Replications){
      results_df = NULL
        for (z in 1:Ztot){
          
          sink(file = paste0(temp_wd_p_vec[p],"/sink-p", p, ".txt"))
          
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
            
            list_stacked = lpa.mi.src::get_stacked_data(z = z, 
                                            list_imputed = list_imputed, 
                                            methods_list = methods_list, 
                                            pctmiss_vec = pctmiss_vec, 
                                            save_it = TRUE,
                                            p = p, 
                                            rep = rep, 
                                            temp_wd_p_vec = temp_wd_p_vec)
            
            dff_complete = list_complete$dffolderfiles
            dff_observed= list_observed$dffolderfiles
            dff_majority = subset(list_imputed$dffolderfiles, !is.na(m))
            dff_averaged = subset(list_imputed$dffolderfiles, is.na(m))
            dff_stacked = list_stacked$dffolderfiles
            
            
            ### Enumerate and save results ###
            results_z <- LL_z <- NULL
            for (kk in 1:KKmax){
              
 
                if (kk == 1){
                  starts0_k = 0
                  mX_max_k = 2
                  starts_txt_k = "0;"
                }  
                if (kk==2 | kk == 3){
                  starts0_k = 20
                  mX_max_k = 6 
                  starts_txt_k = "20 8;"
                }
                if (kk == 4){
                  starts0_k = 320
                  mX_max_k = 3
                  starts_txt_k = "320 80;"
                }
             
              
                if(kk == pop_params$K_z){
                  pop_params_kk = pop_params;
                } else {
                  pop_params_kk = pop_params;
                  pop_params_kk$K_z = kk
                  pop_params_kk$mu_z = NULL;
                  pop_params_kk$S_z = NULL
                  pop_params_kk$pi_z = NULL;
                }
                
                
          ##### Complete data #####
                    print("")
                    print("")
                    print("")
                    print("Complete data")
                    results0_kk = NULL
                    # Delete old Mplus .inp and .out files
                    file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".inp", full.names = TRUE, recursive= TRUE))
                    file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
                    # Fit the model
                    out_ftc_complete = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                  dff_target = dff_complete, 
                                                  target_wd = paste0(temp_wd_p_vec[p], "/Complete data"), 
                                                  starts0 = starts0_k, mX_max = mX_max_k, pop_params_kk = pop_params_kk, type_imputation = FALSE)
                  
                    if (out_ftc_complete$problem == FALSE){
                      svals_txt_complete<-extract_svals(file = list.files(path = paste0(temp_wd_p_vec[p], "/Complete data"), pattern = ".out"), 
                                                        path =  paste0(temp_wd_p_vec[p], "/Complete data"))
                      model_txt_complete = c("MODEL:\n",svals_txt_complete)
                    } else {
                      model_txt_complete = NULL
                    }
                    results0_kk = lpa.mi.src::extract_ICs_naive(TYPE = "Complete data",DATATYPE = "Complete data", p = p, temp_wd_p_vec = temp_wd_p_vec,  k= kk, z = z, rep = rep , starts_txt = starts_txt_k)
                    results0_kk = transform(results0_kk, Problem = ifelse(kk==1,FALSE,out_ftc_complete$problem))
                    
                    
               
                    
          ##### Observed data #####
                  print("")
                  print("")
                  print("")
                  print("Observed Data")
                  results1_kk = NULL
                  # Delete old Mplus .inp and .out files
                  file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".inp", full.names = TRUE, recursive= TRUE))
                  file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
                  # Set target directory
                  for (pm in 1:length(pctmiss_vec)){
                      idx_pm = which(endsWith(as.character(dff_observed$folders),paste0("pm",pm)))
                      # Fit the model
                      out_ftc_observed = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                             dff_target = dff_observed[idx_pm, ], 
                                                             target_wd = paste0(temp_wd_p_vec[p], "/Observed data/pm",pm),
                                                             starts0 = starts0_k, mX_max = mX_max_k, pop_params_kk = pop_params_kk, type_imputation = FALSE)
                      #  Extract IC estimates 
                      results1_kk_pm = lpa.mi.src::extract_ICs_naive(TYPE = "Observed data", DATATYPE = "Observed data", p = p, temp_wd_p_vec = temp_wd_p_vec,  k= kk, z = z, rep = rep , starts_txt = starts_txt_k, pm = pm)
                      results1_kk_pm = transform(results1_kk_pm, Problem = ifelse(kk==1,FALSE,out_ftc_observed$problem))
                      # Extract IC estimates 
                      results1_kk = rbind(results1_kk, results1_kk_pm)
                  }
                        
                
                
          ##### Averaged ICs #####
                    print("")
                    print("")
                    print("")
                    print("Averaged ICs")
                    results2_kk = NULL  
                  # Delete old Mplus .inp and .out files
                    file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
                    file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
                    svals_list_rr = list(NULL); vv = 0;
                    for(pm in 1:length(pctmiss_vec)){ 
                      for(pva in 1:length(methods_list$procedure)){
                        idx_pm_pva = which(endsWith(as.character(dff_averaged$folders),paste0("pm",pm,"/pva",pva)))
                        out_ftc_averaged = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                               dff_target = dff_averaged[idx_pm_pva, ], 
                                                               target_wd = paste0(temp_wd_p_vec[p], "/Imputed data/pm",pm,"/pva",pva),
                                                               starts0 = starts0_k, mX_max = mX_max_k, pop_params_kk = pop_params_kk, 
                                                               type_imputation = TRUE)
                        #  Extract IC estimates 
                        results2_kk_pm_pva = lpa.mi.src::extract_ICs_naive(TYPE = "Averaged", DATATYPE = "Imputed data", 
                                                                           p = p, temp_wd_p_vec = temp_wd_p_vec,  k= kk, z = z, rep = rep , 
                                                                           starts_txt = starts_txt_k, pm = pm, pva = pva, 
                                                                           type_imputation = TRUE)
                        results2_kk_pm_pva = transform(results2_kk_pm_pva, Problem = ifelse(kk==1,FALSE,out_ftc_observed$problem))
                        results2_kk = rbind(results2_kk, 
                                            results2_kk_pm_pva)
                        
                        # Exract Rubin-Rules starting values
                        wd_pm_pva = paste0(temp_wd_p_vec[p], "/Imputed data/pm",pm,"/pva",pva)
                        vv = vv+1;
                        svals_list_rr[[vv]]<-extract_svals(file = list.files(path = wd_pm_pva, pattern = ".out"), 
                                                    path =  wd_pm_pva)
                        names(svals_list_rr)[[vv]] = paste0("svals_pm",pm,"_pva",pva)
                      }
                    }
               
                    
               ##### Majority Vote ICs #####
                    print("")
                    print("")
                    print("")
                    print("Majority vote ICs")
                    results3_kk = NULL  
                  # Delete old Mplus .inp and .out files
                    file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
                    file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
                    results3_kk = NULL
                    for(pm in 1:length(pctmiss_vec)){ 
                      for(pva in 1:length(methods_list$procedure)){
                        M_pva = methods_list$args[[pva]]$m
                        for (m in 1:M_pva){
                            
                            idx_pm_pva_m = which(endsWith(as.character(dff_majority$folders),paste0("pm",pm,"/pva",pva)) & 
                                                   endsWith(as.character(dff_majority$files), paste0("imp_",m,".dat")))
                            out_ftc_m = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                                   dff_target = dff_majority[idx_pm_pva_m,], 
                                                                   target_wd = paste0(temp_wd_p_vec[p], "/Imputed data/pm",pm,"/pva",pva),
                                                                   starts0 = starts0_k, mX_max = mX_max_k, pop_params_kk = pop_params_kk, 
                                                                   type_imputation = FALSE, Model_txt = model_txt_complete, m = m)
                            
                            # # Fit the models and save the starting values
                            # invisible(runModels(target = paste0(temp_wd_p_vec[p], "/Imputed data/pm",pm,"/pva",pva)))
                            # # Extract IC estimates 
                            # results3_kk_pm_pva_m = rbind(results3_kk_pm_pva_m, 
                            #                             lpa.mi.src::extract_ICs_majority(p = p, temp_wd_p_vec = temp_wd_p_vec,  k= kk, z = z, rep = rep , starts_txt = starts_txt, pm = pm, pva = pva))
                            #)
                        
                        }
                        #  Extract IC estimates 
                        results3_kk_pm_pva = extract_ICs_majority(p = p, temp_wd_p_vec, kk = kk, z = z, rep = rep , starts_txt = starts_txt_k, pm = pm, pva = pva)
                        #(TYPE = "Majority vote", DATATYPE = "Imputed data", 
                        #                                      p = p, temp_wd_p_vec = temp_wd_p_vec,  k= kk, z = z, rep = rep , 
                        #                                      starts_txt = starts_txt_k, pm = pm, pva = pva, 
                        #                                      type_imputation = FALSE)
                        results3_kk_pm_pva = transform(results3_kk_pm_pva, Problem = ifelse(kk==1,FALSE,out_ftc_m$problem))
                        results3_kk = rbind(results3_kk, results3_kk_pm_pva)
                        
                      }
                    }
                    
                    
              ##### Stacked Vote Ics #####
                    print("")
                    print("")
                    print("")
                    print("Stacked  ICs")
                    results4_kk = NULL   
                  # Delete old Mplus .inp and .out files
                    file.remove(list.files(paste0(temp_wd_p_vec[p],"/Stacked data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
                    file.remove(list.files(paste0(temp_wd_p_vec[p],"/Stacked data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
                    bb = 0 
                    for(pm in 1:length(pctmiss_vec)){ 
                      for(pva in 1:length(methods_list$procedure)){

                        idx_pm_pva = which(endsWith(as.character(dff_averaged$folders),paste0("pm",pm,"/pva",pva)))
                        out_ftc_stacked = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                               dff_target = dff_stacked[idx_pm_pva, ], 
                                                               target_wd = paste0(temp_wd_p_vec[p], "/Stacked data/pm",pm,"/pva",pva),
                                                               starts0 = starts0_k, mX_max = mX_max_k, pop_params_kk = pop_params_kk, 
                                                               type_imputation = FALSE, Model_txt = model_txt_complete, 
                                                               weight = TRUE)

                        #  Extract IC estimates 
                        # results4_kk_pm_pva = lpa.mi.src::extract_ICs_naive(TYPE = "Stacked", DATATYPE = "Stacked data", 
                        #                                                    p = p, temp_wd_p_vec = temp_wd_p_vec,  k= kk, z = z, rep = rep , 
                        #                                                    starts_txt = starts_txt_k, pm = pm, pva = pva, 
                        #                                                    type_imputation = FALSE)
                        bb = bb + 1
                        results4_kk_pm_pva = lpa.mi.src::extract_ICs_stacked(p = p, temp_wd_p_vec = temp_wd_p_vec, kk = kk, z = z, rep = rep, 
                                                                             starts_txt = starts_txt_k, mids_pm_pva = list_imputed$obj_call[[bb]], pm = pm, pva = pva)
                                                                             
                        results4_kk_pm_pva = transform(results4_kk_pm_pva, Problem = ifelse(kk==1,FALSE,out_ftc_stacked$problem))
                        results4_kk = rbind(results4_kk, 
                                            results4_kk_pm_pva)
                        
                      }
                    }
                    
          ##### D-statistic ####
                  print("")
                  print("")
                  print("")
                  print("Obtain constrained LL for for calculating D-statistic")
              # Delete old Mplus .inp and .out files
                  file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
                  file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
                  LL_kk = NULL
                  for(pm in 1:length(pctmiss_vec)){ 
                    for(pva in 1:length(methods_list$procedure)){
                     
                      LL_tmp = data.frame(params=NA, Converged=NA, LL_free=NA, LL_fixed=NA, kk = kk, z = z, rep = rep, p = p, pm = pm, pva = pva)
                      
                    
                      idx_pm_pva = which(results2_kk$pm==pm & results2_kk$pva == pva)
                      converged_pm_pva = results2_kk$Converged[idx_pm_pva]
                      LL_tmp$params = results2_kk$Parameters[idx_pm_pva]
                      LL_tmp$Converged =converged_pm_pva
                      
                      if(converged_pm_pva==TRUE){
                          LL_tmp$LL_free = results2_kk$LL[idx_pm_pva]
                          
                          idx_pm_pva = which(names(svals_list_rr) == paste0("svals_pm",pm,"_pva",pva))
                          model_txt_pm_pva = c("MODEL:\n", gsub("*","@",svals_list_rr[[idx_pm_pva]], fixed = TRUE))
                          inp_RRfixed = lpa.mi.src::create_naiveMplus_inpfile(z = z, 
                                                                               out_get_FMM = pop_params_kk, 
                                                                               dffolderfiles = dff_averaged[idx_pm_pva,], 
                                                                               temp_wd_p = temp_wd_p_vec[p], 
                                                                               starts_txt = "0;",
                                                                               type_imputation = TRUE, 
                                                                               Model_txt = model_txt_pm_pva)
                          
                          # Fit the models and save the starting values
                          runModels(target = paste0(temp_wd_p_vec[p], "/Imputed data/pm",pm,"/pva",pva))
                          # Save the fixed log-likelihood
                          out_RRfixed = readModels(target =paste0(temp_wd_p_vec[p], "/Imputed data/pm",pm,"/pva",pva), what = c("summaries"))
                          LL_tmp$LL_fixed = out_RRfixed$summaries$LL_Mean
                      }
                      LL_kk = rbind(LL_kk, LL_tmp)
                    }
                  }
                  
              # Append results
                results_z = rbind(results_z, results0_kk, results1_kk, results2_kk, results3_kk, results4_kk)
                LL_z = rbind(LL_z, LL_kk)
               
            } # End for kk 
            
          sink()
        } #end for z = 1:Ztot
      
      #Save the progress
      progress_df = data.frame(progress = round(100*(rep/Replications),2))
      write.csv(x = progress_df, file = paste0(results_wd,"/P", p, "-progress.csv"), row.names = FALSE)
      
      # Save the result
      save(results_df, file = paste0(results_wd,"/p", p, " rep",rep,".RData"))
    } #end for rep = 1:Replications
} #for each p = 1,...,P

stopImplicitCluster()
save.image("environment.RData")
