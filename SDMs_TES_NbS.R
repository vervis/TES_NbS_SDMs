###################################################################################
## Species distribution modelling using MaxEnt
###################################################################################
setwd("D:/")

#Load (and install) packages
if(!"devtools" %in% installed.packages()) install.packages("devtools")
if(!"blockCV" %in% installed.packages())remotes::install_github("rvalavi/blockCV", 
                                                                dependencies = TRUE)
if(!"covsel" %in% installed.packages()) devtools::install_github("N-SDM/covsel", 
                                                                 dependencies = TRUE)
if(!"rmaxent" %in% installed.packages()) devtools::install_github('johnbaums/rmaxent', 
                                                                  dependencies = TRUE)
if(!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, terra, sf, enmSdmX, covsel, blockCV, ENMeval,
               maps, dismo, rmaxent, radiant.data)

######################################################################
# Maxent modelling
######################################################################
#---------------------------------------------------------
#Presence data organisation
#---------------------------------------------------------

# # Read in summary of occurrence counts
# summ_occ_cleaned_df = read.csv("Output/Summary_occurrences_cleaning.csv") %>%
#   filter(spp_file!="NA_NA.csv")
# 
# # Get list of species with enough occurrences (>20)
# spp_maxent = summ_occ_cleaned_df %>%
#   filter(!is.na(num_occ_sregion) & num_occ_sregion>20)
# 
# # Read in ID raster
# rastNA = rast("Output/rastNA.tif")
# 
# # Read in weighted grid raster
# rast_back = rast("Output/rast_back.tif")
# 
# # Remove duplicate occurrences within raster pixels. Then recheck whether there are enough occurrences for modelling
# for(sp in 1:nrow(spp_maxent)){
#   # Get species occurrence data
#   pres = read.csv(paste0("Output/GBIF_cleaned/",spp_maxent$spp_file[sp]))
#   # Extract IDs for occurrences
#   pres_ids = terra::extract(rastNA, cbind(pres$lon, pres$lat)) 
#   # Remove duplicates
#   pres = pres[!duplicated(pres_ids) & !is.na(pres_ids),] 
#   # Change number of occurrences in study region
#   spp_maxent$num_occ_sregion[sp] = nrow(pres)
# }
# # Get list of species with enough occurrences (>20)
# spp_maxent = spp_maxent %>%
#   filter(!is.na(num_occ_sregion) & num_occ_sregion>20)    
# write.csv(spp_maxent, "Output/spp_maxent.csv", row.names = F)
spp_maxent = read.csv("Output/spp_maxent.csv")

# Get species names (and not file names)
spp_maxent$species = unlist(strsplit(spp_maxent$spp_file, split = ".csv"))

#---------------------------------------------------------
# Some initial setup
#---------------------------------------------------------
# Get prepPara function for setting maxent parameters
# Read in from GitHub
# devtools::source_url("https://raw.githubusercontent.com/shandongfx/nimbios_enm/master/Appendix2_prepPara.R")
# OR read from saved R file
source("Code/prepPara.R")

#RCM names
clim_mods = c("CCCma-CanESM2",
              "CNRM-CERFACS-CNRM-CM5",
              "MOHC-HadGEM2-ES",
              "ICHEC-EC-EARTH",
              "CMPI-M-MPI-ESM-LR",
              "IPSL-IPSL-CM5A-LR",
              "MIROC-MIROC5")

# Environmental predictor set names
env_pred_set_names = c("env_stack_rcp45_bsl",
                       "env_stack_rcp85_bsl",
                       "env_stack_rcp45_gwl1.5",
                       "env_stack_rcp85_gwl1.5",
                       "env_stack_rcp45_gwl2.0",
                       "env_stack_rcp85_gwl2.0",
                       "env_stack_rcp85_gwl4.0")

bsl_stacks = c( paste0("env_stack_", "rcp45_bsl_rcm", c(1,2,3,5)),
                paste0("env_stack_", "rcp85_bsl_rcm", 1:7)
)

# Sequence of regularisation parameters to test in Maxent
reg_params = c(0.5,1,2,3,4)
#https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13639
#https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13639 c(0.5,1,2,3,4)
#https://peerj.com/articles/13337/ seq(1,5,1)

# Maxent features to use
feats = "LQP"

#---------------------------------------------------------
# SDM loop
#---------------------------------------------------------
for(sp in 3500:nrow(spp_maxent)){
  # maxent_run = function(sp){
  
  #---------------------------------------------------------
  #Some initialisation stuff
  #---------------------------------------------------------
  # Get original names of species
  sp_name = spp_maxent$species[sp]
  
  results_list = list() #List for storing results
  
  #Create output directory for Maxent models
  if(!"Output/Maxent_models" %in% list.dirs("Output/")){dir.create("Output/Maxent_models/")}
  output_dir = paste0('Output/Maxent_models/',spp_maxent$species[sp])
  if(!file.exists(output_dir)){ dir.create(output_dir, showWarnings = FALSE, recursive = T) }
  
  #Create name of results output files
  results_folds_file = paste0(output_dir,"/SDM_folds_results.csv")
  results_file = paste0(output_dir,"/SDM_results.csv")
  
  # Get species occurrence data
  pres = read.csv(paste0("Output/GBIF_cleaned/",spp_maxent$spp_file[sp]))
  # Extract IDs for occurrences
  pres_ids = terra::extract(rastNA, cbind(pres$lon, pres$lat)) 
  # Remove duplicates
  pres = pres[!duplicated(pres_ids) & !is.na(pres_ids),] 
  
  #Create sf object from pres
  pres_sf = st_as_sf(x = pres,
                     coords = c("lon", "lat"),
                     crs = crs(rastNA)) 
  
  counter = 0 #Used for adding accuracy results to results file
  
  #----------------------------------------------------------------- 
  # Get background points
  #-----------------------------------------------------------------
  # Create 500 km circular buffers around presence localities for sampling background points
  sf_use_s2(FALSE) #Prevents getting an error
  buffs = sf::st_buffer(pres_sf, dist = 500/1.86/60)
  # Merge all buffer circles into one shapefile layer
  buffs = st_union(buffs)
  buffs = vect(buffs)
  buffs$ID = as.character(1)
  # Convert shapefile to raster
  buffs_rast = rasterize(buffs, y = rastNA, field = "ID", fun = 'sum')
  
  # Mask rast_back to background extent
  rast_back_sp = terra::mask(rast_back, buffs_rast)
  
  #Mask out pixels for which species of interest occurs in
  #Get the cell index for each occurrence record of the species of interest and make a table
  count_pres_rast = table(cellFromXY(rast_back_sp, cbind(pres$lon, pres$lat)))
  #Make cells NA in which presences occur
  rast_back_sp[as.numeric(names(count_pres_rast))] = NA
  
  #Get weighted sample of 30,000 background points from within weighted grid
  set.seed(1)
  backs_xy = as.data.frame(sampleRast(rast_back_sp, n = 50000, prob = T)) #Sample extra because it also samples NA cells
  backs_na = terra::extract(rastNA, cbind(backs_xy$x, backs_xy$y))
  backs_xy = bind_cols(backs_xy, backs_na)
  # summary(backs_xy)
  
  #Remove NAs and extra background points
  backs_xy = na.omit(backs_xy)
  backs_xy = backs_xy[sample(nrow(backs_xy), size = 30000, replace = F),]
  backs_xy = backs_xy %>% dplyr::select(x, y)
  
  #Plot background points (if wanted)
  # ggplot(backs_xy, aes(x, y)) + geom_point()
  
  #Write background points as CSV
  write.csv(backs_xy, paste0(output_dir,"/backs_xy.csv"), row.names = F)
  
  
  #---------------------------------------------------------
  # Climate stacks loop
  #---------------------------------------------------------
  for(b in c(1:9,11)){ #1:length(bsl_stacks)
    #Get rcp
    rcp = str_split_i(bsl_stacks[b], "_", 3)
    #Get rcm
    rcm = str_split_i(bsl_stacks[b], "_", 5) %>%
      str_split_i("rcm", 2) %>% as.numeric()
    
    # Read in RCM stacks and mask out common NA pixels
    env_stack_bsl = rast(paste0("Output/Predictors/",bsl_stacks[b],".tif"))
    env_stack_bsl = terra::mask(env_stack_bsl, rastNA)
    env_stack_1.5 = rast(paste0("Output/Predictors/env_stack_", rcp, "_gwl1.5_rcm", rcm, ".tif"))
    env_stack_2.0 = rast(paste0("Output/Predictors/env_stack_", rcp, "_gwl2.0_rcm", rcm, ".tif"))
    if(rcp=="rcp85" & rcm!=7){
      env_stack_4.0 = rast(paste0("Output/Predictors/env_stack_", rcp, "_gwl4.0_rcm", rcm, ".tif"))
    }
    
    # Mask environmental layers to background extent
    env_stack_bsl_buff = terra::mask(env_stack_bsl, buffs_rast)
    env_stack_1.5_buff = terra::mask(env_stack_1.5, buffs_rast)
    env_stack_2.0_buff = terra::mask(env_stack_2.0, buffs_rast)
    if(rcp=="rcp85" & rcm!=7){
      env_stack_4.0_buff = terra::mask(env_stack_4.0, buffs_rast)
    }
    
    #-----------------------------------------------------------------
    #Join presence and background data together
    #-----------------------------------------------------------------
    pa = rbind(data.frame(pres_abs=1, x=pres$lon, y=pres$lat), 
               data.frame(pres_abs=0, x=backs_xy$x, y=backs_xy$y))
    
    #-----------------------------------------------------------------
    #Add environmental data to PA data
    #-----------------------------------------------------------------
    pa_env = cbind(pa, terra::extract(env_stack_bsl, cbind(pa$x, pa$y)))
    pa_env = na.omit(pa_env)
    
    #Create sf object for use below using blockCV
    pa_sf = st_as_sf(x = pa_env %>% dplyr::select(pres_abs, x, y),
                     coords = c("x", "y"),
                     crs = crs(env_stack_bsl))
    
    #-----------------------------------------------------------------
    #Test for collinearity among environmental layers
    #-----------------------------------------------------------------
    if(b==1){
      # https://github.com/N-SDM/covsel
      cov_data = pa_env %>% dplyr::select(-c(pres_abs, x, y)) #Select only environmental data
      covdata_filter = covsel.filteralgo(covdata = cov_data,
                                         pa = pa_env$pres_abs, #Select presence-background column
                                         corcut = 0.7)
      #Dormann et al. (2013) Collinearity
      # names(covdata_filter) #Show remaining variables
      # names(cov_data)[!names(cov_data) %in% names(covdata_filter)] #Show excluded variables
      
      #Create stack of only uncorrelated variables for buffer extent
      env_stack_bsl_buff_filter = env_stack_bsl_buff[[names(covdata_filter)]] 
      #Create stack of only uncorrelated variables for full extent
      env_stack_bsl_filter = env_stack_bsl[[names(covdata_filter)]] 
    }
    #Remove collinear variables from the pa_env data frame
    pa_env = pa_env %>% dplyr::select(c(pres_abs, x, y, names(env_stack_bsl_filter)))
    
    #----------------------------------------------------------------------
    # blockCV for getting spatially independent blocks for cross validation
    #----------------------------------------------------------------------
    # Use a fixed block size of 100 km
    block_size = 100000 
    if(b==1){
      try({
        num_folds = 5
        sb = cv_spatial(x = pa_sf,
                        column = "pres_abs",
                        r = env_stack_bsl_buff_filter[[1]],
                        k = num_folds,
                        size = block_size,  #size of the blocks
                        selection = "random",
                        iteration = 50,
                        biomod2 = F,
                        report = F,
                        plot = F,
                        seed = 1)
        fold_IDs = sb$folds_ids
      })
    }
    # Assign records to a fold
    if(!missing(sb)){
      pa_env$fold_IDs = fold_IDs
    } 
    
    # Check number of presences and absences in each fold
    (num_pa_folds = pa_env %>%
        group_by(pres_abs, fold_IDs) %>%
        dplyr::filter(pres_abs==1) %>%
        dplyr::summarise(n = n()))
    
    
    # If the cv_spatial fails or if there are too few occurrences per fold, 
    # use jackknife leave-one-out cross validation 
    #(generally < 25 localities) (Pearson et al. 2007; Shcheglovitova and Anderson 2013).
    #https://www.sciencedirect.com/science/article/pii/S0304380013004043
    # Jackknife instead of folds:
    do_jackknife = missing(sb) | nrow(pres)<25
    if(!missing(sb) & b==length(bsl_stacks)){rm(sb)}
    if(do_jackknife){
      pa_env$fold_IDs = pa_env$pres_abs
      pa_env$fold_IDs[pa_env$pres_abs==1] = 1:length(pa_env$fold_IDs[pa_env$pres_abs==1])
      num_folds = length(pa_env$fold_IDs[pa_env$pres_abs==1])
    }
    
    
    #----------------------------------------------------------------------
    # Try different regularisation parameters
    #----------------------------------------------------------------------
    for(reg_param in 1:length(reg_params)){
      try({
        #Arguments using prepPara
        args = prepPara(userfeatures=feats, outputformat='cloglog',
                        betamultiplier=reg_params[reg_param], replicates=1,
                        doclamp=T, removeduplicates=F, responsecurves=T, writeplotdata = F)
        args = c(args, 'noaddsamplestobackground', 'noplots')
        
        # Run a model for each cross-validation fold using selected fold to test the model
        # and the remaining data to train the model
        for(fold in 1:num_folds){
          
          #Get training and test data
          pa_train = pa_env %>% filter(fold_IDs!=fold)
          if(do_jackknife){
            pa_train$fold_IDs[pa_train$pres_abs==0] = fold #All background points are used
            # https://jamiemkass.github.io/ENMeval/reference/partitions.html
          }
          pa_test = pa_env %>% filter(fold_IDs==fold)
          if(do_jackknife){
            pa_test = pa_test %>%
              bind_rows(pa_env %>% filter(fold_IDs==0))
            pa_test$fold_IDs[pa_test$pres_abs==0] = fold
          }
          
          #Maxent
          mx = dismo::maxent(x = pa_train %>% dplyr::select(-c(pres_abs, x, y, fold_IDs)),
                             p = pa_train$pres_abs,
                             path = output_dir,
                             args = args,
                             silent = F)
          
          #Predict model results
          if(!any(is.na(mx@results))){
            prediction = rmaxent::project(mx, newdata = pa_test)$prediction_cloglog
            prediction = data.frame(prediction_cloglog = prediction) %>%
              bind_cols(pa_test %>% dplyr::select(pres_abs, fold_IDs))
            
            #Model evaluation
            e = dismo::evaluate(
              p = prediction %>% dplyr::filter(pres_abs==1) %>% pull(prediction_cloglog),
              a = prediction %>% dplyr::filter(pres_abs==0) %>% pull(prediction_cloglog)
            )
            boyce = evalContBoyce(
              pres = prediction %>% dplyr::filter(pres_abs==1) %>% pull(prediction_cloglog),
              contrast = prediction %>% dplyr::filter(pres_abs==0) %>% pull(prediction_cloglog)
            )
            
            #Summary results
            counter = counter + 1
            save(counter, file="Output/Maxent_models/counter.RData")
            results_list[[counter]] = data.frame(
              species=sp_name,
              bsl_stack=bsl_stacks[b],
              rcp=rcp,
              rcm=clim_mods[rcm],
              fold=fold,
              reg_param=reg_params[reg_param],
              np=e@np,
              na=e@na,
              auc=e@auc,
              boyce=boyce,
              thresh_max_sens_spec = dismo::threshold(e, stat = "spec_sens")) #Thresholding for binary predictions. SPARC: Binary map produced with both 5th percentile training presence and the calculated maximum sensitivity + specificity thresholds
            results_dat = as.data.frame(do.call(rbind, results_list))
            
            #Write results to file
            if(b==1){
              write.csv(results_dat, file = results_folds_file, row.names=F)
            } else {
              write.csv(read_csv(results_folds_file, show_col_types = F) %>% 
                          bind_rows(results_dat), 
                        file = results_folds_file, row.names=F)
            }
          } else {
            counter = counter + 1
          } #End  if(!any(is.na(mx@results)))
          
          print(paste0("sp = ",sp,". ",sp_name,", b = ", b, ", reg_param = ", reg_params[reg_param], ", fold = ", fold))
          
        } #End folds
      })
    } #End loop: different regularisation parameters
    
    #Find best regularisation parameter value
    results_RP = results_dat %>%
      group_by(reg_param, bsl_stack) %>%
      dplyr::summarise(boyce_mean = mean(boyce, na.rm = TRUE),
                       auc_mean = mean(auc, na.rm = TRUE)) 
    results_RP = results_RP[results_RP$bsl_stack==bsl_stacks[b],]
    reg_param_best = which(results_RP$auc_mean==max(results_RP$auc_mean, na.rm = T))[1]
    
    
    # Run Maxent using "optimal" settings
    #----------------------------------------------
    # Get optimal settings and put into prepPara format
    args_best = prepPara(userfeatures = feats, 
                         outputformat='cloglog',
                         betamultiplier = reg_params[reg_param_best], 
                         replicates = 1, doclamp=T, removeduplicates=F, responsecurves=F,
                         writeplotdata = F)
    args_best = c(args_best, 'noaddsamplestobackground', 'noplots')
    #Maxent
    mx_best = dismo::maxent(x = pa_env %>% dplyr::select(-c(pres_abs, x, y, fold_IDs)), 
                            p = pa_env$pres_abs, 
                            path = output_dir,
                            args = args_best,
                            silent = F)
    
    if(!any(is.na(mx_best@results))){
      #Project modelled distribution spatially
      mod_name = paste0(spp_maxent$species[sp],'_',rcp,'_bsl_rcm',rcm)
      pred_rast = dismo::predict(mx_best, env_stack_bsl, 
                                 filename = paste0(output_dir,'/',mod_name,'_fullExtent.tif'),
                                 progress = "", args = 'outputformat=Cloglog', overwrite = T)
      dismo::predict(mx_best, env_stack_bsl_buff, 
                     filename = paste0(output_dir,'/',mod_name,'_buffExtent.tif'),
                     progress = "", args = 'outputformat=Cloglog', overwrite = T)
      
      #Model evaluation
      e = dismo::evaluate(
        p = prediction %>% dplyr::filter(pres_abs==1) %>% pull(prediction_cloglog),
        a = prediction %>% dplyr::filter(pres_abs==0) %>% pull(prediction_cloglog)
      )
      
      
      #-----------------------------------
      #Get response curves data
      #-----------------------------------
      resp_dat = {} #Create empty object to store response curves data
      #Get min and max of each variable and create a vector for each
      mins_maxs = terra::minmax(env_stack_bsl_filter)
      var_names = c(names(env_stack_bsl_filter))
      min_max_dat = {}
      for(v in seq_along(var_names)){
        min_max_dat_var = data.frame(
          var_value = seq(mins_maxs[,var_names[v]][1], 
                          mins_maxs[,var_names[v]][2], length.out = 50))
        names(min_max_dat_var) = var_names[v]
        min_max_dat = min_max_dat %>% bind_cols(min_max_dat_var)
      }
      #Dataframe of modal values for each variable
      modal_dat = min_max_dat %>%
        mutate(across(everything(), \(x) radiant.data::modal(x, na.rm = TRUE)))
      
      for(v in seq_along(var_names)){
        vdat = modal_dat
        vdat[,names(vdat)%in%names(vdat)[v]] = min_max_dat[,v]
        suitability = rmaxent::project(mx_best, newdata = vdat)$prediction_cloglog
        
        resp_dat_var = bind_cols(data.frame(var = rep(var_names[v], nrow(vdat))),
                                 data.frame(var_value = vdat[,v], suitability = suitability))
        resp_dat = resp_dat %>% bind_rows(resp_dat_var)
      }
      write.csv(resp_dat,
                paste0(output_dir,"/Response_curves_data_",rcp,'_bsl_rcm',rcm,".csv"),
                row.names = F)
      
      #-----------------------------------
      # Save some useful information
      #-----------------------------------
      # Lower 95% CI calculation function
      ci95low = function(x){
        mean(x, na.rm = T) - 1.96 * sd(x, na.rm = T) / sqrt(length(x[!is.na(x)]))
      }
      # Get fold results of best regularisation parameter
      results_best_regparam = results_dat %>% 
        dplyr::filter(reg_param %in% reg_params[reg_param_best])
      
      # Get summary of results
      results_summary = data.frame(
        species = sp_name,
        rcp = rcp,
        bsl_stack = bsl_stacks[b],
        num_pres = nrow(pa_env %>% dplyr::filter(pres_abs==1)),
        num_pres_per_fold = paste(num_pa_folds$n, collapse = ";"),
        jackknife = do_jackknife,
        auc = mean(results_best_regparam$auc, na.rm = T),
        auc95 = ci95low(results_best_regparam$auc),
        boyce = mean(results_best_regparam$boyce, na.rm = T),
        boyce95 = ci95low(results_best_regparam$boyce),
        thresh_max_sens_spec = mean(results_best_regparam$thresh_max_sens_spec, na.rm = T)
      )
      #Add environmental predictors used
      if(b==1){
        results_summary = 
          results_summary %>%
          bind_cols(
            setNames(data.frame(
              matrix(ncol = length(names(env_stack_bsl)), nrow = 1)),
              names(env_stack_bsl))
          )
      }
      results_summary[,names(covdata_filter)] = 1
      #Write to CSV
      if(b==1){
        write.csv(results_summary, results_file, row.names = F)
      } else {
        write.csv(read_csv(results_file, show_col_types = F) %>% bind_rows(results_summary), 
                  results_file, row.names = F)
      }
      
      #--------------------------------------------------
      #Get Multivariate Environmental Suitability Surface
      #--------------------------------------------------
      mess_map = mess(x = stack(env_stack_bsl_filter),
                      v = pa_env %>% dplyr::select(-c(pres_abs,x,y,fold_IDs)),
                      filename =  paste0(output_dir,'/MESS_',rcp,'_bsl_rcm',rcm,'.tif'), 
                      overwrite = T)
      
      
      #--------------------------------------------------
      #Project to future climates
      #--------------------------------------------------
      # Project to 1.5 deg future
      mod_name = paste0(spp_maxent$species[sp],'_',rcp,'_1.5deg_rcm',rcm)
      pred_rast1.5 = dismo::predict(mx_best, env_stack_1.5, 
                                    filename = paste0(output_dir,'/',mod_name,'_fullExtent.tif'),
                                    progress = "", args = 'outputformat=Cloglog', overwrite = T)
      dismo::predict(mx_best, env_stack_1.5_buff, 
                     filename = paste0(output_dir,'/',mod_name,'_buffExtent.tif'),
                     progress = "", args = 'outputformat=Cloglog', overwrite = T)
      
      # Project to 2.0 deg future
      mod_name = paste0(spp_maxent$species[sp],'_',rcp,'_2deg_rcm',rcm)
      pred_rast2.0 = dismo::predict(mx_best, env_stack_2.0, 
                                    filename = paste0(output_dir,'/',mod_name,'_fullExtent.tif'),
                                    progress = "", args = 'outputformat=Cloglog', 
                                    overwrite = T)
      dismo::predict(mx_best, env_stack_2.0_buff, 
                     filename = paste0(output_dir,'/',mod_name,'_buffExtent.tif'),
                     progress = "", args = 'outputformat=Cloglog', overwrite = T)
      
      if(rcp=="rcp85" & rcm!=7){
        # Project to 4.0 deg future
        mod_name = paste0(spp_maxent$species[sp],'_',rcp,'_4deg_rcm',rcm)
        dismo::predict(mx_best, env_stack_4.0, 
                       filename = paste0(output_dir,'/',mod_name,'_fullExtent.tif'),
                       progress = "", args = 'outputformat=Cloglog', overwrite = T)
        dismo::predict(mx_best, env_stack_4.0_buff, 
                       filename = paste0(output_dir,'/',mod_name,'_buffExtent.tif'),
                       progress = "", args = 'outputformat=Cloglog', overwrite = T)
      } #End if(rcp=="rcp85" & rcm!=7)
      
    } #End if(!any(is.na(mx_best@results)))
  } #End bsl_stacks loop
} #End sp (species loop)
