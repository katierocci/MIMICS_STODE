###########################################
# MIMICS repeat run function for litterbag simulations 
###########################################

source("functions/MIMICS_sim_litterbag.R")
source("functions/MIMICS_calc_steady_state_pools.R")
source("functions/calc_Tpars.R")
source("functions/RXEQ.R")

MIMrepeat <- function(forcing_df, litBAG, dailyInput, rparams) { #need to remove litBAG if you want same litter at all sites rather than site specific litter
  
  
  # # Set global model parameters
  # Vslope <<- Vslope_default * rparams$Vslope_x[1]
  # Vint <<- Vint_default * rparams$Vint_x[1]
  # Kslope <<- Kslope_default * rparams$Kslope_x[1]
  # Kint <<- Kint_default * rparams$Kint_x[1]
  # Tau_MULT <<- Tau_MULT_default * rparams$Tau_x[1] 
  #  tau_r <<- c(tau_r_default[1], tau_r_default[2] * rparams$Tau_r[1]) #can multiply first part of tau_r to multiply total value or second part to infleunce multiplier on fMET
  #  tau_K <<- c(tau_K_default[1], tau_K_default[2] * rparams$Tau_K[1]) #can multiply first part of tau_r to multiply total value or second part to infleunce multiplier on fMET
  #  CUE <<- CUE_default * rparams$CUE_x[1]
  #  CUE <<- c(CUE_default[1], CUE_default[2], CUE_default[3] * rparams$CUE_x[1], CUE_default[4] * rparams$CUE_x[1])  #add indexing to just get r or K-selected
  #  CUE <<- c(CUE_default[1] * rparams$CUE_r[1], CUE_default[2] * rparams$CUE_r[1], CUE_default[3] * rparams$CUE_k[1], CUE_default[4] * rparams$CUE_k[1]) #seperate for r and K
  #  vMOD <<- vMOD_default * rparams$vMOD_x[1]
   vMOD <<- c(vMOD_default[1] * rparams$vMOD_m[1], vMOD_default[2] * rparams$vMOD_s[1], vMOD_default[3], vMOD_default[4] * rparams$vMOD_m[1], vMOD_default[5] * rparams$vMOD_s[1], vMOD_default[6]) #for just r or K
   #kMOD <<- kMOD_default * rparams$kMOD_x[1]
  # kMOD <<- c(kMOD_default[1] * rparams$kMOD_x[1], kMOD_default[2] * rparams$kMOD_x[1], kMOD_default[3] * rparams$kMOD_x[1], kMOD_default[4], kMOD_default[5], kMOD_default[6]) #for just r or K
  # aV <<- aV_default * rparams$aV_x[1]
 # fmet_p <<- c(fmet_p_default[1], fmet_p_default[2] * rparams$fM_x[1], fmet_p_default[3])
#  beta <<- beta_default * rparams$beta_x[1]
   #beta <<- c(beta_default[1] * rparams$beta_x[1], beta_default[2]) #add indexing for r r K-selected
   beta <<- c(beta_default[1] * rparams$beta_r[1], beta_default[2] * rparams$beta_k[1]) #changing r and K at teh same time
  #below should go through each site with the same litter at all
  #MIMrun <- forcing_df %>% split(1:nrow(forcing_df)) %>% map(~MIMICS_LITBAG(forcing_df = ., litBAG = BAGS_BART, nspin_yrs=2,
  #                                                                         nspin_days=0, litadd_day=10)) %>% bind_rows() #_BART
  
   #STEADY STATE
  #below should vary litter and site in a coupled way (litterbag row coupled with site row)
   # BAGS_input <- split(litBAG, 1:nrow(litBAG))
   # forcing_input <- split(forcing_df, 1:nrow(forcing_df))
   # MIMrun <- map2(forcing_input, BAGS_input, ~MIMICS_LITBAG(.x, .y, nspin_yrs=2, nspin_days=0, litadd_day=10, verbose=T)) %>% bind_rows()
   
   
   #DAILY INPUT
   #below should vary litter and site in a coupled way (litterbag row coupled with site row) with daily climate input
   # MIMrun = data.frame()
   #     for (site in Mic_sites) {
   #       BAGS_input <- filter(litBAG, SITE == site) #changed to all caps for new bag input
   #       forcing_input <- filter(forcing_df, SITE == site)
   #       daily_input <- filter(dailyInput, SITE == site)
   #       BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
   #       MIMrun <- rbind(MIMrun,BO_DI)
   #     }
   

   #DAILY INPUT WITH VARYING SM AND LQ
  MIMrun = data.frame()
  SM = c("mean", "max", "min")
  for (SM_type2 in SM) {
    data_in <- filter(forcing_df, SM_type==SM_type2)
    DailyInput_in <- filter(dailyInput, SM_type==SM_type2)
    LQ = c("LIG_N_sp1", "LIG_N_sp2", "LIG_N_sp3")
    for (bag_type in LQ) {
      BAGS_mean <- filter(litBAG, TYPE==bag_type)
      data_in$LIG_N = data_in[[bag_type]]
        for (site in Mic_sites) {
        BAGS_input <- filter(BAGS_mean, SITE == site) #changed to all caps for new bag input
        forcing_input <- filter(data_in, SITE == site)
        daily_input <- filter(DailyInput_in, SITE == site)
        BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
        MIMrun <- rbind(MIMrun,BO_DI)
        }
    }
  }


  
  #Optional combine MIMout with forcing data
  #MIMrun <- MIMrun %>% left_join(forcing_df %>% select(-SITE), by="ID")
  
  #add run number
  MIMrun$run_num <- rparams$run_num[1]
  
  # return MIMrun
  return(MIMrun)
}


