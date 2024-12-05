###########################################
# MIMICS repeat run function for litterbag simulations 
###########################################

source("functions/MIMICS_sim_litterbag.R")
source("functions/MIMICS_calc_steady_state_pools.R")
source("functions/calc_Tpars.R")
source("functions/RXEQ.R")

MIMrepeat <- function(forcing_df, litBAG, dailyInput, rparams) { #need to remove litBAG if you want same litter at all sites rather than site specific litter
  
  
#Set global model parameters
   tau_r <<- c(tau_r_default[1], tau_r_default[2] * rparams$Tau_r[1]) #can multiply first part of tau_r to multiply total value or second part to infleunce multiplier on fMET
   vMOD <<- c(vMOD_default[1] * rparams$vMOD_m[1], vMOD_default[2] * rparams$vMOD_s[1], vMOD_default[3], vMOD_default[4] * rparams$vMOD_m[1], vMOD_default[5] * rparams$vMOD_s[1], vMOD_default[6]) #for just r or K
   beta <<- beta_default * rparams$beta_x[1]
   

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


  
  #add run number
  MIMrun$run_num <- rparams$run_num[1]
  
  # return MIMrun
  return(MIMrun)
}


