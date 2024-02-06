###########################################
# MIMICS repeat run function for litterbag simulations 
###########################################

source("functions/MIMICS_sim_litterbag.R")
source("functions/MIMICS_calc_steady_state_pools.R")
source("functions/calc_Tpars.R")
source("functions/RXEQ.R")

MIMrepeat <- function(forcing_df, litBAG, rparams) { #need to remove litBAG if you want same litter at all sites rather than site specific litter
  
  
  # # Set global model parameters
  # Vslope <<- Vslope_default * rparams$Vslope_x[1]
  # Vint <<- Vint_default * rparams$Vint_x[1]
  # Kslope <<- Kslope_default * rparams$Kslope_x[1]
  # Kint <<- Kint_default * rparams$Kint_x[1]
   Tau_MULT <<- Tau_MULT_default * rparams$Tau_x[1] #change where Tau_MULT comes in in MIMICS_INC_daily to make this just influence r or K
  #  Tau_MULT.r <<- Tau_MULT_default * rparams$Tau_r[1]
  #  Tau_MULT.k <<- Tau_MULT_default * rparams$Tau_k[1]
  CUE <<- CUE_default * rparams$CUE_x[1]
  #  CUE <<- c(CUE_default[1] * rparams$CUE_x[1], CUE_default[2] * rparams$CUE_x[1], CUE_default[3], CUE_default[4])  #add indexing to just get r or K-selected
  #  CUE <<- c(CUE_default[1] * rparams$CUE_r[1], CUE_default[2] * rparams$CUE_r[1], CUE_default[3] * rparams$CUE_k[1], CUE_default[4] * rparams$CUE_k[1]) #seperate for r and K
   vMOD <<- vMOD_default * rparams$vMOD_x[1]
  # kMOD <<- kMOD_default * rparams$kMOD_x[1]
   aV <<- aV_default * rparams$aV_x[1]
  
  
  #below should go through each site with the same litter at all
  #MIMrun <- forcing_df %>% split(1:nrow(forcing_df)) %>% map(~MIMICS_LITBAG(forcing_df = ., litBAG = BAGS_BART, nspin_yrs=2,
  #                                                                         nspin_days=0, litadd_day=10)) %>% bind_rows() #_BART
  
  #below should vary litter and site in a coupled way (litterbag row coupled with site row)
  BAGS_input <- split(litBAG, 1:nrow(litBAG))
  forcing_input <- split(forcing_df, 1:nrow(forcing_df))
  MIMrun <- map2(forcing_input, BAGS_input, ~MIMICS_LITBAG(.x, .y, nspin_yrs=2, nspin_days=0, litadd_day=10, verbose=T)) %>% bind_rows()
  
  #Optional combine MIMout with forcing data
  #MIMrun <- MIMrun %>% left_join(forcing_df %>% select(-SITE), by="ID")
  
  #add run number
  MIMrun$run_num <- rparams$run_num[1]
  
  # return MIMrun
  return(MIMrun)
}


