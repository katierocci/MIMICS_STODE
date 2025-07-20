#----------------------------------------------------------------
### Ftn to calc MIMICS parameters (take inputs, return Tpars)
#----------------------------------------------------------------

# Load MIMICS parameters
#source("Parameters/MIMICS_parameters_sandbox_20231129.R")

# Assumes:
# -- ANPP is gC, not gDW
# -- Clay is a fraction (0-1)
# -- fWmethod is describes how to calculate fW 0=none 1=corpse, 2=calibrated, 3=moisture scalar from other model
# -- historic is a logical for using historic MAT to modify Vslope & Vint

calc_Tpars_Conly <- function(ANPP, fCLAY, TSOI, MAT=NA, CN, LIG, LIG_N=NA, AlFe=NA,
                             theta_liq=NA, theta_frzn=NA, W_SCALAR=NA, litfall=NA) {
  
  # Set lig:N value if not given
  if (is.na(LIG_N)) {
    LIG_N <- (LIG/100)/(1/(CN/2.5))
  }

  # Set moisture control on kinetics (fW)
  if (fWmethod==0) {
    fW=1
  } else if (fWmethod==1) {

    air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
    fW = (theta_liq^3 * air_filled_porosity^2.5)/0.022600567942709
    fW = max(0.05, fW)
    
  } else if (fWmethod==2) {
    
    air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
    f <- function(x, p1, p2) {x^p1 * (1-x)^p2}
    fW_p3 <- optimize(f, interval=c(0.01,1), p1=fW_p1, p2=fW_p2, maximum = T)$objective
    fW = (theta_liq^fW_p1 * air_filled_porosity^fW_p2)/fW_p3
    fW = max(0.05, fW) 
    
  } else if (fWmethod==3) {
    fW=W_SCALAR
  }
  
  # For historic MAT dependent kinetics
  if (historic==TRUE) {
    Vslope = Vslope + (MAT*0.00104) 
    Vint = Vint - (MAT*0.0228) 
  }
  

  # set fMET
  if(!fixed_fMET){
    #fMET  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * LIG_N)  
    fMET  <- ifelse(LIG_N<65, (fmet_p[1] * (fmet_p[2] - fmet_p[3] * LIG_N)), 0.01) #to account for high LIG_N that cause negtaives
    } else {
    # Fixed Fmet, ONLY for model validation against MIMICS et al. 2015 Sandbox 
    fMET <- 0.3846423
  }
  
  # Calc litter input rate from annual or daily flux then convert units
  if (is.na(litfall)) {
    EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
  } else {
    EST_LIT <- (litfall / 24) * 1e3 / 1e4
  }
  
  # ------------ calculate time varying parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV * fW   #<-- Moisture scalar applied
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  if (tauMethod == 'NPP') {
    Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
    Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
    Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
    beta=c(1,1) # turns off beta
  } else if (tauMethod=='beta') {
    Tau_MOD1 <- 1 # turns off NPP effects on turnover
    beta=beta #[1] # use Kat's density dependent function  
    #add in another parameter that adds 1 to 
  }
  
  Tau_MOD2 <- Tau_MOD[4]                        
  
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
  tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT
  #print(tau)
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  
  if (desorbMethod == 'clay') {
    desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))
  } else if (desorbMethod=='metal') {
    fMETAL <- ifelse(AlFe<6367, AlFe/6366.6, 1)
    desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fMETAL))
  }
                 
  desorb   <- desorb * desorb_MULT
  fPHYS    <- fPHYS * fPHYS_MULT
  
  if (psMethod == 'clay') {
    pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  } else if (psMethod=='metal') {
    pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fMETAL)))  #Scalar for texture effects on SOMp
  }
  
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
  I       <- array(NA, dim=2)             
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  Inputs  <- I
  
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- c(NA,NA,NA,NA,NA,NA)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  Tpars <- list(I = I, VMAX = VMAX, KM = KM, CUE = CUE,
                 fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                 tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                 desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO, 
                 beta=beta)
  
  return(Tpars)
}


#-------------------------------------------------------------------------
## as Inputs## as above, but for CN code
#source("Parameters/MIMICS-CN/set_MIMICS_params_CN.R")

calc_Tpars_CN <- function(TSOI, ANPP, CLAY, CN, LIG, x, fW=1,nUPmod=1) {
  
  ANPP        <- ANPP/2
  fCLAY       <- CLAY/100
  lig_N = (LIG/100)/(1/(CN/2.5))
  fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N) 
  CN_s        <<- (CN-CN_m*fMET)/(1-fMET)
  CN_r = CN_r  * sqrt(cnModNum/fMET)
  CN_K = CN_K * sqrt(cnModNum/fMET)
  
  
  # Calc litter input rate
  EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
  #print(EST_LIT)# gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 
  
  # ------------ calculate time varying parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV 
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
  Tau_MOD2 <- Tau_MOD[4]                        
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
  
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), #*fMET
           tau_K[1]*exp(tau_K[2]*fMET)) #*fMET  
  tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT 
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  
  
  desorb <- desorb * desorb_MULT
  fPHYS <- fPHYS * fPHYS_MULT
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD * fW
  KM       <- Km / k_MOD
    
  I       <- array(NA, dim=3)             
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  I[3]    <- 0  
  Inputs <- I
  
  # initialize pools with small values
  lit     <- array(I[1], dim=2)      
  mic     <- array(I[1], dim=2)   
  som     <- array(I[1], 3) 
  

  LITmin  <- rep(NA, dim=4)
  MICtrn  <- c(NA,NA,NA,NA,NA,NA)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  LITminN   <<- array(NA, dim=4)
  MICtrnN   <<- array(NA, dim=6)
  SOMminN   <<- array(NA, dim=2)
  DEsorbN   <<- array(NA, dim=1)
  OXIDATN   <<- array(NA, dim=1)
  
  DINup     <<- array(NA, dim=2)
  Overflow  <<- array(NA, dim=2)
  Nspill    <<- array(NA, dim=2)
  CNup      <<- array(NA, dim=2)
  upMIC_1   <<-  array(NA, dim=1)
  upMIC_1_N <<-  array(NA, dim=1)
  upMIC_2   <<-  array(NA, dim=1)
  upMIC_2_N <<-  array(NA, dim=1)
  
  Nleak = Nleak * nUPmod
  
  Tpars <- list( Inputs = I, VMAX = VMAX, KM = KM, CUE = CUE, 
               fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
               tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
               desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO,
               LITminN = LITminN, SOMminN = SOMminN, MICtrnN = MICtrnN,
               DEsorbN = DEsorbN, OXIDATN = OXIDATN, densDep=densDep,
               CNup=CNup, DINup=DINup, Nspill=Nspill, Overflow=Overflow, 
               upMIC_1=upMIC_1, upMIC_1_N=upMIC_1_N,
               upMIC_2=upMIC_2, upMIC_2_N=upMIC_2_N,
               NUE=NUE, CN_m=CN_m, CN_s=CN_s, CN_r=CN_r, CN_K=CN_K, Nleak=Nleak)
  
  
  return(Tpars)
}