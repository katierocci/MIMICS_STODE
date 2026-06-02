### Ftn to calc MIMICS parameters (take inputs, return Tpars)
#----------------------------------------------------------------

# Load MIMICS parameters
source("set_MIMICS_params.R")

calc_Tpars <- function(TSOI, ANPP, CLAY, CN, LIG) {
  
  #debug
  # TSOI = 7
  # ANPP = 500
  # CLAY = 40
  # CN = 25
  # LIG = 15
  
  ANPP        <- ANPP/2
  fCLAY       <- CLAY/100
  lig_N = (LIG/100)/(1/(CN/2.5))
  fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N) 
  CN_s        <<- (CN-CN_m*fMET)/(1-fMET)
  
  
  # Calc litter input rate
  EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
  #print(EST_LIT)# gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 
  
  # ------------ caclulate parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV 
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
  Tau_MOD2 <- Tau_MOD[4]                        
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
  
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
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
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
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
  
  Tpars <- list(Inputs = I, VMAX = VMAX, KM = KM, CUE = CUE, 
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

#e.g.
t1 = calc_Tpars(TSOI = 7, ANPP = 500, CLAY = 40, CN =20, LIG = 15)
t2 = calc_Tpars(TSOI = 8, ANPP = 500, CLAY = 40, CN =20, LIG = 15)

str(t1)

t1['VMAX']
t2['VMAX']
