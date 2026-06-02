# MIMICS-CN Rev-MM Eq.

CN_RXEQ <- function(t, y, pars) {
  
  #print(pars['VMAX1'])
  
  with (as.list(c(y, pars)),{
    
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + MIC_1)
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + MIC_1)
    MICtrn[1] = MIC_1^densDep * tau[1] * fPHYS[1] 
    MICtrn[2] = MIC_1^densDep * tau[1] * fCHEM[1] 
    MICtrn[3] = MIC_1^densDep * tau[1] * fAVAI[1] 
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + MIC_1)
    
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + MIC_2) 
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + MIC_2)
    MICtrn[4] = MIC_2^densDep * tau[2] * fPHYS[2]
    MICtrn[5] = MIC_2^densDep * tau[2] * fCHEM[2] 
    MICtrn[6] = MIC_2^densDep * tau[2] * fAVAI[2]     
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + MIC_2)
    
    DEsorb    = SOM_1 * desorb
    
    OXIDAT    = ((MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + MIC_1)) +
                   (MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + MIC_2)))
    
    upMIC_1    = CUE[1]*(LITmin[1] + SOMmin[1]) + CUE[2]*(LITmin[2])
    upMIC_2    = CUE[3]*(LITmin[3] + SOMmin[2]) + CUE[4]*(LITmin[4])   
    
    dLIT_1 = Inputs[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_1 = upMIC_1 - (MICtrn[[1]] + MICtrn[[2]] + MICtrn[[3]])# - Overflow[1]    
    dSOM_1 = Inputs[1]*FI[1] + MICtrn[1] + MICtrn[4] - DEsorb     
    
    dLIT_2 = Inputs[2]*(1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_2 = upMIC_2 - (MICtrn[[4]] + MICtrn[[5]] + MICtrn[[6]])# - Overflow[2]    
    dSOM_2 = Inputs[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
    
    dSOM_3 = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    #------ N cycle ---------
    
    MICr_recip = 1.0 / (MIC_1 + 1e-10)
    LITminN[1] =  LITmin[1]*LIT_1_N/(LIT_1 + 1e-10)
    LITminN[2] =  LITmin[2]*LIT_2_N/(LIT_2 + 1e-10)
    MICtrnN[1] =  MICtrn[1]*MIC_1_N * MICr_recip
    MICtrnN[2] =  MICtrn[2]*MIC_1_N * MICr_recip
    MICtrnN[3] =  MICtrn[3]*MIC_1_N * MICr_recip
    SOMminN[1] =  SOMmin[1]*(SOM_3_N/(SOM_3 + 1e-10))
  
    MICk_recip = 1.0 / (MIC_2 + 1.0e-10)
    LITminN[3] =  LITmin[3]*LIT_1_N/(LIT_1 + 1e-10)
    LITminN[4] =  LITmin[4]*LIT_2_N / (LIT_2 + 1e-10)
    MICtrnN[4] =  MICtrn[4] * MIC_2_N * MICk_recip
    MICtrnN[5] =  MICtrn[5] * MIC_2_N * MICk_recip
    MICtrnN[6] =  MICtrn[6] * MIC_2_N * MICk_recip
    SOMminN[2] =  SOMmin[2]*(SOM_3_N/(SOM_3 + 1e-10))
  
    DEsorbN    =  DEsorb[1]*(SOM_1_N/(SOM_1 + 1e-10))
  
    OXIDATN    =  OXIDAT[1] * (SOM_2_N/(SOM_2 + 1e-10))

  #! Partitions available DIN between microbial pools based on relative biomass
  #! The amount of DIN uptake is not demand driven.  It is a fraction of a large pool. -mdh 7/29/2019

    #! Assume that microbes take up all mineral N that is available. The unneeded portion will be spilled later.
    DINup[1]   = DIN*MIC_1/(MIC_1+MIC_2+ 1e-10) #!!! Emily's code first multiplies by (1-Nleak)
    DINup[2]   = DIN*MIC_2/(MIC_1+MIC_2+ 1e-10) #!!! Emily's code first multiplies by (1-Nleak)
  
    upMIC_1_N  = NUE[1]*(LITminN[1] + SOMminN[1]) + NUE[2]*(LITminN[2]) + DINup[1] #!!! Emily's code only has one NUE value
  
    CNup[1]    = (upMIC_1)/(upMIC_1_N + 1e-10)
    Overflow[1] = (upMIC_1) - (upMIC_1_N)*min(CN_r, CNup[1])
    Nspill[1]   = (upMIC_1_N) - (upMIC_1)/max(CN_r, CNup[1])
  
  # #! Add Overflow_r and Overflow_k to output netCDF file. Units conversion occurs in subroutine mimics_caccum. -mdh 10/12/2020
  # #mimicsflux%Overflow_r(npt) = mimicsflux%Overflow_r(npt) + Overflow_r
  
    upMIC_2    = CUE[3]*(LITmin[3] + SOMmin[2]) + CUE[4]*(LITmin[4])
    upMIC_2_N  = NUE[3]*(LITminN[3] + SOMminN[2]) + NUE[4]*(LITminN[4]) + DINup[2]
  
    CNup[2]    = (upMIC_2)/(upMIC_2_N+1e-10)
    Overflow[2] = (upMIC_2) - (upMIC_2_N)*min(CN_K, CNup[2])
    Nspill[2]   = (upMIC_2_N) - (upMIC_2)/max(CN_K, CNup[2])
    
    #print(Overflow)
  
  # #! Add Overflow_r and Overflow_k to output netCDF file. Units conversion occurs in subroutine mimics_caccum. -mdh 10/12/2020
  # #mimicsflux%Overflow_k(npt) = mimicsflux%Overflow_k(npt) + Overflow_k
  
    dLIT_1_N = Inputs[1]*(1-FI[1])/CN_m - LITminN[1] - LITminN[3]
    dMIC_1_N = upMIC_1_N - (MICtrnN[[1]] + MICtrnN[[2]] + MICtrnN[[3]]) - Nspill[1]
    dSOM_1_N = Inputs[1]*FI[1]/CN_m + MICtrnN[1] + MICtrnN[4] - DEsorbN #!!! Dividing by CN-m is not in fortran code. Possibly using different N input? ("NlitInput")
  
    dLIT_2_N = Inputs[2]*(1-FI[2])/CN_s - LITminN[2] - LITminN[4]
    dMIC_2_N = upMIC_2_N - (MICtrnN[[4]] + MICtrnN[[5]] + MICtrnN[[6]]) - Nspill[2]
    dSOM_2_N = Inputs[2]*FI[2]/CN_s + MICtrnN[2] + MICtrnN[5] - OXIDATN
  
    dSOM_3_N = MICtrnN[3] + MICtrnN[6] + DEsorbN + OXIDATN - SOMminN[1] - SOMminN[2]
  
  #!!! CHECK AGAINST TESTBED CODE, THIS IS WRONG...?
    # dDIN = (1-NUE[1])*(LITminN[1] + SOMminN[1]) + (1-NUE[2])*(LITminN[2]) +
    #        (1-NUE[3])*(LITminN[3] + SOMminN[2]) + (1-NUE[4])*(LITminN[4]) +
    #        Nspill[1] + Nspill[2] - DINup[1] - DINup[2]
  
    dDIN = (1-NUE[1])*(LITminN[1] + LITminN[2] + SOMminN[1]) +  #Inputs from r decomp
           (1-NUE[3])*(LITminN[3] + LITminN[4] + SOMminN[2]) +  #Inputs from K decomp
           Nspill[1] + Nspill[2] - DINup[1] - DINup[2]    #Uptake to microbial pools and spillage
    
    LeachingLoss = Nleak*DIN
    dDIN = dDIN-LeachingLoss #N leaching losses
    
    # account for overflow respiration fluxes
    #KR added this to the code 10/24/25 from CN_RXEQ function in SoilStoich repo
    if (!is.na(Overflow[1]) & Overflow[1] > 1e-10) {
      dMIC_1 = dMIC_1 - Overflow[1]
    }
    if (!is.na(Overflow[2]) & Overflow[2] > 1e-10) {
      dMIC_2 = dMIC_2 - Overflow[2]
    }
    
    list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3, dLIT_1_N, dLIT_2_N, dMIC_1_N, dMIC_2_N, dSOM_1_N, dSOM_2_N, dSOM_3_N, dDIN))
  })
}