library(tidyr)
library(dplyr)

source("Parameters/MIMICS_parameters_sandbox_20231129.R")

#input data
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim_final.csv")
#match input data strucutre
#AGNPP should be in gram dry weight! multiply by 2 here to remedy
#if using soil moisture and not water scalar: we are using VWC and not GWC - assuming a BD of 1g/cm3 makes them equivalent but this assumption may be flawed
data1 <- MSBio %>% mutate(SITE = Site, ANPP = AGNPP_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, lig_N = LIG_N, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG_N, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3, GWC, W_SCALAR, lci_SM_ratio, uci_SM_ratio)
#replacing anomalously high ANPP values for OSBS and TALL based on NEON measurements
MSBio2 <- data1
NEON_GPP <- read.csv("Example_simulations/Data/NEON_GPP.csv")
MSBio2$ANPP[MSBio2$SITE == "OSBS"] <- 547 + 0.18*NEON_GPP[6,2] #using relationship between NEON GPP and ANPP
MSBio2$ANPP[MSBio2$SITE == "TALL"] <- 547 + 0.18*NEON_GPP[9,2] #using relationship between NEON GPP and ANPP 
#filtering for sites with microbial data
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
data_sites <- filter(MSBio2, SITE %in% Mic_sites)
data <- data_sites

#loading daily inputs and replacing TALL data to be more realistic
DailyInput <- read.csv("Example_simulations/Data/DailyInput.csv")
DailyInput$LITFALL[DailyInput$SITE == "TALL"] <- DailyInput$LITFALL[DailyInput$SITE == "TALL"]*0.663
DailyInput$ANPP[DailyInput$SITE == "TALL"] <- sum(DailyInput$LITFALL[DailyInput$SITE == "TALL"])

#formatting data for multiple soil moisture types
data_SM <- rbind(data, data, data)
#below creates water scalar over 1 so maybe need to change all maxes where W_SCALAR over 1 is equal to 1? Mathematically, fine to go over 1....
data_SM <- data_SM %>% mutate(SM_type = c(rep("mean", 7), rep("max", 7), rep("min", 7))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
                               SM_type == "max" ~ W_SCALAR*uci_SM_ratio,
                               SM_type == "min" ~ W_SCALAR*lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)
DailyInput_SM <- rbind(DailyInput, DailyInput, DailyInput)
SM_mult <- data %>% select(SITE, uci_SM_ratio, lci_SM_ratio)
DailyInput_SM <- DailyInput_SM %>% left_join(SM_mult, by="SITE") %>% mutate(SM_type = c(rep("mean", 2561), rep("max", 2561), rep("min", 2561))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
                               SM_type == "max" ~ W_SCALAR *uci_SM_ratio,
                               SM_type == "min" ~ W_SCALAR *lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)

#load in MSBio litter bag chemistry
### changed to fMET calculation in STODE script here!! Note that the two options are only somewhat related but less negatives in STODE equation
MSBio_BAGS <- data %>% select(SITE, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3) %>% pivot_longer(2:4, names_to = "TYPE", values_to = "BAG_LIG_N")
MSBio_BAGS$CALC_MET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * (MSBio_BAGS$BAG_LIG_N))
MSBio_BAGS$CALC_MET[MSBio_BAGS$CALC_MET <0] = 0.01 #setting negatives to small number so 99% structural

BAG_init_size <- 100
BAGS <- MSBio_BAGS %>% select(SITE, TYPE, CALC_MET)
BAGS$BAG_LITm <- ((BAG_init_size * 1e3 / 1e4)/ depth) * BAGS$CALC_MET #g/m2 converted to mg/cm3
BAGS$BAG_LITs <- ((BAG_init_size * 1e3 / 1e4)/ depth) * (1-BAGS$CALC_MET) 
BAGS_mean <- filter(BAGS, TYPE == "LIG_N_sp1")
#initial litter = 0.33 because of unit conversions here

####
#bringing input data together with output data
####
#df <- readRDS("C:/github/MIMICS_MSBio/Cheyenne_HPC/HPC_output/MSBio_MIM_MC_runs-1e+05_20221028_112647_.rds")
#MC_MIMICS1 <- readRDS('/glade/work/krocci/MC_output/MSBio_MC_500_20240425_184743_.rds')
MC_MIMICS2 <- readRDS('MC_output/MSBio_MC_500_20240604_210307_.rds') #/glade/work/krocci
MC_MIMICS3 <- readRDS('MC_output/MSBio_MC_500_20240604_210634_.rds') #/glade/work/krocci
MC_MIMICS4 <- readRDS('MC_output/MSBio_MC_500_20240604_211438_.rds') #/glade/work/krocci
MC_MIMICS <- rbind(MC_MIMICS2, MC_MIMICS3, MC_MIMICS4) #MC_MIMICS1, 
df <- left_join(MC_MIMICS, data, by = "SITE") 
#add unique run numbers when combining multiple derecho runs
df$run_num2 <- df$run_num
df$run_num <- rep(1:1500, each=68985)

df$MIM_CO <- as.numeric(df$MICr)/as.numeric(df$MICk)
df$MIC_SOC <- (df$MICr+df$MICk)/(df$SOMc+df$SOMa+df$MICr+df$MICk)
df$LITBAG_tot <- df$LITBAGm + df$LITBAGs

#logical checks
df_LC <- df %>%
  filter(MIM_CO > 0.01) %>%
  filter(MIM_CO < 100) %>%
  filter(MIC_SOC > 0.0001) %>%
  filter(MIC_SOC < 0.40) 

LC.rn <- df_LC %>% group_by(run_num) %>% mutate(S.LQ.SM = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% 
summarize(uniq.site = length(unique(S.LQ.SM))) %>% filter(uniq.site>62)
df_LC2 <- df_LC %>% filter(run_num %in% LC.rn$run_num)

####
#data for litter mass loss cost function
####
Field_LML <- read.csv("Example_simulations/Data/Litter_decomp_all.csv")
LML_sum2 <- Field_LML  %>% group_by(site, time.point) %>% drop_na(percent.loss.litter) %>% summarize(mean.ML = mean(percent.loss.litter*100),
            n.ML = n(),
            sd.ML = sd(percent.loss.litter*100),
            SE = sd.ML/sqrt(n.ML),
            lci.ML = mean.ML - qt(1 - ((1 - 0.95) / 2), n.ML - 1) * SE,
            uci.ML = mean.ML + qt(1 - ((1 - 0.95) / 2), n.ML - 1) * SE,
            min.ML = min(percent.loss.litter*100),
            max.ML = max(percent.loss.litter*100),
            doy = mean(days_elapsed)) %>% mutate(doy=round(doy, digits=0))
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML, lci.ML, uci.ML, min.ML, max.ML)
###
#prepping data for analysis
###
LITi = 0.1
df_LML <- df_LC2  %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100) %>% #%>% mutate(SITE.rn = paste(df$SITE, df$run_num, sep = "")) 
  mutate(field.day=DAY-314) %>% mutate(SITE.DAY=paste(SITE, field.day, sep=".")) %>% right_join(FieldData, by="SITE.DAY")
DI_means <- DailyInput_SM %>% mutate(SITE.SM = paste(SM_type, SITE,sep = ".")) %>% group_by(SITE.SM) %>% 
  summarise(W_SCALAR_mean=mean(W_SCALAR), MAT_mean=mean(MAT)) %>% select(SITE.SM, W_SCALAR_mean, MAT_mean) # SM_type, 
#initial MICrK
MIC_init <- df %>% filter(DAY == 315) %>% mutate(MICrK.i =  MICr/MICk)%>%
  mutate(SITE.SM.LQ.rn = paste(SITE, SM_Type, Litter_Type, run_num, sep = ".")) %>% select(SITE.SM.LQ.rn, MICrK.i) #SM_Type, 
#bag means
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
df_analysis <- df_LML %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ.rn = paste(SITE, SM_Type, Litter_Type, run_num, sep = ".")) %>% mutate(SITE.SM = paste(SM_Type, SITE, sep = ".")) %>% # SM_Type,  SM_Type, 
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% 
  inner_join(MIC_init, by="SITE.SM.LQ.rn") %>% inner_join(BAGS_LIGN, by="SITE.LQ")

###
#doing initial filtering
###
#logical checks
df_check <- df_analysis 
df_rp <- df_check %>% select(run_num, Tau_r, beta_x, CUE_x, vMOD_m, vMOD_s) %>% group_by(run_num) %>% # Tau_K,
  summarise(Tau_r = mean(Tau_r), beta_x = mean(beta_x), CUE_x = mean(CUE_x), vMOD_m = mean(vMOD_m), vMOD_s = mean(vMOD_s)) #Tau_K = mean(Tau_K), 
#filter for LML
#collect run numbers that fit percent loss from observations at time points 1 and 2 - frequency might need to change to account for 9 reps of each site
df_check$SITE <- as.factor(df_check$SITE)
df_rn1 <- df_check %>% filter(time.point == 1) %>% filter(LIT_PerLoss > min.ML & LIT_PerLoss < max.ML) #highest and lowest values of confidence interval for time point 1
LR.rn1 <- df_rn1 %>% group_by(run_num) %>% mutate(S.LQ.SM = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% 
  summarize(uniq.site = length(unique(S.LQ.SM))) #%>% filter(uniq.site>62) 
df_rn2 <- df_check %>% filter(time.point == 2) %>% filter(LIT_PerLoss > min.ML & LIT_PerLoss < max.ML) #highest and lowest values of confidence interval for time point 2
LR.rn2 <- df_rn2 %>% group_by(run_num) %>% mutate(S.LQ.SM = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% 
  summarize(uniq.site = length(unique(S.LQ.SM))) 
#keep only run numbers that are in both dfs
df_rn1.v <- as.vector(LR.rn1$run_num) #as.numeric(LR.rn1$Var1)
df_rn2.v <- as.vector(LR.rn2$run_num) #as.numeric(LR.rn2$Var1)
#rn_final <- which(df_rn2.v %in% df_rn1.v) #this stopped working...
rn_final <- intersect(df_rn1.v,df_rn2.v)
#filter check df to only have run numbers that fit litter mass loss
df_check2 <- df_check %>% filter(run_num %in% rn_final)
#transforming and scaling data for analysis
df_TS <- df_check2 %>% mutate(log_WS = log(W_SCALAR_mean)) %>%select(run_num, SITE, LIT_PerLoss, log_WS, BAG_LIG_N, MICrK.i) %>% 
  mutate_at(vars(c("log_WS", "BAG_LIG_N", "MICrK.i")), ~(scale(.) %>% as.vector))

saveRDS(object=df_TS, file='MC_output/MSBio_df_TS_Strict2.rds')
saveRDS(object=df_rp, file='MC_output/MSBio_df_RP_Strict2.rds')
