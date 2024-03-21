#library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(gridExtra)
library(car)
library(lmerTest)
library(rwa)

#load site data
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim.csv")
data <- MSBio %>% mutate(SITE = Site, ANPP = AGNPP_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, lig_N = LIG_N, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG, C, N, CN, LIG_N, GWC, W_SCALAR)

#bringing input data together with output data
#df <- readRDS("C:/github/MIMICS_MSBio/Cheyenne_HPC/HPC_output/MSBio_MIM_MC_runs-1e+05_20221028_112647_.rds")
# MC_MIMICS1 <- readRDS('Analysis/MC_output/MSBio_MC_3000_20240221_142952_.rds')
# MC_MIMICS2 <- readRDS('Analysis/MC_output/MSBio_MC_3000_20240221_143057_.rds')
# MC_MIMICS3 <- readRDS('Analysis/MC_output/MSBio_MC_3000_20240221_143328_.rds')
MC_MIMICS1 <- readRDS('MC_output/MSBio_MC_3000_20240221_142952_.rds')
MC_MIMICS2 <- readRDS('MC_output/MSBio_MC_3000_20240221_143057_.rds')
MC_MIMICS3 <- readRDS('MC_output/MSBio_MC_3000_20240221_143328_.rds')
df <- rbind(MC_MIMICS1, MC_MIMICS2, MC_MIMICS3)
df$run_num2 <- df$run_num
df$run_num <- rep(1:9000, each=8030) #unique run number for each day of each site (730*11)
#rand_params <- readRDS('Cheyenne_HPC/July 2023/MIMICS_MSBio_KR/MC_output/MSBio_RP_1e+05_20230728_132810_Tc.rds')
df <- left_join(df, data, by = "SITE") #MC output has to be first for CO2 rows below to be right
#df <- left_join(df, rand_params, by = "run_num") #don't have random parameters )':
#df$CO2_of_tot <- rowSums(df[,11:12])/rowSums(df[,4:12])
df$MIM_CO <- as.numeric(df$MICr)/as.numeric(df$MICk)
df$MIC_SOC <- (df$MICr+df$MICk)/(df$SOMc+df$SOMa+df$MICr+df$MICk)
df$LITBAG_tot <- df$LITBAGm + df$LITBAGs




###
#cost functions
###

#data for litter mass loss cost function
Field_LML <- read.csv("Example_simulations/Data/Litter_decomp_all.csv")
LML_sum2 <- Field_LML  %>% group_by(site, time.point) %>% drop_na(percent.loss.litter) %>% summarize(mean.ML = mean(percent.loss.litter*100),
                                                                                                     n.ML = n(),
                                                                                                     sd.ML = sd(percent.loss.litter*100),
                                                                                                     SE = sd.ML/sqrt(n.ML),
                                                                                                     lci.ML = mean.ML - qt(1 - ((1 - 0.95) / 2), n.ML - 1) * SE,
                                                                                                     uci.ML = mean.ML + qt(1 - ((1 - 0.95) / 2), n.ML - 1) * SE,
                                                                                                     doy = mean(days_elapsed)) %>% mutate(doy=round(doy, digits=0))
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML, sd.ML, n.ML)
LIT_init <- df %>% mutate(SITE.rn = paste(df$SITE, df$run_num, sep = "")) %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% select(SITE.rn, LITi)
df_LML <- df %>% mutate(SITE.rn = paste(df$SITE, df$run_num, sep = "")) %>% left_join(LIT_init, by = "SITE.rn")
df_LML <- df_LML %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
# #data for microbial functional group cost function
# #r:K/C:O
# MSBio_micfg <- read.csv("Example_simulations/Data/MSBio_FuncGroups_AllSites_prelim.csv")
# df <- MSBio_micfg %>% filter(material == "soil") %>% mutate(SITE = site) %>% group_by(SITE)  %>% summarise(mean.rK = mean(r_K), n.rK = n(), sd.rK = sd(r_K), SE = sd.rK/sqrt(n.rK), 
#                                                                                                            lci.rK = mean.rK - qt(1 - ((1 - 0.95) / 2), n.rK - 1) * SE,
#                                                                                                            uci.rK = mean.rK + qt(1 - ((1 - 0.95) / 2), n.rK - 1) * SE) %>%
#   select(SITE, mean.rK, n.rK, sd.rK, lci.rK, uci.rK) %>% inner_join(df, by = "SITE")
# df <- MSBio_micfg %>% filter(material == "soil") %>% mutate(SITE = site) %>% group_by(SITE)  %>% summarise(mean.CO = mean(C_O), n.CO = n(),sd.CO = sd(C_O), SE = sd.CO/sqrt(n.CO), 
#                                                                                                            lci.CO = mean.CO - qt(1 - ((1 - 0.95) / 2), n.CO - 1) * SE,
#                                                                                                            uci.CO = mean.CO + qt(1 - ((1 - 0.95) / 2), n.CO - 1) * SE) %>%
#   select(SITE, mean.CO, n.CO, sd.CO, lci.CO, uci.CO) %>% inner_join(df, by = "SITE")

#estimating cost
#Derek's 2022 paper uses RMSE which only accounts for differences and not errors
#consider using maximum likelihood estimation instead? See Richardson and Hollinger, 2005 & Keenan et al., 2011
#initial cost functions
# df$cost.rK <- abs(df$MIM_CO - df$mean.rK) 
# df$cost.CO <- abs(df$MIM_CO - df$mean.CO)
# boxplot(df$cost.rK)
# #alternate cost functions
# #checking normality and variance assumptions of obs
# hist(MSBio_micfg$r_K) #very skewed, regardless if you cutoff outliers
# hist(MSBio_micfg$C_O) #very skewed, regardless if you cutoff outliers
# hist(MSBio_micfg$perc_decomp_T1) #pretty normal
# hist(MSBio_micfg$perc_decomp_T2) #very normal
#skew means these data don't fit the classic MLE assumptions of normality
#trying Keenan et al 2012 cost function with error
# df$cost.rK.2 <- abs(((df$mean.rK- df$MIM_CO)/ df$sd.rK)/df$n.rK) #a few runs have huge costs
# df$cost.CO.2 <- abs(((df$mean.CO- df$MIM_CO)/ df$sd.CO)/df$n.CO) #a few runs have huge costs
#with litter mass loss also accounted for - sum of above and LML version divided by 2(for number of data streams used, also Keenan et al., 2011)
df_LML$cost.ML <-abs(((df_LML$mean.ML- df_LML$LIT_PerLoss)/ df_LML$sd.ML)/df_LML$n.ML)
boxplot(df_LML$cost.ML) #really good!
df$SITE.rn <- paste(df$SITE, df$run_num, sep = ".")
df <- df_LML %>% mutate(SITE.rn = paste(SITE, run_num, sep = ".")) %>% select(SITE.rn, cost.ML) %>% right_join(df, by="SITE.rn") %>%
  
# df$cost.rK.3 <- (abs(((df$mean.rK- df$MIM_CO)/ df$sd.rK)/df$n.rK) + df$cost.ML)/2
# df$cost.CO.3 <- (abs(((df$mean.CO- df$MIM_CO)/ df$sd.CO)/df$n.CO) + df$cost.ML)/2
# boxplot(df$cost.rK.3)
# boxplot(df$cost.CO.3)
  
saveRDS(df_LML, paste0("MC_output/MSBio_df.LML_", format(Sys.time(), "%Y%m%d_%H%M%S_"),  ".rds"))

#####
#cost as matching relative weight analysis or effect size
#####

#filter for reasonable data
df_analysis <- df_LML %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>%
  filter(time.point==2) %>% select(SITE, LIT_PerLoss,
                                   MICrK, MIC, SOC, run_num, Tau_x, CUE_x, vMOD_x,
                                   kMOD_x, fM_x, TSOI, LIG_N, W_SCALAR)
#logical checks
df_check <- df_analysis %>% filter(MICrK > 0.01) %>%
  filter(MICrK < 100) %>%
  filter(MIC/SOC > 0.0001) %>%
  filter(MIC/SOC < 0.40) 

#OPTION 1: cost as RWA - note when observational values are reasonably close, it may pick the same set of lowest parameter sets so no difference will be seen
#relative weights for each run
MIM_rwa <- data.frame()
for (i in 1:100) {
  df_rwa <- filter(df_check, run_num==i)
  rwa_mod <- rwa(df_rwa, "LIT_PerLoss", c("TSOI", "W_SCALAR", "LIG_N", "MICrK"), applysigns = TRUE, plot = FALSE)
  rwa <- as.data.frame(rwa_mod$result)
  rwa$run_num <- i
  MIM_rwa <- rbind(MIM_rwa,rwa)
}
saveRDS(MIM_rwa, paste0("MC_output/MSBio_RWA_", format(Sys.time(), "%Y%m%d_%H%M%S_"),  ".rds"))
# #select weights that best match observation weights- create cost function for weights
# #create dataframe of obs RWA - here using RWAs generated with C:O from Averill et al
# Variables = c("TSOI", "W_SCALAR", "LIG_N", "MICrK") #matching mdoel names
# #obs_rw = c(11.2, 11.9, -63.2, -13.6) #with W_SCALAR
# #obs_rw = c(19.1, 12.1, -55.8, -13) #with MAP
# obs_rw = c(18.5, 24.7, -50.5, -6.3) #with VWC
# obs_rwa <- data.frame(Variables, obs_rw)
# #note that for next line df_LML$cost.ML comes from higher up code
# MIM_rwa2 <- df_LML %>% filter(time.point == 2) %>% group_by(run_num) %>% summarise(mean_cost_ML = mean(cost.ML))%>% right_join(MIM_rwa, by="run_num")
# RWA_cost <- MIM_rwa2 %>% inner_join(obs_rwa, by="Variables") %>% mutate(rw.dif = Sign.Rescaled.RelWeight - obs_rw) %>% 
#   mutate(cost2 = (abs(rw.dif)+mean_cost_ML)/2) %>% mutate(rel.cost = (cost2/abs(obs_rw))*100)
# #low cost based on top models
# RWA_low.cost <- RWA_cost %>% group_by(run_num) %>% summarise(mean_dif = mean(cost2)) %>% slice_min(mean_dif, n=50)
# #low cost based on relative cost value
# #RWA_low.cost <- RWA_cost %>% group_by(run_num) %>% summarise(mean_dif = mean(rel.cost)) #%>% filter(mean_dif<75)
# rand_params <- df_LML %>% select(run_num, Tau_x, CUE_x, vMOD_x, kMOD_x)
# test <- RWA_low.cost %>% left_join(RWA_cost, by="run_num") %>% left_join(rand_params, by="run_num")
# test2 <- test %>% select(run_num, Tau_x, CUE_x, vMOD_x, kMOD_x) %>% pivot_longer(2:5, names_to = "Multiplier", values_to = "value")
# #ridge plots for parameter mulitpliers
# #ggplot(test, aes(x=fM_x)) + geom_density()
# Params_rwa <- ggplot(test2, aes(x = value, y=Multiplier, group=Multiplier, fill=Multiplier))+
#   geom_density_ridges(scale = 2) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   scale_fill_brewer(palette = "Oranges") +
#   labs(title="Parameter multipliers for 50 lowest cost runs out of 8000") +#,
#   #subtitle="n=200 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# Params_rwa
# #ridge plots for cost
# Cost_rwa <- ggplot(test, aes(x = rel.cost, y=Variables, group=Variables, fill=Variables))+
#   geom_density_ridges(scale = 1, rel_min_height=0.01) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   scale_fill_brewer(palette = "Oranges") +
#   labs(title="Percent relative cost by variable for lowest mean cost by run number") +#,
#   #subtitle="n=200 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# Cost_rwa
# test %>% group_by(Variables) %>% summarise(RW.mean = mean(Sign.Rescaled.RelWeight)) %>% ggplot(aes(x=Variables, y=RW.mean)) + 
#   geom_bar(stat="identity", fill="blue") + coord_flip() + geom_text(aes(label=round(RW.mean, digits=1)), color="red", size=7) +theme_bw(base_size = 16)

#OPTION 2: cost as difference in effect size
#effect sizes for each run
MIM_ES <- data.frame() #varaibles correlated? won't let me run the model, seemingly because singularity causes NAs? There are no NAs in the data and that also provides the error I found
for (i in 1:100) {     #W_SCALAR is strongly correlated (0.9) with LIT_PerLoss but removing the WS does not allow the data to run... not sure what's going on!
  df_ES <- filter(df_check, run_num==i)
  #Obs_ES_mod <- lmer(LIT_PerLoss ~ TSOI+W_SCALAR+LIG_N+MICrK+(1|SITE), data = df_ES)
  #Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) #fixed effects coefficients as effect size
  Obs_ES_mod <- lm(LIT_PerLoss ~ TSOI+LIG_N+MICrK+W_SCALAR, data = df_ES) #+W_SCALAR
  Obs_ES <- as.data.frame(Obs_ES_mod$coefficients) #fixed effects coefficients as effect size
  Obs_ES$Vars <- rownames(Obs_ES)
  colnames(Obs_ES)[1] <- "value"
  Obs_ES <- Obs_ES[-1, ]
  Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
  Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
  Obs_ES$run_num <- i
  MIM_ES <- rbind(MIM_ES, Obs_ES)
}
saveRDS(MIM_ES, paste0("MC_output/MSBio_ES_", format(Sys.time(), "%Y%m%d_%H%M%S_"),  ".rds"))

#select effect sizes that best match observation weights- create cost function for weights
#create dataframe of obs RWA - here using RWAs generated with C:O from Averill et al
# Vars = c("TSOI", "W_SCALAR", "LIG_N", "MICrK") #matching mdoel names
# #obs_ES = c(1.7, 88, -4.4, -5.9) #for observational data with W_SCALAR
# #obs_ES = c(20.9, 0.1, -35.4, -43.6) #for observational data with MAP
# obs_ES = c(28.1, 1.6, -41.5, -28.8) #with VWC
# obs_ES_df <- data.frame(Vars, obs_ES)
# #note that for next line df_LML$cost.ML comes from higher up code
# MIM_ES2 <- df_LML %>% filter(time.point == 2) %>% group_by(run_num) %>% summarise(mean_cost_ML = mean(cost.ML))%>% right_join(MIM_ES, by="run_num")
# ES_cost <- MIM_ES2 %>% inner_join(obs_ES_df, by="Vars") %>% mutate(ES.dif = rel_ES - obs_ES) %>% 
#   mutate(cost2 = (abs(ES.dif)+mean_cost_ML)/2) %>% mutate(rel.cost = (cost2/abs(obs_ES))*100)
# #low cost based on defined number of lowest cost parameter sets 
# ES_low.cost <- ES_cost %>% group_by(run_num) %>% summarise(mean_dif = mean(cost2)) %>% slice_min(mean_dif, n=50)
# #low cost based on relative cost value - doesn't work well because W_SCALAR is so low
# #ES_low.cost <- ES_cost %>% group_by(run_num) %>% summarise(mean_dif = mean(rel.cost)) %>% filter(mean_dif<200)
# test <- ES_low.cost %>% left_join(ES_cost, by="run_num") %>% left_join(rand_params, by="run_num")
# test2 <- test %>% select(run_num, Tau_x, CUE_x, vMOD_x, kMOD_x) %>% pivot_longer(2:5, names_to = "Multiplier", values_to = "value")
# Params_ES <- ggplot(test2, aes(x = value, y=Multiplier, group=Multiplier, fill=Multiplier))+
#   geom_density_ridges(scale = 2) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   scale_fill_brewer(palette = "Oranges") +
#   labs(title="Parameter multipliers for 30 lowest cost runs out of 100") +#,
#   #subtitle="n=200 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# Params_ES
# #ridge plots for cost
# Cost_rwa <- ggplot(test, aes(x = cost2, y=Vars, group=Vars, fill=Vars))+
#   geom_density_ridges(scale = 1, rel_min_height=0.01) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   scale_fill_brewer(palette = "Oranges") +
#   labs(title="Cost of lowest 50 cost p-sets by variable") +#,
#   #subtitle="n=200 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# Cost_rwa
# test %>% group_by(Vars) %>% summarise(ES.mean = mean(rel_ES)) %>% ggplot(aes(x=Vars, y=ES.mean)) + 
#   geom_bar(stat="identity", fill="blue") + coord_flip() + geom_text(aes(label=round(ES.mean, digits=1)), color="red", size=7) +theme_bw(base_size = 16)
