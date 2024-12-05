######################################
#code for analyzing MC runs
######################################

#load libraries
library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(gridExtra)
library(car)
library(lmerTest)
library(rwa)
library(ranger)
library(caTools)
library(vip)
library(scatterplot3d)
library(plotly)
library(lattice)
library(rgl)

####
#if running "MIMICS_LitBag_repeat" directly in R and not on a supercomputer
#if using super computer output (as done in manuscript), skip down to line 98
####

#bringing monte carlo (MC) output together with input data and random parameters
df <- left_join(MC_MIMICS, data, by = "SITE") 
df <- left_join(df,rand_params, by="run_num")
#calculating some common metrics
df$MIM_CO <- as.numeric(df$MICr)/as.numeric(df$MICk)
df$MIC_SOC <- (df$MICr+df$MICk)/(df$SOMc+df$SOMa+df$MICr+df$MICk)
df$LITBAG_tot <- df$LITBAGm + df$LITBAGs


#observational data for comparing to
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


#####
#cost as matching effect size
#####

#prepping data for effect size analysis
LITi = 0.1
df_LML <- df %>% mutate(SITE.rn = paste(SITE, run_num, sep = ""))  %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100) %>% 
  mutate(field.day=DAY-314) %>% mutate(SITE.DAY=paste(SITE, field.day, sep=".")) %>% right_join(FieldData, by="SITE.DAY")
DI_means <- DailyInput_SM %>% mutate(SITE.SM = paste(SM_type, SITE,sep = ".")) %>% group_by(SITE.SM) %>% 
  summarise(W_SCALAR_mean=mean(W_SCALAR)) %>% select(SITE.SM, W_SCALAR_mean) # SM_type, 
#initial MICrK
MIC_init <- df %>% filter(DAY == 315) %>% mutate(MICrK.i =  MICr/MICk)%>%
  mutate(SITE.SM.LQ.rn = paste(SITE, SM_Type, Litter_Type, run_num, sep = ".")) %>% select(SITE.SM.LQ.rn, MICrK.i) #SM_Type, 
#bag means
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
df_analysis <- df_LML %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ.rn = paste(SITE, SM_Type, Litter_Type, run_num, sep = ".")) %>% mutate(SITE.SM = paste(SM_Type, SITE, sep = ".")) %>% # SM_Type,  SM_Type, 
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% 
  inner_join(MIC_init, by="SITE.SM.LQ.rn") %>% inner_join(BAGS_LIGN, by="SITE.LQ")
#logical checks
 df_check <- df_analysis %>% filter(MICrK > 0.01) %>%
  filter(MICrK < 100) %>%
  filter(MICrK.i > 0.01) %>%
  filter(MICrK.i < 100) %>%
  filter(MIC/SOC > 0.0001) %>%
  filter(MIC/SOC < 0.40)
#filter for paramter sets that fit within minimum and maximum litter mass loss in observations 
df_check$SITE <- as.factor(df_check$SITE)
df_rn1 <- df_analysis %>% filter(time.point == 1) %>% filter(LIT_PerLoss > min.ML & LIT_PerLoss < max.ML) #highest and lowest values of confidence interval for time point 1
LR.rn1 <- df_rn1 %>% group_by(run_num) %>% mutate(S.LQ.SM = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% 
  summarize(uniq.site = length(unique(S.LQ.SM))) %>% filter(uniq.site>6) #%>% filter(uniq.site>62) #%>% filter(site.levels>4) #as.data.frame(table(df_rn1$run_num)) %>% filter(Freq > 4)
df_rn2 <- df_check %>% filter(time.point == 2) %>% filter(LIT_PerLoss > min.ML & LIT_PerLoss < max.ML) #highest and lowest values of confidence interval for time point 2
LR.rn2 <- df_rn2 %>% group_by(run_num) %>% mutate(S.LQ.SM = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% 
  summarize(uniq.site = length(unique(S.LQ.SM))) %>% filter(uniq.site>6) #%>% filter(uniq.site>62)  #as.data.frame(table(df_rn2$run_num)) %>% filter(Freq > 4)
#keep only run numbers that are in both dfs
df_rn1.v <- as.vector(LR.rn1$run_num) #as.numeric(LR.rn1$Var1)
df_rn2.v <- as.vector(LR.rn2$run_num) #as.numeric(LR.rn2$Var1)
rn_final <- intersect(df_rn1.v,df_rn2.v)
#filter check df to only have run numbers that fit litter mass loss
df_check2 <- df_check %>% filter(run_num %in% rn_final)
#transforming and scaling data for analysis
df_TS <- df_check2 %>% mutate(log_WS = log(W_SCALAR_mean)) %>%select(run_num, SITE, LIT_PerLoss, log_WS, BAG_LIG_N, MICrK.i) %>% 
  mutate_at(vars(c("log_WS", "BAG_LIG_N", "MICrK.i")), ~(scale(.) %>% as.vector))
df_TS$UniqueID <- paste(df_TS$run_num, df_TS$log_WS,df_TS$BAG_LIG_N, sep = ".")
df_RP <- df_check2 %>% select(run_num, beta_x, Tau_r, vMOD_m, vMOD_s) %>% unique()



#if creating df_TS on a super computer:
df_TS.1 <- readRDS('Analysis/MC_output/MSBio_df_TS_2500.2_1.rds')
df_TS.2 <- readRDS('Analysis/MC_output/MSBio_df_TS_2500.2_2.rds')
df_TS <- rbind(df_TS.1, df_TS.2)
df_TS$UniqueID <- paste(df_TS$run_num, df_TS$log_WS,df_TS$BAG_LIG_N, sep = ".")
df_RP.1 <- readRDS('Analysis/MC_output/MSBio_df_RP_2500.2_1.rds')
df_RP.2 <- readRDS('Analysis/MC_output/MSBio_df_RP_2500.2_2.rds')
df_RP <- rbind(df_RP.1, df_RP.2)


#cost as difference in effect size
rn_final <- unique(df_TS$run_num) 
MIM_ES <- data.frame()
for (i in seq_along(rn_final)) {
  df_ES <- filter(df_TS, run_num==rn_final[i])
  df_ES2 <- df_ES %>% mutate_at(vars(c("log_WS", "BAG_LIG_N", "MICrK.i")), ~(scale(.) %>% as.vector))
  Obs_ES_mod <- lmer(LIT_PerLoss ~ log_WS+BAG_LIG_N+MICrK.i+(1|SITE/UniqueID), data = df_ES2)
  Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) #fixed effects coefficients as effect size
  Obs_ES$Vars <- rownames(Obs_ES)
  colnames(Obs_ES)[1] <- "value"
  Obs_ES <- Obs_ES[-1, ]
  Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
  Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
  Obs_ES$run_num <- rn_final[i]
  MIM_ES <- rbind(MIM_ES, Obs_ES)
}
#create dataframe of observed effect sizes- here using effect sizes generated with C:O from Averill et al
Vars = c("log_WS", "BAG_LIG_N", "MICrK.i") 
obs_ES = c(42.8, -34.2, -22.9) 
obs_ES_df <- data.frame(Vars, obs_ES)
#create filtering metrics
#+/- 10 5000 runs
var_upper = c(52.8, -24.2, -12.9) #make uppers 100 if don't want upper bound
var_lower = c(32.8, -44.2, -32.9) #make lowers -100 if don't want lower bound
tol_table <- data.frame(Vars, var_lower, var_upper)
ES_rn <- MIM_ES %>% inner_join(obs_ES_df, by="Vars") %>% inner_join(tol_table, by="Vars") %>% 
  filter(rel_ES < var_upper & rel_ES > var_lower)
ES_rn2 <- as.data.frame(table(ES_rn$run_num)) %>% filter(Freq > 2) #ensuring all variables are present in "ideal" parameter set
ES_rn.v <- ES_rn2$Var1
#next line will only return viable paramter sets
ES_cost <- MIM_ES %>% inner_join(obs_ES_df, by="Vars") %>% mutate(es.dif = rel_ES - obs_ES) %>% filter(run_num %in% ES_rn.v)
#sensitivity to number of random parameter sets
#rn_filter = c(1:4500)
#ES_cost2 <- ES_cost %>% filter(run_num %in% rn_filter)
#combining effect sizes and parameters
test <- ES_cost %>% left_join(df_RP, by="run_num")
#format for plotting
test2 <- test %>% select(run_num, Tau_r, vMOD_m, vMOD_s, beta_x) %>% pivot_longer(2:5, names_to = "Multiplier", values_to = "value")
#plot paramter distributions
Params_ES <- ggplot(test2, aes(x = value, y=Multiplier, group=Multiplier, fill=Multiplier))+
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Oranges") +
  labs(title="Parameter multipliers") + # for 227 best runs out of 5000 using effect size") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
Params_ES
