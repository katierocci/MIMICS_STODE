library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(gridExtra)
library(car)

##need to figure out how to 
#df <- readRDS("C:/github/MIMICS_MSBio/Cheyenne_HPC/HPC_output/MSBio_MIM_MC_runs-1e+05_20221028_112647_.rds")
#MC_MIMICS <- readRDS('Derecho_HPC/Jan 2024/MC_output/MSBio_MC_10_20240117_104907_.rds')
df <- MC_MIMICS
#rand_params <- readRDS('Cheyenne_HPC/July 2023/MIMICS_MSBio_KR/MC_output/MSBio_RP_1e+05_20230728_132810_Tc.rds')
df <- left_join(MC_MIMICS, data, by = "SITE") #MC output has to be first for CO2 rows below to be right
df <- left_join(df, rand_params, by = "run_num")
#df$CO2_of_tot <- rowSums(df[,11:12])/rowSums(df[,4:12])
df$MIM_CO <- as.numeric(df$MICr)/as.numeric(df$MICk)
df$MIC_SOC <- (df$MICr+df$MICk)/(df$SOMc+df$SOMa+df$MICr+df$MICk)
df$LITBAG_tot <- df$LITBAGm + df$LITBAGs

#looking at parameter effects
df_test <- df %>% filter(MIM_CO < 100)
plot(df_test$beta_r, df_test$MIM_CO) #crazy vales near 1.2
plot(df$beta_k, df$MIM_CO)
plot(df_730$beta_k, df_730$LIT_PerLoss) #general increase with increased beta
plot(df$beta_r, df$beta_k)

###
#cost functions
###
#for each run number calculate initial LIT_tot, make that a column, filter to the last day (730) and then subtract that value from that column for percent loss
df$SITE.rn = paste(df$SITE, df$run_num, sep = "")
LIT_init <- df %>% filter(DAY == 10) %>% mutate(LITi = LITBAG_tot) %>% select(SITE.rn, LITi) 
boxplot(LIT_init$LITi)
df <- df %>% left_join(LIT_init, by = "SITE.rn") 
df_730 <- df %>% filter(DAY==730) %>% mutate(LIT_PerLoss = ((LITi - LITBAG_tot)/LITi)*100) #(sample - recovered)/sample *100
boxplot(df_730$LIT_PerLoss) #looks reasonable
#more accurate cost function! Test this tomorrow!
Field_LML <- read.csv("Example_simulations/Data/Litter_decomp_all.csv")
LML_sum2 <- Field_LML  %>% group_by(site, time.point) %>% drop_na(percent.loss.litter) %>% mutate(percent.loss.C = (percent.loss.litter/2)*100) %>% 
  summarize(mean.ML = mean(percent.loss.C),
            n = n(),
            sd = sd(percent.loss.C/2),
            SE = sd/sqrt(n),
            lci.ML = mean.ML - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
            uci.ML = mean.ML + qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
            doy = mean(days_elapsed)) %>% mutate(doy=round(doy, digits=0))
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML)
LIT_init <- df %>% mutate(SITE.rn = paste(df$SITE, df$run_num, sep = "")) %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% select(SITE.rn, LITi)
df_LML <- df %>% mutate(SITE.rn = paste(df$SITE, df$run_num, sep = "")) %>% left_join(LIT_init, by = "SITE.rn")
df_LML <- df_LML %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
#bringing in empirical comparison data
#r:K/C:O
MSBio_micfg <- read.csv("Example_simulations/Data/MSBio_FuncGroups_AllSites_prelim.csv")
df <- MSBio_micfg %>% filter(material == "soil") %>% mutate(SITE = site) %>% group_by(SITE)  %>% summarise(mean.rK = mean(r_K), n.rK = n(), sd.rK = sd(r_K), SE = sd.rK/sqrt(n.rK), 
                                                                                                           lci.rK = mean.rK - qt(1 - ((1 - 0.95) / 2), n.rK - 1) * SE,
                                                                                                           uci.rK = mean.rK + qt(1 - ((1 - 0.95) / 2), n.rK - 1) * SE) %>%
  select(SITE, mean.rK, n.rK, sd.rK, lci.rK, uci.rK) %>% inner_join(df, by = "SITE")
df <- MSBio_micfg %>% filter(material == "soil") %>% mutate(SITE = site) %>% group_by(SITE)  %>% summarise(mean.CO = mean(C_O), n.CO = n(),sd.CO = sd(C_O), SE = sd.CO/sqrt(n.CO), 
                                                                                                           lci.CO = mean.CO - qt(1 - ((1 - 0.95) / 2), n.CO - 1) * SE,
                                                                                                           uci.CO = mean.CO + qt(1 - ((1 - 0.95) / 2), n.CO - 1) * SE) %>%
  select(SITE, mean.CO, n.CO, sd.CO, lci.CO, uci.CO) %>% inner_join(df, by = "SITE")
#litter mass loss
#compare to field litter mass CARBON (divide by 2) loss - currently imperfect since running model for two years but time point 2 is at various years - may need to customize model output for each site!
Field_LML <- read.csv("Example_simulations/Data/Litter_decomp_all.csv")
df_730 <- Field_LML %>% filter(time.point==2) %>% drop_na(percent.loss.litter) %>% mutate(SITE=site) %>%group_by(SITE)%>% summarize(mean.ML = mean(percent.loss.litter/2),
                                                                                                                         n.ML = n(),
                                                                                                                         sd.ML = sd(percent.loss.litter/2),
                                                                                                                         SE.ML = sd.ML/sqrt(n.ML),
                                                                                                                         lci.ML = mean.ML - qt(1 - ((1 - 0.95) / 2), n.ML - 1) * SE.ML,
                                                                                                                         uci.ML = mean.ML + qt(1 - ((1 - 0.95) / 2), n.ML - 1) * SE.ML,
                                                                                                                         doy = mean(days_elapsed)) %>%
  select(SITE, mean.ML, n.ML, sd.ML, lci.ML, uci.ML) %>% inner_join(df_730, by = "SITE")
#estimating cost
#Derek's 2022 paper uses RMSE which only accounts for differences and not errors
#consider using maximum likelihood estimation instead? See Richardson and Hollinger, 2005 & Keenan et al., 2011
#initial cost functions
df$cost.rK <- abs(df$MIM_CO - df$mean.rK) 
df$cost.CO <- abs(df$MIM_CO - df$mean.CO)
boxplot(df$cost.rK)
#alternate cost functions
#checking normality and variance assumptions of obs
hist(MSBio_micfg$r_K) #very skewed, regardless if you cutoff outliers
hist(MSBio_micfg$C_O) #very skewed, regardless if you cutoff outliers
hist(MSBio_micfg$perc_decomp_T1) #pretty normal
hist(MSBio_micfg$perc_decomp_T2) #very normal
#skew means these data don't fit the classic MLE assumptions of normality
#trying Keenan et al 2012 cost function with error
df$cost.rK.2 <- abs(((df$mean.rK- df$MIM_CO)/ df$sd.rK)/df$n.rK) #a few runs have huge costs
df$cost.CO.2 <- abs(((df$mean.CO- df$MIM_CO)/ df$sd.CO)/df$n.CO) #a few runs have huge costs
#with litter mass loss also accounted for - sum of above and LML version divided by 2(for number of data streams used, also Keenan et al., 2011)
df_730$cost.ML <-abs(((df_730$mean.ML- df_730$LIT_PerLoss)/ df_730$sd.ML)/df_730$n.ML)
boxplot(df_730$cost.ML)
df <- df_730 %>% select(SITE.rn, cost.ML) %>% right_join(df, by="SITE.rn")
df$cost.rK.3 <- (abs(((df$mean.rK- df$MIM_CO)/ df$sd.rK)/df$n.rK) + df$cost.ML)/2
df$cost.CO.3 <- (abs(((df$mean.CO- df$MIM_CO)/ df$sd.CO)/df$n.CO) + df$cost.ML)/2
boxplot(df$cost.rK.3)
boxplot(df$cost.CO.3)

#looking at cost
hist(df$cost.rK)
hist(df$cost.CO)
df %>% filter(cost.rK < 100) %>% ggplot(aes(y=cost.rK)) + geom_boxplot() +ggtitle("Absolute value of the difference\nin predicted and obs") +theme_bw(base_size = 20)
df %>% filter(cost.CO < 100) %>% ggplot(aes(y=cost.CO)) + geom_boxplot() +ggtitle("Absolute value of the difference\nin predicted and obs") +theme_bw(base_size = 20)
df %>% filter(cost.rK.2 < 100) %>% ggplot(aes(y=cost.rK.2)) + geom_boxplot() +ggtitle("Error and sample size weighted\ndifference between predicted & obs") +theme_bw(base_size = 20)
df %>% filter(cost.CO.2 < 100) %>% ggplot(aes(y=cost.CO.2)) + geom_boxplot() +ggtitle("Error and sample size weighted\ndifference between predicted & obs") +theme_bw(base_size = 20)
#cost vs multiplier
df %>% filter(cost.rK < 10) %>% ggplot(aes(x=CUE_x, y=cost.rK, group = CUE_x)) + geom_boxplot() + ggtitle("All sites cost (<10) with obs r:K vs CUE multiplier") #higher CUE generally lower cost #filter(run_num != 5) %>%
df %>%  filter(cost.rK < 10) %>% ggplot(aes(x=Tau_x, y=cost.rK, group = Tau_x)) + geom_boxplot() + ggtitle("All sites cost (<10) with obs r:K vs Tau multiplier") #lowest tau lowest cost
df %>%  filter(cost.CO < 100) %>% ggplot(aes(x=CUE_x, y=cost.CO, group = CUE_x)) + geom_boxplot() #higher CUE generally lower cost
df %>% filter(cost.CO < 100) %>% ggplot(aes(x=Tau_x, y=cost.CO, group = Tau_x)) + geom_boxplot() #lowest tau lowest cost
df %>% filter(run_num == 10) %>% ggplot(aes(x=DAY, y=cost.rK)) + geom_point() #cost doesn't vary over time
#cool seems feasible but might need to do some sort of cleaning to get rid of unrealistic values
df %>% filter(MIM_CO < 100) %>% ggplot(aes(y=MIM_CO, x=CUE_x)) +geom_point() +ggtitle("CUE_x multiplier vs r:K for CUE_x-only MC of 100 runs") #r:K increases as CUE_x increases
df %>% filter(MIM_CO < 100) %>% ggplot(aes(y=MIM_CO, x=CUE_k)) +geom_point() +ggtitle("CUE_k multiplier vs r:K for CUE_r+CUE_k MC of 100 runs")
#library(plotly)
df_lowrK <- df %>% filter(MIM_CO < 100)
plot_ly(df_lowrK, x = ~CUE_r, y = ~MIM_CO, z = ~CUE_k)


### Ftn to find best psets (this is an objective function?)
best_psets <- function(df, best_n) {
  cost_df <- df %>% group_by(run_num) %>% summarize(SITE = unique(SITE), n = n(), cost = sum(cost))
  return(cost_df %>% slice_min(cost, n=best_n))
}

best_data <- function(df, best_n) {
  cost_df <- df %>% group_by(run_num) %>% summarize(SITE = unique(SITE),
                                                    n = n(),
                                                    cost = sum(cost),
                                                    min_MIM_CO = min(MIM_CO),
                                                    max_MIM_CO = max(MIM_CO),
                                                    min_MIC_SOC = min(MIC_SOC),
                                                    max_MIC_SOC = max(MIC_SOC))
  low_cost <- cost_df %>%
    filter(min_MIM_CO > 0.01) %>%
    filter(max_MIM_CO < 100) %>%
    filter(min_MIC_SOC > 0.0001) %>%
    filter(max_MIC_SOC < 0.40) %>%
    slice_min(cost, n=best_n)
  
  out <- df %>% filter(run_num %in% low_cost$run_num)
  return(out)
}

### Collect all the best psets by site
sites <- unique(df$SITE)
df$cost <- df$cost.rK.3
top_df <- NULL
for(i in 1:length(sites)) {
  df1 <- best_data(df = df %>% filter(SITE == sites[i]), best_n = 30)
  top_df <- rbind(top_df, df1)
}

#plot distirbution of lowest cost parameters for each site
#putting in order of MAT
MAT_order <- c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO")
top_df$SITE <- factor(top_df$SITE, levels=MAT_order, labels=MAT_order)

plot(top_df$beta_r, top_df$beta_k)

pCOST <- ggplot(top_df, aes(x = cost, y = SITE, fill=SITE, group=SITE)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Oranges") +
  labs(title="Cost at each site for 30 lowest cost runs out of 100") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pCOST #cost still quite high for 3 out of 4 sites

pbeta.k <- ggplot(top_df, aes(x = beta_r, y = SITE, fill=SITE, group=SITE)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Oranges") +
  labs(title="beta Multiplier for 30 lowest cost runs out of 100") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pbeta.k

pCUE <- ggplot(top_df, aes(x = CUE_x, y = SITE, fill=SITE, group=SITE)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Oranges") +
  labs(title="CUE Multiplier for 30 lowest cost runs out of 100") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pCUE

pTAU <- ggplot(top_df, aes(x = Tau_x, y = SITE, fill=SITE, group=SITE)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Oranges") +
  labs(title="Tau Multiplier for 30 lowest cost runs out of 100") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pTAU

#comparing multipliers with environmental data
#linear models
#top_df_envi <- top_df %>% inner_join(data, by= "SITE")
CUE_mod <- lm(CUE_x ~ TSOI + grav_moisture + LIG_N, data = top_df)
plot(CUE_mod) #seems fine
Anova(CUE_mod, type=3) #all sig
ggplot(top_df, aes(x=grav_moisture, y=CUE_x, fill=SITE)) + geom_boxplot() #no clear relationships with any envi variable
Tau_mod <- lm(Tau_x ~ TSOI + GRAV_MOISTURE + LIG_N, data = top_df_envi)
plot(Tau_mod) #seems fine
Anova(Tau_mod, type=3) #all sig
ggplot(top_df_envi, aes(x=GRAV_MOISTURE, y=Tau_x, fill=SITE)) + geom_boxplot() #no clear relationships with any envi variable
beta_mod <- lm(beta_x ~ TSOI + grav_moisture + LIG_N, data = top_df)
plot(beta_mod) #seems fine
Anova(beta_mod, type=3) #all sig
ggplot(top_df, aes(x=LIG_N, y=beta_x, fill=SITE)) + geom_boxplot() +theme_bw() #beta has a U-shaped relationship with soil moisture, persists when varying CUE as well








#keeping all this for now but may want to get rid of
#alt cost function - trying to match max CO2 prop across WHCs
###Will suggests the alt cost function!
Max_CO2 <- data %>% group_by(SITE) %>% summarise(max=max(CO2C_prop)) #provides max incubation CO2
#below is checking empirical relatiosnhip between micorsite moisture and max CO2, currently for Bartlett
#Max_CO2$GWC <- c(0.623, 0.674, 0.658, 0.47, 0.621, 0.491, 0.544, 0.466, 0.591, 0.656, 0.51, 0.376)
#ggplot(Max_CO2, aes(x=GWC, y=max)) +geom_point() + geom_smooth() #vaguely positive
###cost function
df_cost <- df %>% filter(DAY == 105) %>% group_by(SITE.x, run_num) %>% summarise(max_model= max(CO2_of_tot))
#BELOW NEEDS TO BE MODIFIED DEPENDING ON SITE
#max incubation co2 by site, repeated 1000 time for 1000 MC runs
#BART, GRSM, TREE
df_cost$max_incub <- c(rep(Max_CO2$max[1], 1000), rep(Max_CO2$max[2], 1000), rep(Max_CO2$max[3], 1000), rep(Max_CO2$max[4], 1000),
                       rep(Max_CO2$max[5], 1000),rep(Max_CO2$max[6], 1000), rep(Max_CO2$max[7], 1000), rep(Max_CO2$max[8], 1000),
                       rep(Max_CO2$max[9], 1000), rep(Max_CO2$max[10], 1000),rep(Max_CO2$max[11], 1000), rep(Max_CO2$max[12], 1000))
# #HARV
# df_cost$max_incub <- c(rep(Max_CO2$max[1], 1000), rep(Max_CO2$max[2], 1000), rep(Max_CO2$max[3], 1000), rep(Max_CO2$max[4], 1000), rep(Max_CO2$max[5], 1000),
#                        rep(Max_CO2$max[6], 1000), rep(Max_CO2$max[7], 1000), rep(Max_CO2$max[8], 1000), rep(Max_CO2$max[9], 1000), rep(Max_CO2$max[10], 1000),
#                        rep(Max_CO2$max[11], 1000), rep(Max_CO2$max[12], 1000), rep(Max_CO2$max[13], 1000), rep(Max_CO2$max[14], 1000), rep(Max_CO2$max[15], 1000),
#                        rep(Max_CO2$max[16], 1000), rep(Max_CO2$max[17], 1000), rep(Max_CO2$max[18], 1000), rep(Max_CO2$max[19], 1000), rep(Max_CO2$max[20], 1000),
#                        rep(Max_CO2$max[21], 1000), rep(Max_CO2$max[22], 1000), rep(Max_CO2$max[23], 1000), rep(Max_CO2$max[24], 1000), rep(Max_CO2$max[25], 1000),
#                        rep(Max_CO2$max[26], 1000), rep(Max_CO2$max[27], 1000))
df_cost <-left_join(df_cost, rand_params, by = "run_num")
df_MIM <- df %>% filter(DAY == 105) %>% group_by(SITE.x, run_num) %>% summarise(MIM_CO= mean(MIM_CO), MIC_SOC= mean(MIC_SOC))
df_cost$Site_RN <- paste(df_cost$SITE.x, df_cost$run_num, sep = "_")
df_MIM$Site_RN <- paste(df_MIM$SITE.x, df_MIM$run_num, sep = "_")
df_cost <-left_join(df_cost, df_MIM, by = "Site_RN")
df_cost$cost <- abs(df_cost$max_model - df_cost$max_incub)
df_cost$cost2 <- df_cost$max_model - df_cost$max_incub #always negative for CUE and tau
df_cost <- df_cost %>% drop_na(cost) # still NAs
max(df_cost$cost2)
mean(df_cost$cost2)
#debugging
ggplot(df_cost, aes(x=SITE.x.x, y=cost, group=SITE.x.x)) + geom_boxplot()
#cost is different between sites but distributed the same
ggplot(df_cost, aes(x=CUE_x, y=cost, group=SITE.x.x, color=SITE.x.x)) + geom_point()
#CUE so tightly tied to CO2 that the lowest cost is always the same parameter
#tying vMOD below to see if also strongly constrained like CUE
ggplot(df_cost, aes(x=vMOD_x, y=max_model, group=SITE.x.x, color=SITE.x.x)) + geom_point()
#tau also looks the same...
ggplot(df_cost, aes(x=Tau_x, y=cost, group=SITE.x.x, color=SITE.x.x)) + geom_point()
#trying to create model vs incubation data plots
#try on a subset of points - too big right now
df_100 <- filter(df, run_num < 100)
ggplot(df_100, aes(x=CUE_x, y = CO2_of_tot, color = WHC)) + geom_point()
ggplot(df_100, aes(x=CUE_x, y = CO2_of_tot)) + geom_point() +
  geom_line(aes(color = WHC)) + facet_grid(.~SITE.x)
#plots of multiplier control of max CO2
ggplot(df_cost, aes(x=kMOD_x, y = max_model)) + geom_line()

### Ftn to find best psets (this is an objective function?)

best_psets <- function(df, best_n) {
  cost_df <- df %>% group_by(run_num) %>% summarize(SITE = unique(SITE.x.x), n = n(), cost = sum(cost))
  return(cost_df %>% slice_min(cost, n=best_n))
}

best_data <- function(df, best_n) {
  cost_df <- df %>% group_by(run_num) %>% summarize(SITE = unique(SITE.x.x),
                                                    n = n(),
                                                    cost = sum(cost),
                                                    min_MIM_CO = min(MIM_CO),
                                                    max_MIM_CO = max(MIM_CO))#,
  #min_MIC_SOC = min(MIC_SOC),
  #max_MIC_SOC = max(MIC_SOC))
  low_cost <- cost_df %>%
    filter(min_MIM_CO > 0.01) %>%
    filter(max_MIM_CO < 100) %>%
    #filter(min_MIC_SOC > 0.0001) %>%
    #filter(max_MIC_SOC < 0.40) %>%
    slice_min(cost, n=best_n)
  
  out <- df %>% filter(run_num %in% low_cost$run_num)
  return(out)
}




### Plot best psets
# df$SITE<- df$SITE.x
# BART.2 <- best_psets(df = df %>% filter(SITE == "2"), best_n = 100) #originally BART and 1000 for this and below
# BART.2_data <- best_data(df = df %>% filter(SITE == "2"), best_n = 100)
#
# plot_df <- df %>% filter(SITE == "2") %>%
#                   #filter(MAT == 15) %>%
#                   filter(run_num %in% BART.2$run_num)
#
# #not sure how Derek's data table is set up so not sure this code works for my purposes
# pairs(plot_df[,c(30,28,29, 11,12,6,7,10,9,24:27)],
#       col = alpha("black", 0.2), # Change color
#       pch = 16, # Change shape of points
#       cex=1.5,
#       labels = c("cost", "CO2_frac_of_tot", "LIT", "CO2-r", "CO2-K", "MICr", "MICK", "SOMa", "SOMc", "Vslope", "Vint", "Kslope", "Kint"), # Change labels of diagonal
#       main = "n=30 Lowest Cost Parameter Sets",
#       upper.panel = NULL,
#       breaks=c(0,0.2))
#
#
# ### 3D plot example
# library(plotly)
# plot_ly(plot_df, x = ~Vint_x, y = ~Vslope_x, z = ~Kslope_x, color = ~cost) #, colors = c('#BF382A', '#0C4B8E'))
# plot_ly(plot_df, x = ~Vint_x, y = ~Kint_x, z = ~-cost, color = ~cost)
# plot_ly(top_df, x = ~CUE_x, y = ~Tau_x, z = ~cost, color = ~SITE.x.x)
# plot_ly(df_cost, x = ~CUE_x, y = ~Tau_x, z = ~max_incub, color = ~SITE.x.x)
# #is the lack of difference because the model max is never reaching the incubation max?
#    #So, to match it will always be the highest point?

####
#this is most helpful code
####

### Collect all the best psets by site
df_cost$run_num<- df_cost$run_num.x
sites <- unique(df_cost$SITE.x.x)
top_df <- NULL
for(i in 1:length(sites)) {
  df1 <- best_data(df = df_cost %>% filter(SITE.x.x == sites[i]), best_n = 100)
  top_df <- rbind(top_df, df1)
}
#seems like cost function limits aren't working...
min(top_df$MIC_SOC)
mean(top_df$MIC_SOC)
max(top_df$MIC_SOC) #limited this so it wasn't greater than 40%, which is still pretty high
min(top_df$MIM_CO)
mean(top_df$MIM_CO)
max(top_df$MIM_CO)
#does cutoff make sense when we are starting with 0 soil C? Not sure it does... regardless, all the same across sites

### Ridge plots - had to add group to aesthetics and comment out fill to get to work
#likely because site is number, if make factor would likely work
#data_red<- filter(data, WHC == 60)
# top_df$SITE <- as.character(top_df$SITE.x.x)
# top_df <- inner_join(top_df, data, by = "SITE")
# pCost <- ggplot(top_df, aes(x = (cost/CO2C_prop)*100, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
#   geom_density_ridges(scale = 2) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   #scale_fill_brewer(palette = "Greys") +
#   labs(title="Cumulative %diff from incubation CO2-C target",
#        subtitle="n = 100 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# pCost
# #debugging
# ggplot(top_df, aes(x=CUE_x, y = (cost/CO2C_prop)*100, color = WHC)) + geom_point()
# ggplot(top_df, aes(x=CUE_x, y = (cost/CO2C_prop)*100)) + geom_point() +
#   geom_line(df, aes(color = WHC)) + facet_grid(.~SITE.x)
#
# pVint <- ggplot(top_df, aes(x = Vint_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
#   geom_density_ridges(scale = 2) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   #scale_fill_brewer(palette = "Blues") +
#   labs(title="Vint Multiplier") +#,
#           #subtitle="n=200 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# pVint
#
# pVslope <- ggplot(top_df, aes(x = Vslope_x, y = SITE, fill=SITE, group=SITE)) +
#   geom_density_ridges(scale = 2) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   #scale_fill_brewer(palette = "Blues") +
#   labs(title="Vslope Multiplier") +#,
#        #subtitle="n=200 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# pVslope
#
# pKslope <- ggplot(top_df, aes(x = Kslope_x, y = SITE, fill=SITE, group=SITE)) +
#   geom_density_ridges(scale = 2) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   #scale_fill_brewer(palette = "Oranges") +
#   labs(title="Kslope Multiplier") +#,
#        #subtitle="n=200 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# pKslope
#
# pKint <- ggplot(top_df, aes(x = Kint_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
#   geom_density_ridges(scale = 2) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   theme_ridges() +
#   #scale_fill_brewer(palette = "Oranges") +
#   labs(title="Kint Multiplier") +#,
#        #subtitle="n=200 lowest cost parameter sets") +
#   theme(legend.position = "none") #removes the legend
# pKint


pVmod <- ggplot(top_df, aes(x = vMOD_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="vMOD Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pVmod


pKmod <- ggplot(top_df, aes(x = kMOD_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="kMOD Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pKmod

pCUE <- ggplot(top_df, aes(x = CUE_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="CUE Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pCUE
#all the same for all 1000 parameter multipliers as well so not a top_df issue

pTAU <- ggplot(top_df, aes(x = Tau_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="Tau Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pTAU

#community comp changes with parameters
top_df$SITE <- top_df$SITE.x.x
top_df2 <- top_df %>% inner_join(data, by="SITE")
pCUE.r <- ggplot(top_df2, aes(x = CUE_r, y = MS_GWC, fill=MS_GWC, group=MS_GWC)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="kMOD_r Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pCUE.r

pCUE.k <- ggplot(top_df2, aes(x = CUE_k, y = MS_GWC, fill=MS_GWC, group=MS_GWC)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="kMOD_K Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pCUE.k

#determining the multiplication factor for R vs K
top_df$SITE <- top_df$SITE.x.x
top_df2 <- top_df %>% inner_join(data, by="SITE")
pMF.r <- ggplot(top_df2, aes(x = MF.r_x, y = MS_GWC, fill=MS_GWC, group=MS_GWC)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="MF_r Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pMF.r

pMF.k <- ggplot(top_df2, aes(x = MF.k_x, y = MS_GWC, fill=MS_GWC, group=MS_GWC)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="MF_K Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pMF.k
top_sum <- top_df %>% group_by(SITE.x.x) %>% summarise(mean.r = mean(MF.r_x), stdev.r = sd(MF.r_x), mean.k = mean(MF.k_x), stdev.k = sd(MF.k_x))
top_sum$SITE <- top_sum$SITE.x.x
top_sum.l <- pivot_longer(top_sum, c(2,4), names_to = "Comm", values_to = "mean")
data60 <- data %>% filter(WHC==60)
top_sum2<- inner_join(top_sum.l, data60, by = "SITE")
ggplot(top_sum2, aes(x = MS_GWC, y= mean, color = Comm, group = Comm)) +
  geom_point(size=2) +
  ylab("Mean multiplication factor") + xlab("Microsite GWC (%)") + theme_bw(base_size=16)
write.csv(top_sum, "top_sum_CUE.csv")

ggplot(top_df, aes(x=SITE.x.x, y=cost, group = run_num, color=run_num)) + geom_point()


top_row = ggarrange(pCost, ncol = 1, labels = c("A"))
bottom_row = ggarrange(
  pVint, pVslope,
  pKint, pKslope,
  pVMAX, pKM,
  nrow=3,
  ncol=2,
  labels = c("B", "C", "D", "E"))

png(file = "C:/github/MIMICS_MSBio/Cheyenne_HPC/HPC_output/best_pset_ridge_plot.png", width = 1000, height = 1200, res=108)
ggarrange(top_row, bottom_row,
          nrow=2,
          ncol=1,
          heights = c(1,3))
dev.off()

######
#Katie's code down here - looking for relationship between multipliers and soil moisture
##########
#summarize by mean and median
top_df$SITE <- as.character(top_df$SITE.x.x)
data$SITE <- as.character(data$SITE)
#Vint_sum <- top_df %>% group_by(SITE) %>%
#  summarise(
#                     mean.Vi=mean(Vint_x),
#                     med.Vi = median(Vint_x))
#Vslope_sum <- top_df %>% group_by(SITE) %>%
#  summarise(
#    mean.Vs=mean(Vslope_x),
#    med.Vs = median(Vslope_x))
#Kint_sum <- top_df %>% group_by(SITE) %>%
#  summarise(
#    mean.Ki=mean(Kint_x),
#    med.Ki = median(Kint_x))
#Kslope_sum <- top_df %>% group_by(SITE) %>%
#  summarise(
#    mean.Ks=mean(Kslope_x),
#    med.Ks = median(Kslope_x))
vMOD_sum <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.vMOD=mean(vMOD_x),
    med.vMOD = median(vMOD_x))
kMOD_sum <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.kMOD=mean(kMOD_x),
    med.kMOD = median(kMOD_x))
##all the same for CUE and Tau
CUE_sum <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.CUE=mean(CUE_x),
    med.CUE = median(CUE_x))
Tau_sum <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.tau=mean(Tau_x),
    med.tau = median(Tau_x))
#community comp summaries for parameter multipliers
CUE_sum.r <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.CUE.r=mean(CUE_r),
    med.CUE.r = median(CUE_r))
CUE_sum.k <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.CUE.k=mean(CUE_k),
    med.CUE.k = median(CUE_k))
vMOD_sum.r <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.vMOD.r=mean(vMOD_r),
    med.vMOD.r = median(vMOD_r))
vMOD_sum.k <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.vMOD.k=mean(vMOD_k),
    med.vMOD.k = median(vMOD_k))
#vMOD, kMOD
MC_sum <- vMOD_sum %>% #left_join(Vslope_sum, by = "SITE") %>%
  #left_join(Kint_sum, by = "SITE") %>%
  #left_join(Kslope_sum, by = "SITE") %>%
  #left_join(data, by = "SITE") %>%
  inner_join(kMOD_sum, by = "SITE") %>%
  inner_join(data, by = "SITE") %>%
  pivot_longer(cols=c(mean.vMOD,mean.kMOD), values_to = "mean", names_to = "Parameter") %>%
  pivot_longer(cols=c(med.vMOD,med.kMOD), values_to = "median", names_to = "Parameter2")
#vMOD, CUE or tau
MC_sum <- vMOD_sum %>%
  inner_join(CUE_sum, by = "SITE") %>%
  inner_join(data, by = "SITE") %>%
  pivot_longer(cols=c(mean.vMOD,mean.CUE), values_to = "mean", names_to = "Parameter") %>%
  pivot_longer(cols=c(med.vMOD,med.CUE), values_to = "median", names_to = "Parameter2")
#CUE, tau
MC_sum <- CUE_sum %>% #left_join(Vslope_sum, by = "SITE") %>%
  inner_join(Tau_sum, by = "SITE") %>%
  inner_join(data, by = "SITE") %>%
  pivot_longer(cols=c(mean.CUE,mean.tau), values_to = "mean", names_to = "Parameter") %>%
  pivot_longer(cols=c(med.CUE,med.tau), values_to = "median", names_to = "Parameter2")
#CUE or tau
MC_sum <- CUE_sum %>%
  inner_join(data, by = "SITE") %>%
  pivot_longer(cols=c(mean.CUE), values_to = "mean", names_to = "Parameter") %>%
  pivot_longer(cols=c(med.CUE), values_to = "median", names_to = "Parameter2")
#community change in r vs K
MC_sum <- CUE_sum.r %>%
  inner_join(CUE_sum.k, by = "SITE") %>%
  inner_join(vMOD_sum.r, by = "SITE") %>%
  inner_join(vMOD_sum.k, by = "SITE") %>%
  inner_join(data, by = "SITE") %>%
  pivot_longer(cols=c(mean.CUE.r,mean.CUE.k,mean.vMOD.r,mean.vMOD.k), values_to = "mean", names_to = "Parameter") %>%
  pivot_longer(cols=c(med.CUE.r,med.CUE.k,med.vMOD.r,med.vMOD.k), values_to = "median", names_to = "Parameter2")
#plots
ggplot(MC_sum, aes(x = GWC_site, y=mean, group=Parameter, color=Parameter)) + geom_point() +geom_smooth(method = "lm")
ggplot(MC_sum, aes(x = GWC_site, y=median, group=Parameter2, color=Parameter2)) + geom_point() +geom_smooth(method = "lm")
Tau_sum %>% inner_join(data, by = "SITE") %>% ggplot(aes(x = GWC_site, y=mean.tau)) + geom_point() +geom_smooth(method = "lm")
MC_sum$GWC_site <- MC_sum$GWC_site.x
library(plotly)
plot_ly(MC_sum, x = ~GWC_site, y = ~mean.CUE, z = ~mean.tau) #, colors = c('#BF382A', '#0C4B8E'))
#plot with low moisture removed - pretty similar actually!
ggplot(MC_sum %>% filter(GWC_site > 0.4), aes(x = GWC_site, y=mean, group=Parameter, color=Parameter)) + geom_point() +geom_smooth(method = "lm")
ggplot(MC_sum %>% filter(GWC_site > 0.4), aes(x = GWC_site, y=median, group=Parameter2, color=Parameter2)) + geom_point() +geom_smooth(method = "lm")
#equations
# MC_sum_Vs.mean <- filter(MC_sum, Parameter == "mean.Vs")
# MC_sum_Vs.med <- filter(MC_sum, Parameter2 == "med.Vs")
# Vs.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_Vs.mean)
# summary(Vs.mean_GWC)
# Vs.med_GWC <- lm(median ~ GWC_site, data=MC_sum_Vs.med)
# summary(Vs.med_GWC)
# MC_sum_Ks.mean <- filter(MC_sum, Parameter == "mean.Ks")
# MC_sum_Ks.med <- filter(MC_sum, Parameter2 == "med.Ks")
# Ks.mean_GWC <- lm(GWC ~ mean, data=MC_sum_Ks.mean)
# summary(Ks.mean_GWC)
# Ks.med_GWC <- lm(GWC ~ median, data=MC_sum_Ks.med)
# summary(Ks.med_GWC)
# MC_sum_Vi.mean <- filter(MC_sum, Parameter == "mean.Vi")
# MC_sum_Vi.med <- filter(MC_sum, Parameter2 == "med.Vi")
# Vi.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_Vi.mean)
# summary(Vi.mean_GWC)
# Vi.med_GWC <- lm(median ~ GWC_site, data=MC_sum_Vi.med)
# summary(Vi.med_GWC)
# MC_sum_Ki.mean <- filter(MC_sum, Parameter == "mean.Ki")
# MC_sum_Ki.med <- filter(MC_sum, Parameter2 == "med.Ki")
# Ki.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_Ki.mean)
# summary(Ki.mean_GWC)
# Ki.med_GWC <- lm(median ~ GWC_site, data=MC_sum_Ki.med)
# summary(Ki.med_GWC)
#vMOD and kMOD
MC_sum_vMOD.mean <- filter(MC_sum, Parameter == "mean.vMOD")
MC_sum_vMOD.med <- filter(MC_sum, Parameter2 == "med.vMOD")
vMOD.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_vMOD.mean)
summary(vMOD.mean_GWC)
vMOD.med_GWC <- lm(median ~ GWC_site, data=MC_sum_vMOD.med)
summary(vMOD.med_GWC)
MC_sum_kMOD.mean <- filter(MC_sum, Parameter == "mean.kMOD")
MC_sum_kMOD.med <- filter(MC_sum, Parameter2 == "med.kMOD")
kMOD.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_kMOD.mean)
summary(kMOD.mean_GWC)
kMOD.med_GWC <- lm(median ~ GWC_site, data=MC_sum_kMOD.med)
summary(kMOD.med_GWC)
#CUE and tau
MC_sum_CUE.mean <- filter(MC_sum, Parameter == "mean.CUE")
MC_sum_CUE.med <- filter(MC_sum, Parameter2 == "med.CUE")
CUE.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_CUE.mean)
summary(CUE.mean_GWC)
CUE.med_GWC <- lm(median ~ GWC_site, data=MC_sum_CUE.med)
summary(CUE.med_GWC)
MC_sum_tau.mean <- filter(MC_sum, Parameter == "mean.tau")
MC_sum_tau.med <- filter(MC_sum, Parameter2 == "med.tau")
tau.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_tau.mean)
summary(tau.mean_GWC)
tau.med_GWC <- lm(median ~ GWC_site, data=MC_sum_tau.med)
summary(tau.med_GWC)
#community comp
MC_sum_tau.mean <- filter(MC_sum, Parameter == "mean.Tau.r")
tau.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_tau.mean)
summary(tau.mean_GWC)
#plotting each separately
mean.Vs_plot <- MC_sum %>% filter(Parameter == "mean.Vs") %>% ggplot(aes(x = GWC_site, y=mean, group=Parameter)) +
  geom_point(aes(color = Parameter)) + geom_smooth(method = "lm", aes(color=Parameter))
mean.Vi_plot <- MC_sum %>% filter(Parameter == "mean.Vi") %>% ggplot(aes(x = GWC_site, y=mean, group=Parameter)) +
  geom_point(aes(color = Parameter)) + geom_smooth(method = "lm", aes(color=Parameter))
mean.Ks_plot <- MC_sum %>% filter(Parameter == "mean.Ks") %>% ggplot(aes(x = GWC_site, y=mean, group=Parameter)) +
  geom_point(aes(color = Parameter)) + geom_smooth(method = "lm", aes(color=Parameter))
mean.Ki_plot <- MC_sum %>% filter(Parameter == "mean.Ki") %>% ggplot(aes(x = GWC_site, y=mean, group=Parameter)) +
  geom_point(aes(color = Parameter)) + geom_smooth(method = "lm", aes(color=Parameter))
grid.arrange(mean.Vs_plot, mean.Vi_plot, mean.Ks_plot, mean.Ki_plot, ncol = 2, nrow = 2)
