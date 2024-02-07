LiDET_BAGS <- data.frame(TYPE = c('TRAEf', 'PIREf','THPLf','ACSAf','QUPRf','DRGLf'),
                         BAG_LIG = c(16.2, 19.2, 26.7, 15.9, 23.5, 10.9),
                         BAG_N = c(0.38, 0.59, 0.62, 0.81, 1.03, 1.97), 
                         BAG_CN = c(133.3,92.7, 83.1, 61.8, 50.5, 24.2))

LiDET_BAGS$CALC_N <- (1 / LiDET_BAGS$BAG_CN) / 2.5 * 100   
LiDET_BAGS$CALC_MET <- 0.85 - 0.013 * LiDET_BAGS$BAG_LIG/LiDET_BAGS$CALC_N

BAG_init_size <- 100
BAGS <- LiDET_BAGS %>% select(TYPE, CALC_MET)
BAGS$BAG_LITm <- BAG_init_size * BAGS$CALC_MET  
BAGS$BAG_LITs <- BAG_init_size * (1-BAGS$CALC_MET) 