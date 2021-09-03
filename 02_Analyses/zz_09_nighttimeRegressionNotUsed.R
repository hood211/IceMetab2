# Uses Nighttime regression approach to estimate K 
# JMH
# just a test, not cleaned

# libraries

library(tidyverse)

rd_S9_15 <- read.csv(file = file.path(here::here("01_Data/raw data for metab"),"use_stream9_2015_metabolism copy.csv"))
rd_S9_16 <- read.csv(file = file.path(here::here("01_Data/raw data for metab"),"use_stream9_2016_metabolism copy.csv"))
rd_S9_17 <- read.csv(file = file.path(here::here("01_Data/raw data for metab"),"use_stream9_2017_metabolism copy.csv"))

rd_S9 <- rbind(rd_S9_15, rd_S9_16, rd_S9_17) %>% 
  mutate(Pdt = as.POSIXct(Pdt, format = "%Y-%m-%d %H:%M:%S"),
         Y = as.numeric(strftime(Pdt, format = "%Y")))

ggplot(rd_S9 %>% 
         filter(Y == 2015), aes(y = persat, x = Pdt, color = site)) +
  geom_line() +
  facet_grid(site ~.)

ggplot(rd_S9 %>% 
         filter(Y == 2015), aes(y = lux, x = Pdt, color = site)) +
  geom_line() +
  facet_grid(site ~.)

ggplot(rd_S9, aes(y = temp, x = Pdt, color = site)) +
  geom_line()

ggplot(rd_S9 %>% 
         filter(Y == 2016), aes(y = persat, x = Pdt, color = site)) +
  geom_line() +
  facet_grid(site ~.)

ggplot(rd_S9 %>% 
         filter(Y == 2016), aes(y = lux, x = Pdt, color = site)) +
  geom_line() +
  facet_grid(site ~.)

ggplot(rd_S9 %>% 
         filter(Y == 2017), aes(y = persat, x = Pdt, color = site)) +
  geom_line() +
  facet_grid(site ~.)



osat  <- function(C,bp){
  sato<-(exp(2.00907 + 3.22014 * (log((298.15-C) / (273.15 + C))) + 4.0501 * (log((298.15 - C) / (273.15 + C))) ^ 2 + 4.94457 * (log((298.15 - C) / (273.15 + C))) ^ 3 - 0.256847 * (log((298.15 - C) / (273.15 + C))) ^ 4 + 3.88767 * (log((298.15 - C) / (273.15 + C))) ^ 5)) * 1.4276 * bp / 760
  
  sato
}


# NIGHT TIME REGRESSION FUNCTION; intercept = ER (g O2/m3/d); slope = K (1/d); mg DO/L = g DO/m3
nightreggression <- function(data, bp, time_bw_loggings_d){
  # data <- blah
  # time_bw_loggings_d <- 1/60/24 #in days
  # bp <- mean(blah$mmHg)
  
  temp<-data$C
  oxy<-data$DO
  
  ##moving average on oxy data
  oxyf1<- stats::filter(data$DO,rep(1/3,3), sides=2) # this takes a moving average of a DO meas, the prior meas, and the next meas.
  
  #trime the ends of the oxy data - removes the first and last value
  oxyf2<- oxyf1[c(-1,-length(oxyf1))]
  
  ##calculate delO/delt
  deltaO2<-((oxyf2[-1]-oxyf2[-length(oxyf2)])/time_bw_loggings_d)  # is t2 - t1 / time betwee readings in d
  
  #Trim the first two and last one from the temp data to match the filter oxy data
  temptrim<-temp[c(-2:-1,-length(temp))]
  #calc the dodef
  satdef<-osat(temptrim,bp)-oxyf2[-1]
  
  #calculate regression and plot
  nreg<-lm(deltaO2~satdef)
  plot(satdef,deltaO2)
  abline(nreg)
  
  summary(nreg)
}

blah <- rd_S9 %>% 
  filter(Pdt > as.POSIXct("2016-06-01 22:00:00", format = "%Y-%m-%d %H:%M:%S") &
           Pdt < as.POSIXct("2016-06-02 22:00:00", format = "%Y-%m-%d %H:%M:%S")) %>% 
  filter(site == "DS") %>% 
  filter(lux <1) %>%
  rename(C = temp, DO = DO_mgL) 


mean(blah$Q_DS_m3.s*1000)
plot(blah$lux ~ blah$Pdt)

nightreggression(blah, mean(blah$mmHg), 1/60/24)

# aghhhhh
library(lubridate)
rd_S92_ds <- rd_S9 %>% 
  select(Pdt, site, temp, lux, Q_DS_m3.s, velocity_m.s, mmHg, DO_mgL, depth_m) %>% 
  mutate(Date = as.Date(Pdt),
         Hour = as.numeric(strftime(Pdt, format = "%H"))) %>% 
  mutate(Date2 = as.Date(ifelse(Hour > 22, Date + 1, Date), origin='1970-01-01')) %>% 
  # filter(lux < 0.1) %>% 
  filter(lux == 0) %>% 
  filter(site == "DS") %>% 
  mutate(velocity_m.s = as.numeric(velocity_m.s),
         depth_m = as.numeric(depth_m))%>% 
  filter(!(as.character(Date2) %in% c("2015-07-10", "2016-06-03", "2016-07-09", "2016-07-10",
                                      "2017-05-29", "2017-07-12" )))

rd_S92_dsL3 <- rd_S92_ds %>%
  group_by(Date2) %>%
  summarise(count = n()) %>%
  filter(count <3)
  
NregResults <- matrix(nrow = length(unique(rd_S92$Date2)), ncol = 10)

# pdf("blah.pdf")
for(i in 1:length(unique(rd_S92_ds$Date2))){
  Date2_i = unique(rd_S92_ds$Date2)[i]
  rd_292_ds_i = rd_S92_ds[rd_S92_ds$Date2 == Date2_i,]
  time_bw_loggings_d = 1/60/24 # in days
  
  temp<-rd_292_ds_i$temp
  oxy<-rd_292_ds_i$DO_mgL
  bp <- rd_292_ds_i$mmHg
  depthM <- rd_292_ds_i$depth_m
  
  ##moving average on oxy data
  oxyf1<- stats::filter(oxy,rep(1/3,3), sides=2) # this takes a moving average of a DO meas, the prior meas, and the next meas.
  
  #trime the ends of the oxy data - removes the first and last value
  oxyf2<- oxyf1[c(-1,-length(oxyf1))]
  
  ##calculate delO/delt
  deltaO2_temp<-((oxyf2[-1]-oxyf2[-length(oxyf2)])/time_bw_loggings_d)  # is t2 - t1 / time betwee readings in d
  deltaO2 <- ifelse(deltaO2_temp < -30 | deltaO2_temp > 30, as.numeric("NA"), deltaO2_temp)
  
  #Trim the first two and last one from the temp data to match the filter oxy data
  temptrim<-temp[c(-2:-1,-length(temp))]
  bptrim<-bp[c(-2:-1,-length(bp))]
  
  #calc the dodef
  satdef<-osat(temptrim,bptrim)-oxyf2[-1]
  
  #calculate regression and plot
  nreg<-lm(deltaO2 ~ satdef)
  # plot(satdef,deltaO2)
  # abline(nreg)
  # summary(nreg)
  
  
  ER_mgO2_what_d_i <- coef(nreg)[[1]]
  temp_i <- mean(temp, na.rm = T)
  depth_i <- mean(depthM, na.rm = T)
  K_d_i  <- coef(nreg)[[2]]
  K_d_T20_i <- K_d_i/(1.024^(temp_i - 20))
  R2_i  <- summary(nreg)$r.squared
  Pval_i  <- summary(nreg)$coefficients[2,4]
  Q_Ls_i <- mean(rd_292_ds_i$Q_DS_m3.s, na.rm = T)*1000
  Vel_ms_i <- mean(rd_292_ds_i$velocity_m.s, na.rm = T)
  Results_i <- c(as.character(Date2_i), ER_mgO2_what_d_i, temp_i, depth_i, K_d_i, K_d_T20_i,  R2_i, Pval_i, Q_Ls_i, Vel_ms_i)
  NregResults[i,] <- Results_i
}
# dev.off()
  

NregResults2 <- as.data.frame(NregResults)

names(NregResults2) <- c("Date2", "ER_mgO2_what_d", "TempC", "depthM", "K_d", "K_d_T20", "R2", "Pvals", "Q_Ls", "Vel_ms")

NregResults2 <- NregResults2 %>% 
  mutate(Date2 = as.Date(Date2)) %>% 
  mutate(across(ER_mgO2_what_d:Vel_ms, as.numeric)) %>% 
  mutate(Y = as.character(strftime(Date2, format = "%Y")),
         ER_g02_m2_d = (ER_mgO2_what_d*1000)/depthM,
         DOY = as.numeric(strftime(Date2, format = "%j")))
  
  
ggplot(NregResults2 %>% 
         filter(DOY > 180), aes(y = ER_mgO2_what_d, x = DOY, color = Y)) +
  geom_point()

ggplot(NregResults2 %>% 
         filter(DOY > 180), aes(y = ER_mgO2_what_d, x = Y)) +
  geom_boxplot()

ggplot(NregResults2, aes(y = K_d, x = Date2)) +
    geom_point()
  
ggplot(NregResults2, aes(y = -K_d_T20, x = Q_Ls, color = Y)) +
  geom_point()

summary(lm(log(-K_d_T20) ~ log(Q_Ls), NregResults2 %>% 
                                filter(Pvals < 0.05)))

ggplot(NregResults2 %>% 
         filter(Pvals < 0.05), aes(y = -K_d_T20, x = Q_Ls, color = Y)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  stat_smooth(method = "lm")


ggplot(NregResults2, aes(y = -ER_mgO2_what_d, x = Q_Ls, color = Y)) +
  geom_point()

ggplot(NregResults2, aes(y = log(abs(K_d_T20)), x = log(Q_Ls), color = Y)) +
    geom_point()
  
  ggplot(NregResults2, aes(y = K_d, x = Vel_ms)) +
    geom_point()
  
  
# This is from 08_GasExchangePlots
  geS9 <- ge %>% 
    filter(stream == "st9") %>% 
    select(Date, Q_DS, temp_DS, Vms, K_O2_T20) %>% 
    mutate(meas = "SFSjmh") %>% 
    mutate(Date = as.Date(Date, format = "%m/%d/%y"))
  
  NregResults2sig <- NregResults2 %>% 
    filter(Pvals < 0.1) %>% 
    select(Date2, Q_Ls, TempC, Vel_ms, K_d_T20) %>% 
    mutate(meas = "NightReg",
           K_d_T20 = -K_d_T20)
  
  geS9b <- matrix(nrow = 93, ncol = 6)
  # Benoit data - suspect K_0 is at ambient temp (1/d)
  # K_d_i/(1.024^(temp_i - 20)
  # BOLD apr
  BOLDaprK <- 761/(1.024^(8 - 20)) #temp informed guess
  BOLDAugK <- 419/(1.024^(18 - 20)) #temp informed guess
  geS9b[5,] <- c("2008-08-01", 2.6, as.numeric("NA"), as.numeric("NA"), BOLDAugK, "PropBOLD")  
  geS9b[6,] <- c("2009-04-01", 5.4, as.numeric("NA"), as.numeric("NA"), BOLDaprK, "PropBOLD")  
  geS9b[1:4,] <- as.matrix(geS9)
  geS9b[7:93,] <- as.matrix(NregResults2sig)
  geS9b <- as.data.frame(geS9b)
  names(geS9b) <- c("Date", "Q_Ls", "TempC", "Vel_ms", "K_O2_20", "meas")
  geS9b <- geS9b %>% 
    mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d")) %>% 
    mutate(across(Q_Ls:K_O2_20, as.numeric))
  
  # comparison of different K estimates
  ggplot(data = geS9b, aes(y = K_O2_20, x = Q_Ls, color = meas)) +
    geom_point(data = metBS %>%
                 filter(metab_stream == "st9"), aes(y = median_Kr, x = meanQds*1000), color = "black") +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    stat_smooth(method = "lm", se = F)
    
  
  geS9bGasMeas <- geS9b %>% 
    filter(meas %in% c("PropBOLD", "SFSjmh"))
  
  
  summary(lm(log(K_O2_20) ~ log(Q_Ls),geS9bGasMeas))
  
  
  
  
  