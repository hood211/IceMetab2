# Tables
library(tidyverse)
library(kableExtra)

GPP_N_best <- read.csv(file.path(here::here("06_AICtables"), "03a_GPP_N_best.csv")) 
names(GPP_N_best) <- c("Int", "strTemp", "DOYresXstr", "light", "QresXstr", "dayTemp", "trt", "strTempXdayTemp",
                       "strTempXtrt", "DayTempXtrt", "strTempXDayTempXtrt", "Qres", "DOYres", "str", "DayTempXstr",
                       "trtXstr", "DayTempXtrtXstr", "df", "logLik", "AIC", "delta", "weight")

model.tab <- GPP_N_best %>% 
  select(Int, 
         strTempXDayTempXtrt, strTempXdayTemp, strTempXtrt, DayTempXtrt,
         DayTempXtrtXstr, DayTempXstr, trtXstr, 
         strTemp, dayTemp, trt, str,
         light, 
         DOYresXstr, DOYres, 
         QresXstr, Qres,
         df:weight)

model.tab.d5 <- model.tab %>% 
  filter(delta < 5)

TempInTopD5 <- max(ifelse(is.na(model.tab.d5$strTemp) == TRUE,1,0)) #if true 1 other 0 

model.tab.Temp <- if(TempInTopD5 > 0) {
  model.tab %>% 
    filter(!is.na(model.tab$strTemp)) %>% 
    mutate(delta2 = AIC - min(AIC)) %>% 
    filter(delta2 < 2) %>% 
    select(-delta2)
  
}

model.tab.final <-  rbind(model.tab.d5, model.tab.Temp) %>% 
  mutate(model = paste(
    ifelse(!is.na(strTempXDayTempXtrt),"strTempXDayTempXtrt + ",
           ifelse(is.na(strTempXDayTempXtrt) & !is.na(strTempXtrt) & !is.na(DayTempXtrt),"strTempXtrt + DayTempXtrt +",
                  ifelse(is.na(strTempXDayTempXtrt) & !is.na(strTempXtrt) & is.na(DayTempXtrt), "strTempXtrt + ",
                         ifelse(is.na(strTempXDayTempXtrt) & is.na(strTempXtrt) & !is.na(DayTempXtrt) & is.na(DayTempXtrtXstr),"DayTempXtrt +","")))),
    ifelse(!is.na(DayTempXtrtXstr),"DayTempXtrtXstr + ",
           ifelse(is.na(DayTempXtrtXstr) & !is.na(DayTempXstr) & !is.na(trtXstr),"DayTempXstr + trtXstr +",
                  ifelse(is.na(DayTempXtrtXstr) & !is.na(DayTempXstr) & is.na(trtXstr), "strTempXtrt + ",
                         ifelse(is.na(strTempXDayTempXtrt) & is.na(DayTempXtrtXstr) & !is.na(DayTempXtrt) & is.na(trtXstr),"DayTempXtrt +","")))),
    ifelse(is.na(DayTempXtrtXstr) & is.na(DayTempXstr) & is.na(trtXstr) & !is.na(trt),"trt + ",""),
    ifelse(is.na(strTempXDayTempXtrt) & is.na(strTempXtrt) & is.na(DayTempXtrt) & !is.na(strTemp),"strTemp + ",""),
    ifelse(is.na(DayTempXtrtXstr) & is.na(DayTempXstr) & is.na(strTempXDayTempXtrt) & is.na(DayTempXtrt) & !is.na(dayTemp),"dayTemp + ",""),
    ifelse(!is.na(DOYresXstr),"DOYresXstr + ",""),
    ifelse(!is.na(QresXstr),"QresXstr + ",""),
    ifelse(!is.na(DOYres),"DOYres + ",""),
    ifelse(!is.na(light), "L + ",""),
    ifelse(!is.na(Qres),"Q", "")))


model.tab.final %>% 
  select(df:model) %>% 
  kable()

