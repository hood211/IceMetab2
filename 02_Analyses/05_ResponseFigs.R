#Libraries
library(tidyverse)
library(gamm4)
library(mgcViz)
library(MuMIn)
library(forecast)
library(gridExtra)
library(tidymv)
library(ggpubr)
library(grid)
library(egg)
library(here)


# LOAD BEST MODELS
# going to focus on the temp models
# Ambient
ER_AMB.b <- readRDS(file = file.path(here::here("05_SavedModels"), "02d_ER_AMB_Best.rds")) 
GPP_AMB.t <- readRDS(file = file.path(here::here("05_SavedModels"), "03d_GPP_AMB_Best.rds"))
NEP_AMB.t <- readRDS(file = file.path(here::here("05_SavedModels"), "04d_NEP_AMB_Best.rds")) 

# Nitrogen
ER_N.b <- readRDS(file = file.path(here::here("05_SavedModels"), "02b_ER_N_Best.rds"))
GPP_N.t <- readRDS(file = file.path(here::here("05_SavedModels"), "03a_GPP_N_Best.rds"))
NEP_N.b <- readRDS(file = file.path(here::here("05_SavedModels"), "03a_NEP_N_Best.rds"))

# Phosphorus
ER_P.b <- readRDS(file = file.path(here::here("05_SavedModels"), "02a_ER_P_Best.rds"))
GPP_P.b <- readRDS(file = file.path(here::here("05_SavedModels"), "03a_GPP_P_Best.rds"))
NEP_P.b <- readRDS(file = file.path(here::here("05_SavedModels"), "03a_NEP_P_Best.rds"))


# Load data 
metAll <- readr::read_csv(file = file.path(here::here("01_Data"),"01_MundgedMetabDat.csv")) %>% 
  mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6")),
         treatment = fct_relevel(treatment, c("ambient", "nitrogen", "phosphorus")),
         treatment = fct_recode(treatment, "Ambient" = "ambient", "Nitrogen" = "nitrogen", "Phosphorus" = "phosphorus"),
         Yf = as.factor(Yf),
         streamTreat = as.factor(streamTreat))

metN <- metAll %>% 
  filter(treatment == "Ambient" | treatment == "Nitrogen") %>% 
  droplevels()

metP <- metAll %>% 
  filter(treatment == "Ambient" | treatment == "Phosphorus") %>% 
  droplevels()

metAMB <- metAll %>% 
  filter(treatment == "Ambient") %>% 
  droplevels()


# Predict models
#here's how these covariates change across streams
# median light varies a lot across streams - yikes!
COVARsumStr <- metAll %>% 
  select(stream, Qres, LightPerDay.lcg, invKT.C.StMean) %>% 
  group_by(stream) %>% 
  summarise(across(Qres:invKT.C.StMean, median))

COVARsum <- metAll %>% 
  select(stream, Qres, LightPerDay.lcg ) %>% 
  summarise(across(Qres:LightPerDay.lcg, median))

invKT.C.StMeanSEQ <- unique(round(c(COVARsumStr$invKT.C.StMean[1],
                                    COVARsumStr$invKT.C.StMean[2],
                                    COVARsumStr$invKT.C.StMean[3],
                                    COVARsumStr$invKT.C.StMean[4],
                                    seq(min(metAll$invKT.C.StMean), max(metAll$invKT.C.StMean), by = 0.01)),3))

# Predict models
# Use grand median across all streams so that we are only looking at daily temp, stream identity, and treatment effects WITHOUT COVARS
  # LightSet = COVARsum$LightPerDay.lcg  #grand median
  
# Use stream specific light, because there is so much variation in median light across treams
  LightSet = round(COVARsumStr$LightPerDay.lcg, 3)


# Ambient
  # https://cran.r-project.org/web/packages/tidymv/vignettes/predict-gam.html
P.ER_AMB.b <- predict_gam(ER_AMB.b$gam,  
                          values  = list(LightPerDay.lcg = LightSet,
                                         invKT.C.StMean = invKT.C.StMeanSEQ)) %>% 
                        mutate(met = "ER",
                               treatment = "ambient") %>% 
                        mutate(stream = as.factor(ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[1],3), "st11U",
                                                         ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[2],3), "st18",
                                                                ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[3],3), "st9",
                                                                       ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[4],3),"st6", "NoStream")))))) %>% 
                        mutate(treatment = fct_recode(treatment, Ambient = "ambient")) %>% 
                          filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                                   (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                                   (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                                   stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                        select(met, treatment, stream, invKT.C.StMean, fit:se.fit) %>% 
                        mutate(model = "AmbMod")
                        


P.GPP_AMB.t <- predict_gam(GPP_AMB.t$gam,  
                           values  = list(LightPerDay.lcg = LightSet,
                                          invKT.C.StMean = invKT.C.StMeanSEQ)) %>% 
                          mutate(met = "GPP",
                                 treatment = "ambient") %>% 
                          mutate(stream = as.factor(ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[1],3), "st11U",
                                                           ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[2],3), "st18",
                                                                  ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[3],3), "st9",
                                                                         ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[4],3),"st6", "NoStream")))))) %>% 
                          mutate(treatment = fct_recode(treatment, Ambient = "ambient")) %>% 
                            filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                                     (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                                     (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                                     stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                          select(met, treatment, stream, invKT.C.StMean, fit:se.fit)  %>% 
                          mutate(model = "AmbMod")

P.NEP_AMB.t <- predict_gam(NEP_AMB.t$gam,  
                           values  = list(LightPerDay.lcg = LightSet,
                                          invKT.C.StMean = invKT.C.StMeanSEQ)) %>% 
                            mutate(met = "NEP",
                                   treatment = "ambient") %>% 
                            mutate(stream = as.factor(ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[1],3), "st11U",
                                                             ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[2],3), "st18",
                                                                    ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[3],3), "st9",
                                                                           ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[4],3),"st6", "NoStream")))))) %>% 
                            mutate(treatment = fct_recode(treatment, Ambient = "ambient")) %>% 
                              filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                                       (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                                       (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                                       stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                            select(met, treatment, stream, invKT.C.StMean, fit:se.fit)  %>% 
                            mutate(model = "AmbMod")




# Nitrogen
P.ER_N.b <- predict_gam(ER_N.b$gam,  
                        values  = list(LightPerDay.lcg = LightSet,
                                       invKT.C.StMean = invKT.C.StMeanSEQ)) %>% 
                        mutate(met = "ER") %>% 
                        mutate(stream = as.factor(ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[1],3), "st11U",
                                                         ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[2],3), "st18",
                                                                ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[3],3), "st9",
                                                                       ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[4],3),"st6", "NoStream")))))) %>% 
                        mutate(treatment = fct_recode(treatment, Ambient = "ambient", Nitrogen = "nitrogen")) %>% 
                          filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                                   (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                                   (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                                   stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                        select(met, treatment, stream, invKT.C.StMean, fit:se.fit)  %>% 
                        mutate(model = "NMond")

P.GPP_N.t <- predict_gam(GPP_N.t$gam,  
                    values  = list(LightPerDay.lcg = LightSet,
                                   invKT.C.StMean = invKT.C.StMeanSEQ)) %>% 
                    mutate(met = "GPP") %>% 
                    mutate(stream = as.factor(ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[1],3), "st11U",
                                                     ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[2],3), "st18",
                                                            ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[3],3), "st9",
                                                                   ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[4],3),"st6", "NoStream")))))) %>% 
                    mutate(treatment = fct_recode(treatment, Ambient = "ambient", Nitrogen = "nitrogen")) %>% 
                      filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                               (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                               (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                               stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                    select(met, treatment, stream, invKT.C.StMean, fit:se.fit)   %>% 
                    mutate(model = "NMond")

P.NEP_N.b <- predict_gam(NEP_N.b$gam,  
                         values  = list(LightPerDay.lcg = LightSet,
                                        invKT.C.StMean = invKT.C.StMeanSEQ)) %>% 
                    mutate(met = "NEP") %>% 
                    mutate(stream = as.factor(ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[1],3), "st11U",
                                                     ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[2],3), "st18",
                                                            ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[3],3), "st9",
                                                                   ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[4],3),"st6", "NoStream")))))) %>% 
                    mutate(treatment = fct_recode(treatment, Ambient = "ambient", Nitrogen = "nitrogen")) %>% 
                      filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                               (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                               (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                               stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                    select(met, treatment, stream, invKT.C.StMean, fit:se.fit)   %>% 
                    mutate(model = "NMond")


# Phosphorus
P.ER_P.b <- predict_gam(ER_P.b$gam,  
                       values  = list(LightPerDay.lcg = LightSet)) %>% 
                       mutate(met = "ER") %>% 
                       mutate(treatment = fct_recode(treatment, Ambient = "ambient", Phosphorus = "phosphorus")) %>% 
                       mutate(invKT.C.StMean = as.numeric(ifelse(stream == "st11U", COVARsumStr$invKT.C.StMean[1],
                                                                 ifelse(stream == "st18", COVARsumStr$invKT.C.StMean[2],
                                                                        ifelse(stream == "st9", COVARsumStr$invKT.C.StMean[3],
                                                                               ifelse(stream == "st6", COVARsumStr$invKT.C.StMean[4], "blah")))))) %>% 
                        filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                                 (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                                 (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                                 stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                       select(met, treatment, stream, invKT.C.StMean, fit:se.fit)   %>% 
                       mutate(model = "PMond")

P.GPP_P.b <- predict_gam(GPP_P.b$gam,  
                         values  = list(LightPerDay.lcg = LightSet)) %>% 
                        mutate(met = "GPP") %>% 
                        mutate(treatment = fct_recode(treatment, Ambient = "ambient", Phosphorus = "phosphorus")) %>% 
                        mutate(invKT.C.StMean = as.numeric(ifelse(stream == "st11U", COVARsumStr$invKT.C.StMean[1],
                                                                  ifelse(stream == "st18", COVARsumStr$invKT.C.StMean[2],
                                                                         ifelse(stream == "st9", COVARsumStr$invKT.C.StMean[3],
                                                                                ifelse(stream == "st6", COVARsumStr$invKT.C.StMean[4], "blah")))))) %>% 
                          filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                                   (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                                   (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                                   stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                        select(met, treatment, stream, invKT.C.StMean, fit:se.fit)   %>% 
                        mutate(model = "PMond")

P.NEP_P.b <- predict_gam(NEP_P.b$gam,  
                        values  = list(LightPerDay.lcg = LightSet)) %>% 
                        mutate(met = "NEP") %>% 
                        mutate(treatment = fct_recode(treatment, Ambient = "ambient", Phosphorus = "phosphorus")) %>% 
                        mutate(invKT.C.StMean = as.numeric(ifelse(stream == "st11U", COVARsumStr$invKT.C.StMean[1],
                                                                  ifelse(stream == "st18", COVARsumStr$invKT.C.StMean[2],
                                                                         ifelse(stream == "st9", COVARsumStr$invKT.C.StMean[3],
                                                                                ifelse(stream == "st6", COVARsumStr$invKT.C.StMean[4], "blah")))))) %>% 
                          filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
                                   (stream == "st18" & LightPerDay.lcg == LightSet[2])|
                                   (stream == "st9" & LightPerDay.lcg == LightSet[3])|
                                   stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
                        select(met, treatment, stream, invKT.C.StMean, fit:se.fit)   %>% 
                        mutate(model = "PMond")


# MAKE DATAFRAME FOR RESPONSE FIGS
ResponseFigdf <- rbind(P.ER_AMB.b, P.GPP_AMB.t, P.NEP_AMB.t,
                       P.ER_N.b, P.GPP_N.t, P.NEP_N.b, 
                       P.ER_P.b, P.GPP_P.b, P.NEP_P.b) %>% 
  mutate(stream = fct_recode(stream, S11 = "st11U", S18 = "st18", S9 = "st9", S6 = "st6"),
         met = as.factor(met),
         model = as.factor(model)) %>% 
  #back transform ER and GPP
  mutate(fit = ifelse(met != "NEP", exp(fit), fit),
         se.fit = ifelse(met != "NEP", exp(se.fit), se.fit))



# RESPONSE FIGURES
metAllres <- metAll %>% 
    select(Pdt, stream, treatment, lGPP, lER, NEP, invKT.C.StMean, LightPerDay.lcg) %>% 
   pivot_longer(cols = c(lGPP, lER, NEP), names_to = "met", values_to = "values") %>% 
  # ER and GPP are logged in both datasets
  mutate(met = fct_recode(met, ER = "lER", GPP = "lGPP")) 

plot(metAll$meanTemp ~ metAll$invKT.C.StMean)
TempBackTransform <- lm(meanTemp ~ invKT.C.StMean, metAll); summary(TempBackTransform)

metAllresP <- rbind(metAllres %>% 
                      filter(treatment == "Ambient") %>% 
                      mutate(model = "AmbMod"),
                    metAllres %>% 
                      filter(treatment == "Ambient" | treatment == "Nitrogen") %>% 
                      mutate(model = "NMond"),
                    metAllres %>% 
                      filter(treatment == "Ambient" | treatment == "Phosphorus") %>% 
                      mutate(model = "PMond")) %>% 
                mutate(TempC = coef(TempBackTransform)[[1]] + coef(TempBackTransform)[[2]] * invKT.C.StMean)

ResponseFigdf2 <- ResponseFigdf %>% 
  mutate(TempC = coef(TempBackTransform)[[1]] + coef(TempBackTransform)[[2]] * invKT.C.StMean)


# updated 28 July 21
png(file.path(here("03_Plots"), "FigS1_ResponsePlot_20210728.png"), units = "in", height = 7, width = 8, res = 300)
ggplot(data = ResponseFigdf2 %>%
         filter(stream != "NoStream") %>% 
         mutate(met = fct_relevel(met, c("GPP", "ER", "NEP")),
                treatment = fct_relevel(treatment, c("Ambient", "Phosphorus", "Nitrogen")),
                model = fct_recode(model, "Ambient" = "AmbMod", "Nitrogen" = "NMond", "Phosphorus" = "PMond"),
                model = fct_relevel(model, c("Ambient", "Phosphorus", "Nitrogen"))) %>% 
         rename(Treatment = treatment)) +
  # This doesn't loook great, but reassuring that they are inline
  # geom_point(data = metAllresP, aes(y = values, x = -invKT.C.StMean, color = treatment), size = 0.25, position = "jitter", alpha = 25/100) +
  geom_errorbar(aes(ymin = fit - se.fit, ymax = fit +se.fit, x = TempC, color = Treatment), width = 0.05) +
  geom_point(aes(y = fit, x = TempC, fill = Treatment), shape = 21, size = 5, color = "black") +
  # stat_smooth(data = ResponseFigdf2 %>% 
  #               filter(stream != "NoStream"), aes(y = fit, x = TempC, color = treatment), size = 0.5, alpha = 40/100) +
  facet_grid(met ~ model, scales = "free_y", labeller = labeller(model = label_wrap_gen(5)))  +
  scale_x_continuous(limits = c(9,16), breaks = c(10, 12, 14, 16)) +
  # scale_color_manual(values = c("white", "white", "white")) +
  scale_color_manual(values = c("grey60", "steelblue1", "gold1")) +
  scale_fill_manual(values = c("grey60", "steelblue1", "gold1")) +
  xlab("Mean stream temperature (°C)") +
  ylab(expression(paste("Response (g C ",m^-2," ", d^-1,")"))) +
  theme(legend.position = "top",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.width=unit(1,"cm"),
        legend.box = "horizontal",
        legend.background = element_rect(color = "black"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "transparent", color = "transparent", size = 2))
dev.off()


# EFFECT SIZE PLOTS
# make N trt dataframe
EffectPlotNtrt <- ResponseFigdf2 %>% 
                  filter(model == "NMond") %>% 
                  pivot_wider(id_cols = c(met, stream, TempC), names_from = treatment, values_from = c(fit,se.fit)) %>% 
                  mutate(ES = ifelse(met == "NEP", fit_Nitrogen - fit_Ambient, ((fit_Nitrogen/fit_Ambient)-1)*100),
                         trt = "Nitrogen") %>% 
                  mutate(met = as.character(met))

# est error on effect sizes
          # 
          EffectPlotNtrt2 <- EffectPlotNtrt %>% 
            mutate(ES.lowCI = as.numeric("NA"),
                   ES.upperCI = as.numeric("NA"))
          
# use random draws from se to estimate error
          for(i in 1:dim(EffectPlotNtrt)[1]) {
            # i = 1
            Row_i <- EffectPlotNtrt[i,]
            met_i <- as.character(EffectPlotNtrt[i,"met"])
            Amb_est <- rnorm(500, 
                             ifelse(met_i == "NEP",Row_i$fit_Ambient, log(Row_i$fit_Ambient)),
                             ifelse(met_i == "NEP",Row_i$se.fit_Ambient, log(Row_i$se.fit_Ambient)))
            Trt_est <- rnorm(500, 
                             ifelse(met_i == "NEP",Row_i$fit_Nitrogen, log(Row_i$fit_Nitrogen)),
                             ifelse(met_i == "NEP",Row_i$se.fit_Nitrogen, log(Row_i$se.fit_Nitrogen)))
            ESdf <- as.data.frame(as.matrix(cbind(Amb_est, Trt_est))) %>% 
              mutate(met = met_i,
                     ES = ifelse(met == "NEP", Trt_est - Amb_est, (exp(Trt_est)/exp(Amb_est)-1)*100)) 
            EffectPlotNtrt2[i,]$ES.lowCI <- quantile(ESdf$ES, probs = 0.1)
            EffectPlotNtrt2[i,]$ES.upperCI <- quantile(ESdf$ES, probs = 0.9)
          }

        ggplot(EffectPlotNtrt2, aes(y = ES, x = TempC)) +
          geom_point() +
          geom_errorbar(aes(ymin = ES.lowCI, ymax = ES.upperCI, x = TempC)) +
          facet_wrap(vars(met))


EffectPlotPtrt <- ResponseFigdf2 %>% 
                    filter(model == "PMond") %>% 
                    pivot_wider(id_cols = c(met, stream, TempC), names_from = treatment, values_from = c(fit,se.fit)) %>% 
                    mutate(ES = ifelse(met == "NEP", fit_Phosphorus - fit_Ambient,
                                       ((fit_Phosphorus/fit_Ambient)-1)*100),
                           trt = "Phosphorus")%>% 
                  mutate(met = as.character(met))

        EffectPlotPtrt2 <- EffectPlotPtrt %>% 
          mutate(ES.lowCI = as.numeric("NA"),
                 ES.upperCI = as.numeric("NA"))
        
        # use random draws from se to estimate error       
        for(i in 1:dim(EffectPlotPtrt2)[1]) {
          # i = 1
          Row_i <- EffectPlotPtrt2[i,]
          met_i <- as.character(EffectPlotPtrt2[i,"met"])
          Amb_est <- rnorm(500, 
                           ifelse(met_i == "NEP",Row_i$fit_Ambient, log(Row_i$fit_Ambient)),
                           ifelse(met_i == "NEP",Row_i$se.fit_Ambient, log(Row_i$se.fit_Ambient)))
          Trt_est <- rnorm(500, 
                           ifelse(met_i == "NEP",Row_i$fit_Phosphorus, log(Row_i$fit_Phosphorus)),
                           ifelse(met_i == "NEP",Row_i$se.fit_Phosphorus, log(Row_i$se.fit_Phosphorus)))
          ESdf <- as.data.frame(as.matrix(cbind(Amb_est, Trt_est))) %>% 
            mutate(met = met_i,
                   ES = ifelse(met == "NEP", Trt_est - Amb_est, ((exp(Trt_est)/exp(Amb_est)-1))*100)) 
          EffectPlotPtrt2[i,]$ES.lowCI <- quantile(ESdf$ES, probs = 0.1)
          EffectPlotPtrt2[i,]$ES.upperCI <- quantile(ESdf$ES, probs = 0.9)
        }
        
        ggplot(EffectPlotPtrt2, aes(y = ES, x = TempC)) +
          geom_point() +
          geom_errorbar(aes(ymin = ES.lowCI, ymax = ES.upperCI, x = TempC)) +
          facet_wrap(vars(met))

#####################
        # Figure 2
#####################

EffectPlotdf <- rbind(EffectPlotNtrt2 %>% 
                        select(met:TempC,trt,ES, ES.lowCI, ES.upperCI),
                      EffectPlotPtrt2 %>% 
                        select(met:TempC,trt, ES, ES.lowCI, ES.upperCI)) %>% 
                mutate(met = fct_relevel(met, c("GPP", "ER", "NEP")),
                       trt = fct_relevel(trt, c("Phosphorus", "Nitrogen")))


# scale_color_manual(values = c("grey60", "gold1", "steelblue1")) +

Plt.ESgpp <- ggplot(EffectPlotdf %>% 
                      filter(met == "GPP"), aes(y = ES, x = TempC, fill = trt)) +
  geom_errorbar(aes(ymin = ES.lowCI, ymax = ES.upperCI, x = TempC)) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = c("steelblue1", "gold1")) +
  scale_x_continuous(limits = c(9,16), breaks = c(10, 12, 14, 16)) +
  facet_grid(met ~ trt) +
  ylim(-50,800)+
  # ylab(expression(atop("Percent increase", "during enrichment"))) +
  ylab("Percent increase") +
  # geom_hline(yintercept = 1, color = "black", linetype = "dashed")+
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "transparent", color = "transparent", size = 2)) 


Plt.ES1er <- ggplot(EffectPlotdf %>% 
                      filter(met == "ER"), aes(y = ES, x = TempC, fill = trt)) +
  geom_errorbar(aes(ymin = ES.lowCI, ymax = ES.upperCI, x = TempC)) +            
  geom_point(size = 5, shape = 21) +
  facet_grid(met ~ trt) +
  ylim(-50,800) +
  ylab(NULL) +
  # geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  scale_fill_manual(values = c("steelblue1", "gold1")) +
  scale_x_continuous(limits = c(9,16), breaks = c(10, 12, 14, 16)) +
  # ylab(expression(atop("Percent increase", "during enrichment"))) +
  ylab("Percent increase") +
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 22, face = "bold"),
        strip.background.y = element_rect(fill = "transparent", color = "transparent", size = 2)) 



Plt.ESnep <- ggplot(EffectPlotdf %>% 
                      filter(met == "NEP"), aes(y = ES, x = TempC, fill = trt)) +
  geom_errorbar(aes(ymin = ES.lowCI, ymax = ES.upperCI, x = TempC)) +
  geom_point(size = 5, shape = 21) +
  facet_grid(met ~ trt) +
  scale_fill_manual(values = c("steelblue1", "gold1")) +
  scale_x_continuous(limits = c(9,16), breaks = c(10, 12, 14, 16)) +
  ylim(-8,1)+
  ylab(expression(atop("Change in NEP", paste("(g C ", m^-2," ", d^-1,")")))) +
  # ylab(expression(paste("Change in NEP (g C", m^-2, " ", d^-1,")"))) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 22, face = "bold"),
        strip.background.y = element_rect(fill = "transparent", color = "transparent", size = 2)) 


Fig2a <- ggarrange(Plt.ESgpp, Plt.ES1er, Plt.ESnep, nrow = 3, ncol = 1)
Fig2 <- annotate_figure(Fig2a, 
                        # left = text_grob("Effect size", rot = 90, size = 22),
                        bottom = text_grob("Mean stream temperature (°C)", size = 18))

png(file.path(here("03_Plots"), "Fig2_EffectSizePlot_20210728.png"), units = "in", height = 7, width = 8, res = 300)
Fig2
dev.off()

# save.image(file.path(here("04_SavedImages"), "05_ResponseFigs_Rdat"))
# # 
# load(file.path(here("04_SavedImages"), "05_ResponseFigs_Rdat"))






# FIGURING OUT AE WITH INTERACTIONS
NEP_P.t <- readRDS(file = file.path(here::here("05_SavedModels"), "03a_NEP_P_BestTemp.rds"))
# AE under amb is just the coef for invKt, while under treatment is that coef + the invKt*treat interaction
summary(NEP_P.t$gam)


# AVERAGE INCREASES FOR PAPER
EffectPlotdf %>% 
  filter(met == "NEP" & trt == "Phosphorus") %>% 
  group_by(met) %>% 
  summarise(ES = mean(ES))

# How much did NEP increase with temp under amb and +N
ResponseFigdf2 %>% 
  filter(treatment == "Ambient" & met == "NEP" & model == "NMond")















# LIGHT RELATIONSHIPS - NOT USED
ggplot(metAll, aes(y = LightPerDay.lcg, x = julian, color = treatment)) +
  geom_point()+
  facet_wrap(vars(stream))

plot(metAll$LightPerDay ~ metAll$meanTemp)

summary(lm(LightPerDay ~ meanTemp+stream, metAll))

ggplot(metAll, aes(y = LightPerDay, x = meanTemp)) +
  geom_point() +
  facet_wrap(vars(stream))

ggplot(metAll, aes(y = LightPerDay, x = meanQds)) +
  geom_point() +
  facet_wrap(vars(stream))

P.ER_N.bL <- predict_gam(ER_N.b$gam,  
                        values  = list(LightPerDay.lcg = seq(min(metAll$LightPerDay.lcg), max(metAll$LightPerDay.lcg), length = 100),
                                       invKT.C.StMean = invKT.C.StMeanSEQ)) %>% 
  mutate(met = "ER") %>% 
  mutate(stream = as.factor(ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[1],3), "st11U",
                                   ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[2],3), "st18",
                                          ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[3],3), "st9",
                                                 ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[4],3),"st6", "NoStream")))))) %>% 
  mutate(treatment = fct_recode(treatment, Ambient = "ambient", Nitrogen = "nitrogen")) %>% 
  # filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
  #          (stream == "st18" & LightPerDay.lcg == LightSet[2])|
  #          (stream == "st9" & LightPerDay.lcg == LightSet[3])|
  #          stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
  select(met, treatment, stream, invKT.C.StMean,LightPerDay.lcg, fit:se.fit)  %>% 
  mutate(model = "NMond")

ggplot(P.ER_N.bL, aes(y = fit, x = LightPerDay.lcg, color = treatment)) +
  geom_point() +
  facet_wrap(model ~ stream) +
  geom_point(data = metAll %>% 
               filter(treatment != "Phosphorus"), aes(y = lER, x = LightPerDay.lcg, color = treatment))


P.GPP_N.tL <- predict_gam(GPP_N.t$gam,  
                         values  = list(LightPerDay.lcg = seq(min(metAll$LightPerDay.lcg), max(metAll$LightPerDay.lcg), length = 100),
                                        invKT.C.StMean = invKT.C.StMeanSEQ)) %>% 
  mutate(met = "ER") %>% 
  mutate(stream = as.factor(ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[1],3), "st11U",
                                   ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[2],3), "st18",
                                          ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[3],3), "st9",
                                                 ifelse(invKT.C.StMean == round(COVARsumStr$invKT.C.StMean[4],3),"st6", "NoStream")))))) %>% 
  mutate(treatment = fct_recode(treatment, Ambient = "ambient", Nitrogen = "nitrogen")) %>% 
  # filter((stream == "st11U" & LightPerDay.lcg == LightSet[1])|
  #          (stream == "st18" & LightPerDay.lcg == LightSet[2])|
  #          (stream == "st9" & LightPerDay.lcg == LightSet[3])|
  #          stream == "st6" & LightPerDay.lcg == LightSet[4]) %>% 
  select(met, treatment, stream, invKT.C.StMean,LightPerDay.lcg, fit:se.fit)  %>% 
  mutate(model = "NMond")

ggplot(P.GPP_N.tL, aes(y = fit, x = LightPerDay.lcg, color = treatment)) +
  geom_point() +
  facet_wrap(model ~ stream) +
  geom_point(data = metAll %>% 
               filter(treatment != "Nitrogen"), aes(y = lER, x = LightPerDay.lcg, color = treatment))


gv.ER_AMB <- getViz(ER_AMB.b$gam)
plot( gv.ER_AMB) + 
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

