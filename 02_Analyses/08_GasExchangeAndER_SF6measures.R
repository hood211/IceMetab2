
###################
# LIBRARIES
###################
library(tidyverse)
library(ggpubr)

###################
# DATA
###################

# summary of all gas exchange measurements in Hengill valley and beyond
ge2 <-  readr::read_csv(file = file.path(here::here("01_Data"),"GasExchangeEstimatesAndQttEct.csv")) 

# focus on streams with metabolism data
ge <- ge2 %>% 
  filter(stream %in% c("st6", "st9", "st11U")) %>% 
  mutate(K_O2_T20.log = log(K_O2_T20),
         Q_DS.log = log(Q_DS))

###################
# RELATIONSHIPS BETWEEN K AND Q
###################

# S6
S6lm <- lm(K_O2_T20.log ~ Q_DS.log, data = ge %>% 
             filter(stream == "st6"))
summary(S6lm)

# S9
S9lm <- lm(K_O2_T20.log ~ Q_DS.log, data = ge %>% 
             filter(stream == "st9"))
summary(S9lm)

# S11 - no relationship
S11Ulm <- lm(K_O2_T20.log ~ Q_DS.log, data = ge %>% 
               filter(stream == "st11U"))
summary(S11Ulm)


# MIN AND MAX Q AND GE
S6minQ <- exp(min(ge[ge$stream == "st6",]$Q_DS.log))
S6maxQ <- exp(max(ge[ge$stream == "st6",]$Q_DS.log))
S6minGE <- exp(3.202 + 1.097*log(S6minQ))
S6maxGE <- exp(3.202 + 1.097*log(S6maxQ))

S9minQ <- exp(min(ge[ge$stream == "st9",]$Q_DS.log))
S9maxQ <- exp(max(ge[ge$stream == "st9",]$Q_DS.log))
S9minGE <- exp(4.7003 + 1.6980*log(S9minQ))
S9maxGE <- exp(4.7003 + 1.6980*log(S9maxQ))

S11UminQ <- exp(min(ge[ge$stream == "st11U",]$Q_DS.log))
S11UmaxQ <- exp(max(ge[ge$stream == "st11U",]$Q_DS.log))
S11UminGE <- exp(5.56939 + 0.02863*log(S11UminQ))
S11UmaxGE <- exp(5.56939 + 0.02863*log(S11UmaxQ))


# Plots for supplemental
p1a <- ggplot(ge %>% 
                mutate(stream = fct_recode(stream, "S11" = "st11U",
                                           "S6" = "st6",
                                           "S9" = "st9"),
                       stream = fct_relevel(stream, c("S11", "S9", "S6"))), aes(y = K_O2_T20, x  = Q_DS, color = stream)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("forestgreen", "blue", "purple")) +
  scale_x_log10(limits = c(1,30)) +
  scale_y_log10(limits = c(100,700)) +
  geom_segment(aes(x = S6minQ, xend = S6maxQ, y =S6minGE, yend = S6maxGE), color = "purple") +
  geom_segment(aes(x =S9minQ , xend = S9maxQ, y = S9minGE, yend = S9maxGE), color = "blue") +
  geom_segment(aes(x =S11UminQ , xend = S11UmaxQ, y = S11UminGE, yend = S11UmaxGE), linetype = "dashed", color = "forestgreen") +
  xlab(expression(paste("Discharge (L ",s^-1, ")"))) +
  ylab(expression(paste(O[2]," exchange rate at 20°C (",d^-1,")"))) +
  theme(legend.position = "top",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.width=unit(1,"cm"),
        legend.box = "horizontal",
        legend.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "transparent", color = "transparent", size = 2)) +
  annotate("text", x = 1, y = 650, label = "a", size = 12, fontface = "bold")



# for comparing GE measurement range with met measurement range
met <- readr::read_csv(file = file.path(here::here("01_Data"),"01_MundgedMetabDat.csv")) %>% 
  mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6")),
         treatment = fct_relevel(treatment, c("ambient", "nitrogen", "phosphorus")),
         Yf = as.factor(Yf),
         streamTreat = as.factor(streamTreat))


# S6 looks pretty good
summary(met[met$stream == "st6",]$meanQds)*1000
summary(ge[ge$stream == "st6",]$Q_DS)

# S9 - get's close to the interquartile range, but not the max
# maybe a bit of a problem in the N and P year
summary(met[met$stream == "st9",]$meanQds)*1000
summary(ge[ge$stream == "st9",]$Q_DS)
ggplot(met %>% 
         filter(stream == "st9"), aes(y = meanQds*1000, x = julian, color = treatment)) +
  geom_point() +
  geom_abline(intercept = 2.9,slope= 0)

#st11 U - looks good
summary(met[met$stream == "st11U",]$meanQds)*1000
summary(ge[ge$stream == "st11U",]$Q_DS)

ggplot(met %>% 
         filter(stream == "st11U"), aes(y = meanQds*1000, x = julian, color = treatment)) +
  geom_point() +
  geom_abline(intercept = 10.9,slope= 0)

metQ <- met %>% 
  filter(stream %in% c("st11U", "st9", "st6")) %>% 
  select(stream, Q_DS = "meanQds") %>% 
  mutate(Q_DS = Q_DS * 1000,
         meas = "met")

geQ <- ge %>% 
  select(stream, Q_DS) %>% 
  mutate(meas = "GE")

Qcomp <- rbind(metQ, geQ)

p1b <- ggplot() +
  geom_boxplot(data= Qcomp %>% 
                 filter(meas == "met") %>% 
                 mutate(stream = fct_recode(stream, "S11" = "st11U",
                                            "S6" = "st6",
                                            "S9" = "st9"),
                        stream = fct_relevel(stream, c("S11", "S9", "S6"))), aes(y = Q_DS, x = stream, fill = stream)) +
  scale_fill_manual(values = c("forestgreen", "blue", "purple")) +
  geom_point(data= Qcomp %>% 
               filter(meas == "GE") %>% 
               mutate(stream = fct_recode(stream, "S11" = "st11U",
                                          "S6" = "st6",
                                          "S9" = "st9")), aes(y = Q_DS, x = stream), shape = 22, fill = "grey", size = 4) +
  ylab(expression(paste("Discharge (L ",s^-1,")"))) +
  xlab(NULL) +
  ylim(0,30)+
  theme(legend.position = "top",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.width=unit(1,"cm"),
        legend.box = "horizontal",
        legend.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "transparent", color = "transparent", size = 2)) +
  annotate("text", x = 0.6, y = 29, label = "b", size = 12, fontface = "bold")

###################
# GET PRIORS AND POSTERIORS
###################

s6 <- read.csv(file = file.path(here::here("01_Data/model_output_4_Hoodie"),"McmcPriorsAndPosteriors_st6 copy.csv")) %>% 
  filter(X != "10000") %>% 
  select(-X, -analysis_date) %>% 
  pivot_wider(id_cols = c(metab_date, metab_stream), names_from = parameter, values_from = PriorMean:U95per) %>% 
  mutate(metab_date = as.POSIXct(metab_date, format = "%m/%d/%y"),
         Y = as.character(strftime(metab_date, format = "%Y")))

s9 <- read.csv(file = file.path(here::here("01_Data/model_output_4_Hoodie"),"McmcPriorsAndPosteriors_st9 copy.csv")) %>% 
  select(-X, -analysis_date) %>% 
  pivot_wider(id_cols = c(metab_date, metab_stream), names_from = parameter, values_from = PriorMean:U95per) %>% 
  mutate(metab_date = as.POSIXct(metab_date, format = "%m/%d/%y"),
         Y = as.character(strftime(metab_date, format = "%Y")))

s11 <- read.csv(file = file.path(here::here("01_Data/model_output_4_Hoodie"),"McmcPriorsandPosteriors_st11 copy.csv")) %>% 
  select(-X, -analysis_date) %>% 
  pivot_wider(id_cols = c(metab_date, metab_stream), names_from = parameter, values_from = PriorMean:U95per) %>% 
  mutate(metab_date = as.POSIXct(metab_date, format = "%Y-%m-%d"),
         Y = as.character(strftime(metab_date, format = "%Y")))

s18 <- read.csv(file = file.path(here::here("01_Data/model_output_4_Hoodie"),"McmcPriorsandPosteriors_st18 copy.csv")) %>% 
  select(-X, -analysis_date) %>% 
  pivot_wider(id_cols = c(metab_date, metab_stream), names_from = parameter, values_from = PriorMean:U95per)%>% 
  mutate(metab_date = as.POSIXct(metab_date, format = "%Y-%m-%d"),
         Y = as.character(strftime(metab_date, format = "%Y")))

metBay <- rbind(s6, s9, s11, s18) %>% 
  mutate(metab_date = as.POSIXct(metab_date, format = "%m/%d/%y"),
         Y = as.character(strftime(metab_date, format = "%Y")),
         doy = as.numeric(strftime(metab_date, format = "%j")))

metSum <- readr::read_csv(file = file.path(here::here("01_Data"),"01_MundgedMetabDat.csv")) %>% 
  mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6")),
         treatment = fct_relevel(treatment, c("ambient", "nitrogen", "phosphorus")),
         treatment = fct_recode(treatment, "Ambient" = "ambient", "Nitrogen" = "nitrogen", "Phosphorus" = "phosphorus"),
         Yf = as.factor(Yf),
         streamTreat = as.factor(streamTreat))

# join bayes output with summary
# note: this only retains the data used in the analysis
metBS <- metBay %>% 
  right_join(metSum, by  = c(metab_date = "Pdt", metab_stream = "stream"))

# are there duplicates
metT <- metBS %>% 
  mutate(ID = paste0(metab_stream, metab_date))

# nope
metT[duplicated(metT$ID),]

ggplot(metBS %>% 
         filter(metab_stream == "st9"), aes(y = meanQds, x = julian, color = treatment)) +
  geom_point() +
  geom_hline(yintercept = 0.003)+
  geom_hline(yintercept = 0.006)

# what's relationship between prior and posterior?
summary(lm(median_Kr ~ PriorMean_Kr, data = metBS))


###################
# POSTERIOR V. PRIOR K'S
# FIG S8
###################

pdf(file = file.path(here::here("03_Plots"),"FigS8.pdf"), height = 5, width = 6)
ggplot(metBS %>% 
         mutate(metab_stream = fct_recode(metab_stream, "S11" = "st11U",
                                          "S6" = "st6",
                                          "S9" = "st9",
                                          "S18" = "st18"),
                metab_stream = fct_relevel(metab_stream, c("S11", "S18", "S9", "S6"))), aes(x = PriorMean_Kr, y = median_Kr, fill = metab_stream)) +
  geom_point(size = 3, shape = 21) +
  geom_abline(intercept = 0, slope = 1) +
  ylab(expression(paste("Median posterior ",italic(k[O2]), " (",h^-1,")"))) +
  xlab(expression(paste("Mean prior ",italic(k[O2]), " (",h^-1,")"))) +
  xlim(0,800) +
  ylim(0,800) +
  theme(legend.position = "top",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.width=unit(1,"cm"),
        legend.box = "horizontal",
        legend.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "transparent", color = "transparent", size = 2)) +
  annotate("text", x = 600, y = 125, label = "y = -6.7 + 1.0x", size = 6) +
  annotate("text", x = 600, y = 70, label = expression(paste(R^2," = 0.99; P < 0.001")), size = 6) 
dev.off()

###################
# change in posteriors over time
# Fig S9
###################
pdf(file = file.path(here::here("03_Plots"),"FigS9.pdf"), height = 5, width = 6)
ggplot(metBS %>% 
         mutate(metab_stream = fct_recode(metab_stream, S11 = "st11U", S18 = "st18", S6 = "st6", S9 = "st9")), aes(x = julian, y = median_Kr, fill = treatment)) +
  geom_line() +
  geom_point(size = 3, shape = 21) +
  facet_wrap(vars(metab_stream))+
  ylab(expression(paste("Median posterior ", italic(k[O2]), " (",h^-1,")"))) +
  xlab("Day of year") +
  ylim(0,800) +
  theme(legend.position = "top",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.width=unit(1,"cm"),
        legend.box = "horizontal",
        legend.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "transparent", color = "transparent", size = 2))
dev.off()

###################
# does er vary with K?
###################

# Exploritory plots
ggplot(metBS, aes(y = -median_err, x = median_Kr, color = Yf)) +
  geom_point() +
  facet_wrap(vars(metab_stream))

ggplot(metBS, aes(y = -median_err, x = meanQds*1000, color = Yf)) +
  geom_point() +
  facet_wrap(vars(metab_stream), scales = "free_x")

# good news here is that in S9 ER doesn't increase with Q after ~3 (where K plateaus)
ggplot(metBS, aes(y = -median_err, x = meanQds*1000, color = Yf)) +
  geom_point() +
  facet_wrap(vars(metab_stream), scales = "free_x")

# Plotted with other figs in 08_GasExchangePlots

###################
# POSTERIOR T_AMB V. Q
# NOT USED
###################

p1c <- ggplot(metBS%>% 
                mutate(metab_stream = fct_recode(metab_stream, "S11" = "st11U",
                                                 "S6" = "st6",
                                                 "S9" = "st9")), aes(y = median_Kr, x = meanQds*1000, fill = Yf)) +
  geom_point(shape = 21) +
  ylab(expression(paste(O[2]," exchange rate at ",T[amb]," (",d^-1,")"))) +
  xlab(expression(paste("Discharge (L ",s^-1,")"))) +
  facet_wrap(vars(metab_stream), scales = "free_x")+
  theme(legend.position = "top",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.width=unit(1,"cm"),
        legend.box = "horizontal",
        legend.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "transparent", color = "transparent", size = 2)) +
  annotate("text", x = 1, y = 650, label = "a", size = 12, fontface = "bold")



pdf(file = file.path(here::here("03_Plots"),"FigS6.pdf"), height = 8, width = 5)
ggarrange(p1a, p1b,nrow = 2, ncol = 1)
dev.off()

###################
# SAVE IMAGE
###################
save.image(file.path(here::here("04_SavedImages"), "08_GasExchangeAndER_rdat"))
