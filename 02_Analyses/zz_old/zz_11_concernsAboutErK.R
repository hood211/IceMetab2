library(tidyverse)
library(here)

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

# change in posteriors over time
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

# does er vary with K?
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
# did Wyatt convert K to stream temp?
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


ggplot(metBS, aes(y = median_Kr, x = doy.x, color = Yf)) +
  geom_point() +
  facet_wrap(vars(metab_stream))

ggplot(metBay, aes(y = median_err, x = doy, color = Y)) +
  geom_point() +
  facet_wrap(vars(metab_stream)) +
  geom_vline(xintercept = 183)

ggplot(met, aes(y = -median_err, x = median_BS_GPP, color = Yf)) +
  geom_point() +
  facet_wrap(vars(metab_stream))  +
  geom_abline(intercept = 0, slope = 1)

