library(tidyverse)
library(ggpubr)

ge2 <-  readr::read_csv(file = file.path(here::here("01_Data"),"GasExchangeEstimatesAndQttEct.csv")) 
ge <- ge2 %>% 
        filter(stream %in% c("st6", "st9", "st11U")) %>% 
        mutate(K_O2_T20.log = log(K_O2_T20),
               Q_DS.log = log(Q_DS))



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


ggplot(ge2, aes(y =K_O2_T20, x = Vms, label = stream, color = stream)) +
  geom_text()+
  scale_x_log10() +
  scale_y_log10()


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

# see zz_concernsAboutErK for p1c

pdf(file = file.path(here::here("03_Plots"),"FigS6.pdf"), height = 8, width = 5)
ggarrange(p1a, p1b,nrow = 2, ncol = 1)
dev.off()



