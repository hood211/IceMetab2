library(tidyverse)
library(here)
library(MuMIn)
library(gridExtra)
library(GGally)

met <- readr::read_csv(file = file.path(here::here("01_Data"),"01_MundgedMetabDat.csv")) %>% 
  mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6")),
         treatment = fct_relevel(treatment, c("ambient", "nitrogen", "phosphorus")),
         Yf = as.factor(Yf),
         streamTreat = as.factor(streamTreat)) 



ggpairs(met %>% 
          filter(stream == "st11U") %>% 
          select(LightPerDay, meanTemp, meanQds))

ggpairs(met %>% 
          filter(stream == "st18") %>% 
          select(LightPerDay, meanTemp, meanQds))

ggpairs(met %>% 
          filter(stream == "st9") %>% 
          select(LightPerDay, meanTemp, meanQds))

ggpairs(met %>% 
          filter(stream == "st6") %>% 
          select(LightPerDay, meanTemp, meanQds))




# mean Temp - all different
ggplot(met, aes(y = meanTemp, x = stream, color = treatment)) +
  geom_boxplot()

kruskal.test(meanTemp ~ treatment, data = met %>% 
                                                filter(stream == "st11U"))
kruskal.test(meanTemp ~ treatment, data = met %>% 
               filter(stream == "st18"))
kruskal.test(meanTemp ~ treatment, data = met %>% 
               filter(stream == "st9"))
kruskal.test(meanTemp ~ treatment, data = met %>% 
               filter(stream == "st6"))

# mean Q - all different
ggplot(met, aes(y = meanQds, x = stream, color = treatment)) +
  geom_boxplot()

kruskal.test(meanQds ~ treatment, data = met %>% 
               filter(stream == "st11U"))
kruskal.test(meanQds ~ treatment, data = met %>% 
               filter(stream == "st18"))
kruskal.test(meanQds ~ treatment, data = met %>% 
               filter(stream == "st9"))
kruskal.test(meanQds ~ treatment, data = met %>% 
               filter(stream == "st6"))

# LightPerDay - all different
ggplot(met, aes(y = LightPerDay, x = stream, color = treatment)) +
  geom_boxplot()
kruskal.test(LightPerDay ~ treatment, data = met %>% 
               filter(stream == "st11U"))
kruskal.test(LightPerDay ~ treatment, data = met %>% 
               filter(stream == "st18"))
kruskal.test(LightPerDay ~ treatment, data = met %>% 
               filter(stream == "st9"))
kruskal.test(LightPerDay ~ treatment, data = met %>% 
               filter(stream == "st6"))

ggplot(met, aes(y = lGPP, x = lER, color = stream)) +
  geom_point()


met %>% 
  group_by(stream, treatment) %>% 
  summarise(across(c(LightPerDay, meanTemp, meanQds), median)) %>% 
  mutate(meanQds = meanQds * 1000)
met %>% 
  group_by(stream) %>% 
  summarise(across(c(LightPerDay, meanTemp, meanQds), median)) %>% 
  mutate(meanQds = meanQds * 1000)

S11temp <- 10.3-8.89
S18temp <- 11.2-9.88
S9temp <- 13-11.9
S6temp <- 16-15.2

S11q <- (8.39-3.47)
S18q <- (25.3-14.6)
S9q <- (4.85-2.64)
S6q <- (14.7-10.3)
median(S11q, S18q, S9q, S6q)

S11q <- (8.39-3.47)/4.07
S18q <- (25.3-14.6)/20.3
S9q <- (4.85-2.64)/2.96
S6q <- (14.7-10.3)/13.3

# does temp and Q differ among years?
DayTemp.aov1 <- aov(meanTemp ~ stream * treatment, data = met)
DayTemp.aov2 <- aov(meanTemp ~ stream + treatment, data = met)
DayTemp.aov3 <- aov(meanTemp ~ stream, data = met)
DayTemp.aov4 <- aov(meanTemp ~ treatment, data = met)

model.sel(DayTemp.aov1, DayTemp.aov2, DayTemp.aov3, DayTemp.aov4)
summary(DayTemp.aov1)
plot(DayTemp.aov1)
ggplot()