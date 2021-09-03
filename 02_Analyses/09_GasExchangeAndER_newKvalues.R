# Comparison of effect of different K priors on ER
# JMH

#################
#  libraries
#################
library(tidyverse)

#################
# import data
#################

met <- read.csv(file = file.path(here::here("01_Data"),"ERestimatesWdiffKvQrelationships4Sup.csv"), header = T) %>% 
  mutate(metab_date = as.POSIXct(metab_date, format = "%m/%d/%y"),
         treatment = ifelse(metab_date < "2017-01-01", "Phosphorus", "Nitrogen"),
         method2 = ifelse(method == "Hood", "SF6", "SF6propane"),
         ER_C = (ER_post_median * (0.85) * (12/32)),
         treatment = fct_relevel(treatment, c("Phosphorus", "Nitrogen")))

str(met)

# Data from actual met models
s9Bay <- read.csv(file = file.path(here::here("01_Data/model_output_4_Hoodie"),"McmcPriorsAndPosteriors_st9 copy.csv")) %>% 
  select(-X, -analysis_date) %>% 
  pivot_wider(id_cols = c(metab_date, metab_stream), names_from = parameter, values_from = PriorMean:U95per) %>% 
  mutate(metab_date = as.POSIXct(metab_date, format = "%m/%d/%y"),
         Y = as.character(strftime(metab_date, format = "%Y")),
         ER_C = median_err * (0.85) * (12/32))

# combine with summary to get Q
s9Sum <- readr::read_csv(file = file.path(here::here("01_Data"),"01_MundgedMetabDat.csv")) %>% 
  mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6")),
         treatment = fct_relevel(treatment, c("ambient", "nitrogen", "phosphorus")),
         treatment = fct_recode(treatment, "Ambient" = "ambient", "Nitrogen" = "nitrogen", "Phosphorus" = "phosphorus"),
         Yf = as.factor(Yf),
         streamTreat = as.factor(streamTreat)) %>% 
  filter(stream == "st9")

# join bayes output with summary
# note: this only retains the data used in the analysis
metBS <- s9Bay %>% 
  right_join(s9Sum, by  = c(metab_date = "Pdt", metab_stream = "stream")) %>% 
  mutate(treatment = fct_relevel(treatment, c("Ambient", "Phosphorus", "Nitrogen"))) %>% 
  filter(julian > 183)

ggplot(metBS, aes(y = PriorMean_Kr, x = meanQds*1000)) +
  geom_point() +
  facet_wrap(vars(metab_stream)) +
  xlim(0,3)

summary(metBS[metBS$meanQds > 0.004, ]$PriorMean_Kr)


metW <- met %>% 
  select(-method) %>% 
  pivot_wider(id_cols = c(metab_date, metab_stream, treatment), names_from = method2, values_from = c(Q_m3.s:K_U95per,ER_C)) %>% 
  mutate(metab_date = as.Date(metab_date),
         Y = strftime(metab_date, format = "%Y"),
         treatment = fct_relevel(treatment, c("Phosphorus", "Nitrogen")))

# how much greater was new model at Q > 3
metW %>% 
  filter(Q_m3.s_SF6 > 0.0031) %>% 
  # summarize(across(c(ER_post_median_SF6, ER_post_median_SF6propane), mean)) %>% 
  mutate(ER_post_median_SF6propane_C = ER_post_median_SF6propane* (0.85) * (12/32),
         ER_post_median_SF6_C = ER_post_median_SF6* (0.85) * (12/32),
    HowMuchGreater = ER_post_median_SF6propane_C/ER_post_median_SF6_C -1,
    HowMuchLess = ER_post_median_SF6_C/ER_post_median_SF6propane_C -1) %>% 
  select(metab_date, treatment, Q_m3.s_SF6, ER_post_median_SF6propane_C, ER_post_median_SF6_C, HowMuchGreater, HowMuchLess)

# Sup plot
pdf(file.path(here::here("03_Plots"),"FigS7.pdf"), height = 6, width = 10)
ggplot() + 
  geom_point(data = metBS, aes(y = ER_C, x = meanQds*1000, fill  = treatment), shape = 21, size = 3, alpha = 33/100) +
  geom_segment(data = poopW, aes(x = Q_m3.s_SF6*1000, xend = Q_m3.s_SF6propane*1000, 
                                 y = ER_post_median_SF6 * (0.85) * (12/32), yend = ER_post_median_SF6propane * (0.85) * (12/32))) +
  geom_point(data = poop, aes(x = Q_m3.s*1000, y = ER_C, shape = method2, fill = treatment), size = 5) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = c("grey60", "steelblue1", "gold1")) +
  xlab(expression(paste("Discharge (L ",s^-1,")"))) +
  ylab(expression(paste("Ecosystem respiration (g C ", m^-2," ", d^-1,")"))) +
  ylim(-15,0) +
  xlim(0,10) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    axis.line = element_line(color = "black", size = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.key.width=unit(1,"cm"),
    # legend.box = "horizontal",
    legend.background = element_rect(color = "transparent"),
    panel.background = element_rect(fill = "transparent", size = 2),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_line(color = "transparent"),
    strip.text = element_text(size = 22, face = "bold"),
    strip.background = element_rect(fill = "transparent", color = "transparent", size = 2)) +
  annotate("text", x = 6, y = -2.5, label = "Ambient", size = 7, fontface = "bold", hjust = 1) +
  annotate("text", x = 6, y = -3.5, label = "Phosphorus", size = 7, fontface = "bold", hjust = 1) +
  annotate("text", x = 6, y = -4.5, label = "Nitrogen", size = 7, fontface = "bold", hjust = 1) +
  annotate("text", x = 6.7, y = -0.5, label = "All ER", size = 7, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 6.7, y = -1.3, label = "data (SF6)", size = 7, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 8.0, y = -1, label = "SF6", size = 7, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 9.25, y = -0.5, label = "SF6 +", size = 7, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 9.25, y = -1.3, label = "Propane", size = 7, fontface = "bold", hjust = 0.5) +
  annotate("point", x = 6.7, y = -2.5, size = 3, fill = "grey60", shape = 21, alpha = 33/100) +
  annotate("point", x = 6.7, y = -3.5, size = 3, fill = "steelblue1", shape = 21, alpha = 33/100) +
  annotate("point", x = 6.7, y = -4.5, size = 3, fill = "gold1", shape = 21, alpha = 33/100)  +
  # annotate("point", x = 8.5, y = -2, size = 5, fill = "grey60", shape = 21) +
  annotate("point", x = 8, y = -3.5, size = 5, fill = "steelblue1", shape = 21) +
  annotate("point", x = 8, y = -4.5, size = 5, fill = "gold1", shape = 21)  +
  # annotate("point", x = 9.75, y = -2, size = 5, fill = "grey60", shape = 22) +
  annotate("point", x = 9.25, y = -3.5, size = 5, fill = "steelblue1", shape = 22) +
  annotate("point", x = 9.25, y = -4.5, size = 5, fill = "gold1", shape = 22) +
  geom_rect(aes(xmin = 3.9, xmax = 9.99, ymin = -5.1, ymax = -0.001), fill = "transparent", color = "black") +
  geom_vline(xintercept = 3.1)
dev.off()

###################
# SAVE IMAGE
###################
save.image(file.path(here::here("04_SavedImages"), "09_GasExchangeAndER_newKvalue_rdat"))
