library(tidyverse)
library(gamm4)
# library(mgcViz)
# library(MuMIn)
# library(forecast)
# library(gridExtra)
library(tidymv)
library(ggpubr)
library(grid)
library(egg)
library(here)
library(scico)


load(file.path(here("04_SavedImages"), "05_ResponseFigsRdat"))

P.GPP_N.t2 <- P.GPP_N.t %>% 
  pivot_wider(id_cols = c(stream, met, invKT.C.StMean, Tanom.iktCs), names_from = treatment, values_from = c(fit,se.fit))
  
P.GPP_N.t %>% 
  filter(duplicated(stream, met, invKT.C.StMean, Tanom.iktCs))

EffectPlotNtrt <- rbind(P.GPP_N.t %>% 
                           select(treatment, invKT.C.StMean, Tanom.iktCs, stream, fit, se.fit) %>% 
                           mutate(metab = "GPP"), 
                         P.ER_N.b %>% 
                           select(treatment, invKT.C.StMean, Tanom.iktCs, stream, fit, se.fit) %>% 
                           mutate(metab = "ER")) %>% 
          pivot_wider(id_cols = c(stream, metab, invKT.C.StMean, Tanom.iktCs),names_from = treatment, values_from = c(fit,se.fit)) %>% 
          mutate(EffectSize = ifelse(metab == "NEP", fit_Nitrogen - fit_Ambient,
                                     exp(fit_Nitrogen)/exp(fit_Ambient))) %>% 
          select(stream, metab, invKT.C.StMean, Tanom.iktCs, "fit_amb" = "fit_Ambient", "fit_nut" = "fit_Nitrogen",
                 'se.fit_amb' = "se.fit_Ambient", "se.fit_nut" = "se.fit_Nitrogen",EffectSize) %>% 
          mutate(trt = "N")

# Effect sizes for Wyatt
summary(EffectPlotNtrt[EffectPlotNtrt$trt == "N" & EffectPlotNtrt$metab == "GPP",]$EffectSize) # 1.7
summary(EffectPlotNtrt[EffectPlotNtrt$trt == "N" & EffectPlotNtrt$metab == "ER",]$EffectSize) # 5.5

EffectPlotNtrtNEP <- rbind(P.NEP_N.b %>% 
                          select(treatment, invKT.C.StMean, Tanom.iktCs, stream, fit, se.fit) %>% 
                          mutate(metab = "NEP")) %>% 
  pivot_wider(id_cols = c(stream, metab, invKT.C.StMean, Tanom.iktCs),names_from = treatment, values_from = c(fit,se.fit)) %>% 
   mutate(EffectSize = ifelse(metab == "NEP", fit_Nitrogen - fit_Ambient,
                             exp(fit_Nitrogen)/exp(fit_Ambient))) %>%
  select(stream, metab, invKT.C.StMean, Tanom.iktCs, "fit_amb" = "fit_Ambient", "fit_nut" = "fit_Nitrogen",
         'se.fit_amb' = "se.fit_Ambient", "se.fit_nut" = "se.fit_Nitrogen",EffectSize) %>% 
  mutate(trt = "N")

EffectPlotPtrt <- rbind(P.GPP_P.b %>% 
                          select(treatment, invKT.C.StMean,  Tanom.iktCs, stream, fit, se.fit) %>% 
                          mutate(metab = "GPP"),
                        P.ER_P.b %>% 
                          select(treatment, invKT.C.StMean,  Tanom.iktCs, stream, fit, se.fit) %>% 
                          mutate(metab = "ER")) %>% 
  pivot_wider(id_cols = c(stream, metab, invKT.C.StMean, Tanom.iktCs),names_from = treatment, values_from = c(fit,se.fit)) %>% 
            mutate(EffectSize = ifelse(metab == "NEP", fit_Phosphorus - fit_Ambient,
                                       exp(fit_Phosphorus)/exp(fit_Ambient))) %>% 
            select(stream, metab, invKT.C.StMean, Tanom.iktCs, "fit_amb" = "fit_Ambient", "fit_nut" = "fit_Phosphorus",
                   'se.fit_amb' = "se.fit_Ambient", "se.fit_nut" = "se.fit_Phosphorus",EffectSize)%>% 
        mutate(trt = "P")

EffectPlotPtrtNEP <- rbind(P.NEP_P.b %>% 
                          select(treatment,  invKT.C.StMean, Tanom.iktCs, stream, fit, se.fit) %>% 
                          mutate(metab = "NEP")) %>% 
  pivot_wider(id_cols = c(stream, metab, invKT.C.StMean, Tanom.iktCs),names_from = treatment, values_from = c(fit,se.fit)) %>% 
  mutate(EffectSize = ifelse(metab == "NEP", fit_Phosphorus - fit_Ambient,
                             exp(fit_Phosphorus)/exp(fit_Ambient))) %>% 
  select(stream, metab, invKT.C.StMean, Tanom.iktCs, "fit_amb" = "fit_Ambient", "fit_nut" = "fit_Phosphorus",
         'se.fit_amb' = "se.fit_Ambient", "se.fit_nut" = "se.fit_Phosphorus",EffectSize)%>% 
  mutate(trt = "P")

EffectPlot.df <- rbind(EffectPlotNtrt, EffectPlotPtrt) 
EffectPlotNEP.df <- rbind(EffectPlotNtrtNEP, EffectPlotPtrtNEP)



EFSPa <- ggplot(EffectPlot.df %>% 
                  filter(trt == "N" & metab == "GPP") %>% 
                  mutate(metab = fct_recode(metab, "GPP (N treatment)" = "GPP")), aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = EffectSize)) +
  geom_raster() +
  # scale_fill_distiller(palette = "Greens", limits = c(0,10), direction = 1) +
  scale_fill_scico(palette = "batlow", limits = c(0,8), direction = 1) +
  ylab(expression(atop("Mean summer temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  xlab(expression(atop("Mean daily temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  labs(fill = "Effect size (ES)") +
  facet_grid(. ~ metab) +
  theme(legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

EFSPb <- ggplot(EffectPlot.df %>% 
                  filter(trt == "N" & metab == "ER")%>% 
                  mutate(metab = fct_recode(metab, "ER (N treatment)" = "ER")), aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = EffectSize)) +
  geom_raster() +
  # scale_fill_distiller(palette = "Reds", limits = c(0,10), direction = 1) +
  scale_fill_scico(palette = "batlow", limits = c(0,8), direction = 1) +
  ylab(expression(atop("Mean summer temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  xlab(expression(atop("Mean daily temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  facet_grid(. ~ metab) +
  labs(fill = "Effect size (ES)")+
  theme(legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

EFSPc <- ggplot(EffectPlotNEP.df %>% 
                  filter(trt == "N" & metab == "NEP")%>% 
                  mutate(metab = fct_recode(metab, "NEP (N treatment)" = "NEP")), aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = EffectSize)) +
  geom_raster() +
  # scale_fill_distiller(palette = "Purples", limits = c(-25,0)) +
  scale_fill_distiller(palette = "Oranges", limits = c(-8.5,0)) +
  ylab(expression(atop("Mean summer temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  xlab(expression(atop("Mean daily temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  facet_grid(. ~ metab) +
  labs(fill = "Effect size (ES)")+
  theme(legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

EFSPd <- ggplot(EffectPlot.df %>% 
                  filter(trt == "P" & metab == "GPP") %>% 
                  mutate(metab = fct_recode(metab, "GPP (P treatment)" = "GPP")), aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = EffectSize)) +
  geom_raster() +
  # scale_fill_distiller(palette = "Greens", limits = c(0,10), direction = 1) +
  scale_fill_scico(palette = "batlow", limits = c(0,8), direction = 1) +
  ylab(expression(atop("Mean summer temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  xlab(expression(atop("Mean daily temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  labs(fill = "Effect size (ES)") +
  facet_grid(. ~ metab) +
  theme(legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
  
  

EFSPe <- ggplot(EffectPlot.df %>% 
                  filter(trt == "P" & metab == "ER")%>% 
                  mutate(metab = fct_recode(metab, "ER (P treatment)" = "ER")), aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = EffectSize)) +
  geom_raster() +
  # scale_fill_distiller(palette = "Reds", limits = c(0,10), direction = 1) +
  scale_fill_scico(palette = "batlow", limits = c(0,8), direction = 1) +
  ylab(expression(atop("Mean summer temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  xlab(expression(atop("Mean daily temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  facet_grid(. ~ metab) +
  labs(fill = "Effect size (ES)")+
  theme(legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))



EFSPf <- ggplot(EffectPlotNEP.df %>% 
                 filter(trt == "P" & metab == "NEP")%>% 
                 mutate(metab = fct_recode(metab, "NEP (P treatment)" = "NEP")), aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = EffectSize)) +
  geom_raster() +
  # scale_fill_distiller(palette = "Purples", limits = c(-25,0)) +
  scale_fill_distiller(palette = "Oranges", limits = c(-8.5,0)) +
  ylab(expression(atop("Mean summer temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  xlab(expression(atop("Mean daily temperature", paste(frac(1,paste("k",T[amb]))," - ",frac(1,paste("k",T[mean])))))) +
  facet_grid(. ~ metab) +
  labs(fill = "Effect size (ES)")+
  theme(legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))


Fig1a.g <- ggplotGrob(EFSPa)
Fig1b.g <- ggplotGrob(EFSPb)
Fig1c.g <- ggplotGrob(EFSPc)
Fig1d.g <- ggplotGrob(EFSPd)
Fig1e.g <- ggplotGrob(EFSPe)
Fig1f.g <- ggplotGrob(EFSPf)

Fig1a.gtf <- gtable_frame(Fig1a.g,  width = unit(1, "null"), height = unit(1, "null"))
Fig1b.gtf <- gtable_frame(Fig1b.g,  width = unit(1, "null"), height = unit(1, "null"))
Fig1c.gtf <- gtable_frame(Fig1c.g,  width = unit(1, "null"), height = unit(1, "null"))
Fig1d.gtf <- gtable_frame(Fig1d.g,  width = unit(1, "null"), height = unit(1, "null"))
Fig1e.gtf <- gtable_frame(Fig1e.g,  width = unit(1, "null"), height = unit(1, "null"))
Fig1f.gtf <- gtable_frame(Fig1f.g,  width = unit(1, "null"), height = unit(1, "null"))


Fig1ad.gtf <- gtable_frame(gtable_rbind(Fig1a.gtf, Fig1d.gtf),  width = unit(1, "null"), height = unit(2, "null"))
Fig1be.gtf <- gtable_frame(gtable_rbind(Fig1b.gtf, Fig1e.gtf),  width = unit(1, "null"), height = unit(2, "null"))
Fig1cf.gtf <- gtable_frame(gtable_rbind(Fig1c.gtf, Fig1f.gtf),  width = unit(1, "null"), height = unit(2, "null"))


Fig1.gtf <- gtable_frame(gtable_cbind(Fig1ad.gtf, Fig1be.gtf, Fig1cf.gtf),  width = unit(3, "null"), height = unit(2, "null"))

# need to let grid.draw complete before moving on
png(file.path(here("03_Plots"), "EffectsPlot3.png"), units = "in", height = 12, width = 20, res = 300)
grid.newpage()
grid.draw(Fig1.gtf)
grid.text("(a)", x = unit(0.01,"npc"), y = unit(0.94,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("(b)", x = unit(0.36,"npc"), y = unit(0.94,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("(c)", x = unit(0.69,"npc"), y = unit(0.94,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("(d)", x = unit(0.01,"npc"), y = unit(0.425,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("(e)", x = unit(0.36,"npc"), y = unit(0.425,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("(f)", x = unit(0.69,"npc"), y = unit(0.425,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text(expression(paste(ES[GPP]," = ", frac("treatment","ambient"))), 
          x = unit(0.15,"npc"), y = unit(0.8,"npc"), gp=gpar(fontsize = 19, fontface = 2, col = "white"))
grid.text(expression(paste(ES[ER]," = ", frac("treatment","ambient"))), 
          x = unit(0.48,"npc"), y = unit(0.8,"npc"), gp=gpar(fontsize = 19, fontface = 2, col = "white"))
grid.text(expression(paste(ES[NEP]," = ", "treatment - ambient")), 
          x = unit(0.85,"npc"), y = unit(0.8,"npc"), gp=gpar(fontsize = 19, fontface = 2, col = "white"))
dev.off()


# Look at temp interactions
S11temp <- unique(metAll[,c("stream","invKT.C.StMean")])[[3,2]]
S18temp <- unique(metAll[,c("stream","invKT.C.StMean")])[[4,2]]
S9temp <- -0.364392602 #unique(metAll[,c("stream","invKT.C.StMean")])[[2,2]]
S6temp <- unique(metAll[,c("stream","invKT.C.StMean")])[[1,2]]

plot_GPPamb <- ggplot(P.GPP_AMB.t %>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3))) %>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), aes(y = fit, x = -Tanom.iktCs, color = stream)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("GPP Ambient")

plot_ERamb <- ggplot(P.ER_AMB.b %>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3)))%>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), aes(y = fit, x = -Tanom.iktCs, color = stream)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("ER Ambient")


plot_NEPamb <- ggplot(P.NEP_AMB.t %>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3)))%>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), aes(y = fit, x = -Tanom.iktCs, color = stream)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("NEP ambient")

# Phosphorus
plot_GPP_P <- ggplot(P.GPP_P.b %>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3)))%>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), 
       aes(y = fit, x = -Tanom.iktCs, color = stream, linetype = treatment)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("GPP Phosphorus")


plot_ER_P <- ggplot(P.ER_P.b %>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3)))%>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), 
       aes(y = fit, x = -Tanom.iktCs, color = stream, linetype = treatment)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("ER Phosphorus")

plot_NEP_P <- ggplot(P.NEP_P.b %>% 
         mutate(stream2 = ifelse(round(invKT.C.StMean,3) == round(S11temp,3), "S11",
                            ifelse(round(invKT.C.StMean,3) == round(S18temp,3), "S18",
                             ifelse(round(invKT.C.StMean,3) == round(S9temp,3), "S9",
                              ifelse(round(invKT.C.StMean,3) == round(S6temp,3), "S6","blah")))),
                stream2 = as.factor(stream2)) %>% 
         filter(stream2 != "blah") %>% 
         mutate(stream2 = fct_relevel(stream2, c("S11", "S18", "S9", "S6"))), 
       aes(y = fit, x = -Tanom.iktCs, color = stream2, linetype = treatment)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("NEP Phosphorus")

# N
plot_GPP_N <- ggplot(P.GPP_N.t %>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3)))%>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), 
       aes(y = fit, x = -Tanom.iktCs, color = stream, linetype = treatment)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("GPP Nitrogen")

plot_ER_N <- ggplot(P.ER_N.b %>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3)))%>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), 
       aes(y = fit, x = -Tanom.iktCs, color = stream, linetype = treatment)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(2.5,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2))+
  ggtitle("ER Nitrogen")

plot_NEP_N <- ggplot(P.NEP_N.b %>% 
         mutate(stream2 = ifelse(round(invKT.C.StMean,3) == round(S11temp,3), "S11",
                                 ifelse(round(invKT.C.StMean,3) == round(S18temp,3), "S18",
                                        ifelse(round(invKT.C.StMean,3) == round(S9temp,3), "S9",
                                               ifelse(round(invKT.C.StMean,3) == round(S6temp,3), "S6","blah")))),
                stream2 = as.factor(stream2)) %>% 
         filter(stream2 != "blah") %>% 
         mutate(stream2 = fct_relevel(stream2, c("S11", "S18", "S9", "S6"))), 
       aes(y = fit, x = -Tanom.iktCs, color = stream2, linetype = treatment)) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "right",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.width=unit(1,"cm"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 2),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 26, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 2))+
  ggtitle("NEP Nitrogen")

png(file.path(here("03_Plots"), "InteractionsPlot.png"), units = "in", height = 15, width = 25, res = 300)
ggarrange(plot_GPPamb, plot_ERamb, plot_NEPamb,
          plot_GPP_N, plot_ER_N, plot_NEP_N,
          plot_GPP_P, plot_ER_P, plot_NEP_P, nrow = 3, ncol = 3)
dev.off()

# Effect size line plots
plotES_GPP_N <- ggplot(EffectPlot.df %>% 
         filter(trt == "N" & metab == "GPP") %>% 
         mutate(metab = fct_recode(metab, "GPP (N treatment)" = "GPP"))%>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3))) %>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), aes(y = EffectSize, x = -Tanom.iktCs, color = stream)) +
  geom_line(size = 1.25) +
  ylim(0,8) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.width=unit(2.5,"cm"),
    panel.background = element_rect(fill = "transparent", color = "black", size = 2),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_line(color = "transparent"),
    strip.text = element_text(size = 26, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("GPP Nitrogen")

plotES_ER_N <- ggplot(EffectPlot.df %>% 
         filter(trt == "N" & metab == "ER") %>% 
         mutate(metab = fct_recode(metab, "ER (N treatment)" = "ER"))%>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3))) %>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), aes(y = EffectSize, x = -Tanom.iktCs, color = stream)) +
  ylim(0,8) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.width=unit(2.5,"cm"),
    panel.background = element_rect(fill = "transparent", color = "black", size = 2),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_line(color = "transparent"),
    strip.text = element_text(size = 26, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("ER Nitrogen")

plotES_NEP_N <- ggplot(EffectPlotNEP.df %>% 
         filter(trt == "N" & metab == "NEP") %>% 
         mutate(metab = fct_recode(metab, "NEP (N treatment)" = "NEP"))%>% 
         mutate(stream2 = ifelse(round(invKT.C.StMean,3) == round(S11temp,3), "S11",
                                 ifelse(round(invKT.C.StMean,3) == round(S18temp,3), "S18",
                                        ifelse(round(invKT.C.StMean,3) == round(S9temp,3), "S9",
                                               ifelse(round(invKT.C.StMean,3) == round(S6temp,3), "S6","blah")))),
                stream2 = as.factor(stream2)) %>% 
         filter(stream2 != "blah") %>% 
         mutate(stream2 = fct_relevel(stream2, c("S11", "S18", "S9", "S6"))), aes(y = EffectSize, x = -Tanom.iktCs, color = stream2)) +
  ylim(-8,0) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.width=unit(2.5,"cm"),
    panel.background = element_rect(fill = "transparent", color = "black", size = 2),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_line(color = "transparent"),
    strip.text = element_text(size = 26, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("NEP Nitrogen")


plotES_GPP_P <- ggplot(EffectPlot.df %>% 
         filter(trt == "P" & metab == "GPP") %>% 
         mutate(metab = fct_recode(metab, "GPP (P treatment)" = "GPP"))%>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3))) %>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), aes(y = EffectSize, x = -Tanom.iktCs, color = stream)) +
  ylim(0,8) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.width=unit(2.5,"cm"),
    panel.background = element_rect(fill = "transparent", color = "black", size = 2),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_line(color = "transparent"),
    strip.text = element_text(size = 26, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("GPP Phosphorus")

plotES_ER_P <- ggplot(EffectPlot.df %>% 
         filter(trt == "P" & metab == "ER") %>% 
         mutate(metab = fct_recode(metab, "ER (P treatment)" = "ER"))%>% 
         filter((stream == "st11U" & round(invKT.C.StMean,3) == round(S11temp,3)) |
                  (stream == "st18" & round(invKT.C.StMean,3) == round(S18temp,3))|
                  (stream == "st9" & round(invKT.C.StMean,3) == round(S9temp,3))|
                  (stream == "st6" & round(invKT.C.StMean,3) == round(S6temp,3))) %>% 
         mutate(stream = fct_relevel(stream, c("st11U", "st18", "st9", "st6"))), aes(y = EffectSize, x = -Tanom.iktCs, color = stream)) +
  ylim(0,8) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.width=unit(2.5,"cm"),
    panel.background = element_rect(fill = "transparent", color = "black", size = 2),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_line(color = "transparent"),
    strip.text = element_text(size = 26, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("ER Phosphorus")


plotES_NEP_P <- ggplot(EffectPlotNEP.df %>% 
         filter(trt == "P" & metab == "NEP") %>% 
         mutate(metab = fct_recode(metab, "NEP (P treatment)" = "NEP"))%>% 
         mutate(stream2 = ifelse(round(invKT.C.StMean,3) == round(S11temp,3), "S11",
                                 ifelse(round(invKT.C.StMean,3) == round(S18temp,3), "S18",
                                        ifelse(round(invKT.C.StMean,3) == round(S9temp,3), "S9",
                                               ifelse(round(invKT.C.StMean,3) == round(S6temp,3), "S6","blah")))),
                stream2 = as.factor(stream2)) %>% 
         filter(stream2 != "blah") %>% 
         mutate(stream2 = fct_relevel(stream2, c("S11", "S18", "S9", "S6"))), aes(y = EffectSize, x = -Tanom.iktCs, color = stream2)) +
  ylim(-8,0) +
  geom_line(size = 1.25) +
  scale_color_manual(values = c("blue", "lightblue", "pink", "red")) +
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red")) +
  theme(#legend.position = "top",
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.width=unit(2.5,"cm"),
    panel.background = element_rect(fill = "transparent", color = "black", size = 2),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_line(color = "transparent"),
    strip.text = element_text(size = 26, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black", size = 2)) +
  ggtitle("NEP Phosphorus")

png(file.path(here("03_Plots"), "EffectSizeLinePlot.png"), units = "in", height = 10, width = 15, res = 300)
ggarrange(plotES_GPP_N, plotES_ER_N, plotES_NEP_N,
          plotES_GPP_P, plotES_ER_P, plotES_NEP_P, nrow = 2, ncol = 3)
dev.off()
