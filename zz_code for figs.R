#########################################
# ER_all
#########################################



# summary Plots
Focusing on this
```{r}
summary(g.er.BEST1$gam)
```

## Predictions
```{r}
pdf(file.path(here::here("03_Plots"), "02a_ERpredictionXdoy.pdf"), height = 4, width = 8)
ggplot() +
  geom_point(data = met, aes(y = lER, x = doy, color = treatment)) +
  geom_line(data = met, aes(y = ERgamPred, x = doy, color = treatment)) +
  geom_ribbon(data = met, aes(ymin = ERgamPred - ERgamPred.se,
                              ymax = ERgamPred + ERgamPred.se, x = doy, fill = treatment), alpha = 0.5) +
  facet_grid(stream ~ treatment)
dev.off()
```



Use model to predict log ER, with Qres, doy and stream set to "centered" values
```{r}
gamPred2 <- predict_gam(g.er.BEST1$gam, exclude_terms = c(
  "s(Qres):streamst11U",
  "s(Qres):streamst18",
  "s(Qres):streamst6",
  "s(Qres):streamst9",
  "s(doy):streamst11U",
  "s(doy):streamst18",
  "s(doy):streamst6",
  "s(doy):streamst9"),
  values  = list(
    Qres = 0,
    doy = 200,
    stream = "st6")) 

# make the temperature axes sensical
gamPred3 <- gamPred2 %>% 
  mutate(Tanom.iktCsFlip = -1 * Tanom.iktCs,
         invKT.C.StMeanFlip = -1 *invKT.C.StMean)

```

**NOTE: TEMPERATURE AXIS ARE FLIPPED SO LARGE #'S ARE WARM**
During amb year, ER increases with daily temp and with annual temp across streams. Duing the N and P years, this switches. We see ER increasing with annual temp again, but decreasing with daily temp. Cool

```{r}
ErPredSur <- ggplot() +
  geom_raster(data = gamPred2, aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = fit)) +
  facet_wrap(vars(treatment)) +
  scale_fill_distiller(palette = "Spectral") +
  geom_point(data = met, aes(y = -invKT.C.StMean, x = -Tanom.iktCs, size = lER), 
             position = position_jitter(width = 0.05, height = 0.025), shape = 21, fill = "transparent", alpha= 0.5) +
  scale_size_binned(n.breaks = 4, nice.breaks = TRUE) +
  ylab("Mean summer temp (-1 * centered 1/kT)") +
  xlab("Mean daily temp (-1 * centered 1/kT)")
```

```{r}
gamPred2w <- gamPred2 %>% 
  pivot_wider(id_cols = Tanom.iktCs:doy, names_from = treatment, values_from = c(fit,se.fit)) %>% 
  mutate(N_amb = exp(fit_nitrogen)/exp(fit_ambient),
         P_amb = exp(fit_phosphorus)/exp(fit_ambient))

gamPred2w2 <- gamPred2w %>% 
  select(Tanom.iktCs:doy, fit_ambient:P_amb) %>% 
  pivot_longer(cols = N_amb:P_amb, names_to = "treatment", values_to = "effectSize") %>% 
  mutate(treatment = fct_recode(treatment, Nitrogen = "N_amb", Phosphorus = "P_amb"))
```


```{r}
pdf(file.path(here::here("03_Plots"), "02b_ERpredictionSurfaces.pdf"), height = 4, width = 8)
ErPredSur
dev.off()
```

```{r}
pdf(file.path(here::here("03_Plots"), "02c_ERdifferencePlots2.pdf"), height = 5, width = 4)
ggplot(gamPred2w2, aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = effectSize)) +
  geom_raster() +
  scale_fill_distiller(palette = "Spectral") +
  ylab("Mean summer temp (-1 * centered 1/kT)") +
  xlab("Mean daily temp (-1 * centered 1/kT)") +
  facet_wrap(vars(treatment), nrow = 2, ncol = 1)
dev.off()
```


#########################################
# ER_N
#########################################


# summary Plots
Focusing on this
```{r}
summary(g.er.BEST1$gam)
```

## Predictions
```{r}
pdf(file.path(here::here("03_Plots"), "02b_ERpredXdoy_Nonly.pdf"), height = 4, width = 8)
ggplot() +
  geom_point(data = met, aes(y = lER, x = doy, color = treatment)) +
  geom_line(data = met, aes(y = ERgamPredNonly, x = doy, color = treatment)) +
  geom_ribbon(data = met, aes(ymin = ERgamPredNonly - ERgamPredNonly.se,
                              ymax = ERgamPredNonly + ERgamPredNonly.se, x = doy, fill = treatment), alpha = 0.5) +
  facet_grid(stream ~ treatment)
dev.off()
```



Use model to predict log ER, with Qres, doy and stream set to "centered" values
```{r}
gamPred2 <- predict_gam(g.er.BEST1$gam, exclude_terms = c(
  "s(Qres):streamst11U",
  "s(Qres):streamst18",
  "s(Qres):streamst6",
  "s(Qres):streamst9",
  "s(DOYres)"),
  values  = list(
    Qres = 0,
    DOYres = 8.1, # this is the median
    stream = "st6")) 



```

**NOTE: TEMPERATURE AXIS ARE FLIPPED SO LARGE #'S ARE WARM**
During amb year, ER increases with daily temp and with annual temp across streams. Duing the N and P years, this switches. We see ER increasing with annual temp again, but decreasing with daily temp. Cool

```{r}
#note inversing temp axis here
ErPredSur <- ggplot() +
  geom_raster(data = gamPred2, aes(y = invKT.C.StMean, x = Tanom.iktCs, fill = fit)) +
  facet_wrap(vars(treatment)) +
  scale_fill_distiller(palette = "Spectral") +
  geom_point(data = met, aes(y = invKT.C.StMean, x = Tanom.iktCs, size = lER), 
             position = position_jitter(width = 0.05, height = 0.025), shape = 21, fill = "transparent", alpha= 0.5) +
  scale_size_binned(n.breaks = 4, nice.breaks = TRUE) +
  ylab("Mean summer temp (-1 * centered 1/kT)") +
  xlab("Mean daily temp (-1 * centered 1/kT)")
```


```{r}
gamPred2w <- gamPred2 %>% 
  pivot_wider(id_cols = invKT.C.StMean:stream, names_from = treatment, values_from = c(fit,se.fit)) %>% 
  mutate(N_amb = exp(fit_nitrogen)/exp(fit_ambient))

gamPred2w2 <- gamPred2w %>% 
  # select(Tanom.iktCs:doy, fit_ambient:P_amb) %>% 
  pivot_longer(cols = N_amb, names_to = "treatment", values_to = "effectSize") %>% 
  mutate(treatment = fct_recode(treatment, Nitrogen = "N_amb"))
```


```{r}
pdf(file.path(here::here("03_Plots"), "02b_ERpredSurNonly.pdf"), height = 4, width = 8)
ErPredSur
dev.off()
```

```{r}
pdf(file.path(here::here("03_Plots"), "02b_ERdifPlotsNonly.pdf"), height = 5, width = 4)
ggplot() +
  geom_raster(data = gamPred2w2, aes(y = -invKT.C.StMean, x = -Tanom.iktCs, fill = effectSize)) +
  geom_point(data = met, aes(y = -invKT.C.StMean, x = -Tanom.iktCs, size = lER), alpha = 30/100)+
  scale_fill_distiller(palette = "Spectral") +
  ylab("Mean summer temp (-1 * centered 1/kT)") +
  xlab("Mean daily temp (-1 * centered 1/kT)")
dev.off()
```
#########################################
# ER_P
#########################################



# summary Plots
Focusing on this
```{r}
summary(g.er.BEST1$gam)
```


## Predictions
```{r}
pdf(file.path(here::here("03_Plots"), "02c_Ponly_ERpredXdoyPonly.pdf"), height = 4, width = 8)
ggplot() +
  geom_point(data = met, aes(y = lER, x = doy, color = treatment)) +
  geom_line(data = met, aes(y = ERgamPred, x = doy, color = treatment)) +
  geom_ribbon(data = met, aes(ymin = ERgamPred - ERgamPred.se,
                              ymax = ERgamPred + ERgamPred.se, x = doy, fill = treatment), alpha = 0.5) +
  facet_grid(stream ~ treatment)
dev.off()
```



Use model to predict log ER, with Qres, doy and stream set to "centered" values
```{r}
gamPred2 <- predict_gam(g.er.BEST1$gam,exclude_terms = c(
  "s(Qres):streamst11U",
  "s(Qres):streamst18",
  "s(Qres):streamst6",
  "s(Qres):streamst9",
  "s(DOYres):streamst11U",
  "s(DOYres):streamst18",
  "s(DOYres):streamst6",
  "s(DOYres):streamst9"),
  values  = list(
    Qres = 0,
    DOYres = 6.1)) # median in this dataset



```

```{r}
gamPred2w <- gamPred2 %>% 
  pivot_wider(id_cols = c(stream:Tanom.iktCs,Qres:DOYres), names_from = treatment, values_from = c(fit,se.fit)) %>% 
  mutate(P_amb = exp(fit_phosphorus)/exp(fit_ambient))

# gamPred2w2 <- gamPred2w %>% 
#   select(Tanom.iktCs:doy, fit_ambient:P_amb) %>% 
#   pivot_longer(cols = N_amb:P_amb, names_to = "treatment", values_to = "effectSize") %>% 
#   mutate(treatment = fct_recode(treatment, Nitrogen = "N_amb", Phosphorus = "P_amb"))
```

**NOTE: TEMPERATURE AXIS ARE FLIPPED SO LARGE #'S ARE WARM**
During amb year, ER increases with daily temp and with annual temp across streams. Duing the N and P years, this switches. We see ER increasing with annual temp again, but decreasing with daily temp. Cool

```{r}
ErPredSur <- ggplot() +
  geom_raster(data = gamPred2, aes(y = stream, x = -Tanom.iktCs, fill = fit)) +
  scale_fill_distiller(palette = "Spectral") +
  geom_point(data = met, aes(y = stream, x = -Tanom.iktCs, size = lER), 
             position = position_jitter(width = 0, height = 0.075), shape = 21, fill = "transparent", alpha= 0.5) +
  scale_size_binned(n.breaks = 4, nice.breaks = TRUE) +
  facet_wrap(vars(treatment)) +
  ylab("stream") +
  xlab("Mean daily temp (-1 * centered 1/kT)")
```


```{r}
pdf(file.path(here::here("03_Plots"), "02c_ERpredSur_Ponly.pdf"), height = 4, width = 8)
ErPredSur
dev.off()
```

```{r}
pdf(file.path(here::here("03_Plots"), "02c_ERdifPlotsPonly.pdf"), height = 5, width = 4)
ggplot() +
  geom_raster(data = gamPred2w, aes(y = stream, x = -Tanom.iktCs, fill = P_amb)) +
  geom_point(data = met, aes(y = stream, x = -Tanom.iktCs, size = lER), 
             position = position_jitter(width = 0, height = 0.075), shape = 21, fill = "transparent", alpha= 0.5) +
  scale_fill_distiller(palette = "Spectral") +
  xlab("Mean daily temp (-1 * centered 1/kT)") +
  ylab("stream") 
dev.off()
```


#########################################
# GPP_All
#########################################
# These aren't too useful for interactions
```{r echo=FALSE}
GPPquickPlot1 <- plot(g.gpp.BEST1.gv, allTerms = 1) +
  l_points(color = "pink") + l_fitLine(linetype = 3) + l_fitContour() + 
      l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + theme_get() + labs(title = NULL)

  print(GPPquickPlot1, pages = 1)
```

```{r echo=FALSE}
GPPquickPlot2 <- plot(g.gpp.BEST2.gv, allTerms = 1) +
  l_points(color = "pink") + l_fitLine(linetype = 3) + l_fitContour() + 
      l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + theme_get() + labs(title = NULL)

  print(GPPquickPlot2, pages = 1)
```

```{r echo=FALSE}
GPPquickPlot3 <- plot(g.gpp.BEST3.gv, allTerms = 1) +
  l_points(color = "pink") + l_fitLine(linetype = 3) + l_fitContour() + 
      l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + theme_get() + labs(title = NULL)

  print(GPPquickPlot3, pages = 1)
```



Predict for future ploting
```{r echo= FALSE}
met$GPPgamPred <- predict.gam(g.gppBEST3$gam, met, type = "response")
met$GPPgamPred.se <- predict.gam(g.gppBEST3$gam, met, type = "response", se.fit = TRUE)$se.fit
```


#########################################
# GPP_N
#########################################

# These aren't too useful for interactions
```{r echo=FALSE}
GPPquickPlot1 <- plot(g.gpp.BEST1.gv, allTerms = 1) +
  l_points(color = "pink") + l_fitLine(linetype = 3) + l_fitContour() + 
      l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + theme_get() + labs(title = NULL)

  print(GPPquickPlot1, pages = 1)
```

```{r echo=FALSE}
GPPquickPlot2 <- plot(g.gpp.BEST2.gv, allTerms = 1) +
  l_points(color = "pink") + l_fitLine(linetype = 3) + l_fitContour() + 
      l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + theme_get() + labs(title = NULL)

  print(GPPquickPlot2, pages = 1)
```

```{r echo=FALSE}
GPPquickPlot3 <- plot(g.gpp.BEST3.gv, allTerms = 1) +
  l_points(color = "pink") + l_fitLine(linetype = 3) + l_fitContour() + 
      l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + theme_get() + labs(title = NULL)

  print(GPPquickPlot3, pages = 1)
```

