summary(g.er.BEST1$gam)
library(tidymv)

gamPred <- predict_gam(g.er.BEST1$gam, exclude_terms = c(
                                            "s(Qres):streamst11U",
                                              "s(Qres):streamst18",
                                              "s(Qres):streamst6",
                                              "s(Qres):streamst9",
                                              "s(doy):streamst11U",
                                            "s(doy):streamst18",
                                            "s(doy):streamst6",
                                              "s(doy):streamst9"),
                                        values  = list(
                                          "s(Qres):streamst11U" = 0,
                                          "s(Qres):streamst18" = 0,
                                          "s(Qres):streamst6" = 0,
                                          "s(Qres):streamst9" = 0,
                                          "s(doy):streamst11U" = 200,
                                          "s(doy):streamst18" = 200,
                                          "s(doy):streamst6" = 200,
                                          "s(doy):streamst9" = 200))

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

gamPred2w <- gamPred2 %>% 
  pivot_wider(id_cols = Tanom.iktCs:doy, names_from = treatment, values_from = c(fit,se.fit)) %>% 
  mutate(N_amb = exp(fit_nitrogen) - exp(fit_ambient),
         P_amb = exp(fit_phosphorus) - exp(fit_ambient))

ggplot(gamPred2w, aes(y = invKT.C.StMean, x = Tanom.iktCs, fill = N_amb)) +
  geom_raster()

ggplot(gamPred2w, aes(y = invKT.C.StMean, x = Tanom.iktCs, fill = P_amb)) +
  geom_raster()


ggplot(met, aes(y = lER , x = Qres, color = stream)) +
  geom_point() +
  facet_wrap(vars(treatment))

ggplot(gamPred2, aes(y = invKT.C.StMean , x = Tanom.iktCs, color = exp(fit))) +
  geom_point() +
  facet_wrap(vars(treatment))

# These are turned around so high is warm
ggplot(gamPred2, aes(y = fit , x = -invKT.C.StMean , color = -Tanom.iktCs)) +
  geom_point() +
  # geom_point(data = met, aes(y = lER , x = -invKT.C.StMean), color = "pink") +
  facet_wrap(vars(treatment)) +
  ylab("log ER") +
  xlab("mean summer stream temp (-1/kT)")




M4f.c2aPred <-  predict_gam(M4f.c2a$gam, exclude_terms = c("s(stream)", "ArrheniusC"), 
                            values = list(stream = NULL, ArrheniusC = -0.00150)) 


names(gamPred) <- paste0("gam.", names(gamPred))
gamPred2 <- as.data.frame(cbind(met, gamPred))


ggplot(gamPred2) +
  geom_point(aes(y = lER, x = treatment), color = "red", position = "jitter") +
  geom_point(aes(y = gam.treatment, x = treatment)) 


ggplot(gamPred2) +
  geom_point(aes(y = gam.treatment:invKT.C.StMean:Tanom.iktCs, x = invKT.C.StMean, color= treatment))


ggplot(gamPred2) +
  geom_point(aes(y = lER, x = Tanom.iktCs, color = treatment))

ggplot(gamPred2) +
  geom_point(aes(y = lER,x  = Tanom.iktCs, shape = treatment, color = invKT.C.StMean), position = "jitter")

ggplot(gamPred2) +
  geom_point(aes(y = lER, ))