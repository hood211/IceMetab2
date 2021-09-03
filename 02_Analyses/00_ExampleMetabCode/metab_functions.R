####################################
# CALCULATES O2 SATURATION, given temperature and BP (mmHg). From Garcia and Gordon 1992 L&O
osat  <- function(temp,bp){
  sato<-(exp(2.00907 + 3.22014 * (log((298.15-temp) / (273.15 + temp))) + 4.0501 * (log((298.15 - temp) / (273.15 + temp))) ^ 2 + 4.94457 * (log((298.15 - temp) / (273.15 + temp))) ^ 3 - 0.256847 * (log((298.15 - temp) / (273.15 + temp))) ^ 4 + 3.88767 * (log((298.15 - temp) / (273.15 + temp))) ^ 5)) * 1.4276 * bp / 760
  
  sato
}

# osat = function(C,bp){ #Alex's function
# sato = 1.834355-0.305725*C+0.016474*bp
# sato
# }

# CONVERT K600 to value at specific temp by Schmidt number scaling #ASK JIM
Kcor  <- function(temp, meanK_O2_T20)	{
  meanK_O2_T20 * 1.024^(temp - 20)
}

#####################################



#####  TWOSTATIONTUNINGGW  ##############################
#########################################################
# CREATES A PLOT FOR EXPLORING PRIORS
twostationtuningGW <-function(pmax, alpha, ER, todays_master, K, z,  bp, tt, up, down, Qgw, DOgw) {
  oxydf = twostationsatdfGW(pmax, alpha, ER, todays_master, K, z, tt, up, down, Qgw, DOgw)
  #modeledO2 = twosbayesplot2(oxydf = oxydf)
  #PredictActual = predictedActualPlotFun(oxydf = oxydf)
  #GPPlightPlot = GPPlightPlotFun(df = oxydf)
  #TuningPlots = grid.arrange(modeledO2, PredictActual, GPPlightPlot)
  #TuningPlots
}
#########################################################
#########################################################


######  TWOSTATIONSATDFGW  #################################################################################################################################################
#########################################################################################################################################################################

# CREATES DATA.FRAME WITH MEASURED AND MODELED DO DATA
twostationsatdfGW <-function(pmax, alpha, ER, todays_master, K, z, tt, up, down, Qgw, DOgw) {
  
  #number of 5 min readings bewteen up and down probe corresponding to travel time tt  #ASK JIM IF THIS IS STILL CORRECT
  lag<-round(mean(todays_master[todays_master$site== "US",]$tt_d, na.rm = TRUE)/time_bw_loggings_d) 
  
  ##trim the ends of the oxy and temp data by the lag so that oxydown[1] is the value that is the travel time later than oxy up.  The below calls are designed to work with our data structure exactly as I sent it to you.
  updata<-todays_master[todays_master$site== "US",]
  tempup<-updata$temp[1:(length(updata$temp)-lag)] # trim the end by the lag
  oxyup<-updata$DO_mgL[1:(length(updata$temp)-lag)]
  BPup<-updata$mmHg[1:(length(updata$temp)-lag)] #BP at T0 from US logger, mmHg
  lightup<-updata$lux[1:(length(updata$temp)-lag)]
  depthup<-updata$depth_m[1:(length(updata$temp)-lag)]
  TTup<-updata$tt_d[1:(length(updata$temp)-lag)]
  Areaup<-updata$area_m2[1:(length(updata$temp)-lag)]
  
  downdata<-todays_master[todays_master$site== "DS",]
  tempdown<-downdata$temp[(1+lag):length(downdata$temp)]
  oxydown<-downdata$DO_mgL[(1+lag):length(downdata$temp)]
  BPdown<-downdata$mmHg[(1+lag):length(downdata$temp)]  #BP at T0 + tt from DS logger, mmHg
  lightdown<-downdata$lux[(1+lag):length(downdata$temp)]
  depthdown<-downdata$depth_m[(1+lag):length(downdata$temp)]  
  TTdown<-downdata$tt_d[(1+lag):length(downdata$temp)]
  Areadown<-downdata$area_m2[(1+lag):length(downdata$temp)]
  timedown<-downdata$Pdt[(1+lag):length(downdata$temp)]
  
  bp <- (BPup + BPdown)/2
  light <- (lightup + lightdown)/2
  z = (depthup + depthdown)/2
  tt = (TTup + TTdown)/2
  A = (Areaup + Areadown)/2
  
  #adding in a new column to receive the metab data in function below; it's really 'mass balanced' DO downstream
  metab<-numeric(length(oxyup))
  
  for (i in 1:length(oxyup))  {metab[i] <- (oxyup[i] + 
                            (pmax*tanh(alpha*light[i]/pmax)/z[i]) + 
                            ((ER/z[i])*tt[i]) + 
                            (Kcor(tempup[i], K)*tt[i])*((osat(tempup[i], bp[i]) + osat(tempdown[i], bp[i]) - oxyup[i])/2) +
                            ((Qgw*tt[i])/(A[i]*z[i]) * (DOgw - (oxyup[i]/2))))/
                            (1 + (Qgw*tt[i])/(A[i]*z[i]*2)+((Kcor(tempup[i],K)*tt[i])/2))
    }
  
  #adding in a new column to receive the GPP_diff data in function below			
  GPP_diff =  numeric(length(oxyup)) #this is GPP per unit travel time (g DO/m2/tt) calculated by mass balance
  
  for (i in 1:length(oxyup))  {
    GPP_diff[i] = (oxydown[i] * (1 + Qgw*tt[i]/(A[i]*z[i]*2) + Kcor(tempup[i],K)*tt[i]/2) -
                     (oxyup[i] +
                        (ER)*tt[i]/z[i]+
                            (Kcor(tempup[i],K))*tt[i]*((osat(tempup[i],bp[i])-oxyup[i]+ osat(tempdown[i],bp[i])))/2 + 
                                Qgw*tt[i]/(A[i]*z[i])*(DOgw - (oxyup[i]/2)))) *z[i]
  }
  
  #adding in a new column to receive the GPP_tt data in function below				
  GPP_tt <- numeric(length(oxyup))
  
  #this little chunk is basically the predicted GPP with light, given parameters we entered above for pmax and alpha
  for(i in 1:length(oxyup)) {GPP_tt[i] <- pmax * tanh(alpha*light[i] / pmax)} #this is predicted GPP per unit travel time g DO/m2/tt 
  
  
  AnalTime = Sys.time() # tz is UTC
  
  TintLight = (light*60*60*24) * time_bw_loggings_d #now it's Lux per logging time interval 
  
  oxysat<-osat(tempdown,bp)
  oxymodel<-data.frame (timedown, metab, oxydown, oxyup, oxysat, tempdown, tempup, light, bp, GPP_tt, GPP_diff, tt, z, A, Qgw, DOgw, AnalTime, TintLight, alpha, pmax)

  oxymodel
}

####################################################################################################################################################################
####################################################################################################################################################################


#######  TWOSTATIONSATMETPOSTGW  #############################################################################################################################################################
####################################################################################################################################################################
# PERFORMS THE MCMC
twostationsatmetpostGW <- function(start, todays_master, Kmean, Ksd, nbatch, scale, pmaxmean, pmaxsd, alphamean, alphasd, ERmean, ERsd, DOgwmean, DOgwsd, Qgw) {
  
  #number of 5 min readings bewteen up and down probe corresponding to travel time tt
  lag<-round(mean(todays_master[todays_master$site== "US",]$tt_d, na.rm = TRUE)/time_bw_loggings_d)  ###### What is this number?  bob had 0.0037222
  
  ##trim the ends of the oxy and temp data by the lag so that oxydown[1] is the value that is the travel time later than oxy up.  The below calls are designed to work with our data structure exactly as I sent it to you.
  updata<-todays_master[todays_master$site== "US",]
  tempup<-updata$temp[1:(length(updata$temp)-lag)] # trim the end by the lag
  oxyup<-updata$DO_mgL[1:(length(updata$temp)-lag)]
  BPup<-updata$mmHg[1:(length(updata$temp)-lag)] #BP at T0 from US logger, mmHg
  lightup<-updata$lux[1:(length(updata$temp)-lag)]
  depthup<-updata$depth_m[1:(length(updata$temp)-lag)]
  TTup<-updata$tt_d[1:(length(updata$temp)-lag)]
  Areaup<-updata$area_m2[1:(length(updata$temp)-lag)]
  
  downdata<-todays_master[todays_master$site== "DS",]
  tempdown<-downdata$temp[(1+lag):length(downdata$temp)]
  oxydown<-downdata$DO_mgL[(1+lag):length(downdata$temp)]
  BPdown<-downdata$mmHg[(1+lag):length(downdata$temp)]  #BP at T0 + tt from DS logger, mmHg
  lightdown<-downdata$lux[(1+lag):length(downdata$temp)]
  depthdown<-downdata$depth_m[(1+lag):length(downdata$temp)] 
  TTdown<-downdata$tt_d[(1+lag):length(downdata$temp)]
  Areadown<-downdata$area_m2[(1+lag):length(downdata$temp)]
  
  bp <- (BPup + BPdown)/2
  light <- (lightup + lightdown)/2
  z = (depthup + depthdown)/2
  tt = (TTup + TTdown)/2
  A = (Areaup + Areadown)/2
  
  
  #peform MCMC
  met.post<-metrop(tspostsatGW, initial=start, nbatch=nbatch, scale=scale, tempup=tempup, tempdown=tempdown, oxyup=oxyup, oxydown=oxydown,  z=z, bp=bp, light=light, tt=tt, Kmean=Kmean, Ksd=Ksd, pmaxmean= pmaxmean, pmaxsd= pmaxsd, alphamean= alphamean, alphasd= alphasd, ERmean= ERmean, ERsd= ERsd, DOgwmean = DOgwmean, DOgwsd = DOgwsd, Qgw = Qgw, A = A)
  met.post
  
}   ##end of function


########  TSPOSTSATGW  ############################################################################################################################################################
##################################################################################################################################################################################

# CALCULATES PRIOR AND POSTERIOR LOG LIKLIHOODS
tspostsatGW <- function(MET, tempup, tempdown, oxyup, oxydown, light, bp, tt, z, Kmean, Ksd, pmaxmean, pmaxsd, alphamean, alphasd, ERmean, ERsd, DOgwmean, DOgwsd, Qgw, A){
  
  ##assign the parameters we solve for to easy to undertsand values
  pmax<-MET[1]  # units mg/m2/travel time divide by tt to get
  alpha <- MET[2]
  ER<-MET[3]
  K<-MET[4]
  DOgw <- MET[5]
  sigma<-exp(MET[6]) #always model the log of variance so that one does not get a negative standard deviation
  
  lag<-round(mean(todays_master[todays_master$site== "US",]$tt_d, na.rm = TRUE)/time_bw_loggings_d) 
  
  metab<-numeric(length(todays_master$DO_mgL)) #create empty vector
  
  #below is equation 8 in the derivation, solving for downstream O2 at at each 5 min interval.  It references other functions:  Kcor converts K600 to KO2 for a given temperature.  osat calculates oxygen saturation for a given temp and bp.  --pmax has units mg DO/m2/tt; alpha has units mg DO/m2/tt/light and is the initial slope of the P - I relationship; new GW term based on Hall and Tank.
  for (i in 1:length(oxyup))  {metab[i]<-(oxyup[i]+ ((pmax/1e4)*tanh((alpha/1e7)*light[i] / (pmax/1e4)))/z[i] + (ER/1e2)*tt[i]/z[i]+(Kcor(tempup[i],K))*tt[i]*((osat(tempup[i],bp[i])-oxyup[i]+osat(tempdown[i],bp[i])))/2 + Qgw*tt[i]/(A[i]*z[i])*((DOgw/1e2) - oxyup[i]/2))/(1+ Qgw*tt[i]/(A[i]*z[i]*2) + Kcor(tempup[i],K)*tt[i]/2) }
  
  #likelhood is below.  dnorm caculates the probablity density of a normal distribution
  loglik<-sum(dnorm(oxydown, metab, sigma, log=TRUE))
  
  ##priors, 
  prior<-log(dnorm(pmax,mean= pmaxmean, sd = pmaxsd))+log(dnorm(alpha,mean= alphamean, sd=alphasd))+log(dnorm(ER,mean= ERmean,sd=ERsd))+log(dnorm(K,mean=Kmean,sd=Ksd)) +log(dnorm(DOgw,mean=DOgwmean,sd=DOgwsd))	 	 
  
  loglik+prior #because we deal wih log probabilities, we add them to calculate the posterior
}  #end of function

##################################################################################################################################################################################
##################################################################################################################################################################################









########   TWOSTATIONOUTSUM   ###################################
#################################################################

# OUTPUTS SUMMARY OF THE MCMC
twostationoutsum <-  function(met.post, burnin, nbatch) {
  accept<-met.post$accept
  accept4df <- c(as.numeric("NA"),as.numeric("NA"), as.numeric("NA"), accept, as.numeric("NA"))
  pmaxr<-c(pmaxmean, pmaxsd,quantile(met.post$batch[(burnin:nbatch),1], c(0.025, 0.5, 0.975))) /1e4
  alphar <- c(alphamean, alphasd,quantile(met.post$batch[(burnin:nbatch),2], c(0.025, 0.5, 0.975))) /1e7
  err<-c(ERmean, ERsd, quantile(met.post$batch[(burnin:nbatch),3], c(0.025, 0.5, 0.975))) / 1e2
  Kr<-c(Kmean, Ksd, quantile(met.post$batch[(burnin:nbatch),4], c(0.025, 0.5, 0.975))) 
  DOgw <-c(DOgwmean, DOgwsd, quantile(met.post$batch[(burnin:nbatch),5], c(0.025, 0.5, 0.975))) / 1e2
  sr<- c(as.numeric("NA"),as.numeric("NA"),quantile(met.post$batch[(burnin:nbatch),6], c(0.025, 0.5, 0.975)))
  outsum <- rbind(accept4df, pmaxr, alphar, err, Kr, DOgw, sr)
  outsumdf <- data.frame(outsum)
  names(outsumdf) <- c("PriorMean", "PriorSD", "L95per","median", "U95per")
  outsumdf
}		

################################################################
################################################################

################################################################################################################################
################################################################################################################################
# CALCULATES MEDIAN GPP AND 95% CIs
twostationGPPests3 <- function(met.post, burnin, nbatch, oxydf, tt) {
  pmax_median <- quantile(met.post$batch[(burnin:nbatch), 1],0.5)
  alpha_median <- quantile(met.post$batch[(burnin:nbatch), 2],0.5)
  
  # fuction calculates daily GPP
  dailyGPPfun <- function(pmax, alpha, light, tt) {
    light <- oxydf$light
    tt = oxydf$tt
    GPP_int <- (pmax/1e4)*tanh((alpha/1e7)* light/(pmax/1e4)) / tt * time_bw_loggings_d # GPP for each interval
    dailyGPP <- sum(GPP_int) #units mgO2/m2/d
    dailyGPP
  }
  
  # make bootstrap df
  dailyGPPboot <- numeric(length(10000))
  
  # bootstrap
  for(i in 1:10000){
    RowNums = seq(1,dim(met.post$batch[(burnin:nbatch),])[1]) # choosing a line because these are not modeled seperatly
    ranRow = sample(RowNums, size = 1, replace = TRUE)
    pmax <- met.post$batch[ranRow, 1]
    alpha <- met.post$batch[ranRow, 2]
    dailyGPP_i <- dailyGPPfun(pmax, alpha, light, tt = tt)
    dailyGPPboot[i] <- dailyGPP_i
  }
  
  dailyGPPbymed <- dailyGPPfun(pmax = pmax_median, alpha = alpha_median, tt = tt)
  dailyGPPbymeddf <- data.frame(as.numeric("NA"), as.numeric("NA"), as.numeric("NA"), dailyGPPbymed, as.numeric("NA"), row.names = "bymed_GPP")
  names(dailyGPPbymeddf) <- c("PriorMean", "PriorSD", "L95per", "median", "U95per")
  
  
  GPPr <- quantile(dailyGPPboot, c(0.025, 0.5, 0.975))
  GPPrdfBS<- data.frame(as.numeric("NA"), as.numeric("NA"), GPPr[1], GPPr[2], GPPr[3], row.names = "BS_GPP")
  names(GPPrdfBS) <- c("PriorMean", "PriorSD", "L95per", "median", "U95per")
  GPPrdf <- rbind(dailyGPPbymeddf, GPPrdfBS)
  
  #dailyGPPhist <<- hist(dailyGPPboot)
  GPPrdf
}

################################################################################################################################
################################################################################################################################



###################################################################################################################################################################
#####    PLOTS    #################################################################################################################################################
###################################################################################################################################################################

# PLOTS MODELED AND MEASURED DO
twosbayesplot2<-function(oxydf) {
  
  oxyplot <- ggplot(oxydf, aes(y = metab, x = timedown)) + 
    geom_line(color = "red", size = 1.5) +
    geom_line(data = oxydf, aes(y = oxydown, x =timedown), color = "black", size = 0.75) +
    geom_point(data = oxydf, aes(y = oxydown, x =timedown), color = "black", fill = "red", shape = 21) +
    geom_point(data = oxydf, aes(y = oxyup, x =timedown), color = "blue") +
    geom_line(data = oxydf, aes(y = oxyup, x =timedown), color = "blue", size = 0.75) +
    xlab("Time") +
    ylab("Dissolved oxygen (mg/L)") 
  
  oxyplot
}

#########################################

# SCATTERPLOT OF PREDICTED AND MODELED DO
predictedActualPlotFun <- function(oxydf) {
  
  predlm <- lm(oxydown ~ metab, oxydf)
  
  intercept = round(coef(predlm)[1], digits = 2)
  slope = round(coef(predlm)[2], digits = 2)
  R2 = round(summary(predlm)$r.squared, digits = 2)
  
  predActPlot <- ggplot(oxydf, aes(y = oxydown, x = metab, fill = log10(light+1))) +
    geom_point(size = 3, color = "black", shape = 21) +
    geom_abline(intercept = 0, slope = 1) +
    stat_smooth(method = "lm") +
    xlab("Modeled downstream DO (mg/L)") +
    ylab("Measured downstream DO (mg/L)") +
    scale_fill_gradient(low = "black", high = "yellow") +
    annotate("text", x = quantile(oxydf$metab, p = 0.3), y = quantile(oxydf$oxydown, p = 0.95), label = intercept) +
    annotate("text", x = quantile(oxydf$metab, p = 0.3), y = quantile(oxydf$oxydown, p = 0.85), label = slope) +
    annotate("text", x = quantile(oxydf$metab, p = 0.3), y = quantile(oxydf$oxydown, p = 0.75), label = R2) +
    annotate("text", x = quantile(oxydf$metab, p = 0.1), y = quantile(oxydf$oxydown, p = 0.95), label = "Intercept = ") +
    annotate("text", x = quantile(oxydf$metab, p = 0.1), y = quantile(oxydf$oxydown, p = 0.85), label = "Slope = ") +
    annotate("text", x = quantile(oxydf$metab, p = 0.1), y = quantile(oxydf$oxydown, p = 0.75), label = "R2 = ")  +
    theme(legend.position = c(0.9,0.3)) 
  
  
  predActPlot
}


#########################################

# PLOTS P-I CURVE FOR EACH INTERVAL CALCULATED BY BALANCE AND FROM P-I CURVE MODEL
GPPlightPlotFun = function(df) {
  
  df$GPP_h = df$GPP_tt/(df$tt)   #GPP predicted from model, in days
  df$GPP_diff_h = df$GPP_diff/(df$tt) #GPP calculated from "mass balance"
  
  PredActuallm = lm(GPP_diff_h ~ GPP_h, df)
  R2 = round(summary(PredActuallm)$r.squared, digits = 2)
  
  GPPlightPlot = ggplot(df, aes(y = GPP_h, x = light)) +
    geom_point(data = df, aes(y = GPP_diff_h, x = light), shape = 21, fill = "light green", size = 3, alpha = 80/100) + #coord_cartesian(xlim = c(0, 375), ylim = c(-5, 3)) +
    geom_line(color= "red", size = 2) +
    ylab(expression(paste("GPP (mg ",O[2], " ",m^-2," ",d^-1,")"))) +
    xlab(expression(paste("LUX"))) +
    annotate("text", x = quantile(df$light, p = 0.96), y = quantile(df$GPP_diff_h, p = 0.4), label = "R2 = ") +
    annotate("text", x = quantile(df$light, p = 0.97), y = quantile(df$GPP_diff_h, p = 0.4), label = R2) +
    annotate("text", x = quantile(df$light, p = 0.96), y = quantile(df$GPP_diff_h, p = 0.5), label = "Mass balance GPP v. modeled GPP") 
  
  GPPlightPlot
}

#########################################

################################################################

# MCMC TIME SERIES PLOT
twostationoutplot <- function(met.post) {
  plot(ts(met.post$batch))
}

################################################################

#########################################

# NIGHT TIME REGRESSION FUNCTION; intercept = ER (g O2/m3/d); slope = K (1/d); mg DO/L = g DO/m3
nightreggression <- function(data, bp){
  temp<-data$C
  oxy<-data$DO
  
  ##moving average on oxy data
  oxyf1<-filter(data$DO,rep(1/3,3), sides=2) # this takes a moving average of a DO meas, the prior meas, and the next meas.
  
  #trime the ends of the oxy data - removes the first and last value
  oxyf2<- oxyf1[c(-1,-length(oxyf1))]
  
  ##calculate delO/delt
  deltaO2<-((oxyf2[-1]-oxyf2[-length(oxyf2)])/time_bw_loggings_d)  # is t2 - t1 / time betwee readings in d
  
  #Trim the first two and last one from the temp data to match the filter oxy data
  temptrim<-temp[c(-2:-1,-length(temp))]
  #calc the dodef
  satdef<-osat(temptrim,bp)-oxyf2[-1]
  
  #calculate regression and plot
  nreg<-lm(deltaO2~satdef)
  plot(satdef,deltaO2)
  abline(nreg)
  
  summary(nreg)
}


#########################################

# BOWDEN METHOD for estimating ER when there is 24 hours of light
BowdenMethodFun <- function(oxydata, tt, z, Kmean, Qgw, DOgw, A){
  oxydown = oxydata$oxydown
  oxyup = oxydata$oxyup
  tempdown = oxydata$tempdown
  tempup = oxydata$tempup
  bp = oxydata$bp
  light = oxydata$light
  
  satdown = osat(tempdown,bp)
  satup = osat(tempup,bp)
  
  
  metabolism = (oxydown - oxyup)/tt*z - Kmean*z* ((satup + satdown)/2 - (oxydown + oxyup)/2 ) - Qgw/A*(DOgw - (oxydown + oxyup)/2) #eq. 7 from Hall and Tank 2005 with GW term set to zero
  
  LightMetabDF = data.frame(light, metabolism)
  LightMetabDFlowLight = LightMetabDF[LightMetabDF$light <= 50,]
  
  BowdenLm = lm(metabolism ~ light, LightMetabDFlowLight)
  estER = round(coef(BowdenLm)[1], digits = 2)
  
  BowdenPlot <- ggplot(LightMetabDFlowLight, aes(y = metabolism, x = light)) +
    geom_point(shape = 21, fill = "green", size = 2) +
    xlab("PAR") +
    ylab("Ecosystem metabolism") +
    annotate("text", x = quantile(LightMetabDFlowLight$light, p = 0.8), y = quantile(LightMetabDFlowLight$metabolism, p = 0.975), label = estER) +
    annotate("text", x = quantile(LightMetabDFlowLight$light, p = 0.4), y = quantile(LightMetabDFlowLight$metabolism, p = 0.975), label = "Estimated ER =") +
    xlim(0, max(LightMetabDFlowLight$light)*1.25) +
    stat_smooth(method = "lm")
  BowdenPlot
}

#########################################


###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################

