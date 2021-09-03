#Cross, Hood et al. 'Eutrophication intensifies the effects of warming on metabolic balance of running water ecosystems'
#This is the code used to estimate metabolism (in g O2 m-2 d-1) for each day of the study. This code can be altered to run
#multiple dates in one run by adding a repeat feature. Here, the code is set to demonstrate a single day
#(July 21, 2015) in stream 6; 

# SET SYSTEM TIME
  Sys.setenv(TZ = "UTC")

# LOAD LIBRARIES
  library(tidyverse)
	library(chron)
	library(mcmc)
	library(zoo)
	library(plyr)
	library(reshape2)
	library(gridExtra)
	library(scales)
  library(lubridate)

# LOAD FUNCTIONS
	source(file = "metab_functions.R")

	################
	#load data file for a stream and year
	masterfile <- read.csv("stream6_2017_metabolism.csv")

   
	#format date; choose appropriate line below, or alter if you have a different format
  masterfile$Pdt <- as.POSIXct(masterfile$Pdt, format = "%Y-%m-%d %H:%M", tz = "UTC")
  #masterfile$Pdt <- as.POSIXct(masterfile$Pdt, format = "%m/%d/%y %H:%M", tz = "UTC")
  
  #add julian day column
	masterfile$julian <- yday(masterfile$Pdt)
	
	#The following stream days were NOT run because discharge exceeded the 95th percentile:
	#(st11) 2015, julian 194, 207; 2016, julian 214; 2017, julian 195, 196, 199, 200, 201, 220, 222)  
	#(st6) 2015, julian 222; 2017, julian 222
  #(st9) 2015, julian 200, 222; 2016, julian 225; 2017, julian 195, 199, 222
	#(st18) 2015, julian 222

  #enter julian date 
	startJulian <- 202
  
  #select day
	todays_master <- filter(masterfile, julian == startJulian)
	
	#format a few data types if needed
	todays_master$depth_m <- as.numeric(as.character(todays_master$depth_m))
	todays_master$velocity_m.s <- as.numeric(as.character(todays_master$velocity_m.s))
	todays_master$tt_s <- as.numeric(as.character(todays_master$tt_s))
	
	#calculate Q-groundwater
	Qgw = ifelse(mean(todays_master$Q_DS_m3.s - (todays_master$Q_DS_m3.s*todays_master$Q_US_div_DS))>0, mean(todays_master$Q_DS_m3.s - (todays_master$Q_DS_m3.s*todays_master$Q_US_div_DS)), 0) *60*60*24
	
 	#set stream length
	length = mean(todays_master$length_m)

  #get a few other units right - days, meters, and mmHg (except for L)
	todays_master$Q_DS_LperD = todays_master$Q_DS_m3.s*60*60*24 # discharge in liters per day
	todays_master$tt_d = todays_master$tt_s/60/60/24	 #travel time in units days
	mean_tt <- as.numeric(mean(todays_master$tt_d))
 	time_bw_loggings_d <- 1/60/24 # units days 1 mins/ 60 min per h / 24 h per d

	
# INITIAL PLOTS TO EXAMINE THE DATA
	O2plot <- ggplot(todays_master, aes(y = DO_mgL, x = Pdt, color = site)) +
		geom_ribbon(aes(ymin = 8, ymax = (lux/80000)+8), fill = "yellow", show.legend = FALSE) +	
		geom_line(size = 3) +
	  scale_x_datetime(labels = date_format("%H:%M:%S"), date_breaks = "6 hours") +
    theme(legend.justification = c(0,0), legend.position=c(0,0), legend.title = element_blank())  
	  
	tempplot <- ggplot(todays_master, aes(y = temp, x = Pdt, color = site)) +	
		geom_line(size = 3)+
	  scale_x_datetime(labels = date_format("%H:%M"), date_breaks = "6 hours") +
		theme(legend.justification = c(0,1), legend.position=c(0,1), legend.title = element_blank())  
    
	o2satplot <- ggplot(todays_master, aes(y = persat, x = Pdt, color = site)) +	
		geom_line(size = 3)+
		scale_x_datetime(labels = date_format("%H:%M"), date_breaks = "6 hours") +
		theme(legend.justification = c(0,1), legend.position=c(0,1), legend.title = element_blank())  

	lightplot <- ggplot(todays_master[todays_master$site == "US",], aes(y = lux, x = Pdt)) +
		geom_line(size = 1.25, color = "blue") + 
	  scale_x_datetime(labels = date_format("%H:%M"), date_breaks = "6 hours") 
		
	hydroplot = ggplot(todays_master[todays_master$site == "US",], aes(y = Q_DS_m3.s, x = Pdt))+
		geom_line(size = 1.25, color = "blue") +
	 scale_x_datetime(labels = date_format("%H:%M"), date_breaks = "6 hours") 

	TTplot = ggplot(todays_master[todays_master$site == "US",], aes(y = tt_s, x = Pdt))+
		geom_line(size = 1.25, color = "blue") +
		scale_x_datetime(labels = date_format("%H:%M"), date_breaks = "6 hours") 

	Aplot = ggplot(todays_master[todays_master$site == "US",], aes(y = area_m2, x = Pdt))+
		geom_line(size = 1.25, color = "blue") +
		scale_x_datetime(labels = date_format("%H:%M"), date_breaks = "6 hours") 

	 # Set groundwater DO concentrations - this was altered based on the date and the groundwater (see priors data)
	 	 DOgwmean =osat(mean(todays_master$temp, na.rm = TRUE), mean(todays_master$mmHg, na.rm = TRUE))  *  0.9 *1e2   # 0.635 *1e2 # average %sat in springs, and scaling  # 7.7 mg/L
		 DOgwsd = osat(mean(todays_master$temp, na.rm = TRUE), mean(todays_master$mmHg, na.rm = TRUE)) *0.35 *1e2 # SD for %sat in springs		 
		 
   # Define K prior and SD
		Kmean = if((mean(todays_master$Q_DS_m3.s)*31772 + 3.5836) > 629.88) {
	    Kmean = 629.88
	   } else if((mean(todays_master$Q_DS_m3.s)*31772 + 3.5836) < 179.46) {
	     Kmean = 179.46
	   } else {
	     Kmean = (mean(todays_master$Q_DS_m3.s)*31772 + 3.5836)
	   }

	 	   Ksd= 0.25 * Kmean

   #PRIORS TUNING ER, pmax, and alpha (these are dummy values to start)
 	   ERTunning = 0/1e2 
 	   pmaxTuning = 1/1e2  
 	   alphaTunning = 1/1e7 
 	   DOgwTunning = DOgwmean/1e2
 	   
		twostationtuningGW(pmax = pmaxTuning, alpha = alphaTunning , ER = ERTunning, todays_master = todays_master, K = Kmean, z = z, bp = bp, tt = tt, up = "up", down = "down", Qgw = Qgw, DOgw = DOgwTunning)
		poop <- twostationsatdfGW(pmax = pmaxTuning, alpha = alphaTunning , ER = ERTunning, todays_master = todays_master, K = Kmean, z = z, tt = tt, up = "up", down = "down", Qgw = Qgw, DOgw = DOgwTunning)

    #starting to set the priors for ER; this determines where ER should be, based on fact that GPP should go to zero when there is no light	
		pooplowlight <- filter(poop, light < 1500)
		mean_GPP_diff_low_light <- mean(pooplowlight$GPP_diff) #mean change in DO at night after accounting for gas exchange and groundwater
		ER_CHECK <- mean_GPP_diff_low_light/mean(pooplowlight$tt)
		tuned_ER <- -(ER_CHECK)/-0.01
		
		#AGAIN PRIORS TUNING ER, pmax, and alpha (incorporating the 'tuned ER')
		ERTunning = tuned_ER/1e2 
		pmaxTuning = 1/1e4  
		alphaTunning = 1/1e7 
		DOgwTunning = DOgwmean/1e2
	
		twostationtuningGW(pmax = pmaxTuning, alpha = alphaTunning , ER = ERTunning, todays_master = todays_master, K = Kmean, z = z, bp = bp, tt = tt, up = "up", down = "down", Qgw = Qgw, DOgw = DOgwTunning)
		poop <- twostationsatdfGW(pmax = pmaxTuning, alpha = alphaTunning , ER = ERTunning, todays_master = todays_master, K = Kmean, z = z, tt = tt, up = "up", down = "down", Qgw = Qgw, DOgw = DOgwTunning)
			
		#nls to determine alpha and pmax parameters 
		poop_for_nls <- select(poop, light, GPP_diff, tt)
		poop_for_nls$GPP_green <- poop_for_nls$GPP_diff/poop_for_nls$tt
		poop_for_nls_data <- select(poop_for_nls, light, GPP_green)
		light <- poop_for_nls_data$light
		GPP <- poop_for_nls_data$GPP_green
		plot(light, GPP)
		model <- nls(GPP ~ A * tanh(B * light / A), start = list(A = 0.05, B = 0.0001158)) #I tweaked A to work for later dates
		lines(light, predict(model), col = 3)
		
		params <- coef(model)
		tuned_pmax <- params[1]
		tuned_pmax <- as.numeric(tuned_pmax)
		tuned_alpha <- params[2]
		tuned_alpha <- as.numeric(tuned_alpha)
		
		#AND. .. AGAIN PRIORS TUNING ER, pmax, and alpha (now incorporating the 'tuned' pmax and alpha)
		ERTunning = tuned_ER/1e2 
		pmaxTuning = tuned_pmax * mean_tt  
		alphaTunning = tuned_alpha * mean_tt 
		DOgwTunning = DOgwmean/1e2

		twostationtuningGW(pmax = pmaxTuning, alpha = alphaTunning , ER = ERTunning, todays_master = todays_master, K = Kmean, z = z, bp = bp, tt = tt, up = "up", down = "down", Qgw = Qgw, DOgw = DOgwTunning)
		poop <- twostationsatdfGW(pmax = pmaxTuning, alpha = alphaTunning , ER = ERTunning, todays_master = todays_master, K = Kmean, z = z, tt = tt, up = "up", down = "down", Qgw = Qgw, DOgw = DOgwTunning)
		
				
   # this sets ER, pmax, and alpha priors - should be equal to the final tuning values above, multiplied by the denominator
	 ERmean = ERTunning * 1e2    #once K is set do crash run and then tune ER prior using night-time DO - number is ER * 100
	 ERsd = -ERmean * 0.25 #-ERmean * 0.25		#constraining to 25% of mean
	 
	 pmaxmean = pmaxTuning * 1e4  #once ER is set, tune pmax prior using midday DSDO values
	 pmaxsd = pmaxmean * 0.25				#constraining to 25% of mean
	
	 alphamean = alphaTunning * 1e7  #once pmax is set tune alpha prior using AM rise in DS DO
	 alphasd = 	alphamean * 0.25			#constraining to 25% of mean
	 

# BATCH SIZE AND BURN IN for the MCMC; a batch of 200,000 takes about an hour on a basic laptop
		 nbatch = 200000
		 burnin = 2000

 # SCALING VECTOR - get acceptance rate for each parameter on the same scale.
	 scale_vec = c(pmaxmean/2, alphamean/2, -ERmean/2, Kmean/2, DOgwmean/2, 1)
	 scale_mult = 0.08 #0.12  # DECREASE to increase acceptance rate
	 scaleVecMod = scale_vec * scale_mult 
 
 # THE MCMC RUN WITH MODELED GROUNDWATER DO
            metpostfileGW <-twostationsatmetpostGW(start=c(pmaxmean, alphamean, ERmean, Kmean, DOgwmean, 2), todays_master,  Kmean= Kmean, Ksd=Ksd, nbatch=nbatch, scale= scaleVecMod, pmaxmean= pmaxmean, pmaxsd= pmaxsd, alphamean= alphamean, alphasd= alphasd, ERmean= ERmean, ERsd= ERsd, DOgwmean = DOgwmean, DOgwsd = DOgwsd, Qgw = Qgw) #careful the name of this feeds into all of the other functions.
 
# Code for analyzing MCMC
	#summary df from the MCMC run
	    twostationoutsumdf <-suppressWarnings(twostationoutsum(met.post = metpostfileGW,burnin = burnin, nbatch = nbatch))
	    
	#timeseries plot from MCMC run
	#    twostationoutplot(met.post = metpostfileGW)
	 
	#autocorrelation plot
	#	  acf(metpostfileGW$batch, lag.max = 1000)
	 
	#creates a dataframe with DO data, temp, light, modeled (metab) DS data, and GPP_tt
		modeledoxydf <- twostationsatdfGW(pmax= twostationoutsumdf[2,4], alpha = twostationoutsumdf[3,4], ER= twostationoutsumdf[4,4], todays_master, z= z, tt= tt, K= twostationoutsumdf[5,4], Qgw = Qgw, DOgw = twostationoutsumdf[6,4])	
	  
	 # plots to assess model fit 
	 O2fitPlot = twosbayesplot2(oxydf = modeledoxydf)
	 PredActualPlot = predictedActualPlotFun(oxydf = modeledoxydf)
	 PIcurvePlot = GPPlightPlotFun(df = modeledoxydf)

	 	#median GPP and 95 CIs from boostrap and GPP estimated from median pmax and alpha
		twostationGPPestdf <- suppressWarnings(twostationGPPests3(met.post = metpostfileGW, burnin  = burnin, nbatch = nbatch, oxydf = modeledoxydf, tt = tt))
	
	# add GPP to summary datafframe
		twostationoutsumdf2 <- rbind(twostationoutsumdf, twostationGPPestdf)
		
# EXPORT ACTUAL AND MODELED DATA.FRAME
		year <- "2015"
		stream_code <- "st9"
		date_code <- as.Date(startJulian, origin = "2014-12-31")
		
		modeledoxydfW = modeledoxydf
		modeledoxydfW$metab_data = as.character(date_code)
		modeledoxydfW$stream = as.character(stream_code)
		modeledoxydfW$AnalTime = as.character(modeledoxydf$AnalTime)
		modeledoxydfW$timedown = as.character(modeledoxydf$timedown)
		write.csv(modeledoxydfW, paste("ActualModeledDO","_",year, "_",stream_code,"_",startJulian,".csv", sep = ""))

# EXPORTING MODEL RESULTS AND PRIORS
		metab_date <- rep(as.character(date_code), length = dim(twostationoutsumdf2)[1])
		metab_stream <- rep(as.character(stream_code), length = dim(twostationoutsumdf2)[1])
		analysis_date <- rep(as.character(Sys.time()), length = dim(twostationoutsumdf2)[1])
		twostationoutsumdffinal <- cbind( metab_date, metab_stream, analysis_date, twostationoutsumdf2); 
		rowNumbers = seq(1,9, length = 9)
		twostationoutsumdffinal = data.frame(as.character(rownames(twostationoutsumdffinal)),twostationoutsumdffinal, row.names = rowNumbers)
		names(twostationoutsumdffinal)[1] = "parameter"
		write.csv(twostationoutsumdffinal, file = paste("McmcPriorsAndPosteriors","_",year, "_",stream_code,"_",startJulian,".csv", sep = ""))


# EXPORT SUMMARY DATA ROW 
		LightPerDay = sum(modeledoxydf$TintLight, na.rm = TRUE) #units are now umol/m2/d
		meanTemp = mean(rowMeans(data.frame(modeledoxydf $tempdown,modeledoxydf $tempup), na.rm = TRUE))
		meanQds = mean(todays_master$Q_DS_m3.s)
		meanDepthM = mean(todays_master$depth_m)
		meanWidthM = mean(todays_master$width_m)

	  ModelParamDF = data.frame(as.character(date_code) , as.character(stream_code), LightPerDay, meanTemp, meanQds, meanDepthM, meanWidthM, twostationoutsumdffinal[9,8], twostationoutsumdffinal[4,8], twostationoutsumdffinal[2,8], twostationoutsumdffinal[3,8] ,as.character(Sys.time()))
		names(ModelParamDF) = c("date", "stream", "LightPerDay", "meanTemp", "meanQds", "meanDepthM", "meanWidthM", "GPP_mgO2m2d", "ER_mgO2m2d", "pmax", "alpha", "AnalysisTime")	
	
		write.csv(ModelParamDF, paste("MetabSummaryDF","_",year, "_",stream_code,"_",startJulian,".csv", sep = ""))  #I need this file
		
		
		# CREATE SUMMARY PDF OF PLOTS
		pdf(paste(date_code, "J_",startJulian,"st18",".pdf"), paper = "us")
		#title page		
		plot(0:10, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
		text(5, 9, paste("Stream = ", stream_code), cex = 3)
		text(5, 7.5, paste("Julian = ", startJulian), cex = 3)
		text(5, 6, paste("Date = ", date_code), cex = 3)
		text(5,3, paste("Analysis date = ", Sys.time()), cex = 1, col = "red")
		grid.arrange(O2plot, tempplot, o2satplot, lightplot, hydroplot, TTplot, Aplot, nrow = 3, ncol = 3)
		twostationoutplot(met.post = metpostfileGW)
		acf(metpostfileGW$batch, lag.max = 1000)
		grid.arrange(O2fitPlot, PredActualPlot, PIcurvePlot,nrow = 3, ncol = 1)
		#plot(dailyGPPhist)	 
		dev.off()	
		
		
