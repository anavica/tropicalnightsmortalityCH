#############################################################################################
#### Trends in Tropical Nights and their Effects on Mortality in Switzerland across 50 years
#### @Vanessa Rippstein
#### March 2023
#############################################################################################

library(ggplot2); library(colorRamps); library(RColorBrewer);
library(data.table); library(dplyr); library(lubridate); library(reshape2)
library(tidyverse); library(gnm); library(dlnm);

#----------
#Load data:
#----------

#Time Series of tropical nights and number of deaths:
list_80_cant_TN <- readRDS("data/list_80_cant_TN.rds")

########################
# TIME SERIES ANALYSIS #
########################
#---------------------------------------
#crossbasis SPECIFICATION - ADJUSTED BY TEMPERATURE
#---------------------------------------
# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "bs"
vardegree <- 2 
varper <- c(50,90)
# SPECIFICATION OF THE LAG FUNCTION
lag <- 3 

crTN_cant_TN <- list()
crtemp_cant <- list()
mmt_cant <- list()
crtemptemp_c <- list()
mmttemp_c <- list()
names(list_80_cant_TN)


#---------------------------------------
# ANALYSIS 
#---------------------------------------
for (i in 1:length(list_80_cant_TN)){
  
  #import data 
  data <- list_80_cant_TN[[i]]
  
  #additional definitions for ERC and crossbasis 
  varknots <- quantile(data$daily_mean, varper/100, na.rm=T)
  boundknt <- range(data$daily_mean)
  
  
  
  ##crossbasis temperature
  cbtemp <- crossbasis(data$daily_mean, lag=3,
                       argvar=list(fun=varfun,degree=vardegree, knots=varknots, 
                                   Boundary.knots=boundknt),
                       arglag=list(fun="integer"),group=data$year)
  
  data$stratum <- with(data, factor(year),factor(month),factor(dow),factor(District_name))
  ind <- tapply(data$deaths,data$stratum,sum)[data$stratum]
  
  #crossbasis Tropical night
  cbTN <- crossbasis(data$TN.min,lag=3,argvar=list(fun="strata",breaks=1),
                     arglag=list(fun="integer"),group=data$year)
  
  
  
  ## models
  mod0 <- gnm(deaths ~ cbTN+cbtemp, eliminate = stratum,
              family=quasipoisson(),data=data, no.action="na.exclude",
              subset=ind>0)
  modtemp <- gnm(deaths~cbtemp, eliminate=stratum, 
                 family=quasipoisson(),data=data, no.action="na.exclude",
                 subset=ind>0)
  crTN <- crossreduce(cbTN,mod0,cen=0,at=1)
  crtemp0 <- crossreduce(cbtemp,mod0,cen=15)
  crtemptemp <- crossreduce(cbtemp,modtemp,cen=15)
  pct20 <- round(quantile(data$daily_mean, 0.20))
  pct90 <- round(quantile(data$daily_mean, 0.90))
  
  ##set MMT
  mmtrange0 <- crtemp0$predvar[which(crtemp0$predvar==pct20):which(crtemp0$predvar==pct90)] 
  mmtrangetemp <- crtemptemp$predvar[which(crtemptemp$predvar==pct20):which(crtemptemp$predvar==pct90)]
  
  mmt0 <- mmtrange0[which.min(crtemp0$fit[which(crtemp0$predvar%in%mmtrange0)])] #minimum mortality temperautre, temperature with the lowest mortality.
  mmttemp <- mmtrangetemp[which.min(crtemptemp$fit[which(crtemptemp$predvar%in%mmtrangetemp)])]
  
  plot(crtemp0)
  plot(crtemptemp)
  
  
  crtemp0 <- crossreduce(cbtemp,mod0,cen=mmt0,at=quantile(data$daily_mean,0.99))
  crtemptemp <- crossreduce(cbtemp,modtemp,cen=mmttemp,at=quantile(data$daily_mean,0.99))
  
  ## store different models
  crTN_cant_TN[[i]] <- crTN
  crtemp_cant[[i]] <- crtemp0
  mmt_cant[[i]] <- mmt0
  crtemptemp_c[[i]] <- crtemptemp
  mmttemp_c[[i]] <- mmttemp
}

##assign canton names
names(crTN_cant_TN) <- names(list_80_cant_TN)
names(crtemp_cant) <- names(list_80_cant_TN)
names(mmt_cant) <-names(list_80_cant_TN)
names(crtemptemp_c) <- names(list_80_cant_TN)
names(mmttemp_c) <- names(list_80_cant_TN)


## save 
saveRDS(crTN_cant_TN, paste0(homedir2,"/results/crTN_cant_TN.rds"))
saveRDS(crtemp_cant, paste0(homedir2, "/results/crtemp_cant.rds"))
saveRDS(mmt_cant, paste0(homedir2, "/results/mmt_cant.rds"))








#-------------------------------
#crossbasis WITHOUT daily Tmean:
#-------------------------------
crTN_cant_TN_without <- list()

for (i in 1:length(list_80_cant_TN)){
  data <- list_80_cant_TN[[i]]
  data$stratum <- with(data, factor(year),factor(month),factor(dow),factor(District_name))
  ind <- tapply(data$deaths,data$stratum,sum)[data$stratum]
  cbTN <- crossbasis(data$TN.min,lag=3,argvar=list(fun="strata",breaks=1),
                     arglag=list(fun="integer"),group=data$year)
  mod <- gnm(deaths ~ cbTN, eliminate = stratum,
             family=quasipoisson(),data=data, no.action="na.exclude",
             subset=ind>0)
  crTN <- crossreduce(cbTN,mod,cen=0,at=1)
  crTN_cant_TN_without[[i]] <- crTN
}
names(crTN_cant_TN_without) <- names(list_80_cant_TN)
crTN_cant_TN_without$Zürich
crTN_cant_TN_without







#----------------------------------------------
#creating a summary matrix with all information: 
#----------------------------------------------
results_cant_plot <- tibble(NA)
for (i in 1:length(list_80_cant_TN)){
  results_cant_plot[i,1] <- names(list_80_cant_TN)[i]
  results_cant_plot[i,2] <- round(crTN_cant_TN_without[[i]]$RRfit,2)
  results_cant_plot[i,3] <- round(crTN_cant_TN_without[[i]]$RRlow,2)
  results_cant_plot[i,4] <- round(crTN_cant_TN_without[[i]]$RRhigh,2)
  results_cant_plot[i,5] <- round(crTN_cant_TN[[i]]$RRfit,2)
  results_cant_plot[i,6] <- round(crTN_cant_TN[[i]]$RRlow,2)
  results_cant_plot[i,7] <- round(crTN_cant_TN[[i]]$RRhigh,2)
  results_cant_plot[i,8] <- sum(list_80_cant_TN[[i]]$TN.min)
  results_cant_plot[i,9] <- sum(list_80_cant_TN[[i]]$deaths)
}
colnames(results_cant_plot) <- c("Canton","RRfit", "RRlow", "RRhigh","RRfit tmean","RRlow tmean","RRhigh tmean","TN","deaths")

results_cant_plot<-results_cant_plot[-c(4,5,6,12,17),]
results_cant_plot$group <- factor(c(1,0,1,0,1,0,0,0,0,0,0,1,0,0))










#----------------------------------------------
# RESULTS ####################################
#----------------------------------------------


#plot results --> not adjusted for temperature
RR <- ggplot(results_cant_plot, aes(x=reorder(Canton,desc(Canton)),y=RRfit,ymin=RRlow, ymax=RRhigh))+
  geom_pointrange(color="grey55", size=0.8) +
  geom_hline(yintercept=1,lty=2,color="grey55") +
  coord_flip() +
  xlab("") + ylab("RR (95% CI)") + theme_light()+
  ggtitle("Relative Risk") +
  scale_color_manual(values=c("#999999", "darkorange")) + theme(
    axis.text.y = element_text(size=14,color="black"),
    #axis.text.y = element_blank(),
    axis.text.x = element_text(size=14,color="black"),
    axis.title = element_text(size=14,color="black"),
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA), 
    legend.background = element_rect(fill = "transparent", color = NA),
  )
RR


## plot results when adjusting for temperature 
RRadj <- ggplot(results_cant_plot, aes(x=reorder(Canton,desc(Canton)),y=c(`RRfit tmean`),ymin=c(`RRlow tmean`), ymax=c(`RRhigh tmean`)))+ #,color=group
  #geom_pointrange(show.legend=FALSE,size=0.8) + ylim(0.2,2.3) +
  geom_hline(yintercept=1,lty=2,color="grey55") +
  geom_point(size=3, shape=c(19),stroke=0.5)+
  geom_errorbar(aes(ymin=`RRlow tmean`,ymax=`RRhigh tmean`), width=0.2, cex=0.5)+
  coord_flip() +
  xlab("") + ylab("RR (95% CI)") + theme_light()+
  ggtitle("Realtive Risk of Mortality associated with Tropical Nights in Switzerland (1980-2018)")+
  scale_fill_manual(values="black")+
  scale_color_manual(values="black")+
  theme(
    title = element_text(size=9,color="black"),
    axis.text.y = element_text(size=9,color="black"),
    #axis.text.y = element_blank(),
    axis.text.x = element_text(size=9,color="black"),
    axis.title = element_text(size=9,color="black"),
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA), 
    legend.background = element_rect(fill = "transparent", color = NA),
  )
RRadj





