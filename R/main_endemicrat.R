############################################################################
## Compartmental and Epidemiologic modeling of SDAV data
## This code corresponds to the paper:
# Sarah Mullin, Brent Vander Wyk, Jennifer L Asher, Susan R Compton, Heather G Allore, Caroline J Zeiss, 
# Modeling pandemic to endemic patterns of SARS-CoV-2 transmission using parameters estimated from animal model data, 
# PNAS Nexus, Volume 1, Issue 3, July 2022, pgac096, https://doi.org/10.1093/pnasnexus/pgac096

############################################################################

# install packages
library(tidyverse)
library(dplyr)
library(stringr)
library(lubridate)
library(readr)
library(EpiModel)
library(ggplot2)
library(utils)

setwd("~/Desktop/SarsCOV2_transmission/R")
source('utilities.R')

# Set parameters: Rat
# beta_High
bdc<-0.859
# beta_Low
bic<-0.268
# delta
e.dur<-3.477477
# alpha: shed virus at a rate alpha
alpha<- 0.1019225
# phi: decay rate
phi<-c(0.2, 0.4, 0.6, 0.8)
# psi
i.dur<-2.477477
# m
m_plus<-145.1224
# p
p<-0.7027027
# e_1
e_1<-0.5345345
# e2
  e_2<-c(0.7, 0.8, 0.9)
# v_1
v_1<-c(0.0025,0.005,0.01)
# v_2
v_2<-35
# m_v
m_v<-m_plus

### Model 1
Model1 <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Population size
    num <- s.num + e.num + i.num + rplus.num 
    
    # Effective contact rate 
    lambda <- ((bdc * i.num) + (bic* f.num))/num
    
    dS <- -lambda*s.num + (1/m_plus)*rplus.num + (1-p)*(1/i.dur)*i.num
    dE <- lambda*s.num - (1/e.dur)*e.num
    dI <- (1/e.dur)*e.num - (1/i.dur)*i.num 
    dR_plus <- p*(1/i.dur)*i.num-((1/m_plus)*rplus.num)
    dF <- alpha * i.num - phi * f.num
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dE, dI, dR_plus, dF, 
           se.flow = lambda * s.num,
           ei.flow = (1/e.dur) * e.num,
           irplus.flow = p*(1/i.dur) * i.num,
           isminus.flow = (1-p)*(1/i.dur) * i.num,
           rpluss.flow = (1/m_plus)*rplus.num,
           if.flow = alpha * i.num),
         num = num,
         i.prev = i.num / num,
         ei.prev = (e.num + i.num)/num)
  })}

### Model 1 with Table 1 parameters (rat)
param <- param.dcm(bdc = as.numeric(bdc), bic=as.numeric(bic), e.dur = e.dur, i.dur = i.dur, m_plus=m_plus, p=p,phi = phi, alpha=alpha)
init <- init.dcm(s.num = 1000, e.num = 0, i.num = 19, rplus.num = 0, f.num = 0,
                 se.flow = 0, ei.flow = 0, irplus.flow = 0,isminus.flow=0,if.flow=19, rpluss.flow=0)
# number of time steps is 1825 days with time at each day
control <- control.dcm(nsteps = 1825, dt = 1, new.mod =Model1)

# model
mod <- dcm(param, init, control)
cbp1 <- c( "#E69F00", "#56B4E9", "#009E73",
           "#0072B2", "#D55E00", "#CC79A7","#F0E442", "#999999")
par(mfrow = c(1, 2))

# find medians of all compartment curves
mod2<-addmediansci_model1(mod)
par(mar=c(7,5,5,2),cex.main=1.8, cex.axis=1, cex.lab=1.6)

# plot of medians of all possible parameter combinations
# 365 days
plot.dcm2(mod2, y = c("s.num", "e.num","i.num", "rplus.num", "f.num"),
          run = 1, cex=10, xlim = c(0, 365), grid = TRUE, legend='n', col=cbp1, leg.cex=1,ylim=c(0,1000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)'), main='365 Days (Rat)', xlab='Days',cex.lab=5, cex.axis=5, cex.main=5)
# 1825 days
plot.dcm2(mod2, y = c("s.num", "e.num","i.num", "rplus.num", "f.num"),
          run = 1, xlim = c(0, 1825), grid = TRUE, legend='n', col=cbp1, leg.cex=1.2,ylim=c(0,1000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)'), main='1825 Days (Rat)')
# choose only one run
plot(mod, y = c("s.num", "e.num","i.num", "rplus.num", "f.num"),
     run = 1, xlim = c(400,1000), grid = TRUE, legend='f', col=cbp1, leg.cex=1.2,ylim=c(0,1000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)'), main='1825 Days')
# Plot Prevalence/Incidence
par(mfrow = c(1,2))
plot(mod, y = "i.num",alpha = 1, main = "Model 1 Disease Prevalence", legend='full', leg.cex=1.1, xlab='Days', ylab='Individuals')
plot(mod, y = "se.flow", col = "Greens", alpha = 0.8, leg.cex=1.1,main = "Model 1 Disease Incidence", legend='full',  xlab='Days', ylab='Individuals')

### Model 2

Model2 <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    # Population size
    num <- s.num + e.num + i.num + rplus.num + v1.num + v2.num
    
    # Effective contact rate
    lambda <- ((bdc * i.num) + (bic* f.num))/num
    
    dS <- -lambda*s.num + (1/m_plus)*rplus.num + (1-p)*(1/i.dur)*i.num-v_1*s.num+(1/m_v)*v2.num
    dE <- lambda*(s.num+(1-e_1)*v1.num+(1-e_2)*v2.num) - (1/e.dur)*e.num
    dI <- (1/e.dur)*e.num - (1/i.dur)*i.num 
    dR_plus <- p*(1/i.dur)*i.num-((1/m_plus)*rplus.num)
    dF <- alpha * i.num - phi * f.num
    dV1<- v_1*s.num-(1/v_2)*v1.num -lambda*(1-e_1)*v1.num
    dV2<- (1/v_2)*v1.num -lambda*(1-e_2)*v2.num-(1/m_v)*v2.num
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dE, dI, dR_plus, dF, dV1,dV2,
           se.flow = lambda * s.num,
           ei.flow = (1/e.dur) * e.num,
           irplus.flow = p*(1/i.dur) * i.num,
           isminus.flow = (1-p)*(1/i.dur) * i.num,
           rpluss.flow = (1/m_plus)*rplus.num,
           if.flow = alpha * i.num,
           v1v2.flow = (1/v_2)*v1.num,
           sv1.flow = v_1*s.num,
           v2e.flow = lambda*(1-e_2)*v2.num,
           v1e.flow = lambda*(1-e_1)*v1.num,
           v2s.flow = (1/m_v)*v2.num
    ),
    num = num,
    i.prev = i.num / num,
    ei.prev = (e.num + i.num)/num)
  })}

### Model 2 with Table 1 parameters (rat)
b=expand.grid(bdc = bdc, bic=bic, e.dur = e.dur, i.dur = i.dur, m_plus=m_plus,p=0.7027027,phi=0.2,alpha=0.1019225,v_1=v_1, v_2=v_2, e_1=e_1, e_2=e_2, m_v=m_v )

param <- param.dcm(bdc = b$bdc, bic=b$bic, e.dur = b$e.dur, i.dur = b$i.dur, m_plus=b$m_plus, p=b$p,phi = b$phi, alpha=b$alpha,v_1=b$v_1, v_2=b$v_2 ,e_1=b$e_1, e_2=b$e_2 , m_v=b$m_v)
init <- init.dcm(s.num = 1000, e.num = 0, i.num = 19, rplus.num = 0, f.num = 0,v1.num=0, v2.num=0,
                 se.flow = 0, ei.flow = 0, irplus.flow = 0,isminus.flow=0,if.flow=19, rpluss.flow=0, v1v2.flow=0, sv1.flow=0, v2e.flow=0, v1e.flow=0, v2s.flow=0)
control <- control.dcm(nsteps = 1825, dt = 1, new.mod =Model2)
mod <- dcm(param, init, control)

# find medians of all compartment curves
par(mfrow = c(1, 2))
mod2<-addmediansci_model2(mod)
par(mfrow = c(1,2), mar = c(7,7,7,2), cex.lab=1.6,cex.axis=1.2, cex.main=1.8)

# plot of 325 days
plot.dcm2(mod2, y = c("s.num", "e.num","i.num", "rplus.num", "f.num","v1.num", "v2.num"),
          run = 1, cex=10, xlim = c(0, 365), grid = TRUE, legend='n', col=cbp1, leg.cex=1,ylim=c(0,1000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)','1-dose vaccinated', '2-dose vaccinated'), main='365 Days (Rat)', xlab='Days',cex.lab=5, cex.axis=5, cex.main=5)
# plot of 1825 days
plot.dcm2(mod2, y = c("s.num", "e.num","i.num", "rplus.num", "f.num","v1.num", "v2.num"),
          run = 1, xlim = c(0, 1825), grid = TRUE, legend='f', col=cbp1, leg.cex=1.2,ylim=c(0,1000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)','1-dose vaccinated', '2-dose vaccinated'), main='1825 Days (Rat)')

# Human extrapolated

### Model 1 with Table 2 parameters (human)

a=expand.grid(bdc = c(0.05,0.1, 0.2, 0.4, 0.8), bic=c(0.01,0.125, 0.25), e.dur = c(3,4,5,6), i.dur = c(10,17), m_plus=c(90, 145.12, 180, 365),p=0.7027027,phi=c(0.2),alpha=0.1019225 )
# remove some implausible parameters (e.g. beta_High< beta_Low)
a<-a[!((a$bdc==0.05 & a$i.dur==10)| (a$bdc==0.1 & a$i.dur==10)|(a$bdc==0.1 & a$i.dur==10)|(a$bdc==0.05 & a$bic==0.125)|(a$bdc==0.1 & a$bic==0.125)|(a$bdc==0.1 & a$bic==0.25)| (a$bdc==0.05 & a$bic==0.25) ),]
# run model
param <- param.dcm(bdc = a$bdc, bic=a$bic, e.dur = a$e.dur, i.dur = a$i.dur, m_plus=a$m_plus, p=a$p,phi = a$phi, alpha=a$alpha)
init <- init.dcm(s.num = 100000, e.num = 0, i.num = 19, rplus.num = 0, f.num = 0,
                 se.flow = 0, ei.flow = 0, irplus.flow = 0,isminus.flow=0,if.flow=19, rpluss.flow=0)
control <- control.dcm(nsteps = 1825, dt = 1, new.mod =Model1)
mod_table2 <- dcm(param, init, control)

mod3<-addmediansci_model1(mod_table2)
# set display options for Rstudio
par(mfrow = c(1,2))
options(scipen=6)

# plot
plot.dcm2(mod3, y = c("s.num", "e.num","i.num", "rplus.num", "f.num"),
          run = 1, cex=10, xlim = c(0, 365), grid = TRUE, legend='n', col=cbp1, leg.cex=1,ylim=c(0,100000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)'), main='365 Days (Human)', xlab='Days',cex.lab=5, cex.axis=5, cex.main=5)
plot.dcm2(mod3, y = c("s.num", "e.num","i.num", "rplus.num", "f.num"),
          run = 1, xlim = c(0, 1835), grid = TRUE, legend='n', col=cbp1, leg.cex=1.2,ylim=c(0,100000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)'), main='1825 Days (Human)', xlab='Days',cex.lab=5, cex.axis=5, cex.main=5)
par(mfrow = c(1,1))
plot(mod_table2, y = "i.num",alpha = 1, main = "Model 1 Disease Prevalence (Human)", legend='n', leg.cex=1.1, xlab='Days', ylab='Individuals', lwd=0.4)


### Model 2 with Table 2 parameters (human)

# note: this is Model 2 with v_1 set to 0 for the first year
Model2_mv_step <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    if (t<=365) {v_1=0}
    else  {v_1=v_vac}
    # Population size
    num <- s.num + e.num + i.num + rplus.num + v1.num + v2.num
    
    # Effective contact rate 
    lambda <- ((bdc * i.num) + (bic* f.num))/num
    
    dS <- -lambda*s.num + (1/m_plus)*rplus.num + (1-p)*(1/i.dur)*i.num-v_1*s.num+(1/m_v)*v2.num
    dE <- lambda*(s.num+(1-e_1)*v1.num+(1-e_2)*v2.num) - (1/e.dur)*e.num
    dI <- (1/e.dur)*e.num - (1/i.dur)*i.num 
    dR_plus <- p*(1/i.dur)*i.num-((1/m_plus)*rplus.num)
    dF <- alpha * i.num - phi * f.num
    dV1<- v_1*s.num-(1/v_2)*v1.num -lambda*(1-e_1)*v1.num
    dV2<- (1/v_2)*v1.num -lambda*(1-e_2)*v2.num-(1/m_v)*v2.num
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dE, dI, dR_plus, dF, dV1,dV2,
           se.flow = lambda * s.num,
           ei.flow = (1/e.dur) * e.num,
           irplus.flow = p*(1/i.dur) * i.num,
           isminus.flow = (1-p)*(1/i.dur) * i.num,
           rpluss.flow = (1/m_plus)*rplus.num,
           if.flow = alpha * i.num,
           v1v2.flow = (1/v_2)*v1.num,
           sv1.flow = v_1*s.num,
           v2e.flow = lambda*(1-e_2)*v2.num,
           v1e.flow = lambda*(1-e_1)*v1.num,
           v2s.flow = (1/m_v)*v2.num
    ),
    num = num,
    i.prev = i.num / num,
    ei.prev = (e.num + i.num)/num)
  })}

# define parameters
c=expand.grid(bdc = c(0.05,0.1, 0.2, 0.4, 0.8), bic=c(0.01,0.125, 0.25), e.dur = c(3,4,5,6), i.dur = c(10,17), m_plus=c(90, 145.12, 180, 365),p=0.7027027,phi=c(0.2),alpha=0.1019225, v_vac=v_1, v_2=v_2, e_1=e_1, e_2=e_2, m_v=c(90, 145.12, 365)  )
c<-c[!((c$bdc==0.05 & c$i.dur==10)| (c$bdc==0.1 & c$i.dur==10)|(c$bdc==0.1 & c$i.dur==10)|(c$bdc==0.05 & c$bic==0.125)|(c$bdc==0.1 & c$bic==0.125)|(c$bdc==0.1 & c$bic==0.25)| (c$bdc==0.05 & c$bic==0.25) ),]

param1 <- param.dcm(bdc = c$bdc, bic=c$bic, e.dur = c$e.dur, i.dur = c$i.dur, m_plus=c$m_plus, p=c$p,phi = c$phi, alpha=c$alpha,v_vac=c$v_vac, v_2=c$v_2 ,e_1=c$e_1, e_2=c$e_2 , m_v=c$m_v)
init1 <- init.dcm(s.num = 100000, e.num = 0, i.num = 19, rplus.num = 0, f.num = 0,v1.num=0, v2.num=0,
                  se.flow = 0, ei.flow = 0, irplus.flow = 0,isminus.flow=0,if.flow=19, rpluss.flow=0, v1v2.flow=0, sv1.flow=0, v2e.flow=0, v1e.flow=0, v2s.flow=0)
control1 <- control.dcm(nsteps = 1825, dt = 1, new.mod =Model2_mv_step)
mod_table_2 <- dcm(param1, init1, control1)

# add medians to compartment curves
mod3<-addmediansci_model2(mod_table_2)

# plot compartment curves
plot.dcm2(mod3, y = c("s.num", "e.num","i.num", "rplus.num", "f.num","v1.num", "v2.num"),
          run = 1, cex=10, xlim = c(0, 365), grid = TRUE, legend='n', col=cbp1, leg.cex=1,ylim=c(0,100000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)','1-dose vaccinated', '2-dose vaccinated'), main='365 Days', xlab='Days',cex.lab=8, cex.axis=8, cex.main=8)
plot.dcm2(mod3, y = c("s.num", "e.num","i.num", "rplus.num", "f.num","v1.num", "v2.num"),
          run = 1, xlim = c(0, 1835), grid = TRUE, legend='n', col=cbp1, leg.cex=1.2,ylim=c(0,100000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)','1-dose vaccinated', '2-dose vaccinated'), main='Wild Type (Human)')
par(mfrow = c(1,2))

# plot incidence/prevalence
plot(mod_table_2, y = "i.num",alpha = 1, main = "Model 2 Disease Prevalence", legend='full', leg.cex=1.1, xlab='Days', ylab='Individuals')
plot(mod_table_2, y = "se.flow", col = "Greens", alpha = 0.8, leg.cex=1.1,main = "Model 2 Disease Incidence", legend='full',  xlab='Days', ylab='Individuals')

### introducing variants with Model 2
Model2_mv_step_variant <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    if (t<=365) {v_1=0}
    else  {v_1=v_vac}
    if (t<=180) {bdc=beta+1*(0.05)
    m_plus=90}
    else if (t>180 & t<=360) {bdc=beta+2*(0.05)
    m_plus=90}
    else if (t>360 & t<=540) {bdc=beta+3*(0.05)
    m_plus=90}
    else if (t>540 & t<=720) {bdc=beta+4*(0.05)
    m_plus=90}
    else if (t>720 & t<=900) {bdc=beta+5*(0.05)
    m_plus=90}
    else if (t>900 & t<=1080) {bdc=beta+6*(0.05)
    m_plus=90} 
    else if (t>1080 & t<1260) {bdc=beta+7*(0.05)
    m_plus=90}
    else if (t>1260 & t<1440) {bdc=beta+8*(0.05)
    m_plus=90}
    else if (t>1440 & t<1620) {bdc=beta+9*(0.05)
    m_plus=90}
    else{bdc=beta+10*(0.05) 
    m_plus=90}
    bic=bdc/3.2
    
    #print(e_1)
    # Population size
    num <- s.num + e.num + i.num + rplus.num + v1.num + v2.num
    
    # Effective contact rate 
    lambda <- ((bdc * i.num) + (bic* f.num))/num
    
    dS <- -lambda*s.num + (1/m_plus)*rplus.num + (1-p)*(1/i.dur)*i.num-v_1*s.num+(1/m_v)*v2.num
    dE <- lambda*(s.num+(1-e_1)*v1.num+(1-e_2)*v2.num) - (1/e.dur)*e.num
    dI <- (1/e.dur)*e.num - (1/i.dur)*i.num 
    dR_plus <- p*(1/i.dur)*i.num-((1/m_plus)*rplus.num)
    dF <- alpha * i.num - phi * f.num
    dV1<- v_1*s.num-(1/v_2)*v1.num -lambda*(1-e_1)*v1.num
    dV2<- (1/v_2)*v1.num -lambda*(1-e_2)*v2.num-(1/m_v)*v2.num
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dE, dI, dR_plus, dF, dV1,dV2,
           se.flow = lambda * s.num,
           ei.flow = (1/e.dur) * e.num,
           irplus.flow = p*(1/i.dur) * i.num,
           isminus.flow = (1-p)*(1/i.dur) * i.num,
           rpluss.flow = (1/m_plus)*rplus.num,
           if.flow = alpha * i.num,
           v1v2.flow = (1/v_2)*v1.num,
           sv1.flow = v_1*s.num,
           v2e.flow = lambda*(1-e_2)*v2.num,
           v1e.flow = lambda*(1-e_1)*v1.num,
           v2s.flow = (1/m_v)*v2.num
    ),
    num = num,
    i.prev = i.num / num,
    ei.prev = (e.num + i.num)/num)
  })}

# get parameter options
d=expand.grid(beta=0.2156559, e.dur =4.5 , i.dur = 10, p=p,phi = 0.2, alpha=alpha, v_vac=v_1, v_2=v_2, e_1=e_1, e_2=e_2, m_v=c(90, 145.12, 365)  )

param <- param.dcm(beta = d$beta, e.dur = d$e.dur, i.dur = d$i.dur, p=d$p,phi = d$phi, alpha=d$alpha,v_vac=d$v_vac, v_2=d$v_2 ,e_1=d$e_1, e_2=d$e_2 , m_v=d$m_v)
init <- init.dcm(s.num = 100000, e.num = 0, i.num = 19, rplus.num = 0, f.num = 0,v1.num=0, v2.num=0,
                 se.flow = 0, ei.flow = 0, irplus.flow = 0,isminus.flow=0,if.flow=19, rpluss.flow=0, v1v2.flow=0, sv1.flow=0, v2e.flow=0, v1e.flow=0, v2s.flow=0)
control <- control.dcm(nsteps = 1825, dt = 1, new.mod =Model2_mv_step_variant)

# run model
mod_table22 <- dcm(param, init, control)

# add medians to compartment curves 
mod32<-addmediansci_model2(mod_table22)

# plot compartmental curves for adding variants in a step-wise fashion
plot.dcm2(mod32, y = c("s.num", "e.num","i.num", "rplus.num", "f.num","v1.num", "v2.num"),
          run = 1, xlim = c(0, 1825), grid = TRUE, legend='n', col=cbp1, leg.cex=1.2,ylim=c(0,100000),ylab='Individuals', y_names=c('susceptible', 'exposed', 'infected', 'seropositive recovered', 'low-risk (fomite)','1-dose vaccinated', '2-dose vaccinated'), main='Variants')

