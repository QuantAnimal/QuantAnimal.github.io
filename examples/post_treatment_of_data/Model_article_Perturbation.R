######################################
#Masoomeh Taghipoor - MoSAR INRAE - Version March 2020
#Charaterization of animal response to perturbations: how mathematical models adapted recent progress in  monitoring technologies?
#######################################

rm(list=ls(all=TRUE)) #To remove the hisory 
dev.off() #To close all graphs
library("deSolve")
library(ggplot2)
setwd("C:/Users/Utilisateur/seadrive_root/Masoomeh/Mes bibliothèques/2019_ResiliencePaper/2021_Robustness_model/Model_R")
load(file = "DataPert.Rdata")
########################################
#Example 1 

  y=function(t){ #Ideal trajectory
    yl=ymax-ymax*exp(0.1*(-t))
    return(yl)
  }

ymax=6 #plateau
tspan=seq(0,130,1)


plot(tspan,y(tspan), type="l")
points(tspan,data$obs)

  dalpha <- function(t, state, parameters) {  #ODE funccion of perturbations
    with(as.list(c(state, parameters)), {
      tb1=10
      te1=25
      tb2=75
      te2=85
      # rate of change 
      da1=k11*(-a1)*(t<te1 &t>tb1)+ k21*(1-a1)*(t>te1&t<tb2)
      da2=k12*(-a2)*(t<te2 &t>tb2)+ k22*(1-a2)*(t>te2)
          # return the result
      list(c(da1,da2))
    }) # end with(as.list ...
  }
  
    state <- c(a1 = 1,a2=1) #init var
    parameters <- c(k11 = .05, k12 = .1, k21=0.15,k22=0.05)
    alpha= ode(y = state, times = tspan, func = dalpha, parms = parameters) #alpha1 et 2: 1st and 2nd pertubations

alpha=as.data.frame(alpha)

 data.test=as.data.frame(cbind(alpha,y(tspan),data$obs,(alpha[2]+alpha[3])/2*y(tspan)))
 colnames(data.test)=c("time","alpha1","alpha2","perf","obs","pert")
# alpha1 et 2: 1st and 2nd pertubations
# perf: animal ideal trajectory of performance
# pert: actual perturbed performance 

col=c(1,2,3,4)

ggplot(data = data.test, aes(x = tspan)) + 
  geom_line(aes( y = perf, color = col[1]))+  
  geom_line(aes( y = pert, color = col[2]))+
  geom_point(aes( y = obs, color = col[3]))

#==================================================
#                      Example 2
#==================================================
# install.packages("fda")
rm(list=ls(all=TRUE)) #To remove the hisory 
dev.off() #To close all graphs
library("fda")
library("graphics")

load("data_FDA.Rdata")

plot(data$time,data$obs)
write.table(data$obs, file="data-FDA", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
tspan=seq(0,130,1)

order = 6 
pen.what = 4
lambda = 10^4
range = c(data$time[1],data$time[length(data$time)])
knot.vec = seq(1, length(data$time),2)
break.vec = data$time[knot.vec]

# B-spline regression
    basis.func = create.bspline.basis(rangeval = range, norder = order, breaks = break.vec)
    param.func = fdPar(fdobj = basis.func, Lfdobj = pen.what,lambda = lambda)
    smooth.res = smooth.basis(argvals = data$time,data$obs,param.func)

res.fda = smooth.res$fd
eval.days =data$time
  
    res.estim   = eval.fd(eval.days,res.fda)   #
    vel.estim = eval.fd(eval.days,res.fda,1) #
    acc.estim = eval.fd(eval.days,res.fda,2) # 

data$res = res.estim
data$deriv=vel.estim
data$acc=acc.estim

    par(mfrow=c(2,1))
    plot(eval.days,res.estim,type="l" ,
    ylab=c("Indicator of performance"), xlab=c("Time"), lwd=2,cex.lab=1.5)
    points(data$time,data$obs, type="p")
    points(data$time, data$perf,type='l', col="blue", lty=1)
    plot(eval.days,vel.estim,type='l',cex.lab=1.5, ylab=c("First derivative"), xlab=c("Time"), lwd=2)
    abline(0,0, col="blue", lty=2)
# abline(v=day.min,lty=3)
#============


