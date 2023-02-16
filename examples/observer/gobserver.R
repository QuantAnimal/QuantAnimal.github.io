# Example of an state observer to estimate the perturbed component of the body weight trajectory
# Code performed by Rafael Muñoz-Tamayo (INRAE, France), rafael.munoz-tamayo@inrae.fr, 2021

library("deSolve")
library(ggplot2)
library(ggExtra)
library(gridExtra)


setwd("C:/Users/Utilisateur/seadrive_root/Masoomeh/Mes bibliothèques/2019_ResiliencePaper/2021_Robustness_model/Model_R/20211112_ExampleObserver_R_RMT")
setwd("D:/2019_ResiliencePaper/2021_Robustness_model/Model_R/2212_Casestudy3_RMT")
# Loading the measured noisy data [time (d) animal body weight (kg)]
ynoise  <- read.table("ynoise.txt");
DN  <- data.frame(ynoise);
colnames(DN)<- c("time","BW");
tN <- DN$time
yN <- DN$BW

# Loading the measured filtered data [time (d) animal body weight (kg)]
ysmooth <- read.table("ysmooth.txt");
DS  <- data.frame(ysmooth);
colnames(DS)<- c("time","BW");
tS <- DS$time
yS <- DS$BW

# Loading the perturbation factor used to generate the simulated data 

perturbfactor<- read.table("perturbfactor.txt");
DP  <- data.frame(perturbfactor);
colnames(DP)<- c("time","Fi");
tP <- DP$time
yP <- DP$Fi

# Dynamic model
tspan=tS; 

dy<- function(t, state,parameters) {  #ODE function observer
  with(as.list(c(state,parameters)), {
    
    # model parameters 
    # p1 = 0.05
    # p2 = 0.02
    
    # observer parameters
    #w1 =  4.0
    #w2 =  0.5
    
    #ysmooth <- read.table("ysmooth.txt");
    #DS  <- data.frame(ysmooth);
    #colnames(DS)<- c("time","BW");
    #tS <- DS$time
    #yS <- DS$BW
    
    ydata = approx(tS, yS, xout = t)$y # interpolation
    
    dy1=-ydata*y2 + ydata*p1*exp(-p2*t) + w1*(ydata-y1)
    dy2=-w2*y1*w1*(ydata-y1) 
    # return the result
    list(c(dy1,dy2))
  }) # end with(as.list ...
}

state <- c(y1 = 3,y2=0) #initial conditions
parameters <- c(p1 = 0.05, p2 = 0.02, w1 = 4.0, w2 = 0.5)
yout= ode(y = state, times = tspan, func = dy,parms = parameters) # solving the ODE

Dout  <- data.frame(yout);
colnames(Dout)<- c("time","BWobs","Fiobs");
tout <- Dout$time # Time
y1 <- Dout$BWobs  # Estimated body weight
y2 <- Dout$Fiobs  # Unknown perturbation function 



ggplot(DN, aes(tN, yN)) + geom_point(colour = 'black', size = 2.5)  + labs(title="A", x="Time (d)",  y= "Body weight (kg)") +
  geom_line(data = Dout, aes(tout,y1), size = 1, linetype = 1, colour = "blue") 

ggplot(DP, aes(tP, yP)) + geom_point(colour = 'black', size = 2.5)  + labs(title="B", x="Time (d)", y=expression(phi~ ", perturbation factor")) +
  geom_line(data = Dout, aes(tout,y2), size = 1, linetype = 1, colour = "blue" )


#tiff("PlotObserver.tiff", width = 14, height = 8.7, units = 'cm', res = 3000, compression = "lzw")

#plot1 = ggplot(DN, aes(tN, yN)) + geom_point(colour = 'black', size = 2.5)  + labs(title="A", x="Time (d)",  y= "Body weight (kg)") +
#  geom_line(data = Dout, aes(tout,y1), size = 1, linetype = 1, colour = "blue") 

#plot2 = ggplot(DP, aes(tP, yP)) + geom_point(colour = 'black', size = 2.5)  + labs(title="B", x="Time (d)", y=expression(phi~ ", perturbation factor")) +
#  geom_line(data = Dout, aes(tout,y2), size = 1, linetype = 1, colour = "blue" )

#grid.arrange(plot1,plot2,nrow=1, ncol=2); # 

"dev.off()

  
