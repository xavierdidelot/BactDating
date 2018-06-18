library(BactDating)
library(ape)
rm(list=setdiff(ls(),'ind'))
if (!exists('ind')) ind=1
set.seed(ind)
dates=seq(2000,2010,0.1)

if (ind<=100) {
  #Code for simulation with varying rate
  alpha=5
  phy=simcoaltree(dates,alpha=alpha)
  mu=ind/100*10
  sigma=0
} else if (ind<=200) {
  #Code for simulation with varying alpha
  alpha=(ind-100)/100*10
  phy=simcoaltree(dates,alpha=alpha)
  mu=5
  sigma=0
} else {
  #Code for simulation with varying sigma
  alpha=5
  mu=5
  sigma=(ind-200)/100*10
  phy=simcoaltree(dates,alpha=alpha)
}

obsphy=simobsphy(phy,mu=mu,sigma=sigma,model='relaxedgamma')
obsphy=unroot(obsphy)
res=bactdate(obsphy,dates,showProgress=T,nbIts=1e3)#1e6
save.image(sprintf('run%d.RData',ind))
