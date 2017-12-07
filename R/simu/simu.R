library(CreDating)
library(ape)
rm(list=setdiff(ls(),'ind'))
if (!exists('ind')) ind=1
set.seed(ind)
dates=seq(2000,2010,0.1)

if (ind<=100) {
  #Code for simulation with varying rate
  alpha=10
  phy=simcoaltree(dates,alpha=alpha)
  rate=ind/100*20
} else {
  #Code for simulation with varying alpha
  alpha=(ind-100)/100*20
  phy=simcoaltree(dates,alpha=alpha)
  rate=10
}

obsphy=simobsphy(phy,rate=rate)
obsphy=unroot(obsphy)
res=credate(obsphy,dates,showProgress=T,nbIts=1000000)
save.image(sprintf('run%d.RData',ind))
