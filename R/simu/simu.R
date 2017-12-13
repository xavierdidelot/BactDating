library(CreDating)
library(ape)
rm(list=setdiff(ls(),'ind'))
if (!exists('ind')) ind=1
set.seed(ind)
dates=seq(2000,2010,0.1)

if (ind<=100) {
  #Code for simulation with varying rate
  alpha=5
  phy=simcoaltree(dates,alpha=alpha)
  rate=ind/100*10
} else {
  #Code for simulation with varying alpha
  alpha=(ind-100)/100*10
  phy=simcoaltree(dates,alpha=alpha)
  rate=5
}

obsphy=simobsphy(phy,rate=rate)
obsphy=unroot(obsphy)
res=credate(obsphy,dates,showProgress=T,nbIts=1000000)
save.image(sprintf('run%d.RData',ind))
