library(CreDating)
library(ape)
rm(list=setdiff(ls(),'ind'))
if (!exists('ind')) ind=1
set.seed(ind)
dates=seq(2000,2010,0.1)
phy=simcoaltree(dates,neg=10)
rate=runif(1)*20
obsphy=simobsphy(phy,rate=rate)
obsphy=unroot(obsphy)
res=credate(obsphy,dates,showProgress=T,nbIts=100000)
save.image(sprintf('run%d.RData',ind))
