library(CreDating)
library(ape)
rm(list=setdiff(ls(),'ind'))
if (!exists('ind')) ind=1
set.seed(0)
dates=seq(2000,2010,0.1)
phy=simcoaltree(dates,alpha=10)
set.seed(ind)
rate=ind/100*20
obsphy=simobsphy(phy,rate=rate)
obsphy=unroot(obsphy)
res=credate(obsphy,dates,showProgress=T,nbIts=1000000)
save.image(sprintf('run%d.RData',ind))
