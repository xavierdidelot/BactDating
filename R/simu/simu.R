library(CreDating)
library(ape)
rm(list=setdiff(ls(),'ind'))
if (!exists('ind')) ind=1

if (T) {
#Code for simulation run6
set.seed(0)
dates=seq(2000,2010,0.1)
alpha=10
phy=simcoaltree(dates,alpha=alpha)
set.seed(ind)
rate=ind/100*20
obsphy=simobsphy(phy,rate=rate)
obsphy=unroot(obsphy)
res=credate(obsphy,dates,showProgress=T,nbIts=1000000)
save.image(sprintf('run%d.RData',ind))
} else {
#Code for simulation run7
set.seed(ind)
dates=seq(2000,2010,0.1)
alpha=ind/100*20
phy=simcoaltree(dates,alpha=alpha)
rate=10
obsphy=simobsphy(phy,rate=rate)
obsphy=unroot(obsphy)
res=credate(obsphy,dates,showProgress=T,nbIts=1000000)
save.image(sprintf('run%d.RData',ind))
}
