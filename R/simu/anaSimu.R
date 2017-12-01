rm(list=ls())
dir='~/simuCreDating'

if (Sys.info()["nodename"]=='xavierdidelot') {
  #SERVER SIDE
  library(CreDating)
  store=matrix(NA,100,8)
  for (ind in 1:100) {
    tryCatch({
      myind=ind
      load(sprintf('%s/run%d.RData',dir,ind))
      ind=myind
      store[ind,1]=rate
      v=res$record[,'rate']
      v=sort(v[floor(length(v)/2):length(v)])
      store[ind,c(2,3,4)]=v[floor(length(v)*c(0.025,0.5,0.975))]
      store[ind,5]=alpha
      v=res$record[,'alpha']
      v=sort(v[floor(length(v)/2):length(v)])
      store[ind,c(6,7,8)]=v[floor(length(v)*c(0.025,0.5,0.975))]
    },error=function(e) {},warning=function(w) {})
  }
  rm('res')
  save.image(sprintf('%s/all.RData',dir))
} else {

  #LAPTOP SIDE
  system(sprintf('scp ubuntu@137.205.69.116:%s/all.RData /tmp',dir))
  load('/tmp/all.RData')
  pdf('/tmp/fig.pdf',6,6)
  plot(store[,1],store[,3],pch=16,xlab='Correct',ylab='Inferred')
  segments(store[,1],store[,2],store[,1],store[,4])
  dev.off()
  system('open /tmp/fig.pdf')
}
