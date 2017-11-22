rm(list=ls())
dir='~/simuCreDating'

if (Sys.info()["nodename"]=='xavierdidelot') {
  #SERVER SIDE
  library(CreDating)
  xs=c()
  ys1=c()
  ys2=c()
  ys3=c()
  for (ind in 1:100) {
    tryCatch({
      myind=ind
      load(sprintf('%s/run%d.RData',dir,ind))
      ind=myind
      xs=c(xs,rate)
      v=res$record[,'rate']
      v=sort(v[floor(length(v)/2):length(v)])
      ys1=c(ys1,v[floor(length(v)*0.025)])
      ys2=c(ys2,v[floor(length(v)*0.500)])
      ys3=c(ys3,v[floor(length(v)*0.975)])
    },error=function(e) {},warning=function(w) {})
  }
  rm('res')
  save.image(sprintf('%s/all.RData',dir))
} else {

  #LAPTOP SIDE
  system(sprintf('scp ubuntu@137.205.69.116:%s/all.RData /tmp',dir))
  load('/tmp/all.RData')
  pdf('/tmp/fig.pdf',6,6)
  plot(xs,ys2,pch=16,xlab='Correct',ylab='Inferred')
  segments(xs,ys1,xs,ys3)
  dev.off()
  system('open /tmp/fig.pdf')
}
