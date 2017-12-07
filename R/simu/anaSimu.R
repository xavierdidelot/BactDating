rm(list=ls())
dir='~/simuCreDating'

if (Sys.info()["nodename"]=='xavierdidelot') {
  #SERVER SIDE
  library(CreDating)
  store=matrix(NA,200,8)
  for (ind in 1:200) {
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
  allstore=store
  store=allstore[1:100,]
  pdf('/tmp/fig.pdf',4.5,6)
  par(fig=c(0.2,0.9,0.2,0.7), new=F,mar=c(0,0,0,0))
  plot(store[,1],store[,3],xlim=c(0,20),ylim=c(0,20),pch=16,xlab='',ylab='',axes=F)
  box()
  axis(1,at=c(0,5,10,15,20))
  axis(2,at=c(0,5,10,15,20))
  segments(store[,1],store[,2],store[,1],store[,4])
  lines(c(0,20),c(0,20))
  par(xpd=NA)
  text(10,-5,expression(paste('Correct rate')))
  text(-5,10,expression(paste('Inferred rate')),srt=90)

  par(fig=c(0.2,0.9,0.7,0.9), new=TRUE,mar=c(0,0,0,0))
  plot(store[,1],store[,7],xlim=c(0,20),ylim=c(6,14),pch=16,xlab='',ylab='',axes=F)
  box()
  axis(1,at=c(0,5,10,15,20),labels =NA)
  axis(2,at=c(8,12))
  axis(2,at=c(10))
  segments(store[,1],store[,6],store[,1],store[,8])
  lines(c(0,20),c(10,10))
  par(xpd=NA)
  text(-5,10,expression(paste('Inferred alpha')),srt=90)
  dev.off()

  store=allstore[101:200,]
  pdf('/tmp/fig2.pdf',4.5,6)
  par(fig=c(0.2,0.9,0.2,0.7), new=F,mar=c(0,0,0,0))
  plot(store[,5],store[,7],xlim=c(0,20),ylim=c(0,20),pch=16,xlab='',ylab='',axes=F)
  box()
  axis(1,at=c(0,5,10,15,20))
  axis(2,at=c(0,5,10,15,20))
  segments(store[,5],store[,6],store[,5],store[,8])
  lines(c(0,20),c(0,20))
  par(xpd=NA)
  text(10,-5,expression(paste('Correct alpha')))
  text(-5,10,expression(paste('Inferred alpha')),srt=90)

  par(fig=c(0.2,0.9,0.7,0.9), new=TRUE,mar=c(0,0,0,0))
  plot(store[,5],store[,3],xlim=c(0,20),ylim=c(6,14),pch=16,xlab='',ylab='',axes=F)
  box()
  axis(1,at=c(0,5,10,15,20),labels =NA)
  axis(2,at=c(8,12))
  axis(2,at=c(10))
  segments(store[,5],store[,2],store[,5],store[,4])
  lines(c(0,20),c(10,10))
  par(xpd=NA)
  text(-5,10,expression(paste('Inferred rate')),srt=90)
  dev.off()
  system('open /tmp/fig2.pdf')
}
