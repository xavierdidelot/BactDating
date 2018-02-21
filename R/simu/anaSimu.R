rm(list=ls())
dir='~/simuCreDating'

if (Sys.info()["nodename"]=='server1b') {
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
#  system(sprintf('scp ubuntu@137.205.69.104:%s/all.RData /tmp',dir))
  load('/tmp/all.RData')
  allstore=store
  store=allstore[1:100,]
  pdf('/tmp/fig.pdf',8.5,7)
  par(fig=c(0.15,0.45,0.1,0.5), new=F,mar=c(0,0,0,0),xpd=T)
  plot(store[,1],store[,3],xlim=c(0,10),ylim=c(0,11),pch=16,xlab='',ylab='',axes=F)
  box()
  axis(1,at=c(0,2,4,6,8,10))
  axis(2,at=c(0,2,4,6,8,10))
  segments(store[,1],store[,2],store[,1],store[,4])
  lines(c(0,10),c(0,10))
  par(xpd=NA)
  text(5,-2.5,expression(paste('Correct rate')))
  text(-2.5,5,expression(paste('Inferred rate')),srt=90)

  par(fig=c(0.15,0.45,0.5,0.9), new=TRUE,mar=c(0,0,0,0),xpd=T)
  plot(store[,1],store[,7],xlim=c(0,10),ylim=c(0,11),pch=16,xlab='',ylab='',axes=F)
  box()
  axis(1,at=c(0,2,4,6,8,10),labels =NA)
  axis(2,at=c(0,2,4,6,8,10))
  segments(store[,1],store[,6],store[,1],store[,8])
  lines(c(0,10),c(5,5))
  par(xpd=NA)
  text(-2.5,5,expression(paste('Inferred alpha')),srt=90)
  text(-2.5,12,'A',cex=1.5)
  #dev.off()

  store=allstore[101:200,]
  #pdf('/tmp/fig2.pdf',4.5,6)

  par(fig=c(0.6,0.9,0.1,0.5), new=T,mar=c(0,0,0,0),xpd=T)
  plot(store[,5],store[,3],xlim=c(0,10),ylim=c(0,11),pch=16,xlab='',ylab='',axes=F)
  box()
  axis(1,at=c(0,2,4,6,8,10))
  axis(2,at=c(0,2,4,6,8,10))
  segments(store[,5],store[,2],store[,5],store[,4])
  lines(c(0,10),c(5,5))
  par(xpd=NA)
  text(-2.5,5,expression(paste('Inferred rate')),srt=90)
  text(5,-2.5,expression(paste('Correct alpha')))

  par(fig=c(0.6,0.9,0.5,0.9), new=TRUE,mar=c(0,0,0,0),xpd=T)
  plot(store[,5],store[,7],xlim=c(0,10),ylim=c(0,11),pch=16,xlab='',ylab='',axes=F)
  box()
  axis(1,at=c(0,2,4,6,8,10),labels =NA)
  axis(2,at=c(0,2,4,6,8,10))
  segments(store[,5],store[,6],store[,5],store[,8])
  lines(c(0,10),c(0,10))
  par(xpd=NA)
  text(-2.5,5,expression(paste('Inferred alpha')),srt=90)
  text(-2.5,12,'B',cex=1.5)

  dev.off()
  system('open /tmp/fig.pdf')
}
