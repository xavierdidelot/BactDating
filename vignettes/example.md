Example of using CreDating
================
Xavier Didelot
2017-09-23

Initialisation
--------------

``` r
library(CreDating)
library(ape)
set.seed(0)
```

Data
----

We consider a random timed tree with 20 leaves, with the root in 2000.

``` r
nsam <- 20
phy <- rtree(nsam, br = rexp)
plot(phy)
axisPhylo(root.time = 2000,backward = F)
```

![](example_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-14-1.png)

We can compute the age of the leaves.

``` r
dates=rep(2000,nsam)
  for (i in 1:nsam) {
    w=i
    while (1) {
      r=which(phy$edge[,2]==w)
      if (length(r)==0) break
      dates[i]=dates[i]+phy$edge.length[r]
      w=phy$edge[r,1]
    }
  }
```

On each branch we observe a number of substitutions which is distributed *P**o**i**s**s**o**n*(*r**l*) where *l* is the branch length and *r* = 10 per year is the real substitution rate.

``` r
obsphy=phy
obsphy$edge.length=unlist(lapply(obsphy$edge.length*20,function (x) rpois(1,x)))
plot(obsphy)
axisPhylo(backward = F)
```

![](example_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-16-1.png)

``` r
res=credating(obsphy,dates,nbIts = 1000,rate=10,useCoalPrior = T,updateRate = 1)
plot(res$tree)
axisPhylo(backward = F)
```

![](example_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-17-1.png)

``` r
  nc=ncol(res$record)
  par(mfrow=c(2,2))
  plot(res$record[,nc-1],main='Likelihood',type='l',xlab='Sampled iterations',ylab='')
  plot(res$record[,nc-2],main='Date of root',type='l',xlab='Sampled iterations',ylab='')
  plot(res$record[,nc],main='Rate',type='l',xlab='Sampled iterations',ylab='')
  par(mfrow=c(1,1))
```

![](example_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-18-1.png)
