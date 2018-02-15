#' Rescaling function using dual scale
#' @param tree Timed tree
#' @param minDate Date before which dates are to be compressed
#' @param squash Interval of time to give to all dates prior to minDate
#' @return Rescaled timed tree
#' @export
dualScale = function(tree,minDate=1990,squash=5)
{
  dates=allDates(tree)
  ma=max(dates)
  mi=min(dates)
  for (i in 1:length(dates))
    if (dates[i]<minDate)
      dates[i]=minDate-(minDate-dates[i])*squash/(minDate-mi)
  tree$root.time=minDate-squash
  for (i in 1:nrow(tree$edge)) tree$edge.length[i]=dates[tree$edge[i,2]]-dates[tree$edge[i,1]]
  return(tree)
}

#' Plot tree rescaled using dual scale
#' @param tree Timed tree
#' @param minDate Date before which dates are to be compressed
#' @param squash Interval of time to give to all dates prior to minDate
#' @param ... Remaining parameters are passed on to plot.phylo
#' @return Rescaled timed tree
#' @export
plotDualScale = function(tree,minDate=1990,squash=5,...)
{
  tree=dualScale(tree,minDate,squash)
  plot(tree,...)
  at=axTicks(1)
  labels=at+minDate-squash
  w=which(labels>=minDate)
  at=at[w]
  labels=labels[w]
  axis(1,at=at,labels=labels)
  return(tree)
}
