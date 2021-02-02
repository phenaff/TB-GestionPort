# Utility functions for plots

# log-log plots

drawlogaxis <- function(side,range)
{
  par(tck=0.02)
  #	d <- log(range,10)
  d <- range
  mlog <- floor(min(d))
  Mlog <- ceiling(max(d))
  SeqLog <- c(mlog:Mlog)
  Nlog <- (Mlog-mlog)+1
  axis(side,at=SeqLog,labels=10^SeqLog)
  ats <- log(seq(from=2,to=9,by=1),10)
  mod <- NULL
  for(i in SeqLog)
  {
    mod <- c(mod,rep(i,length(ats)))
  }
  ats <- rep(ats,Nlog)
  ats <- ats+mod
  par(tck=0.02/3)
  axis(side,at=ats,labels=NA)
}

logplot <- function(x,y,log='xy',...,forceylim=c(0,0),forcexlim=c(0,0))
{
  par(tck=0.02)
  xlg <- FALSE
  ylg <- FALSE
  if('x'%in%strsplit(log,'')[[1]]){x <- log(x,10);xlg=TRUE}
  if('y'%in%strsplit(log,'')[[1]]){y <- log(y,10);ylg=TRUE}
  yl <- ifelse(forceylim==c(0,0),range(y),forceylim)
  xl <- ifelse(forcexlim==c(0,0),range(x),forcexlim)
  plot(x,y,...,axes=FALSE,ylim=yl,xlim=xl)
  if(xlg){drawlogaxis(1,xl)}else{axis(1,at=pretty(xl),labels=pretty(xl))}
  if(ylg){drawlogaxis(2,yl)}else{axis(2,at=pretty(yl),labels=pretty(yl))}
  box()
}

add.lines <- function(x,y,log='xy',...)
{
  xlg <- FALSE
  ylg <- FALSE
  if('x'%in%strsplit(log,'')[[1]]){x <- log(x,10);xlg=TRUE}
  if('y'%in%strsplit(log,'')[[1]]){y <- log(y,10);ylg=TRUE}
  lines(x,y,...)
  
}

add.points <- function(x,y,log='xy',...)
{
  xlg <- FALSE
  ylg <- FALSE
  if('x'%in%strsplit(log,'')[[1]]){x <- log(x,10);xlg=TRUE}
  if('y'%in%strsplit(log,'')[[1]]){y <- log(y,10);ylg=TRUE}
  points(x,y,...)
  
}
