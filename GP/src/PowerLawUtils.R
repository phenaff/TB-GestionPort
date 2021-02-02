#######################################
## Utility function for power-law model
## and analysis of vol clustering
#######################################

library(poweRlaw)
library(fields)
library(latex2exp)

############################
## Measure of vol clustering
############################

## ref: Conditional Probability as a measure of 
## Vol Clustering, Chen at el.
## difference between s.d. of last and first quartiles.
## should be >0

vol.clustering.cor <- function(R, NBIN=10) {
  tryCatch({
    Rp <- abs(lag(R, trim=TRUE))
    q <- quantile(as.vector(Rp), probs=seq(0,1,length.out=(NBIN+1)))
    z <- stats.bin(Rp, tail(R,-1), breaks=q)
    z$stats["Std.Dev.", NBIN] - z$stats["Std.Dev.", 1]
  },
  error=function(e) {NA}
  )                        
}

## ACF for one lag

vol.acf <- function(R, lg) {
  cov(R, lag(R,lg), use='pairwise')
}

## Empirical density: remove bins with zero observations

emp.density <-function(r, n=20) {
  tmp <- hist(r,n,plot=FALSE)
  kndx <- tmp$density>0
  list(x=tmp$mids[kndx], y=tmp$density[kndx])
}

## Conditional density plot
## empirical density of r conditional on rp
## fig 1(a) of Chen

plot.conditional.r <- function(ticker=NULL, r=NULL, bins=NULL, bw=NULL, nb.hist=20, 
                               NBIN=5, do.legend=T) {
  if(is.null(bins)) {
    if(is.null(r)) r <- act.series(ticker)
    bins <- make.bins(r, NBIN)
  }
  
  cl <- rainbow(NBIN)
  
  # z <- emp.density(r[indx==1], nb.hist)
  z <- emp.density(bins$r[bins$indx==1], nb.hist)
  x.limits <- c(quantile(r, .005), quantile(r, .995))
  
  logplot(z$x, z$y , log='y', col=cl[1], type='p', lwd=3, pch=20,
          
          xlab='r', ylab='P(r|rp)', forcexlim=x.limits)
  title(paste(ticker, 'Daily Return'), cex.main=.6, line=-2, adj=.1)
  
  for(k in seq(2,NBIN)) {
    # z <- emp.density(r[indx==k], nb.hist)
    z <- emp.density(bins$r[bins$indx==k], nb.hist)
    add.lines(z$x, z$y, col=cl[k], log='y',lwd=2, type='p', pch=20)
  }
  
  if(do.legend) legend('topright', paste('Q', seq(NBIN), sep=''), lty=rep(1,NBIN),
                       lwd=rep(2,NBIN), col=cl)
}

plot.scaling.factor <- function(plot.type='plot', ticker=NULL, r=NULL, NBIN=5, ...) {
  if(is.null(r)) r <- act.series(ticker)
  
  tmp <- make.bins(r, NBIN)
  sd.r <- tmp$sd.r
  mean.rp <- tmp$mean.cond
  
  if(plot.type=='plot') {
    logplot(mean.rp, sd.r, xlab=TeX('$r_p$'), ylab=TeX('$\\sigma(r_p)$'),
            log="xy", ...)
  } else {
    add.lines(mean.rp, sd.r, 
              log="xy", ...)
  }
}

make.bins <- function(r, NBIN=5, model='abs', ...) {
  
  if(model == 'abs'){
    # chen's model
    tmp <- cbind(r, lag(r))[-1,]
    r <- tmp[,1]
    rp <- tmp[,2]
    # the conditioning variable
    r.cond <- abs(rp)
  } else if(model == 'sigma.recursive') {
    r <- r-mean(r)
    nu <- power.gamma(r)
    vol <- sigma.recursive(r, nu, ...)
    indx <- !is.nan(vol)
    r <- r[indx[-length(indx)]]
    r.cond <- vol[indx[-length(indx)]]
  } else if(model == 'sigma.centered') {
    r <- r-mean(r)
    nu <- power.gamma(r)
    vol <- sigma.centered(r, nu, ...)
    indx <- !is.na(vol)
    r <- r[indx]
    r.cond <- vol[indx]
  }
  
  # breaks of NBIN bins of equal size
  
  breaks <- quantile(as.vector(r.cond), 
                     probs=seq(0,1,length.out=(NBIN+1)))
  
  # assign each observation to a bin
  indx <- rep(NA, length(r))
  indx[r.cond <= breaks[2]] <- 1
  for (k in 2:NBIN) {
    indx[r.cond <= breaks[k + 1] & r.cond > breaks[k]] <- k
  }
  
  # mean and s.d. of r in each bin
  mean.r <- sapply(seq(NBIN),
                   function(k) mean(r[indx==k]))
  sd.r <- sapply(seq(NBIN), 
                 function(k) sd(r[indx==k]))
  
  # mean and s.d. of conditioning variable
  mean.cond <- sapply(seq(NBIN),
                      function(k) mean(r.cond[indx==k]))
  sd.cond <- sapply(seq(NBIN), 
                    function(k) sd(r.cond[indx==k]))
  
  if(model == 'abs') {
    # scale r by s.d. of bin
    r.sc = rep(NA, length(r))
    for (k in seq(NBIN)) {
      r.sc[indx==k] <- r[indx==k]/sd.r[k]
    }
  } else if((model == 'sigma.recursive') | (model == 'sigma.centered')) {
    # scale r by day-ahead forecast of vol
    r.sc <- r / r.cond
  }
  
  list(indx=indx, r.sc=r.sc, r=r, breaks=breaks,
       mean.r=mean.r, sd.r=sd.r, mean.cond=mean.cond, sd.cond=sd.cond)
}

plot.cum.cond.pn <- function(ticker=NULL, r=NULL, NBIN=5) {
  # conditional distribution, scaled by sd of rp
  # plot positive and neg cdf on same graph. 
  # One plot per bin
  
  if(is.null(r)) r <- act.series(ticker)
  
  tmp <- make.bins(r, NBIN)
  r.sc <- tmp$r.sc
  indx <- tmp$indx
  
  # plot each bin
  
  colors <- rainbow(3)
  power.fit <- vector("list", 3)
  x.val <- vector("list", 3)
  for (k in seq(NBIN)) {
    x <- r.sc[indx==k]
    x.val[[1]] <- x[x>0]
    x.val[[2]] <- -x[x<0]
    x.val[[3]] <- abs(x[abs(x)>0])
    
    for(j in seq(3)) {
      m <- conpl$new(as.vector(x.val[[j]]))
      est <- estimate_xmin(m)
      print(est)
      power.fit[[j]] <- list(xmin=est$xmin, alpha=est$pars)
    }
    
    
    title.part <- sprintf("( %5.3f < rp <= %5.3f)", tmp$breaks[k], tmp$breaks[k+1])
    
    for(j in seq(3)) {
      powerplot(plot.type=ifelse(j==1,'plot','points'),
                x=x.val[[j]],
                sd=sqrt(var(x.val[[j]])),
                alpha=power.fit[[j]]$alpha, 
                xmin=power.fit[[j]]$xmin, 
                length.out = 20, col=colors[j], pch=20, cex=.6,
                xlab='r/w(rp)',ylab='P(r|rp) * w(rp)',
                ticker=ticker,
                main=paste('CDF of', ticker, title.part, sep=' '),
                theo.cdf=T)
    }
    
    # legend
    
    label.part <- c('r>0', 'r<0', 'all')
    labels <- vector(mode="character", 3) 
    for(k2 in seq(3)) {
      labels[k2] <- sprintf("%s % 5.1f% 4.1f", label.part[k2], power.fit[[k2]]$xmin,
                            power.fit[[k2]]$alpha)
    }
    op <- par(family="mono") 
    legend('bottomleft', labels , lty=rep(1,3),
           lwd=rep(2,3), col=colors)
    par(op)
  }
}

plot.cum.cond.r <- function(ticker=NULL, bins=NULL, r=NULL, sgn='all', NBIN=5,
                            do.legend=T, model='abs', params=NULL) {
  # conditional distribution, scaled by sd of rp
  # Fig 1(b) of Chen
  
  if(is.null(bins)) {
    if(is.null(r)) r <- act.series(ticker)
    bins <- make.bins(r, NBIN, model, params)
  }
  
  r.sc <- bins$r.sc
  indx <- bins$indx
  
  # good indices
  good.indx <- switch(sgn,
                      all = (r.sc>0) | (r.sc<0),
                      pos = (r.sc>0),
                      neg = (r.sc<0))
  
  r.s <- abs(r.sc[good.indx])
  indx <- indx[good.indx]
  
  # plot each bin
  
  colors <- rainbow(NBIN)
  power.fit <- vector("list", NBIN)
  for (k in seq(NBIN)) {
    x <- r.s[indx==k]
    m <- conpl$new(as.vector(x))
    est <- estimate_xmin(m)
    power.fit[[k]] <- list(xmin=est$xmin, alpha=est$pars)
    
    title.part <- switch(sgn,
                         all = 'all',
                         pos = 'r>0',
                         neg = 'r<0')
    
    if(model == 'abs') {
      labels = list(x=TeX('$\\frac{r}{\\sigma(r_p)}$'),
                    y=TeX('$P(r|r_p) \\times \\sigma(r_p)$'))
    } else if(model == 'sigma') {
      labels = list(x='r/sigma',y='P(r|sigma) * sigma')      
    }
    powerplot(plot.type=ifelse(k==1,'plot','points'),
              x=r.s[indx==k],
              sd=sd.cond[k],
              alpha=power.fit[[k]]$alpha, 
              xmin=power.fit[[k]]$xmin, 
              length.out = 20, col=colors[k], pch=20, cex=.6,
              xlab=labels$x, ylab=labels$y,
              ticker=ticker,
              main=paste('CDF of ', ticker, ' (', title.part, ')', sep=''),
              theo.cdf=T)
  }
  # legend
  if(do.legend) {
    labels <- vector(mode="character", NBIN)
    for(k in seq(NBIN)) {
      labels[k] <- sprintf("Q%02d% 5.1f% 4.1f", k, power.fit[[k]]$xmin,
                           power.fit[[k]]$alpha)
    }
    op <- par(family="mono") 
    legend('bottomleft', labels , lty=rep(1,NBIN),
           lwd=rep(2,NBIN), col=colors, 
           title=expression(paste(x[plain(min)], ' ', and, ' ', alpha, ' per bin')))
    
    par(op)
  }
}

## Helper functions for plotting power laws
## Fig 1(b) of Chen et al. Conditional Probability...

lseq = function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out=length.out))    
}

get_data_cdf = function(x, lower.tail=TRUE, pad=FALSE){ 
  
  sort.x <- sort(x)
  ecdf <- (1:length(sort.x))/length(sort.x)
  if(lower.tail)
    list(x=sort.x, y=ecdf)
  else
    list(x=sort.x, y=rev(ecdf))
}

# cumulative cdf of power law
dist_cdf <- function(x, alpha, xmin, x_axs) {
  
  list(x=x_axs, y=(x_axs/xmin)^(-alpha+1))
}

# power plot in log 10 scale, with theoretical cdf

powerplot <- function(plot.type='plot',x, sd, alpha, xmin, length.out = 20,
                      theo.cdf=TRUE, ...) {
  res.cdf = get_data_cdf(x, FALSE)
  x_axs = lseq(min(x), max(x), length.out)
  
  if(plot.type=='plot') 
    logplot(res.cdf$x, res.cdf$y, log="xy", 
            forceylim=log(c(1.e-4,10),10),
            ...)
  else
    add.points(res.cdf$x, res.cdf$y, log="xy", ...)
  
  if(theo.cdf) {
    x.trunc <- x[x>xmin]
    scale <- length(x.trunc)/length(x)
    res.m = dist_cdf(x.trunc, alpha, xmin, x_axs)
    add.lines(res.m$x, res.m$y*scale, lwd=2,
              log='xy', ...)
  }
}

