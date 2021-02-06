## build one data frame with all series
## utility functions

library(tseries)
library(timeSeries)
library(plotrix)
library(fields)
library(FNN)
library(emdist)

# display R matrix in text

write_matex2 <- function(x) {
  begin <- "\\begin{bmatrix}"
  end <- "\\end{bmatrix}"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  paste(c(begin, X, end), collapse = "")
}

# find time index closest to a given date

closest.index <- function(ts, date) {
  delta <- abs(as.Date(time(ts)) - date)
  idx <- which.min(delta)
  idx
}

# a plot when some values are off-the-chart

gap.plot.helper <- function(x, y, from, to, ...) {
  gp <- gap.plot(x, y, gap=c(from, to), ...)
  axis.break(2, from, breakcol="snow", style="gap")
  axis.break(2, from*(1+0.02), breakcol="black", style="slash")
  axis.break(4, from*(1+0.02), breakcol="black", style="slash")
  axis(2, at=from)
  gp
}

# stats from local bootstrap
get.stats <- function(ticker) {
  load(file=paste('stats.',ticker,'.rda',sep=''))
  eval(parse(text=paste('res.', ticker, '$stats', sep='')))
}

##
## iterator for simulated series
## not useful: prevents parallel %dopar%
##

iSeries <- function(ticker, count, chunk=1) {
  i <- 1
  nextEl <- function() {
    if((i+chunk-1)>count)
      stop('StopIteration')
    
    sim.r <- sim.series(ticker, i, i+chunk-1)
    i <<- i+chunk
    list(s=sim.r, id=i)
  }
  
  obj <- list(nextElem=nextEl)
  class(obj) <- c('iSeries', 'abstractiter', 'iter')
  obj
}

## bivariate (r_t, r_{t-1}) series
## The default is to normalize the series, since we are interested
## in the higher moments

r.to.v <- function(r, normalize=T) {
  if (normalize) {
    r.sd <- sqrt(var(r))[1]
    r <- (r-mean(r)) / r.sd
  }
  
  r <- as.vector(r) 
  s <- cbind(r, c(NA, head(abs(r),-1)))[-1,]
  colnames(s) <- c('r(t)', '|r(t-1)|')
  s
}

## bivariate probability binning distance
##

make.proBin <- function(v, minRow=NULL) {
  ff <- flowFrame(as.matrix(v))
  if(is.null(minRow)) minRow <- dim(v)[1]*.05
  bins <- proBin(ff, minRow, channels=colnames(v))
  bins
}

## allocate sample to reference binning partition
##

assign.to.proBin <- function(ref.bin, r=NULL, v=NULL) {
  if(is.null(v)) v <- r.to.v(r)
  ff <- flowFrame(as.matrix(v))
  proBins <- binByRef(ref.bin, ff)
  proBins
}

bin.count.proBin <- function(bins) {
  binId <- bins$table$dataIndx[bins$table$left == 0]
  length(binId)
}

## simple KL divergence on data vectors
##
univ.KL <- function(x,y, length.out=20) {
  # range of data
  if(!is.vector(x)) x <- as.vector(x)
  if(!is.vector(y)) y <- as.vector(y)
  
  zmin <- min(c(x,y))
  zmax <- max(c(x,y))
  breaks <- seq(zmin, zmax, length.out=length.out)
  z.x <- stats.bin(x, x, breaks=breaks)
  z.y <- stats.bin(y, y, breaks=breaks)
  p.x <- (1+z.x$stats["N",]) / (length(x)+length(breaks)-1)
  p.y <- (1+z.y$stats["N",]) / (length(y)+length(breaks)-1)
  sum(p.x * log(p.x/p.y))
}

## Chi-2 p value and KL divergence of binned 
## reference and sample data

binned.chi2 <- function(sim.r, ctrlRes, params) {
  
  sampRes = assign.to.proBin(ctrlRes, v=r.to.v(sim.r, normalize=F))
  
  lgs <- params$lgs
  
  binId <- ctrlRes$table$dataIndx[ctrlRes$table$left == 0]
  binCount <- length(binId)
  input = matrix(0, binCount,2)
  for (i in seq_len(binCount)) {
    input[i,1] = nrow(ctrlRes$data[[as.character(binId[i])]])
    input[i,2] = nrow(sampRes[[as.character(i)]])
  }
  
  # Pearson chi-square test
  pearson.chiSq <- chisq.test(input)$p.value
  print(paste('P-X2', pearson.chiSq))
  
  # Roederer chi.2
  counts <- apply(input,2,sum)
  ctrlCount <- counts[1]
  sampCount <- counts[2]
  
  x <- input[,1]/ctrlCount
  y <- input[,2]/sampCount
  
  chiSq <- (x-y)^2 / (x+y)
  pbStat <- (2 * ctrlCount * sampCount * (sum(chiSq))/
               (ctrlCount + sampCount) - (binCount - 1))/(sqrt(2 * (binCount - 1)))
  # KL divergence (crtl=P, sample=Q)
  ctrl.p <- x
  samp.p <- y
  kl.div <- sum(ctrl.p * log(ctrl.p/samp.p))
  list(p.X2=pearson.chiSq, 
       pb.X2=1-pnorm(pbStat), 
       SD=sqrt(var(sim.r)),
       SK=skewness(sim.r), 
       KU=kurtosis(sim.r),
       VC=vol.clustering.cor(sim.r))
}

bin.centers <- function(bins) {
  binId <- bins$table$dataIndx[bins$table$left == 0]
  binCount <- length(binId)
  
  xy <- foreach(j = seq(binCount), .combine=rbind) %do% {
    colMeans(bins$data[[as.character(binId[j])]])
  }
  xy
}

emd.dist <- function(sim.r, pro.bins, ctrl.v, 
                     params=list(bins='hex', 
                                 type='euclidian',
                                 maxIter=1000)) {
  browser()
  sim.v <- r.to.v(sim.r, normalize=T)
  if(params$bins == 'hex') {
    hB.act <- hexBinning(x=ctrl.v[,1], y=ctrl.v[,2])
    hB.sim <- hexBinning(x=sim.v[,1], y=sim.v[,2])
    A <- cbind(hB.act$z, hB.act$xcm, hB.act$ycm)
    B <- cbind(hB.sim$z, hB.sim$xcm, hB.sim$ycm)
  } else if(params$bins == 'Roederer') {
    sim.bins <- make.proBin(sim.v)
    
    # map actual and simulated bins centers
    
    xy.act <- bin.centers(pro.bins)
    xy.sim <- bin.centers(sim.bins)
    
    A <- cbind(rep(1, nrow(xy.act)), xy.act)
    B <- cbind(rep(1, nrow(xy.sim)), xy.sim)
  }
  
  d <- emd(A,B, maxIter=ifelse(is.null(params$maxIter), 1000, params$maxIter))
  d
}

chi2.dist.pearson <- function(sim.r, pro.bins, ctrl.v, params) {
  # case 1: reuse actual returns
  # case 2: comparison of two samples: combine data to compute bins
  
  sim.v <- r.to.v(sim.r, normalize=T)
  if(params$control == 'known') {
    # case 1
    sampRes = assign.to.proBin(pro.bins, v=sim.v)
    binId <- pro.bins$table$dataIndx[pro.bins$table$left == 0]
    binCount <- length(binId)
    input = matrix(0, binCount,2)
    for (i in seq_len(binCount)) {
      input[i,1] = nrow(pro.bins$data[[as.character(binId[i])]])
      input[i,2] = nrow(sampRes[[as.character(i)]])
    }
    res <- sum((input[,1]-input[,2])^2/(input[,1]+input[,2]))
  } else if(params$control == 'all-cells') {
    # case 2
    # combine sample and control, compute bins,
    # assign controls and sample to bins
    ctrl.r <- rnorm(length(sim.r))
    ctrl.v <- r.to.v(ctrl.r)
    pro.bins <- make.proBin(rbind(ctrl.v,sim.v), minRow=2*params$minRow)
    ctrlRes <- assign.to.proBin(pro.bins, v=ctrl.v)
    sampRes = assign.to.proBin(pro.bins, v=sim.v)
    binId <- pro.bins$table$dataIndx[pro.bins$table$left == 0]
    binCount <- length(binId)
    input = matrix(0, binCount,2)
    for (i in seq_len(binCount)) {
      input[i,1] = nrow(ctrlRes[[as.character(i)]])
      input[i,2] = nrow(sampRes[[as.character(i)]])
    }
    
    if(params=='R') {
      res <- chisq.test(input)$statistic      
    } else {
      counts <- apply(input,2,sum)
      browser()
      n <- sum(counts)
      rowSum <- input[,1] + input[,2]
      e = matrix(0, binCount,2)
      for(i in seq(binCount)) {
        for(j in seq(2)) {
          e[i,j] = rowSum[i] * counts[j] / n
        }
      }
      res <- sum((input-e)^2/e)
    }
    
  }
  res
}

chi2.dist <- function(sim.r, pro.bins, ctrl.v, params=list(type='Roederer')) {
  
  if(params$type == 'Pearson') {
    res <- chi2.dist.pearson(sim.r, pro.bins, ctrl.v, params)
  } else {
    # use bins computed on controls only
    sim.v <- r.to.v(sim.r, normalize=T)
    sampRes = assign.to.proBin(pro.bins, v=sim.v)
    binId <- pro.bins$table$dataIndx[pro.bins$table$left == 0]
    binCount <- length(binId)
    input = matrix(0, binCount,2)
    for (i in seq_len(binCount)) {
      input[i,1] = nrow(pro.bins$data[[as.character(binId[i])]])
      input[i,2] = nrow(sampRes[[as.character(i)]])
    }
    
    counts <- apply(input,2,sum)
    ctrlCount <- counts[1]
    sampCount <- counts[2]
    if(params$type == 'Roederer') {    
      # expected 
      e <- input[,1]/ctrlCount
      # observed
      o <- input[,2]/sampCount
      res <- sum((o-e)^2/(o+e)) * 2 * ctrlCount * sampCount / (ctrlCount + sampCount)
      # browser()
    } else if(params$type == 'Baggerly') {
      e <- input[,1] * sampCount / ctrlCount
      o <- input[,2] 
      res <- sum((o-e)^2/(e*(1+sampCount/ctrlCount)))
    }
  }
  res
}

hellinger.dist <- function(ref.r, sim.r, 
                           params=NULL, ...) {
  
  # make a symmetric distance: bins computed on union of 2 series
  
  # browser()
  ref.v <- r.to.v(ref.r)
  sim.v <- r.to.v(sim.r)
  all.v <- rbind(ref.v, sim.v)
  all.bins <- make.proBin(all.v)
  
  sim.bins = assign.to.proBin(all.bins, v=sim.v)
  ref.bins = assign.to.proBin(all.bins, v=ref.v)
  
  # browser()
  binId <- all.bins$table$dataIndx[all.bins$table$left == 0]
  binCount <- length(binId)
  input = matrix(0, binCount,2)
  for (i in seq_len(binCount)) {
    input[i,1] = nrow(ref.bins[[as.character(i)]])
    input[i,2] = nrow(sim.bins[[as.character(i)]])
  }
  
  counts <- apply(input,2,sum)
  
  p <- input[,1]/counts[1]
  q <- input[,2]/counts[2]
  
  sqrt(1/2) * sqrt(sum((sqrt(p) -sqrt(q))^2))
}

##
## proximity measures of return series 
##

proximity.measures <- function(ctrlRes, act.r, act.v, sim.r, params) {
  
  sim.v <- r.to.v(sim.r)
  sampRes = assign.to.proBin(ctrlRes, v=sim.v)
  
  lgs <- params$lgs
  
  kl.1b <- univ.KL(as.vector(act.v[,1]), as.vector(sim.v[,1]), length.out=lgs)
  kl.1t <- KL.divergence(as.matrix(act.v[,1]), as.matrix(sim.v[,1]), k=lgs, algorithm="kd_tree")
  kl.2 <- .5*(KL.divergence(as.matrix(act.v), as.matrix(sim.v), k=lgs, algorithm="kd_tree") +
                KL.divergence(as.matrix(sim.v), as.matrix(act.v), k=lgs, algorithm="kd_tree"))
  
  binId <- ctrlRes$table$dataIndx[ctrlRes$table$left == 0]
  binCount <- length(binId)
  input = matrix(0, binCount,2)
  for (i in seq_len(binCount)) {
    input[i,1] = nrow(ctrlRes$data[[as.character(binId[i])]])
    input[i,2] = nrow(sampRes[[as.character(i)]])
  }
  
  # Pearson chi-square test
  pearson.chiSq <- chisq.test(input)$p.value
  print(paste('P-X2', pearson.chiSq))
  # Roederer chi.2
  counts <- apply(input,2,sum)
  ctrlCount <- counts[1]
  sampCount <- counts[2]
  
  x <- input[,1]/ctrlCount
  y <- input[,2]/sampCount
  
  chiSq <- (x-y)^2 / (x+y)
  pbStat <- (2 * ctrlCount * sampCount * (sum(chiSq))/(ctrlCount + 
                                                         sampCount) - (binCount - 1))/(sqrt(2 * (binCount - 1)))
  # KL divergence (crtl=P, sample=Q)
  ctrl.p <- x
  samp.p <- y
  kl.div <- sum(ctrl.p * log(ctrl.p/samp.p))
  list(p.X2=pearson.chiSq, 
       pb.X2=1-pnorm(pbStat), 
       KL.1BF = exp(-kl.div),
       KL.1B=exp(-kl.1b), 
       KL.1T=exp(-kl.1t[lgs]), 
       KL.2=exp(-kl.2[lgs]), 
       SD=sqrt(var(sim.r)),
       SK=skewness(sim.r), 
       KU=kurtosis(sim.r),
       VC=0) # vol.clustering.cor(sim.r))
}

build_ts <- function(dir_name, nbRows) { 
  flist <- dir(dir_name, 'Something_*')
  
  for (file_name in flist) {
    series_id = as.numeric(unlist(strsplit(file_name, '_'))[2])
    # V4 is price
    df = read.csv(paste(dir_name,'/',file_name,sep=''), 
                  sep=';', header=F)['V4']
    
    if(nrow(df) != nbRows) {
      print(paste('skipping', file_name))
      next    
    }
    
    names(df) <- paste(dir_name, '_', series_id, sep='')
    
    if(!exists("df_final")) { 
      df_final = df}
    else {
      df_final <- cbind(df_final, df)
    }
    
  }
  t <- as.Date(seq(nbRows), origin='2000-01-01')
  timeSeries(df_final, t)
  
}

spearman.cor <- function(ts, n, method='spearman', 
                         confidence.level=.95) {
  
  # compute the spearman aurocorrelation of time series ts
  # and significance level using Student's test
  
  nb.obs = length(ts)
  
  rho = rep(0, n)
  rho.plus = rep(0, n)
  rho.minus = rep(0, n)
  
  for( i in seq(n)) {
    rho[i] = cor(ts, lag(ts,k=i),
                 method=method, use='pairwise')
    
    t <- rho[i] * sqrt((nb.obs-i-2)/(1-rho[i]^2))
    
    # critical value under H0
    
    stderr = sqrt(1.06 / (nb.obs-i-3))
    delta = abs(qnorm((1-confidence.level)/2)) * stderr
    rho.plus[i] = tanh(delta)
    rho.minus[i] = tanh(-delta)
  }
  list(rho=rho, rho.plus=rho.plus, rho.minus=rho.minus)
}  

## measure of volatility clustering (Chen et al.)
## difference in vol of r_t between 1st and last quantile of
## |r_p| as conditioning variable.
## report annualized vol

vol.clustering.cor <- function(R) {
  tryCatch({
    Rp <- abs(lag(R, trim=TRUE))
    q <- quantile(Rp, probs=seq(0,1,length.out=11))
    z <- stats.bin(Rp, tail(R,-1), breaks=q)
    (z$stats["Std.Dev.", 10] - z$stats["Std.Dev.", 1])*sqrt(252)
  },
           error=function(e) {print(e); NA}
  )                        
}


init.folder.name <- function(ticker) {
  paste('./ATOM/', toupper(ticker), '.PA.INIT/', sep='')
}

sim.file.name <- function(i) {
  paste('Something',i, sep='_')
}



sim.series <- function(ticker, j.start, j.end=NULL) {
  init.folder <- init.folder.name(ticker)
  base <- ifelse(toupper(ticker) == 'BNP', 144000, 168000)
  j <- j.start
  j.end <- ifelse(is.null(j.end), j.start, j.end)
  filename <- paste(init.folder, sim.file.name(base+j), sep='')
  tmp <-read.table(filename, header=FALSE, sep=";")
  t <- as.Date(tmp[,2], origin='2000-01-01')
  sim.p <- timeSeries(tmp[,4], t)
  if(j.end>j.start) {
    for(j in seq(j.start+1, j.end)) {
      filename <- paste(init.folder, sim.file.name(base+j), sep='')
      tmp <-read.table(filename, header=FALSE, sep=";")
      t <- as.Date(tmp[,2], origin='2000-01-01')
      tmp <- timeSeries(tmp[,4], t)
      sim.p <- cbind(sim.p, tmp)
    }
    names(sim.p) <- paste(ticker, j.start:j.end, sep='-')
  } else {
    names(sim.p) <- ticker
  }
  returns(sim.p)
}

# confidence intervals for empirical moments,
# used to speed up the filtering and remove extreme cases

in.range <- function(x, x.min, x.max) {
  (x >= x.min) && (x <= x.max)
}

# test range on moments, return boolean vector, 
# size = nb of cols

moments.OK <- function(Returns, stats=NULL, sk.lim=c(-.5,.5), ku.lim=c(1,5)) { 
  sk <- skewness(Returns)
  ku <- kurtosis(Returns)
  if(!is.null(stats)) {
    res = in.range(sk, stats[['lower', 'skewness']],
                   stats[['upper', 'skewness']]) &&
      in.range(ku, stats[['lower', 'kurtosis']],
               stats[['upper', 'kurtosis']])
  }
  else {
    res = in.range(sk, sk.lim[1], sk.lim[2]) &&
      in.range(ku, ku.lim[1], ku.lim[2])
  }
  res
}

# distance based on moments
moments.dist <- function(x, stats, weights=list(sk=1, ku=1)) {
  sk <- stats[['exp.val','skewness']]
  ku <- stats[['exp.val', 'kurtosis']]
  abs((x$SK-sk)/sk)*weights$sk + abs((x$KU-ku)/ku)*weights$ku
}

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

logplot <- function(x,y,log='xy',forceylim=c(0,0),forcexlim=c(0,0), ...)
{
  par(tck=0.02)
  xlg <- FALSE
  ylg <- FALSE
  if('x'%in%strsplit(log,'')[[1]]){x <- log(x,10);xlg=TRUE}
  if('y'%in%strsplit(log,'')[[1]]){y <- log(y,10);ylg=TRUE}
  yl <- ifelse(forceylim==c(0,0),range(y),forceylim)
  xl <- ifelse(forcexlim==c(0,0),range(x),forcexlim)
  plot(x,y,axes=FALSE,ylim=yl,xlim=xl, ...)
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

sim.stats <- function(x=NULL, ticker=NULL, 
                      iFirst=NULL, iLast=NULL,
                      func=NULL, alpha=NULL){
  
  nb <- length(func)
  
  if(is.null(x)) {
    theta.star <- foreach(j = iFirst:iLast, .combine=rbind) %dopar% {
      r <- as.matrix(sim.series(ticker, j))
      
      if((j %% 100) == 0) {print(paste(j, round(r[1,1], 5)))}
      sapply(func, function(foo) apply(r, 2, foo))
    }
  } else {
    nb.col <- ncol(x)
    theta.star <- foreach(j = seq(nb.col), .combine=rbind) %dopar% {
      r <- as.matrix(x[,j])
      
      if((j %% 100) == 0) {print(paste(j, round(r[1,1], 5)))}
      sapply(func, function(foo) apply(r, 2, foo))
    } 
  }
  
  m <- dim(theta.star)[1]
  for(j in seq(nb)) {
    theta.star[,j] <- sort(theta.star[,j])
  }
  cutoff <- floor((alpha/2)*(m+1)) 
  lower <- theta.star[cutoff,]
  upper <- theta.star[m+1-cutoff,]
  theta <- colMeans(theta.star)
  stats <- t(as.matrix(data.frame(theta=theta,
                                  lower=lower,upper=upper)))
  list(stats=stats, theta.star=theta.star)
}


stylized.facts <- function(ts.all, calcs) {
  
  # estimate stylized facts and tail distribution of time series
  # computations done on normalized series
  
  tickers <- colnames(ts.all)
  n <- length(tickers)
  
  res <- foreach(i = icount(n), .combine=rbind) %dopar% {
    print(paste('row: ',i))
    
    ts <- ts.all[,i]
    
    res.1 <- sapply(calcs, function(foo) apply(ts, 2, foo))
    ts <- abs(ts)
    m <- conpl$new(as.vector(ts[ts>0]))
    est <- estimate_xmin(m)
    
    res.2 <- c(est$xmin, est$pars)
    
    c(res.1, res.2)
  }
  
  rownames(res) <- tickers
  colnames(res) <- c(calcs, 'xmin', 'alpha')
  res
}

distance.matrix <- function(ts.all, normalize=T, dist.type='emd.dist') {
  # compute the distance matrix of a list of tickers
  
  tickers <- colnames(ts.all)
  n <- length(tickers)
  
  res <- foreach(i = icount(n-1), .combine=rbind) %dopar% {
    print(paste('row: ',i))
    
    # form the bivariate data set (r_t, |r_{t-1}|)
    act.v <- r.to.v(ts.all[,i], normalize)
    ref.bins <- make.proBin(act.v)
    
    a.row = rep(0,n)
    for(j in seq(i+1, n)) {
      a.row[j] <- do.call(dist.type, 
                          list(ref.r=ts.all[,i], sim.r = ts.all[,j], pro.bins=ref.bins))
    }
    a.row
  }
  res <- rbind(res, rep(0, n))
  
  for(i in seq(n)) {
    for(j in seq(1,i)) {
      res[i,j] <- res[j,i] 
    }
  }
  
  row.names(res) <- tickers
  res
}

chi2.bivariate <- function(r, nc, ns, probs=seq(0,1,.05)) {
  
  ctrl.r <- r[1:nc,] 
  sim.r <- r[(1+nc):(nc+ns),]
  ctrl.v <- r.to.v(ctrl.r, normalize=F)
  sim.v <- r.to.v(sim.r, normalize=F)
  all.v <- rbind(ctrl.v, sim.v)
  all.bins <- make.proBin(all.v)
  
  sim.bins = assign.to.proBin(all.bins, v=sim.v)
  ctrl.bins = assign.to.proBin(all.bins, v=ctrl.v)
  
  binId <- all.bins$table$dataIndx[all.bins$table$left == 0]
  binCount <- length(binId)
  input = matrix(0, binCount,2)
  for (i in seq_len(binCount)) {
    input[i,1] = nrow(ctrl.bins[[as.character(i)]])
    input[i,2] = nrow(sim.bins[[as.character(i)]])
  }
  
  counts <- apply(input,2,sum)
  ctrlCount <- counts[1]
  sampCount <- counts[2]
  
  e <- input[,1]
  o <- input[,2] 
  res <- sum((o-e)^2/(e*(1+sampCount/ctrlCount)))
  
  res
}

# Pearson's chi-2 test of two sample time series
# univariate binning of returns

chi2.contingency <- function(r, nc, ns, probs=seq(0,1,.05)) {
  
  ctrl.v <- r[1:nc] 
  sim.v <- r[(nc+1):(nc+ns)] 
  
  brk <- quantile(ctrl.v, probs=probs)
  brk[1] <- -Inf
  brk[length(brk)] <- Inf
  h.e <- hist(ctrl.v, breaks=brk, plot=F)$count
  h.o <- hist(sim.v, breaks=brk, plot=F)$count
  res <- chisq.test(rbind(h.e, h.o))$statistic
  res
}

# Baggerly's exact chi-2 test of two samples
# univariate binning

chi2.pb <- function(r, nc, ns, probs=seq(0,1,.05), type='contingency') {
  
  ctrl.v <- r[1:nc]
  sim.v <- r[(nc+1):(nc+ns)]
  
  brk <- quantile(ctrl.v, probs=probs)
  brk[1] <- -Inf
  brk[length(brk)] <- Inf
  h.e <- hist(ctrl.v, breaks=brk, plot=F)$count
  h.o <- hist(sim.v, breaks=brk, plot=F)$count
  
  h.e <- h.e * ns / nc
  if(type=='contingency') {
    res <- sum((h.o-h.e)^2/(h.e*(1+ns/nc)))
  } else if(type=='gof') {
    res <- sum((h.o-h.e)^2/h.e)
  }
  res
}

# report on chi-2 experiment
dist.report <- function(dist, nb.bins, confidence, do.print=T, do.plot=F) {
  x.alpha <- qchisq(confidence, df=nb.bins-1)
  
  pct.sig <- sum(dist>x.alpha)/length(dist)
  dist.sorted <- sort(dist)
  x.alpha.hat <- dist.sorted[round(confidence*length(dist))]
  res <- c(mean(dist), var(dist), pct.sig,
           x.alpha.hat)
  
  print(paste('Mean actual: ', round(res[1],2), ' Theory: ', 
              nb.bins-1))
  print(paste('Var actual: ', round(res[2],2), 
              ' Theory: ', round((2*(nb.bins-1)),2)))
  
  print(paste('Alpha actual: ', round(x.alpha.hat,2), 
              ' Theory: ', round(x.alpha,2)))
  
  print(paste('% significant (alpha=5%): ', round(100*res[3],2)))
  
  if(do.plot) {
    dist.cum <- seq(length(dist.sorted))/length(dist.sorted)
    
    title.cdf = expression(paste('cdf of ', chi^2, ' homogeneity test\nfor normal variables'))
    
    plot(dist.sorted, dist.cum, pch=20, col='blue', xlab='Distance', ylab='CDF',
         main=title.cdf)
    
    # compare to theoretical dist
    lines(dist.sorted, pchisq(dist.sorted, df=nb.bins-1), lwd=3, col='red')
  }
  
  res <- matrix(res, 1, 4)
  colnames(res) <- c('mean', 'var', 'alpha.hat', 'critical.value')
  res
}
