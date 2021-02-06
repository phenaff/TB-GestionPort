## utility functions for manipulating time series

# last element
last <- function(x) {
  tail(x,1)
}

# first element
first <- function(x) {
  head(x,1)
}

# series of actual returns

yahoo.fname <- function(folder, ticker) {
  file.path(get.data.folder(), folder, tolower(ticker),'.rda',sep = '')
}

# fetch data from Yahoo! Finance

get.yahoo.series <- function(ticker, uStart, uEnd) {
  p <- get.hist.quote(
    instrument = ticker,
    start = uStart, end = uEnd,
    quote = c("AdjClose"),
    provider = "yahoo", compression = "d"
  )
  timeSeries(p)
}

# get one time serie from data folder 

get.fname <- function(ticker, folder=NULL, extension) {
  if(is.null(folder)) {
    fname <- file.path(get.data.folder(),
                       paste(tolower(ticker), extension, sep = '.'))
  } else {
    fname <- file.path(get.data.folder(), folder,
                       paste(tolower(ticker),extension, sep = '.'))
  } 
  fname
}
  
get.ts <- function(ticker, folder=NULL, returns=T, 
                   dtStart = NULL, dtEnd = NULL) {
    fname <- get.fname(ticker, folder, 'rda')
    
    if (file.exists(fname)) {
      load(fname)
      act.p <- ts
    } else {
      fname <- get.fname(ticker, folder, 'csv')
      tmp <- read.csv(fname,
        header = T, sep = ",", dec =".")
      t <- as.Date(tmp[,1], origin = '2000-01-01')
      act.p <- timeSeries(tmp[,2], t)
    }
    
    names(act.p) <- ticker
    if (returns)
      act.p <- returns(act.p)
    
    if (!is.null(dtStart) | !is.null(dtEnd)) {
      if (is.null(dtStart))
        dtStart <- first(time(act.p))
      if (is.null(dtEnd))
        dtEnd <- last(time(act.p))
      act.p <- window(act.p, start = dtStart, end = dtEnd)
    }
    act.p
  }

# clean up SBF ticker symbols

clean.sbf <- function(ticker) {
  ticker <- gsub("EU:", "", ticker, fixed = T)
  ticker <- gsub(" ", "", ticker, fixed = T)
  paste(ticker, '.PA', sep = '')
}

# List all tickers in a csv file

get.tickers.from.csv <- function(folder) {
  if ((folder %in% c('NASDAQ', 'EuroStoxx50', 'SBF120')) == FALSE)
    stop('Bad folder name')
  
  if (folder == 'NASDAQ') {
    fname <- file.path(get.data.folder(), 'NASDAQ.csv')
    tmp <- read.csv(
      fname, sep = ',', colClasses = 'character',
      header = T, col.names = c(
        'Ticker', 'Name',
        'dummy-1', 'MarketCap', 'dummy-2',
        'IPOYear', paste('dummy-', seq(3,5), sep =
                           '')
      )
    )
    ticker <- tmp$Ticker
  } else if (folder == 'EuroStoxx50') {
    fname <- file.path(get.data.folder(), 'EuroStoxx50.csv')
    tmp <- read.csv(
      fname, sep = ';', colClasses = 'character',
      header = F, col.names = c('ISIN', 'dummy', 'Ticker', 'Name')
    )
    ticker <- tmp$Ticker
  } else if (folder == 'SBF120') {
    fname <- file.path(get.data.folder(), 'sbf120.txt')
    tmp <-
      scan(
        fname, what = list(name = character()), sep = ',', strip.white = F
      )$name
    ticker <- as.vector(sapply(tmp, clean.sbf))
  }
  ticker
}

# load one time series

# get.ts <- function(folder, ticker) {
#   fname <- file.path(get.data.folder(),
#                      folder,
#                      paste(tolower(ticker),'.rda',sep = ''))
#   load(file = fname, verbose=FALSE)
#   ts
# }

# all tickers in a folder

get.tickers <- function(folder) {
  tmp <- list.files(path = file.path(get.data.folder(),folder))
  tickers <- sapply(tmp, function(t)
    gsub('.rda', '', t, fixed = T))
  tickers
}

#' Get all time series in a folder
#'
#' @param folder 
#' @param tickers 
#' @param returns 
#' @param dt.start 
#' @param combine 
#'
#' @return a list of time series if combine=F, a multivariate time series 
#' if combine=T
#' @export
#'
#' @examples
#' 
get.all.ts <-
  function(folder, tickers = NULL, returns = T,
           dt.start = NULL, combine = F) {
    if (is.null(tickers)) {
      tickers <- get.tickers(folder)
    }
    
    if (combine) {
      ts.all <- NULL
    } else {
      ts.all <- list()
    }
    
    for (i in seq_along(tickers)) {
      ticker <- tickers[i]
      tmp <- tryCatch(
        get.ts(ticker, folder, returns),
        error = function(err)
          FALSE, warning = function(warn)
            FALSE
      )
      
      if (is.logical(tmp))
        next()
      
      colnames(tmp) <- ticker
      if (combine) {
        if (!is.null(dt.start)) {
          dt.first <- time(tmp)[1]
          if (dt.first > timeDate(dt.start)) {
            warning(paste('First obs for ', ticker, dt.first, ' ... skipping'), call=FALSE)
            next
          }
          keep <- time(tmp) >= timeDate(dt.start)
          tmp <- tmp[keep,]
        }
        
        if (is.null(ts.all)) {
          ts.all <- tmp
        } else {
          ts.all <- cbind(ts.all, tmp)
        }
      } else {
        # just append to the list
        ts.all <- c(ts.all, list(tmp))
      }
    }
    
    if (combine) {
      if(!is.null(ts.all)) {
        ts.all <- removeNA(ts.all)
      } else {
        stop('No data found', call=FALSE)
      }
    }
    ts.all
  }
