# Clean up workspace -- 
rm(list=ls())

# Load or install necessary packages if necessary
want <- c("quantmod","dplyr","PerformanceAnalytics", "TTR", "magrittr")
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)
lapply(want, function(i) require(i, character.only=TRUE))
rm(want, need)

# Working directories
dir <- list()
dir$root <- dirname(getwd())
dir$source <- paste(dir$root,"/data",sep="")
dir$output <- paste(dir$root,"/figures",sep="")
lapply(dir, function(i) dir.create(i, showWarnings = F))

# stocks pool, 100 stocks from nasdaq
from = "1980-01-01"
# to = "2021-11-30"
threshold <- 0.05

nasdaq <- c("ATVI", "ADBE", "AMD", "ALGN", "GOOG", "AMZN", "AEP", "AMGN",  
            "ADI", "ANSS", "AAPL", "AMAT", "ASML", "TEAM", "ADSK", "ADP", "BIDU",  
            "BIIB", "BKNG", "AVGO", "CDNS", "CDW", "CERN", "CHTR", "CHKP", "CTAS",  
            "CSCO", "CTSH", "CMCSA", "CPRT", "COST", "CRWD", "CSX", "DXCM", "DOCU",  
            "DLTR", "EBAY", "EA", "EXC", "FB", "FAST", "FISV", "FOX", "GILD",  
            "HON", "IDXX", "ILMN", "INCY", "INTC", "INTU", "ISRG", "JD", "KDP",   
            "KLAC", "KHC", "LRCX", "LULU", "MAR", "MRVL", "MTCH", "MELI", "MCHP",  
            "MU", "MSFT", "MRNA", "MDLZ", "MNST", "NTES", "NFLX", "NVDA", "NXPI",  
            "ORLY", "OKTA", "PCAR", "PAYX", "PYPL", "PTON", "PEP", "PDD", "QCOM",  
            "REGN", "ROST", "SGEN", "SIRI", "SWKS", "SPLK", "SBUX", "SNPS", "TMUS",  
            "TSLA", "TXN", "TCOM", "VRSN", "VRSK", "VRTX", "WBA", "WDAY", "XEL",   
            "XLNX", "ZM")

dj30 <- c ("AXP", "AMGN", "AAPL", "BA", "CAT", "CSCO", "CVX", "GS", "HD", "HON", 
           "IBM", "INTC", "JNJ", "KO", "JPM", "MCD", "MMM", "MRK", "MSFT", "NKE",
           "PG", "TRV", "UNH", "CRM", "VZ", "V", "WBA", "WMT", "DIS", "DOW")

sp500 <- c("MMM", "AOS", "ABT", "ABBV", "ABMD", "ACN", "ATVI", "ADM", "ADBE", "AAP", "AMD", "AES", 
           "AFL", "A", "APD", "AKAM", "ALK", "ALB", "ARE", "ALGN", "ALLE", "LNT", "ALL", "GOOGL",
           "GOOG", "MO", "AMZN", "AMCR", "AEE", "AAL", "AEP", "AXP", "AIG", "AMT", "AWK", "AMP",
           "ABC", "AME", "AMGN", "APH", "ADI", "ANSS", "ANTM", "AON", "APA", "AAPL", "AMAT", "APTV", 
           "ANET", "AJG", "AIZ", "T", "ATO", "ADSK", "ADP", "AZO", "AVB", "AVY", "BKR", "BLL", 
           "BAC", "BBWI", "BAX", "BDX", "BRK-B", "BBY", "BIO", "TECH", "BIIB", "BLK", "BK", "BA",   
           "BKNG", "BWA", "BXP", "BSX", "BMY", "AVGO", "BR", "BRO", "BF-B", "CHRW", "CDNS", "CZR",  
           "CPB", "COF", "CAH", "KMX", "CCL", "CARR", "CTLT", "CAT", "CBOE", "CBRE", "CDW", "CE",   
           "CNC", "CNP", "CDAY", "CERN", "CF", "CRL", "SCHW", "CHTR", "CVX", "CMG", "CB", "CHD", 
           "CI", "CINF", "CTAS", "CSCO", "C", "CFG", "CTXS", "CLX", "CME", "CMS", "KO", "CTSH",
           "CL", "CMCSA", "CMA",  "CAG", "COP", "ED", "STZ", "CPRT", "GLW", "CTVA", "COST", "CTRA", 
           "CCI", "CSX", "CMI", "CVS", "DHI", "DHR", "DRI", "DVA", "DE", "DAL", "XRAY", "DVN",  
           "DXCM", "FANG", "DLR", "DFS",  "DISCA", "DISCK", "DISH", "DG", "DLTR", "D", "DPZ", "DOV", 
           "DOW", "DTE", "DUK", "DRE", "DD", "DXC", "EMN", "ETN", "EBAY", "ECL", "EIX", "EW",   
           "EA", "LLY", "EMR", "ENPH", "ETR", "EOG", "EFX", "EQIX", "EQR", "ESS", "EL", "ETSY", 
           "RE", "EVRG", "ES", "EXC", "EXPE", "EXPD", "EXR", "XOM", "FFIV", "FB", "FAST", "FRT", 
           "FDX", "FIS", "FITB", "FRC", "FE", "FISV", "FLT", "FMC", "F", "FTNT", "FTV", "FBHS",
           "FOXA", "FOX", "BEN", "FCX", "GPS", "GRMN", "IT", "GNRC", "GD", "GE", "GIS", "GM",  
           "GPC", "GILD", "GPN", "GL", "GS", "HAL", "HBI", "HAS", "HCA", "PEAK", "HSIC", "HES", 
           "HPE", "HLT", "HOLX", "HD", "HON", "HRL", "HST", "HWM", "HPQ", "HUM", "HBAN", "HII",
           "IBM", "IEX", "IDXX", "INFO", "ITW", "ILMN", "INCY", "IR", "INTC", "ICE", "IFF", "IP",   
           "IPG", "INTU", "ISRG", "IVZ", "IPGP", "IQV", "IRM", "JBHT", "JKHY", "J", "SJM", "JNJ",  
           "JCI", "JPM", "JNRFX", "KSU", "K", "KEY", "KEYS", "KMB", "KIM", "KMI", "KLAC", "KHC", 
           "KR", "LHX", "LH", "LRCX", "LW", "LVS", "LEG", "LDOS", "LEN", "LNC", "LIN", "LYV",  
           "LKQ", "LMT", "L", "LOW", "LUMN", "LYB", "MTB", "MRO", "MPC", "MKTX", "MAR", "MMC",  
           "MLM", "MAS", "MA", "MTCH", "MKC", "MCD", "MCK", "MDT", "MRK", "MET" , "MTD", "MGM",  
           "MCHP", "MU", "MSFT", "MAA", "MRNA", "MHK", "TAP", "MDLZ", "MPWR",  "MNST", "MCO", "MS",   
           "MSI", "MSCI", "NDAQ", "NTAP", "NFLX", "NWL", "NEM", "NWSA", "NWS", "NEE", "NLSN", "NKE",  
           "NI", "NSC", "NTRS", "NOC", "NLOK", "NCLH", "NRG", "NUE", "NVDA", "NVR", "NXPI", "ORLY",
           "OXY", "ODFL", "OMC", "OKE", "ORCL", "OGN", "OTIS", "PCAR", "PKG", "PH", "PAYX", "PAYC",
           "PYPL", "PENN", "PNR", "PBCT", "PEP", "PKI", "PFE", "PM", "PSX", "PNW", "PXD", "PNC",  
           "POOL", "PPG", "PPL", "PFG", "PG", "PGR", "PLD", "PRU", "PTC", "PEG", "PSA", "PHM",  
           "PVH", "QRVO", "QCOM", "PWR", "DGX", "RL", "RJF", "RTX", "O", "REG", "REGN", "RF",
           "RSG", "RMD", "RHI", "ROK", "ROL", "ROP", "ROST", "RCL", "SPGI", "CRM", "SBAC", "SLB", 
           "STX","SEE", "SRE", "NOW", "SHW", "SPG", "SWKS", "SNA", "SO", "LUV", "SWK", "SBUX",
           "STT", "STE", "SYK", "SIVB", "SYF", "SNPS", "SYY", "TMUS", "TROW", "TTWO", "TPR", "TGT",  
           "TEL", "TDY", "TFX", "TER", "TSLA", "TXN", "TXT", "COO", "HIG", "HSY", "MOS", "TRV",  
           "DIS", "TMO", "TJX", "TSCO", "TT", "TDG", "TRMB", "TFC", "TWTR", "TYL", "TSN", "USB",
           "UDR", "ULTA", "UAA", "UA", "UNP", "UAL", "UPS", "URI", "UNH", "UHS", "VLO", "VTR",  
           "VRSN", "VRSK", "VZ", "VRTX", "VFC", "VIAC", "VTRS", "V", "VNO", "VMC", "WRB", "GWW",  
           "WAB", "WBA", "WMT", "WM", "WAT", "WEC", "WFC", "WELL", "WST", "WDC", "WU", "WRK", 
           "WY", "WHR", "WMB", "WLTW", "WYNN", "XEL", "XLNX", "XYL", "YUM", "ZBRA", "ZBH", "ZION", "ZTS")

hsi <- c("0005", "0011", "0388", "0939", "1299", "1398", "2318", "2388", "2628", "3968", 
         "3988", "0002", "0003", "0006", "1038", "0012", "0016", "0017", "0101", "0688", 
         "0823", "0960", "1109", "1113", "1997", "2007", "6098", "0001", "0027", "0066", 
         "0175", "0241", "0267", "0288", "0386", "0669", "0700", "0762", "0857", "0868", 
         "0883", "0941", "0968", "1044", "1093", "1177", "1211", "1810", "1876", "1928",
         "2018", "2020", "2269", "2313", "2319", "2331", "2382", "3690", "6862", "9988")
hsi <- lapply(1:60, function(i){
  hsi[i] <- paste(hsi[i], ".hk", sep="")
})
hsi <- unlist(hsi)

ss <- c("600000", "600028", "600036", "600104", "600309", "600547", "600588", "600745", 
        "600887", "601012", "601138", "601288", "601398", "601668", "601857", "601995", 
        "603501", "600009", "600030", "600048", "600196", "600438", "600570", "600690",
        "600809", "600893", "601066", "601166", "601318", "601601", "601688", "601888",
        "603259", "603986", "600016", "600031", "600050", "600276", "600519", "600585",
        "600703", "600837", "600918", "601088", "601211", "601336", "601628", "601818", 
        "601899", "603288")
ss <- lapply(1:50, function(i){
  ss[i] <- paste(ss[i], ".ss", sep="")
})
ss <- unlist(ss)

# import stock pool
nasdaq_pool <- lapply(1:length(nasdaq), 
                      function(i) nasdaq[i] <- getSymbols(nasdaq[i], from = from
                                                          , auto.assign = FALSE, reload.Symbols = TRUE))
dj30_pool <- lapply(1:length(dj30), 
                    function(i) dj30[i] <- getSymbols(dj30[i], from = from 
                                                      , auto.assign = FALSE, reload.Symbols = TRUE))
sp500_pool <- lapply(1:length(sp500), 
                     function(i) sp500[i] <- getSymbols(sp500[i], from = from 
                                                       , auto.assign = FALSE, reload.Symbols = TRUE))
hsi_pool <- lapply(1:length(hsi), 
                   function(i) hsi[i] <- getSymbols(hsi[i], from = from
                                                    , auto.assign = FALSE, reload.Symbols = TRUE))
ss_pool <- lapply(1:length(ss), 
                  function(i) ss[i] <- getSymbols(ss[i], from = from
                                                  , auto.assign = FALSE, reload.Symbols = TRUE))

# renames
names(nasdaq_pool) <- nasdaq
names(dj30_pool) <- dj30
names(sp500_pool) <- sp500
names(hsi_pool) <- hsi
names(ss_pool) <- ss

# save to the directory /data in csv format
lapply(1:length(nasdaq), function(i){
  stock <- xts(data.frame(nasdaq_pool[i]), order.by = as.Date(rownames(data.frame(nasdaq_pool[i]))))
  write.csv(data.frame(date = index(stock), coredata(stock)), 
            file = paste(dir$source, "/nasdaq/", nasdaq[i], ".csv", sep=""), row.names = F)})

lapply(1:length(dj30), function(i){
  stock <- xts(data.frame(dj30_pool[i]), order.by = as.Date(rownames(data.frame(dj30_pool[i]))))
  write.csv(data.frame(date = index(stock), coredata(stock)), 
            file = paste(dir$source, "/dj/", dj30[i], ".csv", sep=""), row.names = F)})

lapply(1:length(sp500), function(i){
  stock <- xts(data.frame(sp500_pool[i]), order.by = as.Date(rownames(data.frame(sp500_pool[i]))))
  write.csv(data.frame(date = index(stock), coredata(stock)), 
            file = paste(dir$source, "/sp500/", sp500[i], ".csv", sep=""), row.names = F)})

lapply(1:length(hsi), function(i){
  stock <- xts(data.frame(hsi_pool[i]), order.by = as.Date(rownames(data.frame(hsi_pool[i]))))
  write.csv(data.frame(date = index(stock), coredata(stock)), 
            file = paste(dir$source, "/hk/", hsi[i], ".csv", sep=""), row.names = F)})

lapply(1:length(ss), function(i){
  stock <- xts(data.frame(ss_pool[i]), order.by = as.Date(rownames(data.frame(ss_pool[i]))))
  write.csv(data.frame(date = index(stock), coredata(stock)), 
            file = paste(dir$source, "/ss/", ss[i], ".csv", sep=""), row.names = F)})


#######################################################################
# trading strategies 

# basic strategies

# sma / ema / wma / 
strat_sma <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 50) stop("x should be longer than 50.")
  sma10 <- SMA(Cl, n = 10)
  sma20 <- SMA(Cl, n = 20)
  sma50 <- SMA(Cl, n = 50)
  res <- NULL
  # from 51
  for (i in 51:(n-1)) {
    if ((sma10[i - 1] < sma20[i - 1] || sma20[i - 1] < sma50[i - 1])
        && (sma10[i] >= sma20[i] && sma20[i] >= sma50[i])) {
      signal <- "buy"
    } else if ((sma10[i - 1] > sma20[i - 1] || sma20[i - 1] > sma50[i - 1])
               && (sma10[i] <= sma20[i] && sma20[i] <= sma50[i])) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# macd
strat_macd <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 35) stop("x should be longer than 35.")
  # MACD is a function from TTR.
  MACDS <- MACD(Cl, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA")
  macd <- MACDS[, "macd"]
  macdsignal <- MACDS[, "signal"]
  res <- NULL
  # from 34 = (26 - 1 + 9 - 1) + 1
  for (i in 35:(n-1)) {
    if (macd[i - 1] < macdsignal[i - 1] && macd[i] > macdsignal[i]) {
      signal <- "buy"
    } else if (macd[i - 1] > macdsignal[i - 1] && macd[i] < macdsignal[i]) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# rsi
strat_rsi <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (nlen < 14) stop(paste("x should be longer than 14."))
  rsi <- data.frame(RSI(Cl, 14))
  res <- NULL
  for (i in 35:(n-1)) {
    if (TF30(i,i,rsi$rsi)) {
      signal <- "buy"
    } else if (TF70(i,i,rsi$rsi))
    {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# bollinger's bands
strat_bbands <- function(x, n = 20) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- Cl(x)
  Op <- Op(x)
  nlen <- length(Cl)
  if (nlen < n) stop(paste("x should be longer than", n, "."))
  bbands <- BBands(HLC(x), n = n)
  res <- NULL
  for (i in (n + 1):(nlen-1)) {
    if (Cl[i - 1] > bbands$dn[i - 1] && Cl[i] < bbands$dn[i]) {
      signal <- "buy"
    } else if (Cl[i - 1] < bbands$up[i - 1] && Cl[i] > bbands$up[i])
    {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# dmi
strat_dmi <- function(x, n = 14){
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- Cl(x)
  Op <- Op(x)
  nlen <- length(Cl)
  if (nlen < n) stop(paste("x should be longer than", n, "."))
  adx <- ADX(HLC(x))
  DIp <- adx[, "DIp"]  # green positive Direction Index
  DIn <- adx[, "DIn"]  # red negative Direction Index
  res <- NULL
  for (i in (n + 2):(nlen-1)) {
    if (DIp[i - 1] < DIn[i - 1] && DIp[i] > DIn[i]) {
      signal <- "buy"
    } else if (DIp[i-1] > DIn[i - 1] && DIp[i] < DIn[i])
    {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# roc
strat_roc <- function(x, n = 20){
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- Cl(x)
  Op <- Op(x)
  nlen <- length(Cl)
  if (nlen < n) stop(paste("x should be longer than", n, "."))
  roc <- ROC(Cl, n = n)
  res <- NULL
  for (i in (n + 2):(nlen-1)) {
    if ((roc[i-1] < 0) && (roc[i] > 0)) {
      signal <- "buy"
    } else if ((roc[i-1] > 0) && (roc[i] < 0))
    {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# smi / stochastic / stc
strat_smi <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 35) stop("x should be longer than 35.")
  # MACD is a function from TTR.
  smi_result <- SMI(HLC(x), n = 13, nFast = 2, nSlow = 25, nSig = 9)
  smi <- smi_result[, "SMI"]
  smisignal <- smi_result[, "signal"]
  res <- NULL
  # from 34 = (26 - 1 + 9 - 1) + 1
  for (i in 35:(n-1)) {
    if (smi[i - 1] < smisignal[i - 1] && smi[i] > smisignal[i]) {
      signal <- "buy"
    } else if (smi[i - 1] > smisignal[i - 1] && smi[i] < smisignal[i]) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# william index
strat_wpr <- function(x, n = 20){
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- Cl(x)
  Op <- Op(x)
  nlen <- length(Cl)
  if (nlen < n) stop(paste("x should be longer than", n, "."))
  wpr <- WPR(HLC(x), n = n)
  res <- NULL
  for (i in (n + 1):(nlen-1)) {
    if (wpr[i] < 0.2) {
      signal <- "buy"
    } else if (wpr[i] > 0.8)
    {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# sar
# Parabolic Stop And Reverse (SAR)
strat_sar <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  sar <- SAR(cbind(Hi(x), Lo(x)), accel = c(0.02, 0.2))
  res <- NULL
  for (i in 2:(n-1)) {
    if (sar[i-1]>Cl[i-1] && sar[i]<Cl[i]) {
      signal <- "buy"
    } else if (sar[i-1] < Cl[i-1] && sar[i] > Cl[i]) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# cci
strat_cci <- function(x, n = 20, c = 0.015) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  nlen <- length(Cl)
  if (nlen < 20) stop("x should be longer than 20.")
  cci <- CCI(HLC(x), n = n, c = c)
  res <- NULL
  for (i in 20:(nlen-1)) {
    if (is.na(cci[i])){
      next()
    }
    if (cci[i] <= -100) {
      signal <- "buy"
    } else if (cci[i] >= 100) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# obv
strat_obv <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- as.numeric(Cl(x))
  Op <- as.numeric(Op(x))
  n <- length(Cl)
  obv <- data.frame(OBV(Cl, Vo(x)))
  res <- NULL
  for (i in 10:(n-1)) {
    if ((obv[i,] - obv[i-4,]) > 0 && (Cl[i]-Cl[i-4]) > 0) {
      signal <- "buy"
    } else if ((obv[i,] - obv[i-4,]) < 0 && (Cl[i]-Cl[i-4]) < 0) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

#pvt
strat_pvt <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- as.numeric(Cl(x))
  Op <- as.numeric(Op(x))
  Vo <- as.numeric(Vo(x))
  pvt <- ((Cl-lag(Cl,1))/lag(Cl,1))*Vo
  n <- length(Cl)
  res <- NULL
  for (i in 9:(n-1)) {
    if ((Cl[i]-Cl[i-7]) > 0 && (pvt[i]-pvt[i-7])>0) {
      signal <- "buy"
    } else if ((Cl[i]-Cl[i-7]) < 0 && (pvt[i]-pvt[i-7])<0) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

#######################################################################

# multi-strategies

# pvt & stc
strat_pvtstc <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- as.numeric(Cl(x))
  Op <- as.numeric(Op(x))
  Vo <- as.numeric(Vo(x))
  pvt <- ((Cl-lag(Cl,1))/lag(Cl,1))*Vo
  n <- length(Cl)
  smi <- SMI(HLC(x), n = 13, nFast = 2, nSlow = 25, nSig = 9)
  smil <- as.numeric(smi[, 1])
  smisignal <- as.numeric(smi[, 2])
  res <- NULL
  for (i in 35:(n-1)) {
    if ((Cl[i]-Cl[i-7]) > 0 && (pvt[i]-pvt[i-7])>0 && (smi[i - 1] < smisignal[i - 1]) && (smi[i] > smisignal[i])) {
      signal <- "buy"
    } else if ((Cl[i]-Cl[i-7]) < 0 && (pvt[i]-pvt[i-7])<0 && (smi[i - 1] > smisignal[i - 1]) && (smi[i] < smisignal[i])) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# pvt & macd
strat_pvtmacd <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- as.numeric(Cl(x))
  Op <- as.numeric(Op(x))
  Vo <- as.numeric(Vo(x))
  pvt <- ((Cl-lag(Cl,1))/lag(Cl,1))*Vo
  n <- length(Cl)
  MACDS <- MACD(Cl, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA")
  macd <- MACDS[, "macd"]
  macdsignal <- MACDS[, "signal"]
  res <- NULL
  for (i in 35:(n-1)) {
    if ((Cl[i]-Cl[i-7]) > 0 && (pvt[i]-pvt[i-7])>0 && macd[i - 1] < macdsignal[i - 1] && macd[i] > macdsignal[i]) {
      signal <- "buy"
    } else if ((Cl[i]-Cl[i-7]) < 0 && (pvt[i]-pvt[i-7])<0 && macd[i - 1] > macdsignal[i - 1] && macd[i] < macdsignal[i]) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# pvt & sma
strat_pvtsma <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- as.numeric(Cl(x))
  Op <- as.numeric(Op(x))
  Vo <- as.numeric(Vo(x))
  pvt <- ((Cl-lag(Cl,1))/lag(Cl,1))*Vo
  n <- length(Cl)
  if (n < 50) stop("x should be longer than 50.")
  sma10 <- SMA(Cl, n = 10)
  sma20 <- SMA(Cl, n = 20)
  sma50 <- SMA(Cl, n = 50)
  res <- NULL
  for (i in 51:(n-1)) {
    if ((Cl[i]-Cl[i-7]) > 0 && (pvt[i]-pvt[i-7])>0 && (sma10[i - 1] < sma20[i - 1] || sma20[i - 1] < sma50[i - 1])
        && (sma10[i] >= sma20[i] && sma20[i] >= sma50[i])) {
      signal <- "buy"
    } else if ((Cl[i]-Cl[i-7]) < 0 && (pvt[i]-pvt[i-7])<0 && (sma10[i - 1] > sma20[i - 1] || sma20[i - 1] > sma50[i - 1])
               && (sma10[i] <= sma20[i] && sma20[i] <= sma50[i])) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# obv & sma
strat_obvsma <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- as.numeric(Cl(x))
  Op <- as.numeric(Op(x))
  n <- length(Cl)
  obv <- data.frame(OBV(Cl, Vo(x)))
  res <- NULL
  if (n < 50) stop("x should be longer than 50.")
  sma10 <- SMA(Cl, n = 10)
  sma20 <- SMA(Cl, n = 20)
  sma50 <- SMA(Cl, n = 50)
  for (i in 51:(n-1)) {
    if ((obv[i,] - obv[i-4,]) > 0 && (Cl[i]-Cl[i-4]) > 0 && (sma10[i - 1] < sma20[i - 1] || sma20[i - 1] < sma50[i - 1])
        && (sma10[i] >= sma20[i] && sma20[i] >= sma50[i])) {
      signal <- "buy"
    } else if ((obv[i,] - obv[i-4,]) < 0 && (Cl[i]-Cl[i-4]) < 0 &&  (sma10[i - 1] > sma20[i - 1] || sma20[i - 1] > sma50[i - 1])
               && (sma10[i] <= sma20[i] && sma20[i] <= sma50[i])) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# obv & macd
strat_obvmacd <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- as.numeric(Cl(x))
  Op <- as.numeric(Op(x))
  n <- length(Cl)
  obv <- data.frame(OBV(Cl, Vo(x)))
  MACDS <- MACD(Cl, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA")
  macd <- MACDS[, "macd"]
  macdsignal <- MACDS[, "signal"]
  res <- NULL
  for (i in 35:(n-1)) {
    if ((obv[i,] - obv[i-4,]) > 0 && (Cl[i]-Cl[i-4]) > 0 && macd[i - 1] < macdsignal[i - 1] && macd[i] > macdsignal[i]) {
      signal <- "buy"
    } else if ((obv[i,] - obv[i-4,]) < 0 && (Cl[i]-Cl[i-4]) < 0 && macd[i - 1] > macdsignal[i - 1] && macd[i] < macdsignal[i]) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# macd & rsi
macdrsi <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 35) stop("x should be longer than 35.")
  MACDS <- MACD(Cl, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA")
  macd <- MACDS[, "macd"]
  macdsig <- MACDS[, "signal"]
  rsi <- data.frame(RSI(Cl, 14))
  res <- NULL
  for (i in 42:(n-1)) {
    if (macd[i-4]<macdsig[i-4] && macd[i]>macdsig[i] && (TF30(i-4,i,rsi$rsi) || (rsi$rsi[i]<50 && rsi$rsi[i]-rsi$rsi[i-3]>15))) {
      signal <- "buy"
    } ## else if #<# #>#
    else if (macd[i-4]>macdsig[i-4] && macd[i]<macdsig[i] && (TF70(i-4,i,rsi$rsi) || (rsi$rsi[i]>50 && rsi$rsi[i-3]-rsi$rsi[i]>15)) && strat_out(x, i)) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame( date = date[i+1],
                                  option = signal,
                                  price = Op[i+1]))
  }
  return(res)
}

# RSI signal function, change the parameter in necessary 
TF70 <- function(i, j, x){
  for(a in i:j){
    if(x[a] > 65){
      return (1)
    }else {
      return (0)
    }
  }
}


TF30 <- function(i, j, x){
  for(a in i:j){
    if(x[a] < 35){
      return (1)
    }else {
      return (0)
    }
  }
}

# macd & parabolic stop & reverse
macdsar <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 35) stop("x should be longer than 35.")
  MACDS <- MACD(Cl, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA")
  macd <- MACDS[, "macd"]
  macdsig <- MACDS[, "signal"]
  sar <- SAR(cbind(Hi(x), Lo(x)), accel = c(0.02, 0.2))
  res <- NULL
  for (i in 42:(n-1)) {
    if (macd[i-4]<macdsig[i-4] && macd[i]>macdsig[i] && (sar[i-4]>Cl[i-4] && sar[i]<Cl[i])) {
      signal <- "buy"
    } else if (macd[i-4]<macdsig[i-4] && macd[i]>macdsig[i] && (sar[i-4]>Cl[i-4] && sar[i]>Cl[i])) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame( date = date[i+1],
                                  option = signal,
                                  price = Op[i+1]))
  }
  return(res)
}

# smi/stc, & rsi
strat_smirsi <- function(x, n = 14, nFast = 2, nSlow = 25, nSig = 9) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- Cl(x)
  Op <- Op(x)
  nlen <- length(Cl)
  if (nlen < (nSlow+nSig)) stop(paste("x should be longer than", (nSlow+nSig), "."))
  smi <- SMI(HLC(x), n = n, nFast = nFast, nSlow = nSlow, nSig = nSig)
  smil <- smi[, 1]
  smisignal <- smi[, 2]
  rsi <- RSI(Cl, n = n)
  res <- NULL
  for (i in (nSlow+nSig+9):(nlen-1)) {
    if ((rsi[i] < 35)  && (smil[i-9] < smisignal[i-9] && smil[i] > smisignal[i])) {
      signal <- "buy"
    } else if (rsi[i] > 65 && (smil[i-4] > smisignal[i-4] && smil[i] < smisignal[i]))
    {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# stc & macd
strat_macdsmi <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 35) stop("x should be longer than 35.")
  MACDS <- MACD(Cl, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA")
  macd <- MACDS[, "macd"]
  macdsig <- MACDS[, "signal"]
  smi <- SMI(HLC(x), n = 13, nFast = 2, nSlow = 25, nSig = 9)
  smil <- smi[, 1]
  smisignal <- smi[, 2]
  res <- NULL
  for (i in 42:(n-1)) {
    if (macd[i-6]<macdsig[i-6] && macd[i]>macdsig[i] && (smil[i-2] < smisignal[i-2] && smil[i] > smisignal[i])) {
      signal <- "buy"
    } else if (macd[i-6]>macdsig[i-6] && macd[i]<macdsig[i] && (smil[i-2] > smisignal[i-2] && smil[i] < smisignal[i])) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame( date = date[i+1],
                                  option = signal,
                                  price = Op[i+1]))
  }
  return(res)
}

# ema & stc
# sma / ema / wma / 
# bollinger's bands & rsi
strat_bbrsi <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 20) stop("x should be longer than 20.")
  rsi <- RSI(Cl, 14)
  bbands <- BBands(x[,2:4], n = 20)
  res <- NULL
  # from 51
  for (i in 35:(n-1)) {
    if (rsi[i] < 35 && Cl[i - 2] > bbands$dn[i - 2] && Cl[i] < bbands$dn[i]) {
      signal <- "buy"
    } else if (rsi[i] > 65 && (Cl[i - 2] < bbands$up[i - 2] && Cl[i] > bbands$up[i]) && strat_out(x, i)) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# dmi & stc
strat_dmismi <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 35) stop("x should be longer than 35.")
  smi_result <- SMI(HLC(x), n = 13, nFast = 2, nSlow = 25, nSig = 9)
  smi <- smi_result[, "SMI"]
  smisignal <- smi_result[, "signal"]
  adx <- ADX(HLC(x))
  DIp <- adx[, "DIp"]
  DIn <- adx[, "DIn"]
  res <- NULL
  for (i in 40:(n-1)) {
    if ((smi[i - 4] < smisignal[i - 4] && smi[i] > smisignal[i]) && (DIp[i - 4] < DIn[i - 4] && DIp[i] > DIn[i])) {
      signal <- "buy"
    } else if ((smi[i - 4] > smisignal[i - 4] && smi[i] < smisignal[i]) && (DIp[i - 4] > DIn[i - 4] && DIp[i] < DIn[i])) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame(date = date[i+1],
                                 option = signal,
                                 price = Op[i+1]))
  }
  return(res)
}

# dmi & macd/ema   
strat_dmimacd <- function(x) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'") } else {
      date <- index(x)
    }
  Cl <- Cl(x)
  Op <- Op(x)
  n <- length(Cl)
  if (n < 35) stop("x should be longer than 35.")
  MACDS <- MACD(Cl, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA")
  macd <- MACDS[, "macd"]
  macdsig <- MACDS[, "signal"]
  adx <- ADX(HLC(x))
  DIp <- adx[, "DIp"]
  DIn <- adx[, "DIn"]
  res <- NULL
  for (i in 42:(n-1)) {
    if ((macd[i-4]<macdsig[i-4] && macd[i]>macdsig[i]) && (DIp[i - 4] < DIn[i - 4] && DIp[i] > DIn[i])) {
      signal <- "buy"
    } else if ((macd[i-4]>macdsig[i-4] && macd[i]<macdsig[i]) && (DIp[i - 4] > DIn[i - 4] && DIp[i] < DIn[i])) {
      signal <- "sell"
    } else {
      next()
    }
    res <- rbind(res, data.frame( date = date[i+1],
                                  option = signal,
                                  price = Op[i+1]))
  }
  return(res)
}


###################################################
# win rate & return rate function

# perf function
# find the first "buy" signal and search for the first "sell" signal iteratively
perf <- function(record){
  start <- 1
  res <- NULL
  type_trade <- record[, 2]
  while (type_trade[start] != "buy"){
    start <- start + 1
    if(start >= length(record$option)){
      return(res)
    }
  }
  price = record[start, 3]
  res <- rbind(res, data.frame(date = record[start, 1],
                               option = record[start, 2],
                               price = price))
  for(i in start:length(record$option)){
    if(record[i,2]==record[start,2]){
      if(i == start || record[i-1, 2]==record[i, 2]){
        next()
      }
      else{
        price = record[i, 3]
        res <- rbind(res, data.frame(date = record[i, 1],
                                     option = record[i, 2],
                                     price = price))
      }
    } else if(record[i-1,2]==record[i,2]){
      next()
    } else {price = record[i, 3]
    res <- rbind(res, data.frame(date = record[i, 1],
                                 option = record[i, 2],
                                 price = price))
    }
  }
  return(res)
}

# ret function will return transaction details. a transaction consists of a buy and a sell;
# if the last row of record is "buy", which means there is no "sell" matches "buy".
# it will use the latest close price to represent "sell" price
# param, record is the return value of perf function above, stock should be the particular stock in calculating

# *****the below include the current stock value*****
# ret <- function(record, stock){
#   ret <- NULL
#   length <- nrow(record)
#   if(length != 1){
#     for (i in seq(1, length - 1, 2)){
#       buy_price = record[i, 3]
#       sell_price = record[i + 1, 3]
#       ret <- rbind(ret, data.frame(buy_date = record[i, 1],
#                                    sell_date = record[i + 1, 1],
#                                    buy_price = buy_price,
#                                    sell_price = sell_price,
#                                    change = sell_price - buy_price,
#                                    return = sell_price / buy_price - 1))
#     }}
#   # consider the current value of stocks, so here import the stock
#   if(length%%2==1 && record[length, 1]!=index(stock)[length(index(stock))]){
#     buy_price = record[length, 3]
#     sell_price = data.frame(stock[length(stock[,4]),4])[,1]
#     ret <- rbind(ret, data.frame(buy_date = record[length, 1],
#                                  sell_date = as.Date(rownames(data.frame(stock[length(stock[,4]),4]))),
#                                  buy_price = buy_price,
#                                  sell_price = sell_price,
#                                  change = sell_price - buy_price,
#                                  return = sell_price / buy_price - 1))
#   }
#   return(ret)
# }

# ****the below exclude the current value, so the param stock overlooked
ret <- function(record){
  ret <- NULL
  length <- nrow(record)
  if(length != 1){
    for (i in seq(1, length - 1, 2)){
      buy_price = record[i, 3]
      sell_price = record[i + 1, 3]
      return = sell_price / buy_price - 1
      buy_date = record[i, 1]
      sell_date = record[i + 1, 1]
      date_change = as.numeric(sell_date) - as.numeric(buy_date)
      ret <- rbind(ret, data.frame(buy_date = buy_date,
                                   sell_date = sell_date,
                                   buy_price = buy_price,
                                   sell_price = sell_price,
                                   change = sell_price - buy_price,
                                   return = return,
                                   ann_return = (return/date_change)*365))
    }}
  return(ret)
}

winrate_n_arbr <- NULL
winrate_dj_a <- NULL
winrate_h_arbr <- NULL
winrate_s_arbr <- NULL

winrate_n_pvtsma <- NULL
winrate_dj_pvtsma <- NULL
winrate_h_pvtsma <- NULL
winrate_s_pvtsma <- NULL

winrate_n_pvt <- NULL
winrate_dj_pvt <- NULL
winrate_h_pvt <- NULL
winrate_s_pvt <- NULL

winrate_n_pvtstc <- NULL
winrate_dj_pvtstc <- NULL
winrate_h_pvtstc <- NULL
winrate_s_pvtstc <- NULL

winrate_n_pvtmacd <- NULL
winrate_dj_pvtmacd <- NULL
winrate_h_pvtmacd <- NULL
winrate_s_pvtmacd <- NULL

winrate_n_obv <- NULL
winrate_dj_obv <- NULL
winrate_h_obv <- NULL
winrate_s_obv <- NULL

winrate_n_obvmacd <- NULL
winrate_dj_obvmacd <- NULL
winrate_h_obvmacd <- NULL
winrate_s_obvmacd <- NULL

winrate_n_obvsma <- NULL
winrate_dj_obvsma <- NULL
winrate_h_obvsma <- NULL
winrate_s_obvsma <- NULL
# winrate_sp <- NULL
# re <- NULL

# param: stock pool, stock name, win rate table pre-defined, trading strategy
final <- function(stock_pool, stock_name, winrate){
  length <- length(stock_pool)
  num <- length
  # please replace stock_pool to dj30_pool or nasdaq_pool etc.
  for(i in 1:length){
    print(paste(Sys.time(), "is working on stock", i, names(stock_pool[i])))
    stock <- na.omit(data.frame(stock_pool[i]))
    stock <- na.omit(stock[which(stock[, 1] > 0), ])
    stock <- na.omit(stock[which(stock[, 5] > 0), ])
    stock <- xts(stock, order.by = as.Date(rownames(stock)))
    iterim <- ARBR(stock)
    # whether there is no signal or only 1 signal (only 1 "buy" signal)
    if(is.null(iterim) || nrow(perf(iterim))==1 || is.null(perf(iterim))){
      print(paste("the stock: ", stock_name[i], "has no full transaction signals under current selected trading strategy"))
      num <- num - 1
      winrate <- rbind(winrate, data.frame(win_num = 0,
                                           trans_num = 0,
                                           winrate = 0,
                                           ann_return = 0))
      if(i == length){
        winrate <- rbind(winrate, data.frame(
                        win_num = sum(winrate$win_num),
                        trans_num = sum(winrate$trans_num),
                        winrate = sum(winrate$win_num) / sum(winrate$trans_num),
                        ann_return = sum(winrate$ann_return) / num
                      ))
        break
      }
      next()
    }
    re <- ret(perf(iterim))
    # re <- ret(perf(iterim), stock) # consider current value, even there is no sell signal now
    
    win_num <- length(which(re$return > threshold))  #threshold = 0.05
    ann_return <- mean(re$ann_return)
    trans_num <- nrow(re)
    winrate <- rbind(winrate, data.frame(win_num = win_num,
                                         trans_num = trans_num,
                                         winrate = win_num / trans_num,
                                         ann_return = ann_return))
    if(i == length){
      winrate <- rbind(winrate, data.frame(
        win_num = sum(winrate$win_num),
        trans_num = sum(winrate$trans_num),
        winrate = sum(winrate$win_num) / sum(winrate$trans_num),
        ann_return = sum(winrate$ann_return) / num
      ))
    }
  }
  # add name to winrate
  rownames(winrate) <- c(stock_name, "overall")
  return(winrate)
}
# the 5 winrate benchmarks are resulted by macdrsi
winrate_n_arbr <- final(nasdaq_pool, nasdaq, winrate_n_arbr)
winrate_dj_a <- final(dj30_pool, dj30, winrate_dj_a)
winrate_h_arbr <- final(hsi_pool, hsi, winrate_h_arbr)
winrate_s_arbr <- final(ss_pool, ss, winrate_s_arbr)
# winrate_sp <- final(sp500_pool, sp500, winrate_sp)


###################################################
# execute your strategies below
# iterim <- t4
# if(is.null(iterim) || nrow(perf(iterim))==1){
#   print(paste("the stock: ", iterim, "has no full transaction signals under current selected trading strategy"))
# }

test_d <- NULL
test_n <- NULL
test_h <- NULL
test_s <- NULL
test_sp <- NULL

# bb_rsi 
test_d <- final(dj30_pool, dj30, test_d)
test_n <- final(nasdaq_pool, nasdaq, test_n)
test_h <- final(hsi_pool, hsi, test_h)
test_s <- final(ss_pool, ss, test_s)
test_sp <- final(sp500_pool, sp500, test_sp) 

# macdrsi, hk & ss & nasdaq test_n1
test_h1 <- NULL
test_s1 <- NULL
test_n1 <- NULL
test_d1 <- NULL
test_sp1 <- NULL
test_h1 <- final(hsi_pool, hsi, test_h1)
test_s1 <- final(ss_pool, ss, test_s1)
test_d1 <- final(dj30_pool, dj30, test_d1)
test_n1 <- final(nasdaq_pool, nasdaq, test_n1)
test_sp1 <- final(sp500_pool, sp500, test_sp1)

### sell out conditions
strat_out <- function(x, i) {
  if (!any(class(x) %in% c("xts", "zoo"))) {
    stop("The class of input should be 'xts' or 'zoo.'")
  } else {
    date <- index(x)
  }
  Cl <- data.frame(Cl(x))[,1] #close price
  Op <- data.frame(Op(x))[,1] #open price
  Hi <- data.frame(Hi(x))[,1] #high price
  Lo <- data.frame(Lo(x))[,1] #low price
  n <- length(Cl) 
  if (i < 35) stop("x should be longer than 34.")
  # 断头铡刀, 阴线包住三根均线，5，10，30日均线
  sma5 <- SMA(Cl, n = 5)
  sma10 <- SMA(Cl, n = 10)
  sma30 <- SMA(Cl, n = 30)
  count <- rep(0, 11)
  # 阴烛
  if (Cl[i] < Op[i]){
    # 断头铡刀
    if (sma5[i] <= Op[i] && sma10[i] <= Op[i] && sma30[i] <= Op[i]
        && sma5[i] >= Cl[i] && sma10[i] >= Cl[i] && sma30[i] >= Cl[i]) {
      count[1] <- count[1] + 1
    } 
    if (Cl[i-1]>Op[i-1] && (Op[i]-Cl[i])>(Cl[i-1]-Op[i-1])){ # 阴包阳
      count[2] <- count[2] + 1
    } 
    if (Cl[i-1]>Op[i-1] && Op[i] > Cl[i-1]){  #乌云盖顶，前一日为阳烛，今日高开低走
      count[3] <- count[3] + 1
    } 
    if (abs(Cl[i-4]-Cl[i])<1.5 && Op[i-3] > Cl[i-3] && Op[i-3] < min(Cl[i-4], Op[i-4])){ # 相差小于1.5认为股价回复到原来水平
      count[4] <- count[4] + 1
    } 
    if (Cl[i-1] > Op[i-1] && (Cl[i-1]-Op[i-1]) > (Op[i]-Cl[i])){  # 高位孕线，前一日为阳烛，阳烛的范围大于今日阴烛的范围
      count[5] <- count[5] + 1
    } 
    if (Cl[i-2] > Op[i-2] && Cl[i-1] > Op[i-1] && 
               (Hi[i-1]-Lo[i-1])>2*(Cl[i-1]-Op[i-1])){ # 黄昏之星，前两日为阳烛，第二日最高价与最低价的差明显大于收盘价开盘价之差，第三日即今日为阴烛
      count[6] <- count[6] + 1
    } 
    if (Cl[i-2] < Op[i-2] && Cl[i-1] < Op[i-1] && Cl[i-3] > Op[i-3]){  # 三只乌鸦，连着三日阴烛，而前一天为阳烛
      count[7] <- count[7] + 1
    } 
    if (Cl[i-1] > Op[i-1] && Op[i] >= Cl[i-1] && (Hi[i]-Op[i])>2*(Op[i]-Cl[i]) && (Op[i]-Cl[i])>(Cl[i]-Lo[i])){
      count[8] <- count[8] + 1
    }
  }
  # 阳烛
  if (Cl[i] >= Op[i]){
    # 吊颈线, 在上涨的情况下，开盘价相较最低价的差至少是开盘收盘价之差的2倍以上，才有卖出信号
    if (Cl[i] > Cl[i-1] && (Op[i]-Lo[i]) >= 2*(Cl[i]-Op[i])){
      count[9] <- count[9] + 1
    } 
    if (abs(Cl[i-4]-Cl[i])<1.5 && Op[i-3] > Cl[i-3] && Op[i-3] < min(Cl[i-4], Op[i-4])){ # 相差小于1.5认为股价回复到原来水平
      count[10] <- count[10] + 1
    } 
    if ((Hi[i]-Cl[i])>=2*(Cl[i]-Op[i]) && Cl[i-1] > Op[i-1] && Op[i] >= Cl[i-1] && (Cl[i]-Op[i])>(Op[i]-Lo[i])){  #射击之星
      count[11] <- count[11] + 1
    }
  }
  # print(count)
  if (sum(count) >= 2){return(1)}
  else {return(0)}
}