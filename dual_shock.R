# load packages
library("data.table")
library("tidyverse")
library("tseries")
library("quantmod")
library('Rsolnp')
library('msos')
library('tikzDevice')
library('xtable')

# set working directory
setwd("/Users/mac/Desktop/Research/Post-Shock Prediction/")

## Conoco Phillips
getSymbols('COP', from = "2000-01-01")
COP <- as.data.frame(COP)
COP <- COP %>% mutate(Date = rownames(COP))

## S&P 500
getSymbols('^GSPC', from = "1970-01-01")
GSPC <- as.data.frame(GSPC)
GSPC <- GSPC %>% mutate(Date = rownames(GSPC))

## Brent Crude prices
Brent_Crude <- read.csv("https://pkgstore.datahub.io/core/oil-prices/brent-daily_csv/data/d93216330ab2c27aae3d177b2f0f0921/brent-daily_csv.csv") %>%
  rename(Oil_Close = Price)

## WTI Crude prices
WTI_Crude <- read.csv("https://pkgstore.datahub.io/core/oil-prices/wti-daily_csv/data/c414c9d375ec3c8f9d7c276d866fb2a4/wti-daily_csv.csv") %>%
  rename(WTI_Close = Price)

## Gold Price
getSymbols('GC=F', from = "2000-01-01")
gold <- as.data.frame(`GC=F`)
gold <- na.omit(gold)
gold <- gold %>% mutate(Date = rownames(gold))

## Dollar Index
getSymbols('DX-Y.NYB', from = "2000-01-01")
USD <- as.data.frame(`DX-Y.NYB`)
USD <- na.omit(USD)
USD <- USD %>% mutate(Date = rownames(USD))

## 13-Week T-Bill
getSymbols('^IRX', from = "2000-01-01")
TB <- as.data.frame(IRX)
TB <- na.omit(TB)
TB <- TB %>% mutate(Date = rownames(TB))

## Volatility Index
getSymbols('^VIX', from = "2000-01-01")
VIX <- as.data.frame(VIX)
VIX <- na.omit(VIX)
VIX <- VIX %>% mutate(Date = rownames(VIX))


## inflation adjustment
getSymbols("CPIAUCSL", src = 'FRED')
avg.cpi <- apply.yearly(CPIAUCSL, mean)


inflation_adj <- as.numeric(avg.cpi['2020'])/avg.cpi

inflation_adj <- as.data.frame(inflation_adj)
colnames(inflation_adj) <- c("dollars_2020")
inflation_adj <- inflation_adj %>% mutate(year = 1947:2020)


# Data Preparation
COP_close <- COP %>% dplyr::select(COP.Close, Date) %>% rename(COP_Close = COP.Close)
GSPC_close <- GSPC %>% dplyr::select(GSPC.Close, Date) %>% rename(GSPC_Close = GSPC.Close)
USD_close <- USD %>% dplyr::select(`DX-Y.NYB.Close`, Date) %>% rename(USD_Close = `DX-Y.NYB.Close`)
TB_close <- TB %>% dplyr::select(IRX.Close, Date) %>% rename(TB_Close = IRX.Close)
VIX_close <- VIX %>% dplyr::select(VIX.Close, Date) %>% rename(VIX_Close = VIX.Close)

tom <- list(GSPC_close, WTI_Crude, USD_close, TB_close, VIX_close)
for (i in 1:length(tom)) {
  COP_close <- merge(COP_close, tom[[i]])
}

# response
Y <- COP_close$COP_Close[-1]

# data frame
COP_close <- data.frame(COP_close[-nrow(COP_close), ], Y)

#### Monday, March 17th, 2008

## March 17th; 1 day nowcast

# shock effect date
start <- which(COP_close$Date == "2008-03-14")
start_day_20080317 <- as.numeric(1:nrow(COP_close) == start)
COP_close <- COP_close %>% mutate(start_day_20080317 = start_day_20080317)
TS2 <- COP_close[(start - 30):start, ]
# inflation adjustment
TS2[, 2:8] <- TS2[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008]
m_COP_3_17 <- lm(Y ~ COP_Close + start_day_20080317 + GSPC_Close + WTI_Close +
                    USD_Close + TB_Close,
                 data = TS2)
alpha_3_17 <- summary(m_COP_3_17)$coef[3,1:2]
# shock-effects
alpha_3_17

#### 2008 shock effects

# shock effect date
start_09_08_08 <- which(COP_close$Date == "2008-09-08")
start_09_12_08 <- which(COP_close$Date == "2008-09-12")
start_09_26_08 <- which(COP_close$Date == "2008-09-26")
# three shocks
start_day_09_08_08 <- as.numeric(1:nrow(COP_close) %in% start_09_08_08)
start_day_09_12_08 <- as.numeric(1:nrow(COP_close) %in% start_09_12_08)
start_day_09_26_08 <- as.numeric(1:nrow(COP_close) %in% start_09_26_08)
COP_close <- COP_close %>% mutate(start_day_09_08_08 = start_day_09_08_08,
                                  start_day_09_12_08 = start_day_09_12_08,
                                  start_day_09_26_08 = start_day_09_26_08)
# time window
TS3 <- COP_close[which(COP_close$Date == "2008-08-26"):which(COP_close$Date == "2008-09-26"), ]
# adjust for inflation
TS3[, 2:8] <- TS3[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2008]
# AR(1)
m_COP_Sept_08 <- lm(Y ~ COP_Close + start_day_09_08_08 + start_day_09_12_08 +
                      start_day_09_26_08 + GSPC_Close + WTI_Close + USD_Close + TB_Close, data = TS3)
alpha_Sept_08 <- summary(m_COP_Sept_08)$coef[3:5,1:2]
cov2cor(vcov(m_COP_Sept_08)[3:5, 3:5])

# shock-effects
alpha_Sept_08

# test independence

# test independence of a sequence of covariance matrices
lrindcov <- function(cov, cols, v){
  # cols is the sub-columns (in list) you want to test
  # est is the sample covariance matrix
  # v is the degrees of freedom of the sample covariance matrix under H0
  # test
  # H0: Sigmaij = 0; HA: Sigmaij != 0

  # sub sample covariance matrix
  est <- cov[unlist(cols), unlist(cols)]
  # multiple sub
  subest <- lapply(cols, function(x) cov[x, x])
  # dimensions
  q <- length(unlist(cols))
  # dimensions of sub sample covariance matrices
  qs <- sapply(cols, FUN = function(x) length(x))
  # calculate test statistics: 2log(LR)
  teststat <- v * (sum(sapply(subest,
                              FUN = function(x) if (length(x) == 1) log(x) else logdet(x)))
                   - logdet(est))
  # output result
  list(tlogLR = teststat, pvalue = 1 - pchisq(q = teststat, df = prod(qs)))
}

# not significant => Independence
lrindcov(vcov(m_COP_Sept_08)[3:5, 3:5], cols = list(1, 2, 3), v = df.residual(m_COP_Sept_08))

xtable(vcov(m_COP_Sept_08)[3:5, 3:5], digits = 3)

#### Thursday, November 27, 2014

# shock effect date
start <- which(COP_close$Date == "2014-11-26")
start_day_20141127 <- as.numeric(1:nrow(COP_close) == start)
COP_close <- COP_close %>% mutate(start_day_20141127 = start_day_20141127)
# time window
TS4 <- COP_close[(start - 30):start,]
# adjust for inflation
TS4[, 2:8] <- TS4[, 2:8] * inflation_adj$dollars_2020[inflation_adj$year == 2014]
# AR(1)
m_COP_11_27_14 <- lm(Y ~ COP_Close + start_day_20141127 + GSPC_Close + WTI_Close + USD_Close + TB_Close,
                     data = TS4)
alpha_11_27_14 <- summary(m_COP_11_27_14)$coef[3,1:2]
# shock-effects
alpha_11_27_14


#### The March 9th, 2020 shock effect:

# shock effect date
start <- which(COP_close$Date == "2020-03-06")
start_day_20200309 <- as.numeric(1:nrow(COP_close) == start)
COP_close <- COP_close %>% mutate(start_day_20200309 = start_day_20200309)
# time window
TS1 <- COP_close[(start-30):(start), ]

# shock-effect estimate
m_COP_03_09_20 <- lm(Y ~ COP_Close + start_day_20200309 + GSPC_Close + WTI_Close +
                         USD_Close + TB_Close,
                     data = TS1)
alpha_03_09_20 <- summary(m_COP_03_09_20)$coef[3,1:2]
# shock-effects
alpha_03_09_20


## Shock effect estimators
estimates <- rbind(alpha_3_17, alpha_Sept_08, alpha_11_27_14)
estimates[, 2] <- estimates[, 2] ^ 2
colnames(estimates) <- c("alpha_hat", "var")
rownames(estimates) <- c("m2008","s8y2008","s12y2008","s26y2008","y2014")


# adjustment estimator
alpha_adj <- mean(estimates[, 1])

# IVW estimator
weights <- (1 / estimates[,2]) / sum(1 / estimates[, 2])
alpha_IVW <- sum(weights * estimates[, 1])

# weighted adjustment estimator
Tstar.Date <- c("2020-03-05"
                , "2008-03-13"
                , "2008-09-05"
                , "2008-09-11"
                , "2008-09-25"
                ,  "2014-11-25")
Tstar <- sapply(Tstar.Date, function(x) which(COP_close$Date == x))
# X1
X1 <- as.matrix(TS1[nrow(TS1), 3:7])
# X1 <- as.matrix(COP_close[c(Tstar[1], Tstar[1] + 1), c(3, 4)])
# X0
X0 <- c()
for (i in 1:5) {
 X0[[i]] <- as.matrix(COP_close[Tstar[i + 1] + 1, 3:7])
}

# SCM

dat <- scale(rbind(X1, do.call('rbind', X0)), center = T, scale = T)
X1 <- dat[1, , drop = FALSE]
X0 <- c()
for (i in 1:5) {
  X0[[i]] <- dat[i + 1, , drop = FALSE]
}

# Euclidean metric
# objective function

scmm <- function(X1, X0) {
  weightedX0 <- function(W) {
    # W is a vector of weight of the same length of X0
    n <- length(W)
    p <- ncol(X1)
    XW <- matrix(0, nrow = 1, ncol = p)
    for (i in 1:n) {
      XW <- XW + W[i] * X0[[i]]
    }
    norm <- as.numeric(crossprod(matrix(X1 - XW)))
    return(norm)
  }
  # constraint for W
  Wcons <- function(W) sum(W) - 1
  n <- length(X0)
  # optimization
  outs <- solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, eqB = 0, LB = rep(0, n), UB = rep(1, n))

  # output weights
  Wstar <- outs$pars

  return(Wstar)
}

# objective function is not 0; the fit may not be good
Wstar <- scmm(X1 = X1, X0 = X0)

weightedX0 <- function(W) {
  # W is a vector of weight of the same length of X0
  n <- length(W)
  p <- ncol(X1)
  XW <- matrix(0, nrow = 1, ncol = p)
  for (i in 1:n) {
    XW <- XW + W[i] * X0[[i]]
  }
  norm <- as.numeric(crossprod(matrix(X1 - XW)))
  return(norm)
}
# constraint for W
Wcons <- function(W) sum(W) - 1
n <- length(X0)
# optimization
outs <- solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, eqB = 0, LB = rep(0, n), UB = rep(1, n))

# output weights
Wstar <- outs$pars
weightedX0(round(outs$pars, digits = 3))

# weighted adjustment estimator
alpha_wadj <- sum(Wstar * estimates[,1])

# Parametric Bootstrap Estimation

# Set a seed
set.seed(2020)
# Bootstrap replications
B <- 1000
# List of linear models
lmod <- list(m_COP_3_17, m_COP_Sept_08, m_COP_11_27_14)
# List of Data
TS <- list(TS2, TS3, TS4)
# List of T*
Tstar.Date <- c("2020-03-05", "2008-03-13", "2008-09-05", "2008-09-11", "2008-09-25",  "2014-11-25")
# Empty List for storation
alphas <- vector(mode = 'list', length = 3)
# Loop begins
for (b in 1:B) {

  # Vector for storing alpha hats
  alphahatsb <- c()
  # Weights for IVW Estimator
  weights <- c()

  for (i in 1:3) {
    # preparation
    res <- residuals(lmod[[i]])
    dat <- TS[[i]]
    Ti <- nrow(dat)
    coef <- matrix(coef(lmod[[i]]), nrow = 1)

    # BOOTSTRAP
    resb <- sample(res, size = Ti, replace = TRUE)

    # New response
    yi0 <- dat$COP_Close[1]
    yib <- yi0

    # Case by Case
    if (i == 2) {
      # i = 2 has 3 dates
      Tstari <- which(dat$Date %in% Tstar.Date[c(i + 1, i + 2, i + 3)])

      for (t in 1:Ti) {
        datt <- matrix(c(1, yib[t],
                         ifelse(t == Tstari+ 1, yes = 1, no = 0),
                         as.numeric(dat[t, c('GSPC_Close', 'WTI_Close',
                                             'USD_Close',
                                             'TB_Close')])))
        yib <- c(yib, resb[t] + coef %*% datt)
      }

      # Prepare for new data
      yb <- yib[-1]; yblag <- yib[-(Ti + 1)]
      datbi <- data.frame(yblag,
                          ifelse(1:Ti == Tstari[1] + 1, yes = 1, no = 0),
                          ifelse(1:Ti == Tstari[2] + 1, yes = 1, no = 0),
                          ifelse(1:Ti == Tstari[3] + 1, yes = 1, no = 0),
                          dat[,  c('GSPC_Close', 'WTI_Close',
                                   'USD_Close',
                                   'TB_Close')])
      # New colnames
      colnames(datbi) <- c('yblag', 'shock1', 'shock2', 'shock3',  c('GSPC_Close', 'WTI_Close',
                                                                     'USD_Close',
                                                                     'TB_Close'))

      # New Linear Model
      lmodbi <- lm(yb ~ 1 + yblag + shock1 + shock2 + shock3 + GSPC_Close + WTI_Close
                    + USD_Close + TB_Close, dat = datbi)

      # 3 Shock Effects
      alphahatsb <- c(alphahatsb, coef(lmodbi)[3:5])

      # IVW Weights
      weights <- c(weights, 1 / summary(lmodbi)$coef[3:5, 2] ^ 2)
    } else {

      # i = 1 differs from i = 3
      if (i == 1) {
        Tstari <- which(dat$Date == Tstar.Date[i + 1])
      } else {
        Tstari <- which(dat$Date == Tstar.Date[i + 3])
      }


      for (t in 1:Ti) {
        datt <- matrix(c(1, yib[t], ifelse(t == Tstari + 1, yes = 1, no = 0),
                         as.numeric(dat[t, c('GSPC_Close', 'WTI_Close',
                                             'USD_Close',
                                             'TB_Close')])))
        yib <- c(yib, resb[t] + coef %*% datt)
      }

      # Prepare for new data
      yb <- yib[-1]; yblag <- yib[-(Ti + 1)]
      datbi <- data.frame(yblag, ifelse(1:Ti == Tstari + 1, yes = 1, no = 0), dat[, c('GSPC_Close', 'WTI_Close',
                                                                                      'USD_Close',
                                                                                      'TB_Close')])

      # New colnames
      colnames(datbi) <- c('yblag', 'shock', c('GSPC_Close', 'WTI_Close',
                                               'USD_Close',
                                               'TB_Close'))

      # New Linear Model
      lmodbi <- lm(yb ~ 1 + yblag + shock + GSPC_Close + WTI_Close +
                       USD_Close + TB_Close, dat = datbi)

      # Shock Effects
      alphahatsb <- c(alphahatsb, coef(lmodbi)[3])

      # Weights
      weights <- c(weights, 1 / summary(lmodbi)$coef[3, 2] ^ 2)
    }
  }

  # Store Computed Shock-Effects Estimators
  alphas[[1]] <- c(alphas[[1]], mean(alphahatsb))
  alphas[[2]] <- c(alphas[[2]], sum(alphahatsb * weights / sum(weights)))
  alphas[[3]] <- c(alphas[[3]], sum(Wstar * alphahatsb))
}

# Parameters
means <- c()
vars <- c()
for (j in 1:3) {
  means <- c(means, mean(alphas[[j]]))
  vars <- c(vars, var(alphas[[j]]))
}
names(means) <- names(vars) <- c('adj', 'wadj', 'IVW')

# sample version
risk.reduction2 <- function(est, vars) {
  rr.adj <- (est['wadj']) ^ 2 - vars['adj'] - (est['adj'] - est['wadj']) ^ 2
  rr.wadj <- (est['wadj']) ^ 2 - vars['wadj']
  rr.IVW <- (est['wadj']) ^ 2 - vars['IVW'] - (est['IVW'] - est['wadj']) ^ 2
  rest <- c(rr.adj, rr.wadj, rr.IVW)
  names(rest) <- c('adj', 'wadj', 'IVW')
  return(list(usable = ifelse(rest > 0, yes = 1, no = 0), best = which.max(rest), rr = rest))
}
est <- c(alpha_adj, alpha_wadj, alpha_IVW)
names(est) <- c('adj', 'wadj', 'IVW')
risk.reduction2(est = est, vars = vars)

## Non-Parametric Bootstrap Estimation
# Parametric Bootstrap Estimation
# Set a seed
set.seed(2020)
# Bootstrap replications
B <- 1000
# List of linear models
lmod <- list(m_COP_3_17, m_COP_Sept_08, m_COP_11_27_14)
# List of Data
TS <- list(TS2, TS3, TS4)
# List of T*
Tstar.Date <- c("2020-03-05", "2008-03-13", "2008-09-05", "2008-09-11", "2008-09-25",  "2014-11-25")
# Empty List for storation
alphas <- vector(mode = 'list', length = 3)


sample(1:3, replace = T)

# Loop begins
for (b in 1:B) {

  # Vector for storing alpha hats
  alphahatsb <- c()
  # Weights for IVW Estimator
  weights <- c()

  for (i in 1:3) {
    # preparation
    res <- residuals(lmod[[i]])
    dat <- TS[[i]]
    Ti <- nrow(dat)
    coef <- matrix(coef(lmod[[i]]), nrow = 1)

    # BOOTSTRAP
    resb <- sample(res, size = Ti, replace = TRUE)

    # New response
    yi0 <- dat$COP_Close[1]
    yib <- yi0

    # Case by Case
    if (i == 2) {
      # i = 2 has 3 dates
      Tstari <- which(dat$Date %in% Tstar.Date[c(i + 1, i + 2, i + 3)])

      for (t in 1:Ti) {
        datt <- matrix(c(1, yib[t],
                         ifelse(t == Tstari+ 1, yes = 1, no = 0),
                         as.numeric(dat[t, c('GSPC_Close', 'WTI_Close',
                                             'USD_Close',
                                             'TB_Close')])))
        yib <- c(yib, resb[t] + coef %*% datt)
      }

      # Prepare for new data
      yb <- yib[-1]; yblag <- yib[-(Ti + 1)]
      datbi <- data.frame(yblag,
                          ifelse(1:Ti == Tstari[1] + 1, yes = 1, no = 0),
                          ifelse(1:Ti == Tstari[2] + 1, yes = 1, no = 0),
                          ifelse(1:Ti == Tstari[3] + 1, yes = 1, no = 0),
                          dat[,  c('GSPC_Close', 'WTI_Close',
                                   'USD_Close',
                                   'TB_Close')])
      # New colnames
      colnames(datbi) <- c('yblag', 'shock1', 'shock2', 'shock3',  c('GSPC_Close', 'WTI_Close',
                                                                     'USD_Close',
                                                                     'TB_Close'))

      # New Linear Model
      lmodbi <- lm(yb ~ 1 + yblag + shock1 + shock2 + shock3 + GSPC_Close + WTI_Close
                   + USD_Close + TB_Close, dat = datbi)

      # 3 Shock Effects
      alphahatsb <- c(alphahatsb, coef(lmodbi)[3:5])

      # IVW Weights
      weights <- c(weights, 1 / summary(lmodbi)$coef[3:5, 2] ^ 2)
    } else {

      # i = 1 differs from i = 3
      if (i == 1) {
        Tstari <- which(dat$Date == Tstar.Date[i + 1])
      } else {
        Tstari <- which(dat$Date == Tstar.Date[i + 3])
      }


      for (t in 1:Ti) {
        datt <- matrix(c(1, yib[t], ifelse(t == Tstari + 1, yes = 1, no = 0),
                         as.numeric(dat[t, c('GSPC_Close', 'WTI_Close',
                                             'USD_Close',
                                             'TB_Close')])))
        yib <- c(yib, resb[t] + coef %*% datt)
      }

      # Prepare for new data
      yb <- yib[-1]; yblag <- yib[-(Ti + 1)]
      datbi <- data.frame(yblag, ifelse(1:Ti == Tstari + 1, yes = 1, no = 0), dat[, c('GSPC_Close', 'WTI_Close',
                                                                                      'USD_Close',
                                                                                      'TB_Close')])

      # New colnames
      colnames(datbi) <- c('yblag', 'shock', c('GSPC_Close', 'WTI_Close',
                                               'USD_Close',
                                               'TB_Close'))

      # New Linear Model
      lmodbi <- lm(yb ~ 1 + yblag + shock + GSPC_Close + WTI_Close +
                     USD_Close + TB_Close, dat = datbi)

      # Shock Effects
      alphahatsb <- c(alphahatsb, coef(lmodbi)[3])

      # Weights
      weights <- c(weights, 1 / summary(lmodbi)$coef[3, 2] ^ 2)
    }
  }

  # Store Computed Shock-Effects Estimators
  alphas[[1]] <- c(alphas[[1]], mean(alphahatsb))
  alphas[[2]] <- c(alphas[[2]], sum(alphahatsb * weights / sum(weights)))
  alphas[[3]] <- c(alphas[[3]], sum(Wstar * alphahatsb))
}

# Parameters
means <- c()
vars <- c()
for (j in 1:3) {
  means <- c(means, mean(alphas[[j]]))
  vars <- c(vars, var(alphas[[j]]))
}
names(means) <- names(vars) <- c('adj', 'wadj', 'IVW')

# sample version
risk.reduction2 <- function(est, vars) {
  rr.adj <- (est['wadj']) ^ 2 - vars['adj'] - (est['adj'] - est['wadj']) ^ 2
  rr.wadj <- (est['wadj']) ^ 2 - vars['wadj']
  rr.IVW <- (est['wadj']) ^ 2 - vars['IVW'] - (est['IVW'] - est['wadj']) ^ 2
  rest <- c(rr.adj, rr.wadj, rr.IVW)
  names(rest) <- c('adj', 'wadj', 'IVW')
  return(list(usable = ifelse(rest > 0, yes = 1, no = 0), best = which.max(rest), rr = rest))
}
est <- c(alpha_adj, alpha_wadj, alpha_IVW)
names(est) <- c('adj', 'wadj', 'IVW')
risk.reduction2(est = est, vars = vars)


# All Estimators Are Usable !
# Vote for Weighted Adjustment Estimators !

## Additive adjustment

# adjustment estimator
alpha.adj.additive <- mean(estimates[1:4, 1]) + estimates[5, 1]

# IVW estimator
weights <- (1 / estimates[1:4,2]) / sum(1 / estimates[1:4, 2])
alpha.IVW.additive <- sum(weights * estimates[1:4, 1])
alpha.IVW.additive <- as.numeric(alpha.IVW.additive + estimates[5, 1])

# X1
X1 <- as.matrix(TS1[nrow(TS1), 3:7])
# X1 <- as.matrix(COP_close[c(Tstar[1], Tstar[1] + 1), c(3, 4)])
# X0
X0 <- c()
for (i in 1:4) {
  X0[[i]] <- as.matrix(COP_close[Tstar[i + 1] + 1, 3:7])
}

# SCM

dat <- scale(rbind(X1, do.call('rbind', X0)), center = T, scale = T)
X1 <- dat[1, , drop = FALSE]
X0 <- c()
for (i in 1:4) {
  X0[[i]] <- dat[i + 1, , drop = FALSE]
}

n <- 4
# optimization
outs <- solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, eqB = 0, LB = rep(0, n), UB = rep(1, n))
# output weights
Wstar.additive <- outs$pars
# objective function is not 0; the fit may not be good
weightedX0(round(outs$pars, digits = 3))

# Weighted Adjustment Estimator
alpha.wadj.additive <- as.numeric(sum(Wstar.additive * estimates[1:4,1]) + estimates[5, 1])

## Post-shock forecasts
m_COP_03_06_20 <- lm(Y ~ COP_Close + GSPC_Close + WTI_Close + USD_Close + TB_Close + VIX_Close,
                     data = TS1[-nrow(TS1), ])
# Yhat 1
Yhat_nothing <- coef(m_COP_03_06_20) %*% t(as.matrix(cbind(1, as.vector(TS1[nrow(TS1), 2:7]))))
# Yhat 2
Yhat_adj <- c(adj = alpha_adj + Yhat_nothing, IVW = alpha_IVW + Yhat_nothing, wadj = alpha_wadj + Yhat_nothing)
Yhat_additive <- c(adj = alpha.adj.additive + Yhat_nothing, IVW = alpha.IVW.additive + Yhat_nothing,
                   wadj = alpha.wadj.additive + Yhat_nothing)

## doing nothing completely misses the mark
Yhat_nothing - TS1$Y[nrow(TS1)]

## adjustment gets closer
Yhat_adj - TS1$Y[nrow(TS1)]

## additive effect does well
Yhat_additive - TS1$Y[nrow(TS1)]


# plot data
TS1$id <- 1:nrow(TS1)
mat <- cbind(TS1$id[nrow(TS1)], c(Yhat_adj))
colnames(mat) <- c("id", "Yhat_adj")
dat <- as.data.frame(mat)
dat$id <- c(30.4, 31, 31.6)
colnames(Yhat_nothing) <- "y"
# set working directory
setwd('/Users/mac/Desktop/Research/Post-Shock Prediction/')
# plot setting
tikz('fig2.tex', standAlone = TRUE, width = 7, height = 5)
# plot
ggplot(TS1, mapping = aes(x = id, y = Y)) +
  labs(title = "Conoco Phillips Stock Forecasting (2020 January 24th to 2020 March 9th)",
       x = "Day", y = "Closing Stock price (in USD)") +
  geom_point() +
  geom_point(data = dat, aes(x = id, y = Yhat_adj),
             col = c("magenta", "deepskyblue", "indianred2"),
             pch = 2:4, cex = 1.5) +
  geom_point(data = data.frame(x = unique(dat$id), y = Yhat_nothing),
             aes(x = x, y = y), col = c('white', "violet", 'white'), cex = 1.5) +
  geom_line(aes(x = id, y = c(m_COP_03_06_20$fitted.values, Yhat_nothing)),
            col = "violet") +
  annotate("text", x = 1, y = seq(from = 48, to = 40, length.out = 4),
           label = c("$\\hat{y}_{1, T_1^* + 1}^{1}$",
                     "$\\hat{y}_{1, T_1^* + 1}^{1} + \\hat{\\alpha}_{\\rm adj}$",
                     "$\\hat{y}_{1, T_1^* + 1}^{1} + \\hat{\\alpha}_{\\rm IVW}$",
                     "$\\hat{y}_{1, T_1^* + 1}^{1} + \\hat{\\alpha}_{\\rm wadj}$"),
           hjust = 0, size = 5) +
  annotate("point", x = 0, y = seq(from = 48, to = 40, length.out = 4),
           pch = c(16, 2:4),
           color = c("violet", "magenta", "deepskyblue", "indianred2"),
           size = 2) +
  # add margin
  theme(plot.margin = unit(c(.5, .3, .3, .5), "cm")) +
  # no grid
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
# minimal
  theme_minimal() +
  # center title
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
# output
dev.off()

