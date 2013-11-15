
## @knitr unnamed-chunk-1
library(FLa4av2)
library(plyr)
library(Hmisc)
library(multicore)
load("stocks.RData")
source("utilities.R")


## @knitr unnamed-chunk-2
sessionInfo()


## @knitr unnamed-chunk-3
opts_chunk$set(fig.align='center', fig.pos='H', fig.keep="last", cache=TRUE)


## @knitr unnamed-chunk-4
plot(FLStocks(lapply(stks01, "[[", "stock")))


## @knitr unnamed-chunk-5
set.seed(239246)
stks01 <- mclapply(stks01, genObs)


## @knitr unnamed-chunk-6
fmodel <- ~ bs(age, 4) + factor(year)
qmodel <- list(~ factor(age))


## @knitr unnamed-chunk-7
fits01 <- mclapply(stks01, doFits, fmodel = fmodel, qmodel = qmodel)


## @knitr unnamed-chunk-8
doPlots(fits01)


