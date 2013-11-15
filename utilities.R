
# function to extract results
outStats <- function(object) {	### broken for now as fit is no lonker stock object ###
	if(is(object, "try-error")) {
	    object
	} else {
		fit <- object[["fit"]]
		sim <- object[["sim"]]
		Fobs <- c(harvest(sim))
		Nobs <- c(stock.n(sim))
		Cobs <- c(catch.n(sim))
		Iobs <- rnorm(length(Cobs), 0, 1)
		Ffit <- c(harvest(fit))
		Nfit <- c(stock.n(fit))
		Cfit <- c(catch.n(fit))
		Ifit <- rnorm(length(Cobs), 0, 1)
		vars <- rep(c("F","N","C","I"), rep(length(Cobs),4))
		out <- data.frame( expand.grid(age = dimnames(harvest(fit))[[1]], year = dimnames(harvest(fit))[[2]]), 
		                   obs = c(Fobs, Nobs, Cobs, Iobs), 
		                   fit = c(Ffit, Nfit, Cfit, Ifit), 
		                   var = vars)
		out $ res <- out $ obs - out $ fit
		lst2 <- lapply(split(out, out[,c("age","var")]), transform, resScl = scale(res))
		do.call("rbind", lst2)
	}
}

summStats <- function(object){
	if(is(object, "try-error")) {
	    object
	} else {
		fit <- object[["fit"]]
		sim <- object[["sim"]]
		Fsim <- fbar(sim)
		Ffit <- fbar(fit)
		SSBsim <- ssb(sim)
		SSBfit <- ssb(fit)
		Csim <- catch(sim)
		Cfit <- catch(fit)
		Rsim <- rec(sim)
		Rfit <- rec(fit)
		data.frame(year = as.numeric(as.character(dimnames(Fsim)[[2]])), 
		           src  = c(rep("sim", prod(dim(Fsim))*4), rep("fit", prod(dim(Ffit))*4)), 
		           var  = rep(rep(c("F","S","R","C"), rep(prod(dim(Fsim)), 4)), 2), 
		           val  = c(Fsim, SSBsim, Rsim, Csim, Ffit, SSBfit, Rfit, Cfit))
	}
}

genObs <- function(x) {
  stock <- x $ stock	
	ages <- 1:min(9, range(stock)["max"])
	range(stock)[c("minfbar","maxfbar")] <- c(2, min(5, max(ages)))
	stock <- setPlusGroup(stock, max(ages))
  n <- stock.n(stock)
  z <- harvest(stock) + m(stock)
  logq <- -exp(-exp(0.2 * ages)) - 3 # trawl like catchability
  # observe index in 1st quarter with 10% cv
  index <- FLIndex(index = n * exp(-0.25 * z) * exp(logq + rnorm(prod(dim(n)), 0, .1)))  # 10% cv
  range(index)[c("startf","endf")] <- 0.25
  # observe catch with 10% cv
  catch.n(stock) <- catch.n(stock) * exp(rnorm(prod(dim(n)), 0, .1)) # 10% cv
  catch(stock) <- computeCatch( stock )
  list(stock = stock, index = list(index))
}

doFits <- function(x, fmodel, qmodel, rmodel = ~ factor(year)) {
  msg <- paste0("fitting ", x$stock@name, " in PID: ", Sys.getpid(), "\n")   
  cat(msg)
  fit <- try( a4aFit(fmodel, qmodel, rmodel, x$stock, x$index) )
  if (is(fit, "try-error")) {
    fit 
  } else {
    list(sim = x$stock, fit = fit)
  }
}

doPlots.old <- function(fits) {
  lst <- mclapply(fits, outStats)
  lst2 <- mclapply(fits, summStats)
  nms <- names(lst)
  for(i in 1:length(lst)) {
	  if(!is(lst[[i]], "try-error")) {
      print(xyplot(res ~ fit | var, groups = age, data = lst[[i]], scales=list(relation="free"), layout=c(4,1), main=nms[i], xlab="fitted", ylab="residuals", pch=19, cex=0.5))
      print(xyplot(resScl ~ fit | var, groups = age, data = lst[[i]], scales=list(relation="free"), layout=c(4,1), xlab="fitted", ylab="scaled residuals", pch=19, cex=0.5))
		  print(xyplot(val ~ year | var, groups = src, data = lst2[[i]], scales=list(y=list(relation="free")), type="l", layout=c(4,1)))
	  }
  }
}  


doPlots <- function(fits) 
{
  cols <- 
     c(rgb(215, 48, 39, max = 255), 
       rgb(252, 141, 89, max = 255),
       rgb(254, 224, 144, max = 255),
       rgb(224, 243, 248, max = 255),
       rgb(145, 191, 219, max = 255),
       rgb(69, 117, 180, max = 255))

  for(obj in fits) {
    if(is(obj, "try-error")) next

    x <- obj $ fit
    ages  <- as.numeric(dimnames(stock.n(x)) $ age)
    years <- as.numeric(dimnames(stock.n(x)) $ year)

    par(mfrow = c(1,4))

    # catch
    plotError(years, colSums((stock.wt(obj$sim) * exp(x @ catch.hat))[drop=TRUE]), 0,
              ylab = 'Catch', xlab = "Year", cols = colorRampPalette(cols[5:6])(3))
    lines(years, catch(obj $ sim)[drop=TRUE], col = "red", lwd = 2)

    # fplot
    plotError(years, fbar(x)[,,,,,"mean",drop=TRUE], sqrt(fbar(x)[,,,,,"var",drop=TRUE]), 
              ylab = 'Fbar', xlab = "Year", cols = colorRampPalette(cols[5:6])(3))
    lines(years, fbar(obj $ sim)[drop=TRUE], col = "red", lwd = 2)

    # rec plot
    Xr <- model.matrix(x @ models $ rmodel, list(year = years))
    r <- x @ coefficients[1:length(years)]
    rcov <- x @ covariance[1:length(years),1:length(years)]
    plotError(years, Xr %*% r, sqrt(diag(Xr %*% rcov %*% t(Xr))), FUN = exp,
              ylab = "Rectruitment", xlab = "Year", cols = colorRampPalette(cols[5:6])(3))
    lines(years, rec(obj $ sim)[drop=TRUE], col = "red", lwd = 2)

    # ssb plot
    plotError(years, ssb(x)[,,,,,"mean",drop=TRUE] * 1e-3, sqrt(ssb(x)[,,,,,"var",drop=TRUE]) * 1e-3, 
              ylab = "SSB ('000 tonnes)", xlab = "Year", cols = colorRampPalette(cols[5:6])(3))
    lines(years, ssb(obj $ sim)[drop=TRUE] * 1e-3, col = "red", lwd = 2)

    mtext(name(obj$sim), outer=TRUE)

  }
}  


  panel.3d.levelplot <-
  function(x, y, z, rot.mat, distance, zlim.scaled, at, drape = TRUE, shade = FALSE, ...) 
  {
    zrng <- 0.001 * diff(zlim.scaled) 
    ## vertical range of 2nd surface (>0)     
    z.scaled <- (z - min(z))/diff(range(z))     
    at.scaled <- (at - min(z))/diff(range(z)) 
    new.level <- zlim.scaled[2]    
    new.z <- new.level + zrng * (z.scaled)     
    new.at <- new.level + zrng * (at.scaled)     
    panel.3dwire(x, y, new.z, at = new.at,
                 col = grey(.9),
                 rot.mat = rot.mat, distance = distance,
                 shade = FALSE,
                 drape = TRUE,
                 zlim.scaled = zlim.scaled,
                 alpha = 1, ...)
    dots <- list(...)
    dots $ col.regions <- paste0(colorRampPalette(rev(cols))(100), "AA")
    dots $ background <- "transparent"
    args <- c(dots, list(x=x, y=y, z=z, rot.mat = rot.mat, distance = distance,
                 zlim.scaled = zlim.scaled, at = at, drape = FALSE, shade = FALSE))
    do.call(panel.3dwire, args)
  }
  cols <- 
   c(rgb(215, 48, 39, max = 255), 
    rgb(252, 141, 89, max = 255),
    rgb(254, 224, 144, max = 255),
    rgb(224, 243, 248, max = 255),
    rgb(145, 191, 219, max = 255),
    rgb(69, 117, 180, max = 255))
    
    
    

