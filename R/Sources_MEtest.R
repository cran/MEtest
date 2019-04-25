me.test <- function(W, V, B = 1000, wt = c("Uniform", "Normal"), 
	wt.bd = NULL, wt.prob = 0.99, nGL = 32){
  ##
  W_in <- W[complete.cases(W),]
  V_in <- V[complete.cases(V),]
  
  nx <- nrow(W_in)
  mx <- ncol(W_in)
  ny <- nrow(V_in)
  my <- ncol(V_in)
  
  if(mx < 2L & my >= 2L){ 
    stop(paste0("Not enough replicates for W: mx = ", mx))
  } else if(mx >= 2L & my < 2L) {
    stop(paste0("Not enough replicates for V: my = ", my))
  } else if(mx < 2L & my < 2L){
    stop(paste0("Not enough replicates for both W and V: mx = ", mx, " and my = ", my))
  }
  
  if(nx < 2L & ny >= 2L){ 
    stop(paste0("Not enough data for W: nx = ", nx))
  } else if(nx >= 2L & ny < 2L) {
    stop(paste0("Not enough data for V: ny = ", ny))
  } else if(nx < 2L & ny < 2L){
    stop(paste0("Not enough data for both W and V: nx = ", nx, " and ny = ", ny))
  }
  
  ##
  DNAME <- paste(deparse(substitute(W)), "and", deparse(substitute(V)))
  
  
  ##
  wfun <- match.arg(wt)
  if(is.null(wt.bd)){
    Boot.dist <- Boot.dist(W_in, V_in, nboot = B, wt.bd = wt.bd, wt.prob = wt.prob, nGL = nGL)
    Test.stat <- TS.comp(W = W, V = V, bd = Boot.dist$bd, nGL = nGL)	
  } else{
    Boot.dist <- Boot.dist(W_in, V_in, nboot = B, wt.bd = wt.bd, wt.prob = wt.prob, nGL = nGL)
    Test.stat <- TS.comp(W = W, V = V, bd = wt.bd, nGL = nGL)	
  }
  
	pval.tmp <- NULL
	for(k in 1:2){
		pval.tmp <- c(pval.tmp, sum(Boot.dist$Boot[,k] > Test.stat[k], na.rm = TRUE)/
					length(Boot.dist$Boot[!is.na(Boot.dist$Boot[,k]),k]))
	}
	
	METHOD <- switch(wfun, 
	                 Uniform = "Homogeneity test under Measurement Error with Uniform weight",
	                 Normal = "Homogeneity test under Measurement Error with Normal weight")
	Tstat <- switch(wfun, Uniform = Test.stat[1], Normal = Test.stat[2])
	names(Tstat) <- switch(wfun, Uniform = "T (Uniform)", Normal = "T (Normal)")
	PVAL <- switch(wfun, Uniform = pval.tmp[1], Normal = pval.tmp[2])
	
	re <- list(statistic = Tstat, p.value = PVAL, method = METHOD, data.name = DNAME, 
	           alternative = "True signals come from different distributions",
	           boundary = wt.bd)
	class(re) <- "htest"
	return(re)
}

### Calculating Test statistic via Fortran 90
TS.comp <- function(W, V, bd = c(qnorm(0.005), qnorm(0.995)), nGL = 32){

	nx <- nrow(W); mx <- ncol(W)
	ny <- nrow(V); my <- ncol(V)
	Nx <- nx*mx*(mx-1)/2; Ny <- ny*my*(my-1)/2
	GL1 <- gauss.quad(n = nGL, "legendre")
	GL <- cbind(GL1$nodes, GL1$weights)

	TS.f90 <- .Fortran("TScomp", W = as.double(W), V = as.double(V), nx = as.integer(nx), ny = as.integer(ny), 
					mx = as.integer(mx), my = as.integer(my), GL = as.double(GL), nGL = as.integer(nGL), 
					BDs = as.double(bd), TS = as.double(rep(0,2)))
	return(TS.f90$TS)
}

## Boot.dist gives the bootstrap distribution
Boot.dist <- function(W, V, nboot = 1000, wt.bd = NULL, wt.prob = 0.99, nGL = 32){

	nx <- nrow(W); ny <- nrow(V)
	mx <- ncol(W); my <- ncol(V)
	GL1 <- gauss.quad(n = nGL, "legendre")
	GL <- cbind(GL1$nodes, GL1$weights)

	## 1. Finding optimal Bandwidth for True signal
	h.min.W <- h.min.V <- 0.07; h.max.W <- h.max.V <- 0.25 
	nseq = 300 
	Wseq <- seq(h.min.W, h.max.W, length = nseq)
	Vseq <- seq(h.min.V, h.max.V, length = nseq)
	AMISE.f90 <- .Fortran("AMISE", W = as.double(W), V = as.double(V), nx = as.integer(nx), ny = as.integer(ny), 
						mx = as.integer(mx), my = as.integer(my), GL = as.double(GL), nGL = as.integer(nGL), 
						Wseq = as.double(Wseq), Vseq = as.double(Vseq), nseq = as.integer(nseq), 
						Wobj = as.double(rep(0, nseq)), Vobj = as.double(rep(0, nseq)), 
						varx = as.double(0), vary = as.double(0), Wdiffvar = as.double(0), Vdiffvar = as.double(0))
	h.opt.W <- Wseq[which.min(AMISE.f90$Wobj)]	
	h.opt.V <- Vseq[which.min(AMISE.f90$Vobj)]


	## 2. Drawing from the common distribution
	zsd1 <- zsd2 <- 4 
	Wvar <- AMISE.f90$varx
	Vvar <- AMISE.f90$vary
	Wmean <- mean(as.vector(W))
	Vmean <- mean(as.vector(V))
	zseq <- seq(min(Wmean - zsd1*sqrt(Wvar), Vmean - zsd1*sqrt(Vvar)), 
				max(Wmean + zsd2*sqrt(Wvar), Vmean + zsd2*sqrt(Vvar)), length = nseq)
	Fhat.f90 <- .Fortran("Fhat", W = as.double(W), V = as.double(V), nx = as.integer(nx), ny = as.integer(ny), 
						mx = as.integer(mx), my = as.integer(my), GL = as.double(GL), nGL = as.integer(nGL), 
						zseq = as.double(zseq), nseq = as.integer(nseq), bwx = as.double(h.opt.W),
						bwy = as.double(h.opt.V), Fx = as.double(rep(0, nseq)), Fy = as.double(rep(0, nseq)))
	Fx.mon <- 0; Fy.mon <- 0

	Fx.mon[1:which.min(Fhat.f90$Fx)] <- min(Fhat.f90$Fx)
	Fy.mon[1:which.min(Fhat.f90$Fy)] <- min(Fhat.f90$Fy)
	for(i in (which.min(Fhat.f90$Fx)+1):nseq) Fx.mon[i] <- max(Fhat.f90$Fx[which.min(Fhat.f90$Fx):i], na.rm = TRUE)
	for(i in (which.min(Fhat.f90$Fy)+1):nseq) Fy.mon[i] <- max(Fhat.f90$Fy[which.min(Fhat.f90$Fy):i], na.rm = TRUE)
	F.com <- (nrow(W)*Fx.mon + nrow(V)*Fy.mon)/(nrow(W)+nrow(V))
	F.com[nseq] <- 1
	

	## 3. Drawing from the error distribution
	errsd <- 4
	errseq <- seq(min(-errsd*sqrt(AMISE.f90$Wdiffvar), -errsd*sqrt(AMISE.f90$Vdiffvar)), 
				max(errsd*sqrt(AMISE.f90$Wdiffvar), errsd*sqrt(AMISE.f90$Vdiffvar)), length = nseq)
	Uhat.f90 <- .Fortran("Uhat", W = as.double(W), V = as.double(V), nx = as.integer(nx), 
						ny = as.integer(ny), mx = as.integer(mx), my = as.integer(my), GL = as.double(GL), 
						nGL = as.integer(nGL), errseq = as.double(errseq),	nseq = as.integer(nseq), 
						Ux = as.double(rep(0, nseq)), Uy = as.double(rep(0, nseq)))
	Ux.mon <- 0; Uy.mon <- 0
	for(i in 1:nseq){ 
		Ux.mon[i] <- max(Uhat.f90$Ux[1:i], na.rm = TRUE)
		Uy.mon[i] <- max(Uhat.f90$Uy[1:i], na.rm = TRUE)
	}
	Ux.mon[nseq] <- 1
	Uy.mon[nseq] <- 1
	
	if(is.null(wt.bd)){
	  min.prob <- (1-wt.prob)/2
	  bd <- c( min(max(zseq[Fx.mon <= min.prob]), max(zseq[Fy.mon <= min.prob])),
	           max(max(zseq[Fx.mon <= wt.prob+min.prob]), max(zseq[Fy.mon <= wt.prob+min.prob])))
	} else{
	  bd <- wt.bd
	}
	
	BOOT.sam.f90 <- .Fortran("BootGen", nx = as.integer(nx), ny = as.integer(ny), mx = as.integer(mx), my = as.integer(my), 
						nboot = as.integer(nboot), F = as.double(F.com), Ux = as.double(Ux.mon), Uy = as.double(Uy.mon), 
						zseq = as.double(zseq), errseq = as.double(errseq), nseq = as.integer(nseq), 
						Xb = as.double(rep(0, nx*nboot)), Yb = as.double(rep(0, ny*nboot)), 
						Xerr = as.double(rep(0, nx*mx*nboot)), Yerr = as.double(rep(0, ny*my*nboot)))
						
	BOOT.dist.f90 <- .Fortran("BootDist", nx = as.integer(nx), ny = as.integer(ny), mx = as.integer(mx), 
						my = as.integer(my), nboot = as.integer(nboot), Xb = as.double(BOOT.sam.f90$Xb), 
						Yb = as.double(BOOT.sam.f90$Yb), Xerr = as.double(BOOT.sam.f90$Xerr), 
						Yerr = as.double(BOOT.sam.f90$Yerr), GL = as.double(GL), nGL = as.integer(nGL), 
						BDs = as.double(bd), Bdist = as.double(matrix(0, nrow = nboot, ncol = 2)))
	Boot.dist <- matrix(BOOT.dist.f90$Bdist, nrow = nboot, ncol = 2)
	return(list(Boot = Boot.dist, bd = bd))
}
