library(ncdf4)
library(splines)
library(pbs)
library(quantreg)
# library(doParallel)
library(doMC)
# library(QRM)
# library(plyr)
# library(pryr)
# library(cluster)
# registerDoMC(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))


exeedence_var = function(y, yhat, breaks=30) {
	day_of_year = 1:365
	n = length(yhat) %/% (days_per_year * upper_year)
	res = y - yhat
	posres = res[res > 0]
	postime = rep(day_of_year, n*upper_year)[res > 0]
	h = hist(postime, breaks=breaks, plot=F)
	ti = h$counts[-length(h$counts)]
	T = sd(ti / sum(ti))
	T
}

printf <- function(...) invisible(print(sprintf(...)))

# x appears in blue and y appears in red.
comp.hist <- function(x,y, breaks=20, ...) {
	h.x = hist(x, breaks=breaks, plot=FALSE)
	h.y = hist(y, breaks=breaks, plot=FALSE)
	dots = list(...)
	dx = min(abs(h.x$breaks[2]-h.x$breaks[1]),
	         abs(h.y$breaks[2]-h.y$breaks[1]))

	ylim = range(c(h.x$density,h.y$density))
	u = range(c(h.x$breaks,h.y$breaks))
	breaks = seq(u[1],u[2]+dx,dx)
	if (max(breaks)<max(ylim)) breaks = c(breaks,max(breaks)+dx)
		if ("plot" %in% names(dots) && dots$plot == FALSE) {
			h.x = hist(x,breaks=breaks, plot=FALSE)	
		} else if ("xlim" %in% names(dots)) {
			h.x = hist(x,breaks=breaks,density=15,ylim=ylim,col="black",freq=FALSE, ...)
		} else {
			h.x = hist(x,breaks=breaks,density=15,ylim=ylim,col="black",xlim=u,freq=FALSE, ...)
		}

	# abline(v=median(x), col="blue")
	h.y = hist(y,breaks=breaks,density=15,col="red",add=TRUE,freq=FALSE,angle=-45)
	list(h.x=h.x, h.y=h.y)
}


print.memory.usage = function() {
	size = 0
	for (o in objects(envir=.GlobalEnv)) {
		size = size + object.size(get(o))
		message(o); print(object.size(get(o)), units='auto')
	}
	print(paste("Total memory usage:", size / 1000 / 1000, "Mb"))
}

plot.location.scale = function(data) {
	par(mfrow=c(2,1))
	plot(data$mu)
	plot(data$scale)
}


fit.2d.location.scale = function(t, x, y, q1, q2, days_per_year=365, n_years=250) {	
	n_days = (days_per_year * n_years)
	sample = seq(1, n_days, 1)
	rq.fit = rq(y ~ x*t, tau=c(0.25, 0.75, q1, q2), method="pfn")

	mu1 = matrix(rq.fit$fitted.values[sample,3], ncol=n_years)
	mu2 = matrix(rq.fit$fitted.values[sample,4], ncol=n_years)
	coef1 = rq.fit$coef[,3]
	coef2 = rq.fit$coef[,4]

	scale = mu2 - mu1
	list(mu=mu1, scale=scale, 
		coef_low=rq.fit$coef[,1], coef_high=rq.fit$coef[,2],
		coef1=coef1, coef2=coef2)
}

fit.location.scale = function(x, y, q1, q2, days_per_year=365) {	
	
	rq.fit1 = rq('y ~ .', tau=q1, method="pfn", data=data.frame(x=x, y=y))

	rq.fit2 = rq('y ~ .', tau=q2, method="pfn", data=data.frame(x=x, y=y))
	mu1 = rq.fit1$fitted.values[1:days_per_year]
	mu2 = rq.fit2$fitted.values[1:days_per_year]
	scale = mu2 - mu1
	data.frame(mu=mu1, scale=scale)
}

get.ncdf.timefirst = function(filename, lat_start, lat_count, lon_start, lon_count, t_start, t_count, varname='TREFHT') {
	ncdata <- nc_open(filename)	
		data = ncvar_get(ncdata, 
		varname, 
		start=c(t_start, lon_start, lat_start), 
		count=c(t_count, lon_count, lat_count))
	nc_close(ncdata)
	data
}
get.ncdf.light = function(filename, lat_start, lat_count, lon_start, lon_count, t_start, t_count) {
	ncdata <- nc_open(filename)	
		data = ncvar_get(ncdata, 
		'TREFHT', 
		start=c(lon_start, lat_start, t_start), 
		count=c(lon_count, lat_count, t_count))
	nc_close(ncdata)
	data
}

get.ncdf = function(filename, lat_min, lat_max, lon_min, lon_max) {
	ncdata <- nc_open(filename)	
	lons = ncvar_get(ncdata, 'lon') - 180
	lats = ncvar_get(ncdata, 'lat')
	lon_idx = which((lons < lon_max) & (lons > lon_min))
	lat_idx = which((lats < lat_max) & (lats > lat_min))
	data = ncvar_get(ncdata, 
		'TREFHT', 
		start=c(lon_idx[1], lat_idx[1], 1), 
		count=c(length(lon_idx), length(lat_idx), -1))
	nc_close(ncdata)
	data
}

