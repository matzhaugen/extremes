source("helper.R")
library(doMC)
library(quantreg)
library(ncdf4)
a <- read.table("outBoot100_nb")

# n_cores = Sys.getenv("SLURM_NTASKS_PER_NODE")

n_cores = 10
path.to.files = "../rsriver/timefirst/"

# ncdata <- nc_open('CCSM3_temperatures/4001_5000_i1400_1870_R01_combine_TREFHTMX_US_only.nc')
ncdata <- nc_open('../rsriver/timefirst/trefht_4200.nc')
# print(ncdata)
# data = ncvar_get(ncdata, 'TREFHT')
lons = ncvar_get(ncdata, 'lon')
lons[lons>180] = lons[lons>180] - 360
lats = ncvar_get(ncdata, 'lat')
nc_close(ncdata)

lat_min = 10
lat_max = 60
lon_min = -130
lon_max = -60
n_years_window = 10
n_years_increment = 10
days_per_year = 365
upper_year = 250
q1 = 0.95
q2 = 0.975
q = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
nq = length(q)
B = 10 # number of bootstrap runs
max_time = days_per_year * upper_year
KtoC = 273.15

restricted_lat = lats[(lats < lat_max) & (lats > lat_min)]
restricted_lon = lons[(lons < lon_max) & (lons > lon_min)]
n_lats = length(restricted_lat)
n_lons = length(restricted_lon)
# n_pixels = n_cores # number of pixels to analyse
n_pixels = n_lats * n_lons

n_files = 50
files = list.files("../rsriver/timefirst/")[1:n_files]

pos_idx=175 # a good pixel to look at
my_idx=c(204, 175, 169)

## Choosing the right set of predictors
res_lon_idx = (pos_idx-1) %% n_lons + 1
res_lat_idx = (pos_idx-1) %/% n_lons + 1
lon_idx = which(lons == restricted_lon[res_lon_idx])
lat_idx = which(lats == restricted_lat[res_lat_idx])
print(paste("Reading idx #:", pos_idx))
print(paste("with co-ordinates #:", signif(lats[lat_idx], digits=3), 
									signif(lons[lon_idx], digits=3)))	

# If you want to plot the locations
pdf("CVlocation.pdf")
map(xlim=c(-180, 0),
	  ylim=c(0, 90), xaxt='n')
map.axes()
points(lons[lon_idx], lats[lat_idx], cex=2, col=1)
dev.off()

# Import all data
one_pixel = lapply(files, function(f) {
	path_f = paste(path.to.files, f, sep="")
	get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, -1)
})

n_samples = dim(X)[1]
nfolds = 10
T_tr = rep(0, nfolds)
T_te = rep(0, nfolds)
idx = sample(rep(seq(nfolds), length=50))

models = list(c(5, 3, 3),
			  c(7, 3, 3),
			  c(10, 3, 3),
			  c(10, 3, 4),
			  c(12, 3, 4),
			  c(15, 3, 4),
			  c(15, 3, 5),
			  c(15, 5, 5),
			  c(18, 5, 5),
			  c(18, 7, 6),
			  c(18, 7, 8))
# get spline bases for day of year (x) and number of years (t)
day_of_year = 1:365
x = rep(day_of_year, n_files*upper_year)
t = rep(c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year))), n_files)
x.main.basis = as.matrix(pbs(x, df=m[1]))
x.int.basis = as.matrix(pbs(x, df=m[2]))
t.basis = ns(t, df=m[3])
X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis)
coef = matrix(0, dim(X)[2], length(q))
rm(t,x, t.basis, x.main.basis, x.int.basis)
gc()

# instead of running the loop below you can load "CV.Rdata"
# load("CV.Rdata")
T = mclapply(models, function(m) { 
	
	T = lapply(1:nfolds, function(k) {
		print(k)
		# Fit location and scale parameters
		tr_data = unlist(lapply(one_pixel[idx!=k], function(ts) {
				ts[1:(max_time)]
				}))
		te_data = unlist(lapply(one_pixel[idx==k], function(ts) {
				ts[1:(max_time)]
				}))

		start.time2 <- Sys.time()
		X_tr = X[1:(days_per_year*upper_year*sum(idx != k)), ]
		X_te = X[1:(days_per_year*upper_year*sum(idx == k)), ]
		coef = rq.fit.pfn(X_tr, y=tr_data, tau=0.95, max.bad.fixup=20)$coefficients
		y_tr = X_tr%*%coef
		y_te = X_te%*%coef
		# Change check.seasonality to change which variable to cross-validate
		c(exeedence_var(tr_data, y_tr, check.seasonality=TRUE), 
			exeedence_var(te_data, y_te, check.seasonality=TRUE))
	})
	T_tr = unlist(lapply(T, function(x) x[1]))
	T_te = unlist(lapply(T, function(x) x[2]))
	cbind(T_tr, T_te)

}, mc.cores = 9)
save(T, file="CV_t_loc175.Rdata")
load("CV.Rdata")
load("CV_t.Rdata") #This cross-validates time

# Plot CV
range_d = c(0, 0.008)
t_mean_tr = as.vector(unlist(lapply(T, function(t) colMeans(t)[1])))
t_mean_te = as.vector(unlist(lapply(T, function(t) colMeans(t)[2])))
t_sd_tr = as.vector(unlist(lapply(T, function(t) apply(t,2, sd)[1])))
t_sd_te = as.vector(unlist(lapply(T, function(t) apply(t,2, sd)[2])))
cex.lab = 1.2
m_lat = signif(lats[lat_idx], digits=3)
m_lon = signif(lons[lon_idx], digits=3)
pdf(paste("CV_t_loc175.pdf", sep=""), width=5, height=5)
# pdf(paste("CV.pdf", sep=""), width=5, height=5)
# par(mfrow=c(1,3))
errbar(1:length(T), t_mean_tr, t_mean_tr+t_sd_tr, t_mean_tr-t_sd_tr, ylim=range_d,
	xlab="Model #", ylab="Block Exeedence Std Dev", cex.lab=cex.lab)
# title("Seasonality")
errbar(1:length(T), t_mean_te, t_mean_te+t_sd_te, t_mean_te-t_sd_te, add=TRUE, col=2)
legend("topright", c("Train", "Test"), col=c(1,2), pch=16)
dev.off()
# points(1:length(T), t_mean_te, col=2)


## Plot exceedence with error bars
day_of_year = 1:365
x = rep(day_of_year, n_files*upper_year)
t = rep(c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year))), n_files)
x.main.basis = as.matrix(pbs(x, df=15))
x.int.basis = as.matrix(pbs(x, df=3))
t.basis = ns(t, df=4)
X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis)
coef = matrix(0, dim(X)[2], length(q))
rm(t,x, t.basis, x.main.basis, x.int.basis)
gc()
k=1
X_tr = X[1:(days_per_year*upper_year*sum(idx != k)), ]
X_te = X[1:(days_per_year*upper_year*sum(idx == k)), ]
n_breaks = 30	
out = mclapply(1:10, function(k) {
	print(k)
	tr_data = unlist(lapply(one_pixel[idx!=k], function(ts) {
				ts[1:(max_time)]
				}))
	te_data = unlist(lapply(one_pixel[idx==k], function(ts) {
					ts[1:(max_time)]
					}))

	start.time2 <- Sys.time()
	coef = rq.fit.pfn(X_tr, y=tr_data, tau=0.95, max.bad.fixup=20)$coefficients
	y_tr = X_tr%*%coef
	y_te = X_te%*%coef

	n_tr = length(y_tr) %/% (days_per_year * upper_year)
	n_te = length(y_te) %/% (days_per_year * upper_year)
	res_tr = tr_data - y_tr
	posres_tr = res_tr[res_tr > 0]
	postime_tr = rep(day_of_year, n_tr*upper_year)[res_tr > 0]
	res_te = te_data - y_te
	q_est_te = length(res_te[res_te > 0])/dim(X_te)[1]
	posres_te = res_te[res_te > 0]
	postime_te = rep(day_of_year, n_te*upper_year)[res_te > 0]
	
	h_tr = hist(postime_tr, breaks=n_breaks, plot=F)
	h_te = hist(postime_te, breaks=n_breaks, plot=F)
	list(h_tr, h_te)
}, mc.cores=10)
density_tr = do.call(rbind, lapply(out, function(item) {
	item[[1]]$density[1:(n_breaks-1)]
}))


load("CVout.Rdata")
h_tr = out[[1]][[1]]
h_te = out[[1]][[2]]
n_breaks=length(h_tr$breaks)
qhat_te = do.call(rbind, lapply(out, function(item) {
	item[[2]]$counts[1:(n_breaks-1)] / (5*250*10)
}))
qhat_tr = h_tr$counts[1:(n_breaks-1)] / (45*250*10)
den_mean = apply(qhat_te, 2, mean)
den_sd = apply(qhat_te, 2, sd)

pdf("Exeedence_loc175.pdf", width=5, height=5)
# plot(h_tr$breaks[1:(n_breaks-1)], qhat_tr, type='n',
# 	xlab="Day of year", ylab="Blocked quantile estimate", main="", col=2, lwd=2,
# 	ylim=c(0.04, 0.06))
errbar(h_te$breaks[1:(n_breaks-1)], den_mean, den_mean + den_sd, den_mean - den_sd, 
	add=FALSE,xlab="Day of year", ylab="Blocked quantile estimate", main="", col=1, lwd=2,
	ylim=c(0.04, 0.06), cex.lab=1.2)
legend('topright', c('Test data'), lty=c(0), pch=c(16), col=c(1), lwd=c(0))
dev.off()

# comp.hist(postime_tr, postime_te, breaks=30, xlab="Days in year",
#  ylab="Exeedence density", main="Exeedence Times of Residuals")
# h_tr = hist(postime_tr, breaks=30)
# h_te = hist(postime_te, breaks=30)
# dev.off()
# ti_tr = h_tr$counts[-length(h_tr$counts)]
# ti_te = h_te$counts[-length(h_te$counts)]
# T_tr = sd(ti_tr)
# T_te = sd(ti_te)

# Plot observations of the first 10 years and the associated fit
k=1
tr_data = unlist(lapply(one_pixel[idx!=k], function(ts) {
				ts[1:(max_time)]
				}))
tr_data_C = tr_data - KtoC
X_tr = X[1:(days_per_year*upper_year*sum(idx != k)), ]
X_te = X[1:(days_per_year*upper_year*sum(idx == k)), ]
coef_low = rq.fit.pfn(X_tr, y=tr_data_C, tau=0.025, max.bad.fixup=20)$coefficients
coef_mid = rq.fit.pfn(X_tr, y=tr_data_C, tau=0.50, max.bad.fixup=20)$coefficients
coef_high = rq.fit.pfn(X_tr, y=tr_data_C, tau=0.975, max.bad.fixup=20)$coefficients
y_tr_low = X_tr%*%coef_low
y_tr_mid = X_tr%*%coef_mid
y_tr_high = X_tr%*%coef_high

n_years_begin = 5
nobs_per_run = days_per_year * upper_year
idx_i = c(outer(1:(days_per_year * n_years_begin), (1:sum(idx != k) - 1) * nobs_per_run, '+'))
idx_est = 1:(days_per_year) +  (days_per_year * floor(n_years_begin/2))  
y_hat_mid = y_tr_mid[idx_est]
y_hat_low = y_tr_low[idx_est]
y_hat_high = y_tr_high[idx_est]
pdf("Estimate.pdf", width=6, height=6)
par(oma=c(0,1,0,0))
ex_x = (idx_i-1) %% days_per_year + 1
ex_y = tr_data_C[idx_i]
random_sample = sample(1:length(ex_x), round(length(ex_x)/50))
plot(ex_x[random_sample], ex_y[random_sample],
	xlab="Day of the year", ylab="Temperature (deg C)", cex.lab=1.5,
	main="",
	ylim=range(ex_y) * 1.15,
	col='gray')
lines(1:days_per_year, y_hat_mid, col=1, lwd=3)
lines(1:days_per_year, y_hat_low, col=1, lty=2, lwd=3)
lines(1:days_per_year, y_hat_high, col=1, lty=2, lwd=3)
legend('topright', c('Median', '2.5/97.5th quantile'), lty=c(1,2), col=1, lwd=c(2,2))
dev.off()
save(out, T, ex_x, ex_y, y_hat_mid, y_hat_low, y_hat_high, file="CVout.Rdata")

# Plot the above plot at the three locations
load("CVout.Rdata")
my_lon = c(-101,-82, -98)
my_lat = c(50.1, 42.7, 35.3)
text.pos.y = c(20,21,31)
n_years_begin = 1
nobs_per_run = days_per_year * upper_year
idx_i = c(outer(1:(days_per_year * n_years_begin), (1:length(files) - 1) * nobs_per_run, '+'))
idx_est = 1:(days_per_year) +  (days_per_year * floor(n_years_begin/2))  
ex_y_3loc = list(); length(ex_y_3loc) = 3
pdf("Estimate_3loc.pdf", width=12, height=4)
par(mfrow=c(1,3), oma=c(0,1.5,0,0))
for (i in 1:length(my_lon)) {
	lon_idx = nearest_idx(lons, my_lon[i])
	lat_idx = nearest_idx(lats, my_lat[i])

	one_pixel = lapply(files, function(f) {
		path_f = paste(path.to.files, f, sep="")
		get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, -1)
	})
	data = unlist(lapply(one_pixel, function(ts) {
					ts[1:(max_time)] - KtoC
					}))
	pos_idx = getPos_idx(my_lat[i], my_lon[i] + 360)
	myquantiles = getq(pos_idx, 'avg')

	x_idx = 1:73
	y_hat_mid = getQuantile(pos_idx, 0.5, 'avg')[[1]][x_idx] 
	y_hat_low = getQuantile(pos_idx, 0.025, 'avg')[[1]][x_idx] 
	y_hat_high = getQuantile(pos_idx, 0.975, 'avg')[[1]][x_idx] 
	

	ex_x = (idx_i - 1) %% days_per_year + 1
	ex_y_3loc[[i]] = data[idx_i]
	random_sample = sample(1:length(ex_x), round(length(ex_x)/50))
	plot(ex_x[random_sample], ex_y_3loc[[i]][random_sample],
		xlab="Day of the year", ylab="Temperature (deg C)", cex.lab=1.5,
		main="",
		ylim=range(ex_y_3loc[[i]]) * 1.15,
		col='gray')
	my_x = seq(1, days_per_year, by=5)
	lines(my_x, y_hat_mid, col=1, lwd=3)
	lines(my_x, y_hat_low, col=1, lty=2, lwd=3)
	lines(my_x, y_hat_high, col=1, lty=2, lwd=3)
	if (i ==1) legend('bottomleft', c('Median', '2.5/97.5th quantile'), lty=c(1,2), col=1, lwd=c(2,2))
	text(10, text.pos.y[i], c('a', 'b', 'c')[i], cex=3)

}
dev.off()	

save(ex_y_3loc, out, T, ex_x, ex_y, y_hat_mid, y_hat_low, y_hat_high, file="CVout.Rdata")


