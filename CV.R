source("helper.R")
library(Hmisc)
# n_cores = Sys.getenv("SLURM_NTASKS_PER_NODE")

n_cores = 10
path.to.files = "../rsriver/timefirst/"

# ncdata <- nc_open('CCSM3_temperatures/4001_5000_i1400_1870_R01_combine_TREFHTMX_US_only.nc')
ncdata <- nc_open('../rsriver/timefirst/trefht_4200.nc')
# print(ncdata)
# data = ncvar_get(ncdata, 'TREFHT')
lons = ncvar_get(ncdata, 'lon') - 180
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


restricted_lat = lats[(lats < lat_max) & (lats > lat_min)]
restricted_lon = lons[(lons < lon_max) & (lons > lon_min)]
n_lats = length(restricted_lat)
n_lons = length(restricted_lon)
# n_pixels = n_cores # number of pixels to analyse
n_pixels = n_lats * n_lons

n_files = 50
files = list.files("../rsriver/timefirst/")[1:n_files]

pos_idx=169 # a good pixel to look at
my_idx=c(204, 175, 169)

## Choosing the right set of predictors
res_lon_idx = (pos_idx-1) %% n_lons + 1
res_lat_idx = (pos_idx-1) %/% n_lons + 1
lon_idx = which(lons == restricted_lon[res_lon_idx])
lat_idx = which(lats == restricted_lat[res_lat_idx])
lon_idx = which(lons == -127.5)
lat_idx = which(signif(lats, 3) == 35.3)
print(paste("Reading idx #:", pos_idx))
print(paste("with co-ordinates #:", signif(lats[lat_idx], digits=3), 
									signif(lons[lon_idx], digits=3)))	

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
			  c(18, 5, 5))
# get spline bases for day of year (x) and number of years (t)

# instead of running the loop below you can load "CV.Rdata"
load("CV.Rdata")
T = mclapply(models, function(m) { 
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
		c(exeedence_var(tr_data, y_tr), exeedence_var(te_data, y_te))
	})
	T_tr = unlist(lapply(T, function(x) x[1]))
	T_te = unlist(lapply(T, function(x) x[2]))
	cbind(T_tr, T_te)

}, mc.cores = 9)
save(T, file="CV.Rdata")

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

t_mean_tr = as.vector(unlist(lapply(T, function(t) colMeans(t)[1])))
t_mean_te = as.vector(unlist(lapply(T, function(t) colMeans(t)[2])))
t_sd_tr = as.vector(unlist(lapply(T, function(t) apply(t,2, sd)[1])))
t_sd_te = as.vector(unlist(lapply(T, function(t) apply(t,2, sd)[2])))

m_lat = signif(lats[lat_idx], digits=3)
m_lon = signif(lons[lon_idx], digits=3)
pdf(paste("CV", m_mat, "_", m_lon, ".pdf"), width=6, height=6)
# par(mfrow=c(1,3))
errbar(1:length(T), t_mean_tr, t_mean_tr+t_sd_tr, t_mean_tr-t_sd_tr, ylim=range(unlist(T)*1.1),
	xlab="Model #", ylab="Block Exeedence Variance")
title("Model selection")
errbar(1:length(T), t_mean_te, t_mean_te+t_sd_te, t_mean_te-t_sd_te, add=TRUE, col=2)
legend("topright", c("Train", "Test"), col=c(1,2), pch=16)
dev.off()
# points(1:length(T), t_mean_te, col=2)



n_tr = length(y_tr) %/% (days_per_year * upper_year)
n_te = length(y_te) %/% (days_per_year * upper_year)
res_tr = tr_data - y_tr
posres_tr = res_tr[res_tr > 0]
postime_tr = rep(day_of_year, n_tr*upper_year)[res_tr > 0]
res_te = te_data - y_te
q_est_te = length(res_te[res_te > 0])/dim(X_te)[1]
posres_te = res_te[res_te > 0]
postime_te = rep(day_of_year, n_te*upper_year)[res_te > 0]
pdf("Exeedence.pdf", width=6, height=6)
comp.hist(postime_tr, postime_te, breaks=30, xlab="Days in year", ylab="Exeedence density", main="Exeedence Times of Residuals")
# h_tr = hist(postime_tr, breaks=30)
# h_te = hist(postime_te, breaks=30)
dev.off()
# ti_tr = h_tr$counts[-length(h_tr$counts)]
# ti_te = h_te$counts[-length(h_te$counts)]
# T_tr = sd(ti_tr)
# T_te = sd(ti_te)


hist(postime_tr, breaks=30)
plot(fitted.values[1:365])
length(postime[postime <= 300]) / n_samples * 300/365
pdf('../figures/p204_exceedences_vs_day.pdf')
hist(postime, breaks=80, xlab="Day of the year")
dev.off()
gc()
printf("fit time: %f", Sys.time() - start.time2)
