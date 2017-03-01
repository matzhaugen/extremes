# setwd("~/Cloud/Research/ExtremeEvents/")
# setwd("r_code")
source("helper.R")
# source("game.R")
# n_cores = Sys.getenv("SLURM_NTASKS_PER_NODE")
n_cores = 14
outer.cores = 4 # number of outer loop cores, used to optmize code.
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
max_time = days_per_year * upper_year


restricted_lat = lats[(lats < lat_max) & (lats > lat_min)]
restricted_lon = lons[(lons < lon_max) & (lons > lon_min)]
n_lats = length(restricted_lat)
n_lons = length(restricted_lon)
# n_pixels = n_cores # number of pixels to analyse
n_pixels = n_lats * n_lons

n_files = 50
files = list.files("../rsriver/timefirst/")[1:n_files]

# get spline bases for day of year (x) and number of years (t)
day_of_year = 1:365
x = rep(day_of_year, n_files*upper_year)
t = rep(c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year))), n_files)
x.basis = as.matrix(pbs(x, df=6))
t.basis = ns(t, df=4)
rm(t,x)
gc()
# Turn the 3d array of lat x lon x time to a 1d list each with the time axis,
# a 1d array. This is becasue it will use less space when distributing each entry
# of the list to a core
start.time.tot <- Sys.time()
# out = lapply(1:(n_lats * n_lons), function(pos_idx) {
out = mclapply(1:n_pixels, function(pos_idx) {
	res_lon_idx = (pos_idx-1) %% n_lons + 1
	res_lat_idx = (pos_idx-1) %/% n_lons + 1
	lon_idx = which(lons == restricted_lon[res_lon_idx])
	lat_idx = which(lats == restricted_lat[res_lat_idx])
	print(paste("Reading idx #:", pos_idx))
	print(paste("with co-ordinates #:", signif(lats[lat_idx], digits=3), 
										signif(lons[lon_idx], digits=3)))
	
	one_pixel = lapply(files, function(f) {
		path_f = paste(path.to.files, f, sep="")
		get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, -1)
		})
	
	# Fit location and scale parameters
	ts_data = unlist(lapply(one_pixel, function(ts) {
				ts[1:(max_time)]
				}))
	rm(one_pixel)
	gc()
	start.time2 <- Sys.time()
	loc.scale = fit.2d.location.scale(t.basis, x.basis, ts_data, q1, q2, 
		days_per_year=days_per_year)
	printf("fit time: %f", Sys.time() - start.time2)	

	# Fit shape and scale parameters, xi and beta.
	# nu = log(xi(1+beta))
	mu = rep(c(loc.scale[[1]]), n_files)
	res = ts_data - mu
	idx = which(res > 0)
	over = ts_data[idx]
	idx_x = (idx - 1) %% days_per_year + 1
	idx_t = ((idx - 1) %/% days_per_year + 1 - 1) %% upper_year + 1
	seas = pbs(idx_x, df=6)
	class(seas) = append(class(seas),"matrix")
	t_year = ns(idx_t, df=4)
	# compute size of clusters of extremes, max 50 consequtive days
	event.clust = cluster(idx)[1:50]

	formula = ~(seas.1+seas.2+seas.3+seas.4+seas.5+seas.6)*
				(t_year.1+t_year.2+t_year.3+t_year.4)
	start.time2 <- Sys.time()
	fit = gamGPDfit(x=data.frame(y=res[idx], seas=seas, t_year=t_year), 
			threshold=0,
			datvar="y", 
			xiFrhs=formula,
			nuFrhs=formula,
			init=c(-0.5, 0.5))
	printf("fit time for shape/scale: %f", Sys.time() - start.time2)
	
	# Extract xi and nu estimates
	x = rep(day_of_year, upper_year)
	t = c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year)))
	mat.x = as.matrix(pbs(x, df=6))
	X = model.matrix(~mat.x*ns(t, df=4))
	xi.hat = matrix(X %*% fit$xiObj$coeff, ncol=upper_year)
	nu.hat = matrix(X %*% fit$nuObj$coeff, ncol=upper_year)
	beta.hat = exp(nu.hat) / (1 + xi.hat)

	# Keep only every 5th sample to save space
	sample_x = seq(1, 365, 5)
	sample_y = seq(1, 250, 5)
	loc.scale$xi = xi.hat[sample_x, sample_y]
	loc.scale$beta = beta.hat[sample_x, sample_y]
	loc.scale$mu = loc.scale$mu[sample_x, sample_y]
	loc.scale$scale = loc.scale$scale[sample_x, sample_y]
	loc.scale$event.clust = event.clust
	loc.scale
}, mc.cores=n_cores)
printf("Total time elapsed: %f", Sys.time() - start.time.tot)
save(restricted_lat, restricted_lon, out, file = paste("out.Rdata"))

## Loop through North America
# start.time <- Sys.time()
# out = mclapply(1:length(restricted_lat), function(lat_idx) {
# 	print(lat_idx)
# 	mclapply(1:length(restricted_lon), function(lon_idx) {
# 		print(lon_idx)
# 		loc.scale = lapply(idx, function(my_idx) { 
# 			one_pixel = unlist(lapply(r, function(x) x[lon_idx, lat_idx, my_idx]))
# 			fit.location.scale(basis, one_pixel, q1, q2, days_per_year=days_per_year)
# 		})
# 		mus = do.call(cbind, lapply(loc.scale, function(x) x$mu))
# 		scales = do.call(cbind, lapply(loc.scale, function(x) x$scale))
# 		data.frame(mus=mus, scales=scales)
# 	}, mc.cores=4)
# }, mc.cores=2)
# end.time <- Sys.time()
# print(end.time - start.time)

# save(restricted_lat, restricted_lon, out, file = paste("out.Rdata"))

# # Debug one pixel

# loc.scale = lapply(idx, function(my_idx) { 
# 	one_pixel = unlist(lapply(r, function(x) x[1, 1, my_idx]))
# 	fit.location.scale(basis, one_pixel, q1, q2, days_per_year=days_per_year)
# })
# mus = do.call(cbind, lapply(loc.scale, function(x) x$mu))
# scales = do.call(cbind, lapply(loc.scale, function(x) x$scale))
# one_pixel_out = data.frame(mus=mus, scales=scales)
# persp(x=seq(1, 365, 5), y=seq(1,250 - 1, 10), z=as.matrix(one_pixel_out[seq(1, 365, 5),26:50]), 
# 			theta=35, phi=25, ylab="Time (years)", xlab="Day of year", zlab="Temp.")
