# setwd("~/Cloud/Research/ExtremeEvents/")
# setwd("r_code")

source("helper.R")
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
q = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975)
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

# get spline bases for day of year (x) and number of years (t)
day_of_year = 1:365
x = rep(day_of_year, n_files*upper_year)
t = rep(c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year))), n_files)
x.int.basis = as.matrix(pbs(x, df=3))
x.main.basis = as.matrix(pbs(x, df=16))
t.basis = ns(t, df=4)
X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis)
coef = matrix(0, dim(X)[2], length(q))
rm(t,x, t.basis, x.main.basis, x.int.basis)
gc()
pos_idx=169 # a good pixel to look at
my_idx=c(204, 175, 169)

start.time.tot <- Sys.time()
out = lapply(my_idx, function(pos_idx) {
# for (pos_idx in 1:4) {
	res_lon_idx = (pos_idx-1) %% n_lons + 1
	res_lat_idx = (pos_idx-1) %/% n_lons + 1
	lon_idx = which(lons == restricted_lon[res_lon_idx])
	lat_idx = which(lats == restricted_lat[res_lat_idx])
	print(paste("Reading idx #:", pos_idx))
	print(paste("with co-ordinates #:", signif(lats[lat_idx], digits=3), 
										signif(lons[lon_idx], digits=3)))	

	# The first run is the point estimate
	rq.fit = mclapply(1:B, function(b) {
		print(b)
		if (b != 1) { 
			files_resampled = sample(files, replace=TRUE)
		} else {
			files_resampled = files
		}
		one_pixel = lapply(files_resampled, function(f) {
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
		for (i in 1:length(q)) {
			coef[,i] = rq.fit.pfn(X, y=ts_data, tau=q[i], max.bad.fixup=20)$coefficients	
			gc()
		}
		printf("fit time: %f", Sys.time() - start.time2)
		coef
	}, mc.cores=n_cores)	
	gc()
	rq.fit
})

printf("Total time elapsed: %f", Sys.time() - start.time.tot)
save(restricted_lat, restricted_lon, out, file = paste("outBoot.Rdata"))

