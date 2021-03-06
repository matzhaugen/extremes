suppressMessages(source("helper.R"))
mc.cores = 10
library(ncdf4)
path.to.files = "../rsriver/timefirst/"
# ncdata <- nc_open('CCSM3_temperatures/4001_5000_i1400_1870_R01_combine_TREFHTMX_US_only.nc')
ncdata <- nc_open('../rsriver/timefirst/trefht_4200.nc')
# print(ncdata)
varname = 'TREFHT'
data = ncvar_get(ncdata, 'TREFHT')
lons = signif(ncvar_get(ncdata, 'lon'), digits=3)
lons[lons>180] = lons[lons>180] - 360
lats = ncvar_get(ncdata, 'lat')
nc_close(ncdata)
KtoC = 273.15
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
B = 100 # number of bootstrap runs
max_time = days_per_year * upper_year


restricted_lat = lats[(lats < lat_max) & (lats > lat_min)]
restricted_lon = lons[(lons < lon_max) & (lons > lon_min)]
n_lats = length(restricted_lat)
n_lons = length(restricted_lon)
# n_pixels = n_cores # number of pixels to analyse
n_pixels = n_lats * n_lons

n_files = 50
files = list.files(path.to.files)[1:n_files]

# get spline bases for day of year (x) and number of years (t)
day_of_year = 1:365
pos_idx=100 # a good pixel to look at
my_idx=c(204, 175, 169)
t_grace=15
d_grace=15
start.time.tot <- Sys.time()
breaks = 100
seasons = list(c(seq(365-31, 365), c(1:58)), #winter
	c(60:(60+91)),  #spring
	c(153:(153+91)), #summer
	c(246:(246+90))) #fall

hist.data = mclapply(1:n_pixels, function(pos_idx) {
# hist.data = mclapply(1:10, function(pos_idx) {
	res_lon_idx = (pos_idx-1) %% n_lons + 1
	res_lat_idx = (pos_idx-1) %/% n_lons + 1
	lon_idx = which(lons == restricted_lon[res_lon_idx])
	lat_idx = which(lats == restricted_lat[res_lat_idx])
	print(paste("Reading idx #:", pos_idx))
	print(paste("with co-ordinates #:", signif(lats[lat_idx], digits=3), 
										signif(lons[lon_idx], digits=3)))	

	one_pixel_i = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, days_per_year * t_grace,
					varname=varname) - KtoC
			}))

	one_pixel_f = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, days_per_year * (upper_year - t_grace), days_per_year * t_grace,
					varname=varname) - KtoC
			}))
	temp_i = array(one_pixel_i, dim=c(days_per_year, t_grace, n_files))
	temp_f = array(one_pixel_f, dim=c(days_per_year, t_grace, n_files))
	temp_range = range(c(temp_i, temp_f))

	# out = lapply(seq(1, 365, d_grace), function(start_idx) {
	# 	idx = 1:d_grace + min(start_idx, 366 - d_grace) - 1
	out = lapply(seasons, function(start_idx) {
		idx = start_idx
		hist_mom = lapply(1:B, function(b) {
			if (b != 1) {
				file_idx = sample(1:n_files, replace=T)
			} else {
				file_idx = 1:n_files
				obs.i = c(temp_i[idx,,file_idx])
				obs.f = c(temp_f[idx,,file_idx])
				original_breaks <<- hist(c(obs.i, obs.f), breaks=breaks, plot=FALSE)$breaks
			}
			obs.i = c(temp_i[idx,,file_idx])
			obs.f = c(temp_f[idx,,file_idx])
			h.i = hist(obs.i, breaks=original_breaks, plot=FALSE)
			h.f = hist(obs.f, breaks=original_breaks, plot=FALSE)
			moments.i = c(mean(obs.i), sd(obs.i), skewness(obs.i))
			moments.f = c(mean(obs.f), sd(obs.f), skewness(obs.f))
			list(h.i = h.i, h.f = h.f, moments.i = moments.i, moments.f = moments.f)
		})
		hist_stats_i = hist_stats(do.call(rbind, lapply(hist_mom, function(el) {el$h.i$counts})))
		hist_stats_f = hist_stats(do.call(rbind, lapply(hist_mom, function(el) {el$h.f$counts})))
		mom_i = do.call(rbind, lapply(hist_mom, function(el) {el$moments.i}))
		mom_f = do.call(rbind, lapply(hist_mom, function(el) {el$moments.f}))
		list(breaks = original_breaks,
			hist_stats_i = hist_stats_i,
			hist_stats_f = hist_stats_f,
			moments.i = hist_stats(mom_i),
			moments.f = hist_stats(mom_f),
			moments.d = hist_stats(mom_f - mom_i))
	})
		
		# h.f = hist(c(temp_f[idx,,]), breaks=breaks, plot=FALSE)
	out
}, mc.cores = mc.cores)

# save(hist.data, file="hist.Rdata")
save(hist.data, file="histSeasons.Rdata")


#test 
print(c(lats[lat_idx], lons[lon_idx]))
path.to.files = "../rsriver/spacefirst/"
files = list.files(path.to.files)
one_pixel_i_s = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.light(path_f, lat_idx, 1, lon_idx, 1, 1, days_per_year * t_grace) - KtoC
			}))
temp_i_s = array(one_pixel_i_s, dim=c(days_per_year, t_grace, n_files))
path.to.files = "../rsriver/timefirst/"
files = list.files(path.to.files)
one_pixel_i_t = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, days_per_year * t_grace) - KtoC
			}))
temp_i_t = array(one_pixel_i_t, dim=c(days_per_year, t_grace, n_files))
hist(temp_i_s[1:10,,])
filename = paste(path.to.files, files[1], sep="")
ncdata <- nc_open(filename)	
data = ncvar_get(ncdata, 
'TREFHT', 
start=c(1, 1, 1), 
count=c(3650, 1, 1)) - KtoC
nc_close(ncdata)
	data
get.ncdf.timefirst(path_f, 96, 1, 48, 1, 1, days_per_year * t_grace) - KtoC
# plot data
library(ggplot2)
load("hist.Rdata")
breaks = seq(-70,50, 1)
pos_idx=204 # a good pixel to look at
blue <- rgb(0, 0, 1, alpha=0.4)
red <- rgb(1, 0, 0, alpha=0.4)
day = 1
time_idx = floor(day / 10) + 1
pos_idx = 204
res_lon_idx = (pos_idx-1) %% n_lons + 1
res_lat_idx = (pos_idx-1) %/% n_lons + 1
lon_idx = which(lons == restricted_lon[res_lon_idx])
lat_idx = which(lats == restricted_lat[res_lat_idx])
xval = (breaks[-1] + breaks[-length(breaks)]) / 2
yval.i = hist.data[[pos_idx]][[time_idx]]$counts.i
yval.f = hist.data[[pos_idx]][[time_idx]]$counts.f
mi = min(which(yval.i != 0))
mi.f = min(which(yval.f != 0))
ma = max(which(yval.f != 0))

initial = data.frame(Temp=breaks[mi:ma], Count=yval.i[mi:ma], Model="Initial")
final = data.frame(Temp=breaks[mi:ma], Count=yval.f[mi:ma], Model="Final")
simple = rbind(initial, final)
base = ggplot(simple, aes(x=Temp, y=Count, fill=type))
print(base + geom_bar(stat="identity", alpha=0.5, width=1, position="identity"))