library(ncdf4)
# Read a subset of the netcdf file based on starting lat/lon/time indeces. 
# Also include the index count for each variable.
get.ncdf.timefirst = function(filename, lat_start, lat_count, lon_start, lon_count, t_start, t_count, varname='TREFHT') {
	ncdata <- nc_open(filename)	
		data = ncvar_get(ncdata, 
		varname, 
		start=c(t_start, lon_start, lat_start), 
		count=c(t_count, lon_count, lat_count))
	nc_close(ncdata)
	data
}

path.to.files = "../rsriver/timefirst/"

ncdata <- nc_open(paste(path.to.files, 'trefht_4200.nc', sep=""))
lons = ncvar_get(ncdata, 'lon')
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
t_grace=10
d_grace=10
start.time.tot <- Sys.time()
breaks = seq(-70,50, 1)

pos_idx = 100 # this corresponds to Lousiana approximately
res_lon_idx = (pos_idx-1) %% n_lons + 1
res_lat_idx = (pos_idx-1) %/% n_lons + 1
lon_idx = which(lons == restricted_lon[res_lon_idx])
lat_idx = which(lats == restricted_lat[res_lat_idx])
print(paste("Reading idx #:", pos_idx))
print(paste("with co-ordinates #:", signif(lats[lat_idx], digits=3), 
									signif(lons[lon_idx], digits=3)))

#test 
print(c(lats[lat_idx], lons[lon_idx]))
path.to.files = "../rsriver/timefirst/"
files = list.files(path.to.files)
one_pixel_i_s = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, days_per_year * t_grace) - KtoC
			}))
temp_i_s = array(one_pixel_i_s, dim=c(days_per_year, t_grace, n_files))
hist(temp_i_s[1:30,,]) 