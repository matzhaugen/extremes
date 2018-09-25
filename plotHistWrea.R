# Run in the r_code subfolder
########################################
# PART 1: Get the data and store it ####
########################################
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
# lons[lons>180] = lons[lons>180] - 360
lats = signif(ncvar_get(ncdata, 'lat'), 3)

nc_close(ncdata)
KtoC = 273.15
lat_min = 0
lat_max = 90
lon_min = 180
lon_max = 360
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
latlon = cbind(rep(restricted_lat, n_lons), 
				rep(restricted_lon, n_lats))
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
seas_x = as.matrix(pbs(rep(rep(day_of_year, t_grace), n_files), df=15))

start.time.tot <- Sys.time()
breaks = 100
seasons = list(c(seq(365-31, 365), c(1:58)), #winter
	c(60:(60+91)),  #spring
	c(153:(153+91)), #summer
	c(246:(246+90))) #fall

# Reanalysis
filename = "../reanalysis/t_19790101_19940101.nc"

ncdata <- nc_open(filename)	
lat = ncvar_get(ncdata, 'latitude')
lon = ncvar_get(ncdata, 'longitude')
time = ncvar_get(ncdata, 'time')
last.t = max(which(time!=0))
time = time[1:last.t]

time.s=as.POSIXct('1979-02-01 00:00',tz='UTC')
time.e=as.POSIXct('1994-01-01 00:00',tz='UTC')
tseq=seq(time.s, time.e, by='6 hours')
times=as.POSIXct(time*60*60, origin='1900-01-01 00:00:0.0', tz='UTC')
t1=which(times==time.s)
t2=which(times==time.e)
lim.lat = c(90,0)
lim.lon = c(180, 359.25)
lat1 = which(lat==lim.lat[1])
lat2 = which(lat==lim.lat[2])
lon1 = which(lon==lim.lon[1])
lon2 = which(lon==lim.lon[2])
dlat = lat2-lat1+1
dlon = lon2-lon1+1
latsub=lat[lat2:lat1]
lonsub = lon[lon1:lon2]
timesub = times[seq(t1, t2, by=1)]
nt = length(timesub) - 1
timesub = times[seq(t1, t2-1, by=4)]
dt=t2-t1+1
rea_data = ncvar_get(ncdata, 't2m', start=c(lon1, lat1, t1),
	count=c(dlon, dlat, dt)) - 273.15

# rea_data_avg = 1/4 * ( rea_data[,length(latsub):1,seq(1, nt, by=4)] + rea_data[,length(latsub):1,seq(2, nt, by=4)] + 
# 	rea_data[,length(latsub):1,seq(3, nt, by=4)] + rea_data[,length(latsub):1,seq(4, nt, by=4)] )
rea_data_avg = 1/4 * ( rea_data[,length(latsub):1,seq(1, nt, by=4)] + rea_data[,length(latsub):1,seq(3, nt, by=4)] + 
	rea_data[,length(latsub):1,seq(2, nt, by=4)] + rea_data[,length(latsub):1,seq(4, nt, by=4)] )

nc_close(ncdata)

times.month = as.numeric(format(timesub, format="%m"))
times.day = as.numeric(format(timesub, format="%d"))
x = as.matrix(pbs(times.day, df=15))
season.vector = list(c(12,1,2), c(3,4,5), c(6,7,8), c(9,10,11))
nseasons=length(season.vector)
rea_moments = list()
rea_moments_detrend = list()


# Pin reanalysis data onto the coarse scale CESM grid
dlat_coarse = 3.75
dlon_coarse = 3.75
fine2coarse_lat = list(); length(fine2coarse_lat) = n_lats
fine2coarse_lon = list(); length(fine2coarse_lon) = n_lons
for (j in 1:(n_lats)) {
	fine2coarse_lat[[j]] = which((latsub < restricted_lat[j]) & (latsub > (restricted_lat[j] - dlat_coarse)))
}
for (j in 1:(n_lons)) {
	fine2coarse_lon[[j]] = which((lonsub < restricted_lon[j]) & (lonsub > (restricted_lon[j] - dlon_coarse)))
}
rea_data_agg = list(); length(rea_data_agg) = length(season.vector)
for (i in 1:length(season.vector)) {
	rea_data_agg[[i]] = list()
	length(rea_data_agg[[i]]) = n_lons
	for (j in 1:n_lons) {
		rea_data_agg[[i]][[j]] = list()
		length(rea_data_agg[[i]][[j]]) = n_lats
	}
}
rea_moments_agg = rea_data_agg
for (i in 1:length(season.vector)) {
	print(i)
	for (j in 1:n_lons) {
		for (k in 1:n_lats) {
			
			rea_data_agg[[i]][[j]][[k]] = c(rea_data_avg[fine2coarse_lon[[j]], 
												fine2coarse_lat[[k]],
												times.month %in% season.vector[[i]] ])
			ts = rea_data_agg[[i]][[j]][[k]]
			rea_moments_agg[[i]][[j]][[k]] = c(mean(ts), sd(ts),  skewness(ts))
		}
	}
}
# Get reanalysis moments
# Original scale, not downscaled currently not in use
# for (i in 1:length(season.vector)){
# 	print(i)
# 	rea_moments[[i]] = apply(rea_data_avg[,,times.month %in% season.vector[[i]]], 
# 		c(1,2), function(x) {
# 		 c(mean(x), sd(x),  skewness(x))
# 	})
	
# }
# Store 3 locations histograms
h_ra = list(list());length(h_ra) = nseasons
h_i = list(list());length(h_i) = nseasons
h_m = list(list());length(h_m) = nseasons
h_f = list(list());length(h_f) = nseasons
season.vector = list(c(12,1,2), c(3,4,5), c(6,7,8), c(9,10,11))
my_lon = c(-101, -82, -98) + 360
my_lat = c(50.1, 42.7, 35.3)
dres=2
nbreaks=200
ra_lon_idx = lapply(my_lon, function(m_lon) {
	which((lonsub <= (m_lon+dres)) & (lonsub >= (m_lon-dres)))
})
ra_lat_idx = lapply(my_lat, function(m_lat) {
	which((latsub <= (m_lat+dres)) & (latsub >= (m_lat-dres)))
})

# Get histogram data from restricted locations
for (i in 1: length(my_lon)) {
	lon_idx = which(signif(lons, 3) == my_lon[i])
	lat_idx = which(signif(lats, 3) == my_lat[i])
	one_pixel_i = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, days_per_year * t_grace,
					varname=varname) - KtoC
			}))

	one_pixel_m = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, days_per_year * t_grace + days_per_year*(1979-1850),
					varname=varname) - KtoC
			}))

	one_pixel_f = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, days_per_year * (upper_year - t_grace), days_per_year * t_grace,
					varname=varname) - KtoC
			}))
	temp_i = array(one_pixel_i, dim=c(days_per_year, t_grace, n_files))
	temp_m = array(one_pixel_m, dim=c(days_per_year, t_grace, n_files))
	temp_f = array(one_pixel_f, dim=c(days_per_year, t_grace, n_files))

	pixels.alltime = rea_data_avg[ra_lon_idx[[i]],ra_lat_idx[[i]], 
		times.month %in% unlist(season.vector) ]
	ah = hist(c(temp_i, temp_f, c(pixels.alltime)), breaks=nbreaks, plot=FALSE)
	for (j in 1:nseasons){
		pixels = rea_data_avg[ra_lon_idx[[i]],ra_lat_idx[[i]], 
			times.month %in% season.vector[[j]]]

		h_ra[[j]][[i]] = hist(pixels, breaks=ah$breaks, plot=FALSE)
		h_i[[j]][[i]]  = hist(temp_i[seasons[[j]],,], breaks=ah$breaks, plot=FALSE)
		h_f[[j]][[i]]  = hist(temp_f[seasons[[j]],,], breaks=ah$breaks, plot=FALSE)
		h_m[[j]][[i]]  = hist(temp_m[seasons[[j]],,], breaks=ah$breaks, plot=FALSE)
	}	
}

plot(h_ra[[3]][[2]]$mids, h_ra[[3]][[2]]$density, type='h', lwd=2, lend=1)

rm(rea_data)
# Compare 1979 with 2100 and reanalysis between 1979-1994
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

	one_pixel_m = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, days_per_year * t_grace + days_per_year * (1979-1850),
					varname=varname) - KtoC
			}))

	one_pixel_f = unlist(lapply(files, function(f) {
				path_f = paste(path.to.files, f, sep="")
				get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, days_per_year * (upper_year - t_grace), days_per_year * t_grace,
					varname=varname) - KtoC
			}))
	deseasoned_i = lm(one_pixel_i~seas_x)$res
	deseasoned_f = lm(one_pixel_i~seas_x)$res
	temp_i = array(one_pixel_i, dim=c(days_per_year, t_grace, n_files))
	temp_m = array(one_pixel_m, dim=c(days_per_year, t_grace, n_files))
	temp_f = array(one_pixel_f, dim=c(days_per_year, t_grace, n_files))

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
			}
			obs.i = c(temp_i[idx,,file_idx])
			obs.f = c(temp_f[idx,,file_idx])
			moments.i = c(mean(obs.i), sd(obs.i), skewness(obs.i))
			moments.f = c(mean(obs.f), sd(obs.f), skewness(obs.f))
			list(moments.i = moments.i, moments.f = moments.f)
		})
		obs.m = c(temp_i[idx,,])
		mom_m = c(mean(obs.m), sd(obs.m), skewness(obs.m))
		mom_i = do.call(rbind, lapply(hist_mom, function(el) {el$moments.i}))
		mom_f = do.call(rbind, lapply(hist_mom, function(el) {el$moments.f}))
		list(moments.i = hist_stats(mom_i),
			moments.f = hist_stats(mom_f),
			moments.d = hist_stats(mom_f - mom_i),
			moments.m = mom_m)
	})
		
		# h.f = hist(c(temp_f[idx,,]), breaks=breaks, plot=FALSE)
	out
}, mc.cores = mc.cores)

# save(hist.data, file="hist.Rdata")
save(hist.data, rea_moments,
	latlon, rea_moments_agg, 
	restricted_lon, restricted_lat, 
	lonsub, latsub, file="histSeasons.Rdata")
save(h_ra, h_i, h_f, h_m, my_lat, my_lon, nseasons, season.vector, file="histodata.Rdata")
################################
# PART 2: Plot the data ########
################################
# Run the below for the North American comaprison between CESM and ERA-interim 
load("histSeasons.Rdata")
load("histodata.Rdata")
n_lons=length(unique(restricted_lon))
n_lats=length(unique(restricted_lat))
my_lon = c(-101, -82, -98) + 360
my_lat = c(50.1, 42.7, 35.3)

hist_mom_winter = lapply(hist.data, function(el) el[[1]])
hist_mom_summer = lapply(hist.data, function(el) el[[3]])
rea_moments_agg_mat = list()

# plot locations
pdf("reanalysis_locations_w_model_pixels.pdf")
map(xlim=c(-180, 0), ylim=c(0,90))
points(x=my_lon-360, y=my_lat, col=2, lty=2, lwd=3)
for (i in 1:length(my_lon)){
	mx = lonsub[ra_lon_idx[[i]]] - 360
	my = latsub[ra_lat_idx[[i]]]
	mmx = matrix(rep(mx, length(my)), length(mx), length(my))
	mmy = matrix(rep(my, length(mx)), length(my), length(mx))
	mlatlon = cbind(c(t(mmy)),c(mmx))
	points(x=mlatlon[,2], y=mlatlon[,1])
}
dev.off()

# season_index = 1 # winter
season_index = 3 # summer
for (i in 1:3) {
	rea_moments_agg_mat[[i]] = as.matrix(do.call(rbind, lapply(rea_moments_agg[[season_index]], function(mlist) lapply(mlist, function(mm) mm[i]))))
}

hist_mom = hist_mom_summer # change this for different seasons
pdf("model_reanalysis_world_JJA_agg.pdf", width=6,height=6)
# hist_mom = hist_mom_winter # change this for different seasons
# pdf("model_reanalysis_world_DJF_agg.pdf", width=6,height=6)
rea_agg = rea_moments_agg_mat
mss=1
lettersize = 1.5
labelsize = 1.7
colbar_fontsize = 1
model_data = do.call(rbind, lapply(hist_mom, function(el) el$moments.m[mss]))
my_levels = make_zlims(c(model_data[,1], unlist(rea_agg[[mss]])), 2)
z = matrix(model_data[,1], n_lons, n_lats)
#top left
par(new = "TRUE",              
    plt = c(0.08,0.48,0.67,0.95),   # using plt instead of mfcol (compare                                  # coordinates in other plots)
    las = 1,                      # orientation of axis labels
    cex.axis = 1,                 # size of axis annotation
    tck = -0.02, 
    oma=c(1,1,2,1))

#Top left plot:
#
# the filled contour - coloured areas
xcoords = restricted_lon
ycoords = restricted_lat
# my.filled.contour3(restricted_lon, restricted_lat, 
# 		z=z, my_levels)
my.filled.contour3(restricted_lon, restricted_lat, z, my_levels)
text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=lettersize)
# par(xpd=NA)
text(x=xcoords[1]-10, y=ycoords[11], adj=2, label="Mean", srt=90, 3,
	cex=labelsize, xpd=NA)
mtext("CESM", 3, line=1, cex=labelsize)
par(new = "TRUE",plt = c(0.52,0.9,0.67,0.95),las = 1,cex.axis = 1)
# top left
my.filled.contour3(restricted_lon, restricted_lat, rea_agg[[mss]], my_levels)
text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=lettersize)
mtext("ERA-interim", 3, line=1, cex=labelsize)
par(new = "TRUE",plt = c(0.91,0.93,0.67,0.95),las = 1,cex.axis = 1)
filled.legend(xcoords,ycoords, z,color.palette = rgb.palette,
	xlab = "",
	ylab = "", levels=my_levels, cex.axis=colbar_fontsize)

#Middle row
mss=2
model_data = do.call(rbind, lapply(hist_mom, function(el) el$moments.i[,mss]))

my_levels = make_zlims(c(model_data[,1], unlist(rea_agg[[mss]])), 2)
z = matrix(model_data[,1], n_lons, n_lats)

wh.palette=colorRampPalette(c('white','red2'),interpolate='spline')


par(new = "TRUE", plt = c(0.08,0.48,0.35,0.63), 
	las = 1, cex.axis = 1, tck = -0.02, oma=c(1,1,2,1))
my.filled.contour3(restricted_lon, restricted_lat, 
		z=z, levels=my_levels, color.palette=wh.palette)
# text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=2)
# par(xpd=NA)
text(x=xcoords[1]-10, y=ycoords[11], adj=2, label="Std. Dev.", srt=90, 3,
	cex=labelsize, xpd=NA)
par(new = "TRUE",plt = c(0.52,0.90,0.35,0.63),las = 1,cex.axis = 1)
# top left
# my.filled.contour3(lonsub, latsub, rea_data, 
# 				levels=my_levels, color = wh.palette)
my.filled.contour3(restricted_lon, restricted_lat, z=rea_agg[[mss]],
 levels=my_levels, color=wh.palette)
# text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=2)
par(new = "TRUE",plt = c(0.91,0.93,0.35,0.63),las = 1,cex.axis = 1)
filled.legend(xcoords,ycoords, z, color.palette=wh.palette,
	xlab = "",
	ylab = "", levels=my_levels, cex.axis=colbar_fontsize)

#Bottom row
mss=3
model_data = do.call(rbind, lapply(hist_mom, function(el) el$moments.i[,mss]))
my_levels = make_zlims(c(model_data[,1], unlist(rea_agg[[mss]])), 2)
z = matrix(model_data[,1], n_lons, n_lats)
# my_levels = make_zlims(c(unlist(rea_agg[[mss]]), c(rea_data)), 1)
# z = rea_agg[[mss]]


par(new = "TRUE", plt = c(0.08,0.48,0.03,0.31), 
	las = 1, cex.axis = 1, tck = -0.02, oma=c(1,1,2,1))
my.filled.contour3(restricted_lon, restricted_lat, 
		z=z, my_levels)
# text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=2)
# par(xpd=NA)
text(x=xcoords[1]-10, y=ycoords[11], adj=2, label="Skewness", srt=90, 3,
	cex=labelsize, xpd=NA)
par(new = "TRUE",plt = c(0.52,0.90,0.03,0.31), las=1, cex.axis=1)
# top left
my.filled.contour3(restricted_lon, restricted_lat, 
		rea_agg[[mss]], my_levels)
# text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=2)
par(new = "TRUE", plt = c(0.91,0.93,0.03,0.31), las=1, cex.axis=1)
filled.legend(xcoords, ycoords, z, color.palette=rgb.palette,
	xlab = "",
	ylab = "", levels=my_levels, cex.axis=colbar_fontsize)
dev.off()

# 3 location plot histograms with reanalysis
pdf('rea_temp_locations.pdf', width=8, height=5)
plot.3locations()
dev.off()

# 3 location plot histograms without reanalysis
pdf('temp_locations.pdf', width=8, height=5)
plot.3locations(with.rea=FALSE)
dev.off()



