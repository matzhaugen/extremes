# This code is intended to swap the time dimension with the space dimensions 
# in order to access the time data faster

library(ncdf4)

n_cores=14
source("helper.R")
varname = 'TREFHTMX'
path.to.output = "../rsriver/mx/timefirst/"
path.to.input= "../rsriver/mx/spacefirst/"
files = list.files(path.to.input)

r = mclapply(files, function(f) {
	printf("Reading and writing file %s", f)
	path_f = paste(path.to.input, f, sep="")
	ncdata <- nc_open(path_f)
	data = ncvar_get(ncdata, varname)
	lons = ncvar_get(ncdata, 'lon')
	lats = ncvar_get(ncdata, 'lat')
	latdim = ncdata$dim[[1]]
	londim = ncdata$dim[[2]]
	timedim = ncdata$dim[[3]]
	timedim$unlim = FALSE
	nc_close(ncdata)

	var3d = ncvar_def(varname, 'K', list(timedim, londim, latdim), 
		longname="Reference height temperature max")

	path_f = paste(path.to.output, f, sep="")
	nc = nc_create(path_f, var3d)
	for( i in 1:length(lons)) {
		print(i)
		for (j in 1:length(lats)){
				temp = ncvar_put(nc, var3d, data[i, j,], start=c(1,i,j), count=c(-1,1,1) )
			}
	}
	nc_close(nc)
}, mc.cores=n_cores)

# Testing
# ncdata = nc_open(paste(path.to.output, f, sep=""))
# print(ncdata)
# nc_close(ncdata)