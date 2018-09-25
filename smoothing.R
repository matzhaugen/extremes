# file: smoothing.R
# Here I attempt to do quantile regression with a smoothing 
# spline in 2 dimensions where one dimension is periodic.
# The penalty function corresponds to the following covariance
# function in x, t, the two dimensions
# K(x,t) = materm(x) * 1 - 30t^2(1-t)^2


source("/project/moyer/mahaugen/r_code/helper.R")
library(binhf)
library(LatticeKrig)
ncdata <- nc_open('/project/moyer/mahaugen/rsriver/timefirst/trefht_4200.nc')
print(ncdata)
lons = signif(ncvar_get(ncdata, 'lon'), 3)
lons[lons>180] = lons[lons>180] - 360
lats = signif(ncvar_get(ncdata, 'lat'), 3)
nc_close(ncdata)

lat_min = 10
lat_max = 60
lon_min = -130
lon_max = -60
restricted_lat = lats[(lats < lat_max) & (lats > lat_min)]
restricted_lon = lons[(lons < lon_max) & (lons > lon_min)]
n_lats = length(restricted_lat)
n_lons = length(restricted_lon)
days_per_year = 365
upper_year = 250
q1 = 0.95
q2 = 0.975
q = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
nq = length(q)
B = 10 # number of bootstrap runs
max_time = days_per_year * upper_year

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
path.to.files = "../rsriver/timefirst/"
files = list.files("../rsriver/timefirst/")
n_files = length(files)
one_pixel = lapply(files, function(f) {
	path_f = paste(path.to.files, f, sep="")
	get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, -1)
})

tr_data = unlist(lapply(one_pixel, function(ts) {
			ts[1:(max_time)]
			}))
day_of_year = 1:365
x = rep(day_of_year, n_files*upper_year)
t = rep(c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year))), n_files)

ni = 250*365
lkout = LatticeKrig(cbind(x,t)[1:ni,], tr_data[1:ni])

library(fields)
n1 = 250
n2 = 365
data = array(data=tr_data, dim=c(365,250,50))
zmean = apply(data, c(1,2), mean)
theta = 1
t = seq(0,100, len=n1)
k1 = Matern(rdist(t, t), smoothness=1.5)
x = seq(0,40, len=n2)
r2 = rdist(x, x)
k2 = 1 - 30*r2^2*(1-r2)^2
svd1 = svd(k1)
svd2 = svd(k2)
d1 = svd1$d
d2 = svd2$d
L = 1 / (d2%*%t(d1) + theta)
T1 = t(svd2$v) %*% zmean %*% svd1$v
A = L * T1 # elementwise multiplication
B = svd2$u %*% A %*% t(svd1$u)
zhat = t(k2) %*% B %*% k1

#extract i^th row from k1 x k2 
i=200;j=250
kij = rep(k2[,i], n1) * c(t(matrix(rep(k1[j,], n2), ncol=n2)))
W1 = t(svd2$v) %*% matrix(kij, n2, n1) %*% svd1$v
lambda = svd2$u %*% (L * W1) %*% t(svd1$u)
pdf("plot1.pdf")
par(mfrow=c(2,2))
image.plot(1:365, 1:250, log(abs(lambda)), main="log-weights for day 200, \nyear 250",
	xlab="", ylab="")
image.plot(1:365, 1:250, zhat, main="Cond. mean estimate", xlab="", ylab="")
plot(zhat[,1], type='l', main="Year 1850", xlab="Day")
plot(zhat[200,], type='l', main="Day 200", xlab="Year")
dev.off()



# That took forever
# Let's try Doug's Krig
time.cov <- function (x1, x2 = NULL, theta = 1, C = NA, marginal = FALSE) 
{
    if (!is.null(x2)) {
        x2 <- x1
    }
    if (is.na(C[1]) & !marginal) {
        m1 = Matern(rdist(x1[,1], x2[,1]), smoothness=1.5)
		r2 = rdist(x1[,2],x2[,2])
		m2 = 1 - 30*r2^2*(1-r2)^2
		return(m1*m2)
    }
    if (!is.na(C[1])) {
        m1 = Matern(rdist(x1[,1], x2[,1]), smoothness=1.5)
		r2 = rdist(x1[,2],x2[,2])
		m2 = 1 - 30*r2^2*(1-r2)^2
		return(m1 * m2 * C[1])
    }
    if (marginal) {
        return(rep(1, nrow(x1)))
    }
}

Krig(x=cbind(x,t)[1:(365*250), ], Y = tr_data[1:(365*250)], 
	cov.function=time.cov)
#that didn't work either
# Let's calculate the spatial empirical covariance
# First look at stationarity assumption

d = round(seq(1,365, by=1))
y = round(seq(1, 250, by=1))
nd = length(d) # Days grid
ny = length(y) # years grid
cov.vs.dist = matrix(0, nd, ny)
origin.idx = c(outer(d[1]:(d[2]-1), (y[1]:(y[2]-1) - 1)*365, "+"))
multi.origin.idx = c(outer(origin.idx, (1:50 - 1)*365*250, "+"))
origin.data = tr_data[multi.origin.idx]
for (i.d in 2:length(d)) {
	print (i.d)
	for (i.y in 2:length(y)) {
		one.run.idx = c(outer(d[i.d-1]:(d[i.d]-1), (y[i.y-1]:(y[i.y]-1) - 1)*365, "+"))
		multi.run.idx = c(outer(one.run.idx, (1:50 - 1)*365*250, "+"))

		limited.data = tr_data[multi.run.idx]
		s = cor(origin.data, limited.data)
		cov.vs.dist[i.d, i.y] = as.numeric(s)
	}
}

image(d, y, cov.vs.dist)

ny = 250
nfiles = 50
nd = 365
std = array(0, dim=c(nd,nd,ny))
data = array(data=tr_data, dim=c(365,250,50))
for (i in 1:ny) {
	print(i)
	timeslice = t(as.matrix(data[,i,]))
	# std[,,i] = cov(timeslice)
	std[,,i] = t(timeslice) %*% (timeslice) 
	stp = t(apply(cbind(0:364, st[[i]]), 1, function(v) {
			shift(v[-1], v[1], dir="left")
		}))
}

ny = 250
nfiles = 50
nd = 365
sty = array(0, dim=c(ny,ny,nd))
data = array(data=tr_data, dim=c(365,250,50))
for (i in 1:nd) {
	print(i)
	timeslice = t(as.matrix(data[i,,]))
	sty[,,i] = t(timeslice) %*% (timeslice) 
	# st[,,i] = cov(timeslice)
	stp = t(apply(cbind(0:249, st[[i]]), 1, function(v) {
			shift(v[-1], v[1], dir="left")
		}))
}
par(mfrow=c(1,2))
image(1:nd, 1:nd, apply(std, c(1,2), mean))
image(1:ny, 1:ny, apply(sty, c(1,2), mean))
