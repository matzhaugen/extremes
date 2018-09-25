setwd("~/Extremes")
a <- read.table("~/Extremes/data/outBoot100_nb")
a_big <- read.table("~/Extremes/data/outbig")
a_mn <- read.table("~/Extremes/data/outBoot100_mn")
a_mx <- read.table("~/Extremes/data/outBoot100_mx")
load("~/Extremes/data/hist.Rdata")
load("~/Extremes/data/CVout.Rdata")
load("~/Extremes/data/histSeasons.Rdata")
load("~/Extremes/data/histodata.Rdata")

# source("~/Extremes/shinyapp2/helper.R")

source("~/Extremes/r_code/helper.R")

# Plot map of quantile differences change from 1850-2100
day = c(1, 186)
year = 241
q1 = c(0.025, 0.25, 0.95)
q2 = c(0.05, 0.75, 0.975)
my_lon = c(-101,-82, -98) + 360
my_lat = c(50.1, 42.7, 35.3)
colormap = colorRampPalette(c("blue", "white", "red"))(28)
q_string = c("Low", "IQR", "High")
season = c("Winter", "Summer")
data=c('min','avg','max')
data=rep('avg',3)
maxabs=c(2, 6.5, 2)
maxabs=c(.5, 6.5, .5)
lettersize = 1.3
labelsize = 1.5
colbar_fontsize = 1


minmax_breaks = make_zlims(c(-0.4, 0.4))
avg_breaks = make_zlims(c(-10.4, 3.44))
n_stds = 3
# Another version of the plots
placement = list(list(c(0.05,0.32,0.55,0.97), c(0.37,0.64,0.55,0.97), c(0.69,0.96,0.55,0.97)),
				list(c(0.05,0.32,0.05,0.47), c(0.37,0.64,0.05,0.47), c(0.69,0.96,0.05,0.47)))
legend.placement = list(c(0.05,0.32,0.02,0.03), c(0.37,0.64,0.02,0.03), c(0.69,0.96,0.02,0.03))

pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/winter_summer_var_avg_relnorm.pdf", width=11, height=7)
plot.winter.summer.avg(rep('avg', 3), limits=c(-.6, .6), normalize=TRUE, relnorm=TRUE)
dev.off()
pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/winter_summer_var_mn.pdf", width=11, height=7)
plot.winter.summer.avg(rep('min', 3), limits=c(-.6, .6), normalize=TRUE, relnorm=TRUE)
dev.off()
pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/winter_summer_var_mx.pdf", width=11, height=7)
plot.winter.summer.avg(rep('max', 3), limits=c(-.6, .6), normalize=TRUE, relnorm=TRUE)
dev.off()
pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/winter_summer_var_avg_big.pdf", width=11, height=7)
plot.winter.summer.avg(rep('big', 3), limits=c(-.6, .6), normalize=TRUE, relnorm=TRUE)
dev.off()


# plot "winter_summer_qchange" - the change in each quantile from 1850-2100
colormap = colorRampPalette(c("white", "black"))(100)
cex.lab=2.5; cex.main=2.7
zlim=19
plot_sig_pts = FALSE
seasonstr = c('Winter', 'Summer')
main = c('Low', 'Median', 'High')
my_quantile = c(0.025, 0.5, 0.975)
data_type=c('min','avg','max')
pdf("~/Box Sync/Matz_ensembles/figures/winter_summer_qchange.pdf", width=11, height=6)
par(mfrow=c(2,3), oma=c(0,3,0,2))
k = 1
for (i in 1:length(day)) {
	for (j in 1:length(my_quantile)) {
		# screen(ind[k])
		latlon = get.latlon(data_type[j])
		mydt = getdt(my_quantile[j], day[i], data_type[j])
		mystats = cbind(mydt$d_ypoint, mydt$d_ylow, mydt$d_yhigh)
		notsig_pts = !((mystats[,2] < 0 & mystats[,3] > 0))
		if (i == 1) {
			par(mar=c(0.5,1,3,1))
			plot.contour(unique(latlon[,2]), unique(latlon[,1]), 
				matrix(mydt$d_ypoint, ncol=length(unique(latlon[,1]))), 
				zlim=19, nlevels=12, col='red', main=main[j], 
				xaxt='n', yaxt='n', cex.main=cex.main, lwd=2, labcex=cex.lab)
		} else {
			par(mar=c(0.5,1,1,1))
			plot.contour(unique(latlon[,2]), unique(latlon[,1]), 
				matrix(mydt$d_ypoint, ncol=length(unique(latlon[,1]))), 
				zlim=19, nlevels=12, col='red', 
				xaxt='n', yaxt='n', lwd=2, labcex=cex.lab)
		}
		text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=3)
		if (j==1) mtext(seasonstr[i], 2, line=1, ylbias=0, cex=cex.lab)
		if (plot_sig_pts) points(latlon[notsig_pts,2], latlon[notsig_pts,1], 
			ch=3, col='gray')
		k = k + 1
	}
}
dev.off()


# Plot mean/std/skewness maps
winter_idx = floor(as.numeric(1) / 15) + 1
summer_idx = floor(as.numeric(186) / 15) + 1
hist_mom_winter = lapply(hist.data, function(el) {el[[winter_idx]]})
hist_mom_summer = lapply(hist.data, function(el) {el[[summer_idx]]})

# v4 of 3-moments plot, one for each season, with colbar_fontsize adjustable
lettersize = 1.5
labelsize = 2
colbar_fontsize = 1

left.windows = list(c(0.08,0.45,0.67,0.95), c(0.08,0.45,0.33,0.61), c(0.08,0.45,0.01,0.29))
right.windows = lapply(left.windows, function(el) el + c(0.47, 0.47, 0, 0))
left.titles = c("Pre-Industrial", "", "")
right.titles = c("Change", "", "")
ylabs = c("Mean", "Std. Dev.", "Skewness")
filenames = c("~/Box Sync/Matz_ensembles/Paper_1/meanStdSkew_Winter_wcrosses.pdf", 
	"~/Box Sync/Matz_ensembles/Paper_1/meanStdSkew_Summer_wcrosses.pdf")
# filenames = c("~/Box Sync/Matz_ensembles/figures/meanStdSkew_Winter_nonchange.pdf", 
# 	"~/Box Sync/Matz_ensembles/figures/meanStdSkew_Summer_nonchange.pdf")
hist_mom = list(hist_mom_winter, hist_mom_summer)
lats = unique(latlon_small[,1])
lons = unique(latlon_small[,2])
n_lats = length(lats)
n_lons = length(lons)
for (j in 1:2){
	pdf(filenames[j], width=7, height=8)
	# png(filenames[j], width=700, height=800)
	for (i in 1:3) {
		mss=i
		z.i = matrix(do.call(rbind, lapply(hist_mom[[j]], function(el) el$moments.i[,mss]))[,1], ncol=n_lats)
		z.d = matrix(do.call(rbind, lapply(hist_mom[[j]], function(el) el$moments.d[,mss]))[,1], ncol=n_lats)
		z.d.low = matrix(do.call(rbind, lapply(hist_mom[[j]], function(el) el$moments.d[,mss]))[,2], ncol=n_lats)
		z.d.high = matrix(do.call(rbind, lapply(hist_mom[[j]], function(el) el$moments.d[,mss]))[,3], ncol=n_lats)
		sig_pts = z.d.low*z.d.high > 0
		plot.mini.panel(lons, lats, z.i, my_lon, my_lat,
			main.window=left.windows[[i]], add=i!=1, main=left.titles[[i]], ylab=ylabs[i], labelsize=labelsize, 
			lettersize=lettersize, colbar_fontsize=colbar_fontsize)
		plot.mini.panel(lons, lats, z.d,  my_lon, my_lat,
			main.window=right.windows[[i]], sig_pts=sig_pts, add=TRUE, main=right.titles[[i]], ylab='', labelsize=labelsize, 
			lettersize=lettersize, colbar_fontsize=colbar_fontsize)
	}
	dev.off()
}

# Check for crossings in the quantile estimates in the study area - takes a minute to run
latlons = outer_concat(la, lo)
crossings = apply(latlons, 1, function(pos) {
	pos_idx = getPos_idx(pos[1], pos[2])
	myquantiles = getq(pos_idx, 'avg')
	getcrossings(myquantiles$point)
})



####################
## New plot ########
####################
pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/dq_vs_day_norm.pdf", width=9, height=6)
par(mfrow=c(3,3), oma=c(3,4,3,4), mar=c(3,2,1,1), xpd=FALSE)
labelsize=1.5
ypos = rep(1,3)
# ypos = rep(2,3)
# ypos = c(0.9, 0.28, .35)
n_lat = length(my_lat)
for (i in 1:n_lat) {
	pos_idx = getPos_idx(round(my_lat[i]), round(my_lon[i]))
	plot_dist_diff(pos_idx, ylim=c(0.6, 1.4), div_by_iqr=FALSE, div_by_initial=TRUE)
	text(380, ypos[i], letters[i], cex=3, xpd=NA, adj=c(0, .5))
}
mtext(text=" Day of the year",side=1,line=1,outer=TRUE, cex=labelsize)
# mtext(text="Temperature (deg C)",side=2,line=1,outer=TRUE, cex=labelsize)
mtext(text="Normalized (Unitless)",side=2,line=1,outer=TRUE, cex=labelsize)
# mtext(text="c                    b                    a",side=4,line=1,outer=TRUE, cex=labelsize*1.2)
# text(380, ypos[i], letters[i], cex=3, xpd=NA, adj=c(0, .5))
# text(x=-1020, y=2.1, adj=2, label="Temperature (Deg C)", srt=90, 3,
# 	cex=2.2, xpd=NA)
mtext(text="Low                                IQR                               High",side=3,line=0,outer=TRUE, cex=labelsize)
dev.off()

#####################
#### NEW PLOT ########
# plot first day for each year when 0.05th quantile is above freezing
####################
i = 2
THRESHOLD = 0
THRESHOLD2 = -2.2
QUANTILES_TO_PLOT = c(0.05, .25, 0.5)
data_type = 'avg'
pos_idx = getPos_idx(round(my_lat[i]), round(my_lon[i]))
m_quantiles = getq(pos_idx, 'avg')
q_mat = matrix(m_quantiles[[1]], 73, 8)
x = rep(day_of_year, n_years)
t = c(t(matrix(rep(1:n_years, days_per_year), ncol=days_per_year)))
x.int.basis = as.matrix(pbs(x, df=3))
x.main.basis = as.matrix(pbs(x, df=15))
t.basis = ns(t, df=4)
X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis)
start_idx = B * np * nq * (pos_idx - 1) + 1
end_idx = start_idx + B * np * nq
# my_pixel = updateData(data_type, start_idx - 1, end_idx)
my_pixel = getData(data_type, start_idx:end_idx)
my_beta = array(my_pixel, dim=c(np, nq, B))
quantiles = array(apply(my_beta[,q %in% QUANTILES_TO_PLOT,], 3, function(bi) {
	bi[1, ] = bi[1, ] - KtoC
	X %*% bi
}), dim=c(dim(X)[1], length(QUANTILES_TO_PLOT), B))
ypoint = apply(quantiles, c(1,2), quantile, 0.5)
cex.lab = 1.5
q_mat = array(c(ypoint), dim = c(365, 250, 3))
when_freeze = rep(0, dim(q_mat)[2])
when_freeze_2 = rep(0, dim(q_mat)[2])
for (j in 1:length(QUANTILES_TO_PLOT)){
	for (i in 1:dim(q_mat)[2]) {
		when_freeze[i] = min(which(q_mat[,i,j] > THRESHOLD))
		when_freeze_2[i] = min(which(q_mat[,i,j] > THRESHOLD2))
	}
	if (j == 1){
		par(oma=c(0,1,0,0))
		plot(1850:2099, when_freeze, type='l',
			ylab='Quantile crossing time (days since Jan 1)', xlab='Year', col=j,
			ylim=c(0, 100), cex.lab=cex.lab)
		lines(1850:2099, when_freeze_2, lty=2, col=j)	
	} else {
		lines(1850:2099, when_freeze, lty=1, col=j)	
		lines(1850:2099, when_freeze_2, lty=2, col=j)	
	}
}
	
ticks = c(31, 59, 90)
tick_labels = c(15, 45, 75)
axis(4, at=ticks, labels=c("", "", ""))
axis(4, at=tick_labels, labels=c("Jan", "Feb", "Mar"), tick=FALSE)


####################
## New plot: 3 locations quantiles  ########
####################
pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/q_vs_day.pdf", width=9, height=6)
par(mfrow=c(3,3), oma=c(3,4,3,4), mar=c(3,2,1,1), xpd=FALSE)
labelsize=1.5
ypos = rep(1,3)
# ypos = rep(2,3)
# ypos = c(0.9, 0.28, .35)
n_lat = length(my_lat)
for (i in 1:n_lat) {
	pos_idx = getPos_idx(round(my_lat[i]), round(my_lon[i]))
	plot_dist_quantiles(pos_idx)
	text(380, ypos[i], letters[i], cex=3, xpd=NA, adj=c(0, .5))
}
mtext(text=" Day of the year",side=1,line=1,outer=TRUE, cex=labelsize)
# mtext(text="Temperature (deg C)",side=2,line=1,outer=TRUE, cex=labelsize)
mtext(text="Temperature (degrees C)",side=2,line=1,outer=TRUE, cex=labelsize)
# mtext(text="c                    b                    a",side=4,line=1,outer=TRUE, cex=labelsize*1.2)
# text(380, ypos[i], letters[i], cex=3, xpd=NA, adj=c(0, .5))
# text(x=-1020, y=2.1, adj=2, label="Temperature (Deg C)", srt=90, 3,
# 	cex=2.2, xpd=NA)
mtext(text="Low                                IQR                               High",side=3,line=0,outer=TRUE, cex=labelsize)
dev.off()

models = do.call(rbind, list(c(5, 3, 3),
			  c(7, 3, 3),
			  c(10, 3, 3),
			  c(10, 3, 4),
			  c(12, 3, 4),
			  c(15, 3, 4),
			  c(15, 3, 5),
			  c(15, 5, 5),
			  c(18, 5, 5)))
colnames(models) = c("Seasonal", "Seasonal-Int.", "Temporal")
rownames(models) = c(1:9)
library(xtable)
xtable(models, digits=0, caption="Degrees of freedom in the spline basis for each 
	independent variable including their interaction term. ")


# Plot a conceptual picture describing why we are looking at 
probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975)
col=c(2,2,3,3,4,4)
x = seq(0,3, length=1000)
pdf('paper/schema.pdf', width=7, height=5)
par(mar=c(6,4,4,2))
plot(x, dweibull(x, 2, scale = 1), type='l', lwd=3, xlab="Value", 
		ylab="Probability density", yaxt="n", xaxt="n", cex.lab=2, ylbias=1, bty='n', mgp=c(1,0,0))
abline(v=0);abline(h=0)
quant = qweibull(probs, 2, scale=1)
for (i in 1:length(probs)) {
	abline(v=quant[i], col=col[i])
}
arrows(quant[1], 0.1, quant[2], 0.1, code=3, col=2, cex=0.5, length=0.05)
text(quant[2], 0.1, labels=TeX('$\\Delta q_{l}$'), adj=c(0,0), cex=1.8)
arrows(quant[3], 0.4, quant[4], 0.4, code=3, col=3, cex=0.5, length=0.05)
text((quant[3] + quant[4]) / 2, 0.4, labels=TeX('IQR'), adj=c(0.5,0), cex=2)
arrows(quant[5], 0.2, quant[6], 0.2, code=3, col=4, cex=0.5, length=0.05)
text(quant[6], 0.2, labels=TeX('$\\Delta q_{h}$'), adj=c(0,0), cex=2)
dev.off()



##############################
##### Cross-validation part #####
##############################
# load("~/Extremes/data/CV.Rdata")
load("~/Extremes/data/CV_t.Rdata") #This cross-validates time

# Plot CV
range_d = c(0, 0.008)
t_mean_tr = as.vector(unlist(lapply(T, function(t) colMeans(t)[1])))
t_mean_te = as.vector(unlist(lapply(T, function(t) colMeans(t)[2])))
t_sd_tr = as.vector(unlist(lapply(T, function(t) apply(t,2, sd)[1])))
t_sd_te = as.vector(unlist(lapply(T, function(t) apply(t,2, sd)[2])))
cex.lab = 1.2
m_lat = signif(lats[lat_idx], digits=3)
m_lon = signif(lons[lon_idx], digits=3)
pdf(paste("~/Box Sync/Matz_ensembles/Paper_1/figures/CV_t.pdf", sep=""), width=5, height=5)
# pdf(paste("~/Box Sync/Matz_ensembles/Paper_1/figures/CV.pdf", sep=""), width=5, height=5)
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
upper_year = 250
n_files = 50
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

# save(out, T, file="CVout.Rdata")
load("~/Extremes/data/CVout.Rdata")
h_tr = out[[1]][[1]]
h_te = out[[1]][[2]]
n_breaks=length(h_tr$breaks)
qhat_te = do.call(rbind, lapply(out, function(item) {
	item[[2]]$counts[1:(n_breaks-1)] / (5*250*10)
}))
qhat_tr = h_tr$counts[1:(n_breaks-1)] / (45*250*10)
den_mean = apply(qhat_te, 2, mean)
den_sd = apply(qhat_te, 2, sd)

pdf("Exeedence.pdf", width=5, height=5)
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

####################################
############# NEW PLOT################
####################################
# 3 Locations quantile fits
load("~/Extremes/data/CVout.Rdata")
my_lon = c(-101,-82, -98) + 360
my_lat = c(50.1, 42.7, 35.3)
n_years_begin = 1
n_files = 50

upper_year = 250
nobs_per_run = days_per_year * upper_year
idx_i = c(outer(1:(days_per_year * n_years_begin), (1:n_files- 1) * nobs_per_run, '+'))
idx_est = 1:(days_per_year) +  (days_per_year * floor(n_years_begin/2))  
# ex_y_3loc = list(); length(ex_y_3loc) = 3
pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/Estimate_3loc.pdf", width=12, height=4)
par(mfrow=c(1,3), oma=c(0,3,0,0))
ylabs = c("Temperature (deg C)", "", "")
text.pos.y = c(20,21,31)
for (i in 1:length(my_lon)) {
	
	# one_pixel = lapply(files, function(f) {
	# 	path_f = paste(path.to.files, f, sep="")
	# 	get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, -1)
	# })
	# data = unlist(lapply(one_pixel, function(ts) {
	# 				ts[1:(max_time)] - KtoC
	# 				}))
	pos_idx = getPos_idx(my_lat[i], my_lon[i])
	myquantiles = getq(pos_idx, 'avg')

	x_idx = 1:73
	y_hat_mid = getQuantile(pos_idx, 0.5, 'avg')[[1]][x_idx] 
	y_hat_low = getQuantile(pos_idx, 0.025, 'avg')[[1]][x_idx] 
	y_hat_high = getQuantile(pos_idx, 0.975, 'avg')[[1]][x_idx] 	

	ex_x = (idx_i - 1) %% days_per_year + 1
	# ex_y_3loc[[i]] = data[idx_i]
	random_sample = sample(1:length(ex_x), round(length(ex_x)/10))
	plot(ex_x[random_sample], ex_y_3loc[[i]][random_sample],
		xlab="Day of the year", ylab="", cex.lab=2,
		main="",
		ylim=range(ex_y_3loc[[i]]) * 1.15,
		col='gray', cex.axis=1.8)
	if (i==1) mtext(ylabs[i], side=2, line=3, cex=1.5)
	my_x = seq(1, days_per_year, by=5)
	lines(my_x, y_hat_mid, col=1, lwd=3)
	lines(my_x, y_hat_low, col=1, lty=2, lwd=3)
	lines(my_x, y_hat_high, col=1, lty=2, lwd=3)
	if (i ==1) legend('bottomleft', c('Median', '2.5/97.5th quantile'), lty=c(1,2), 
		col=1, lwd=c(2,2), cex=1.5)
	text(10, text.pos.y[i], c('a', 'b', 'c')[i], cex=3)

}
dev.off()

################################
############# End CV ###########
################################
################################
# PART 2: Plot the data ########
################################
# Run the below for the North American comaprison between CESM and ERA-interim 
load("~/Extremes/data/histSeasons.Rdata")
load("~/Extremes/data/histodata.Rdata")
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

season_index = 1 # winter
# season_index = 3 # summer
for (i in 1:3) {
	rea_moments_agg_mat[[i]] = as.matrix(do.call(rbind, lapply(rea_moments_agg[[season_index]], function(mlist) lapply(mlist, function(mm) mm[i]))))
}

# hist_mom = hist_mom_summer # change this for different seasons
# pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/model_reanalysis_world_JJA_agg.pdf", width=6,height=6)
hist_mom = hist_mom_winter # change this for different seasons
pdf("~/Box Sync/Matz_ensembles/Paper_1/figures/model_reanalysis_world_DJF_agg.pdf", width=6,height=6)
rea_agg = rea_moments_agg_mat
mss=1
lettersize = 1.2
labelsize = 1.7
colbar_fontsize = 1
model_data = do.call(rbind, lapply(hist_mom, function(el) el$moments.m[mss]))
my_levels = make_zlims(c(model_data[,1], unlist(rea_agg[[mss]])), 2)
n_lon_big = length(unique(latlon_big[,2]))
n_lat_big = length(unique(latlon_big[,1]))
z = matrix(model_data[,1], n_lon_big, n_lat_big)
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
xcoords = unique(latlon_big[,2])
ycoords = unique(latlon_big[,1])
# my.filled.contour3(restricted_lon, restricted_lat, 
# 		z=z, my_levels)
my.filled.contour3(xcoords, ycoords, z, my_levels)
text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=lettersize)
# par(xpd=NA)
text(x=xcoords[1]-10, y=ycoords[11], adj=2, label="Mean", srt=90, 3,
	cex=labelsize, xpd=NA)
mtext("CESM", 3, line=1, cex=labelsize)
par(new = "TRUE",plt = c(0.52,0.9,0.67,0.95),las = 1,cex.axis = 1)
# top left
my.filled.contour3(xcoords, ycoords, rea_agg[[mss]], my_levels)
text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=lettersize)
mtext("ERA-interim", 3, line=1, cex=labelsize)
par(new = "TRUE",plt = c(0.91,0.93,0.67,0.95),las = 1,cex.axis = 1)
filled.legend(xcoords,ycoords, z,color.palette = rgb.palette,
	xlab = "",
	ylab = "", levels=my_levels, cex.axis=colbar_fontsize)

#Middle row
mss=2
model_data = do.call(rbind, lapply(hist_mom, function(el) el$moments.i[,mss]))
# rea_data = rea[mss,,]

my_levels = make_zlims(c(model_data[,1], unlist(rea_agg[[mss]])), 2)
z = matrix(model_data[,1], n_lon_big, n_lat_big)

wh.palette=colorRampPalette(c('white','red2'),interpolate='spline')


par(new = "TRUE", plt = c(0.08,0.48,0.35,0.63), 
	las = 1, cex.axis = 1, tck = -0.02, oma=c(1,1,2,1))
my.filled.contour3(xcoords, ycoords, 
		z=z, levels=my_levels, color.palette=wh.palette)
# text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=2)
# par(xpd=NA)
text(x=xcoords[1]-10, y=ycoords[11], adj=2, label="Std. Dev.", srt=90, 3,
	cex=labelsize, xpd=NA)
par(new = "TRUE",plt = c(0.52,0.90,0.35,0.63),las = 1,cex.axis = 1)
# top left
# my.filled.contour3(lonsub, latsub, rea_data, 
# 				levels=my_levels, color = wh.palette)
my.filled.contour3(xcoords, ycoords, z=rea_agg[[mss]],
 levels=my_levels, color=wh.palette)
# text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=2)
par(new = "TRUE",plt = c(0.91,0.93,0.35,0.63),las = 1,cex.axis = 1)
filled.legend(xcoords,ycoords, z, color.palette=wh.palette,
	xlab = "",
	ylab = "", levels=my_levels, cex.axis=colbar_fontsize)

#Bottom row
mss=3
model_data = do.call(rbind, lapply(hist_mom, function(el) el$moments.i[,mss]))
# rea_data = rea[mss,,]
my_levels = make_zlims(c(model_data[,1], unlist(rea_agg[[mss]])), 2)
z = matrix(model_data[,1], n_lon_big, n_lat_big)
# my_levels = make_zlims(c(unlist(rea_agg[[mss]]), c(rea_data)), 1)
# z = rea_agg[[mss]]


par(new = "TRUE", plt = c(0.08,0.48,0.03,0.31), 
	las = 1, cex.axis = 1, tck = -0.02, oma=c(1,1,2,1))
my.filled.contour3(xcoords, ycoords, 
		z=z, my_levels)
# text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=2)
# par(xpd=NA)
text(x=xcoords[1]-10, y=ycoords[11], adj=2, label="Skewness", srt=90, 3,
	cex=labelsize, xpd=NA)
par(new = "TRUE",plt = c(0.52,0.90,0.03,0.31), las=1, cex.axis=1)
# top left
my.filled.contour3(xcoords, ycoords, 
		rea_agg[[mss]], my_levels)
# text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=2)
par(new = "TRUE", plt = c(0.91,0.93,0.03,0.31), las=1, cex.axis=1)
filled.legend(xcoords, ycoords, z, color.palette=rgb.palette,
	xlab = "",
	ylab = "", levels=my_levels, cex.axis=colbar_fontsize)
dev.off()

# 3 location plot histograms with reanalysis
pdf('~/Box Sync/Matz_ensembles/Paper_1/figures/rea_temp_locations.pdf', width=8, height=5)
plot.3locations()
dev.off()

# 3 location plot histograms without reanalysis
pdf('~/Box Sync/Matz_ensembles/Paper_1/figures/temp_locations.pdf', width=8, height=5)
plot.3locations(with.rea=FALSE)
dev.off()

