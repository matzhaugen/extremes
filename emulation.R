# With this project we attempt to project current temperatures
# into the future using model simulations combined with reanalysis 
# data. The goes is to make a mapping between current and future 
# quantiles. For example, if we want to know the 95th quantile
# temperature in 2080 we first look at the current 95th quantile
# from reanalysis. 
source('emulationHelper.R')

# load('emu.Rdata')
q = c(.1, .18, 0.25, .35, .42, 0.5, .58, .65, 0.75, .82, 0.9)
q_tail = c(0.01, 0.1, .25, 0.5, .75)
q_low = q_tail
q_high = rev(1-q_tail)
q_all = c(q_low*q[1], q, q[length(q)] + q_high*(1-q[length(q)]))
cex.main=1.5
cex.lab=1.5
nq = length(q)
q_norm = c(.1, .5, .9)
n_files = 50
Sys.setenv(TZ='UTC')
# t_grace = 15

# Pick a location (lat, lon)
# my_idx=c(204, 175, 169)
# South Texas
mlat = 30 # -90 to 90
mlon = 360 - 100 # 0 to 360
# Chicago
# mlat = 42 # -90 to 90
# mlon = 360 - 90 # 0 to 360
# Seattle
# mlat = 48 # -90 to 90
# mlon = 360 - 122 # 0 to 360
# LA
# mlat = 34 # -90 to 90
# mlon = 360 - 118.5 # 0 to 360




# How well do we capture a quantile? Let's do cross validation with 5 hold-out
# simulations, 50 years and 90 days hold-out blocks within the set of 5 held out
# simulations.
# model.cv = cv(X, y=y_norm, lambda=1,
# 	train.func=function(x, y, lambda) {
# 		rq.fit.pfn(x=x, y=y, tau=.99, max.bad.fixup=20)$coefficients
# 	},
# 	loss.func=function(y,y.hat) {model_exceedence_loss(y, y.hat, tau=.99)}
# 	)

# loss = array(unlist(model.cv), dim=c(10, 5, 5))
# rf <- colorRampPalette(rev(brewer.pal(9,'Greys')))
# image.plot(sqrt(apply(loss, c(1,2), mean)^2), col=rf(32))
# mean(c(loss))

y = getModelData(mlat, mlon)
y.lens = getLENSData(mlat, mlon)
tseq= timeBasedSeq('18500101::20991231')
tseq = tseq[strftime(tseq, format = "%j")!='366'] # remove leep years
y.xts = xts(matrix(y, 365*250, 50), order.by=tseq)['1920/2099']
p.ryan = dim(y.xts)[2]
y.lens = y.lens['1920/2099']
p.lens = dim(y.lens)[2]
years_i = '1979/2016'

nfolds = 5
foldid = c(t(matrix(rep(seq(nfolds), p.ryan/nfolds), nfolds, p.ryan/nfolds)))[1:p.ryan]
foldid.lens = c(t(matrix(rep(seq(nfolds), p.lens/nfolds), nfolds, p.lens/nfolds)))[1:p.lens]
which.run = 1:10
which.run.lens = 1:8
tseq_f = timeBasedSeq("20590101::20961231::d")
years_f = '2059/2096'
# years_f = '2019/2056'
tseq_f = tseq_f[strftime(tseq_f, format = "%j")!='366'] # remove leep years

# Import reanalysis data 
# and cross-validate to find the best set of predictors

reanalysis_y = getReanalysisData2(mlat, mlon)
rea.y = reanalysis_y[.indexhour(reanalysis_y) %in% 12]
dayofyear = as.numeric(strftime(time(rea.y), format = "%j"))
rea.y = rea.y[dayofyear!=366]
dayofyear = as.numeric(strftime(time(rea.y), format = "%j"))
rea.year = as.numeric(strftime(time(rea.y), format = "%Y"))
index = dayofyear + 365*(rea.year - min(rea.year))
rea.x = getPredictors(n_files=1, df.x=rea.x.df[1], df.t=rea.x.df[2], 
    df.xt=rea.x.df[3], year.range=range(rea.year), get.volc=TRUE, mlat=mlat)
# rea_x = x_one[dayofyear + 365*(year - min(year)), ]
rea_p = dim(rea.x)[2]
n_quantiles = length(q_norm)

# Cross-validate
nbins.days = 4
nbins.years = 5
rea.cv = cv(rea.x, y=rea.y, lambda=1,
	train.func=function(x, y, lambda) {
		rq.fit(x=x, y=as.numeric(y), tau=.9)$coefficients
	},
	loss.func=function(y,y.hat) {exceedence_loss(y, y.hat, nbins.days=4, nbins.years = 5) }
	)
loss = array(unlist(rea.cv), dim=c(nbins.days, nbins.years, 5))
avg_loss = apply(loss, c(1,2), mean)
mean(c(loss))
rea_fit = rq(as.numeric(rea.y)~rea.x-1, tau=q_norm)
summary(rea_fit)
pdf("figures/ReaFit.pdf")
plot.rq.fit(dayofyear, rea.y, rea_fit)
legend("topleft", as.character(c(min(rea.year), max(rea.year))), col=c(2,3), lwd=2)
dev.off()

# Load data from main.R runs
# This file contains the variables qout and out.map, qout contains the quantile fit 
# while out.map is the resultant transformation on some data. It is the former variable
# that is computationally intensive, but the former also takes less storage space. 
# load('qout_simple.Rdata')
# load(paste('qout_', mlat,'_', mlon, '.pdf', sep=''))
# load('qout_SouthTexas.Rdata')
# load('qout_SouthTexas_1533_1033_310_1510.Rdata')

# load('qout_SouthTexas_10-13_1013_310_1010.Rdata')
# load('qout_SouthTexas_log_10-13_1013_310_1010.Rdata')
load('qout_SouthTexas_10-13_10-13_310_1010.Rdata') # Prettu good but too complex
load('qout_SouthTexas_10-13_6-13_310_1010.Rdata') 
load('qout_SouthTexas_1053_653_310_1010.Rdata') 
load('qout_SouthTexas_1063_643_310_1010.Rdata') # The best
load('qout_SouthTexas_1463_643_310_1010.Rdata')
load('qout_SouthTexas_1463_1463_310_1010.Rdata')

# Plot staistics of the two samples' qq-plot 
t.xts = time(y.xts)
mid.range = seq(2019, 2056)
late.range = seq(2059, 2096)
time.month = as.numeric(strftime(t.xts, format = "%m")) - 1
time.year = as.numeric(strftime(t.xts, format = "%Y"))
season = list(c(11, 0, 1), c(2,3,4), c(5, 6, 7), c(8,9,10))
period = list(mid.range, late.range)
qqryan = model_qqstats(map.ryan, 
	# y.xts[time.month %in% c(12,1,2) & time.year %in% late.range,], indexmon=c(12,1,2))
	y.xts[time.year %in% late.range,])
qqlens = model_qqstats(map.lens, 
	y.lens[time.month %in% c(11,0,1) & time.year %in% late.range,], indexmon=c(11,0,1))
qnames = c("Mean Diff. (deg C)", "Median Diff. (deg C)", "Max Diff. (deg C)")
pdf("figures/qqstats.pdf", width=9, height=3)
par(mfrow=c(1,3))
for (j in 1:3){
	comp.histlike(list(qqryan[,j], qqlens[,j]), breaks=8, xlab=qnames[j],
		ylab="Density", cex.lab=1.5)
}
dev.off()
j=1
# y.hat = map.lens[[1]]$rea.y.f[.indexmon(map.lens[[1]]$rea.y.f)%in%season[[j]],1]
# y = y.lens[time.month %in% season[[j]] & time.year %in% late.range,]
y.hat = map.lens[[1]]$rea.y.f[(.indexmon(map.lens[[1]]$rea.y.f))  %in% season[[j]],1]
y = y.lens[(.indexmon(y.lens)) %in% season[[j]] & (.indexyear(y.lens) + 1900) %in% late.range,]
pdf('figures/test.pdf')
marginal.summary.plot(y, y.hat)
dev.off()
pheight = 500
pwidth = 660
# Lens out-sample tests using ryans model
for (k in seq(nfolds)) {		
	png(paste('figures/EmuProjectionRyanOnLens_',k,'.png', sep=''), width=pwidth, height=pheight)
	seasonal.marginal.summary.plot(y.lens[years_f,], map.ryan.lens[[k]]$rea.y.f[,1],
		season, season.names)
	dev.off()
}

# Lens in-sample tests
for (k in seq(nfolds)) {		
	png(paste('figures/EmuProjectionLens_',k,'.png', sep=''), width=pwidth, height=pheight)
	seasonal.marginal.summary.plot(y.lens[years_f,], map.lens[[k]]$rea.y.f[,1],
		season, season.names)
	dev.off()
}

# Ryans in-sample tests
season.names = c('Winter', 'Spring', 'Summer', 'Autumn')
for (k in seq(nfolds)) {		
	png(paste('figures/EmuProjectionRyan_',k,'.png', sep=''), width=pwidth, height=pheight)
	seasonal.marginal.summary.plot(y.xts[years_f,], map.ryan[[k]]$rea.y.f[,1],
		season, season.names)
	dev.off()
}

# Plot all out-of-sample projections on a qq-plot
y.f = y.xts[years_f, ]
rea.y.f = xts(do.call(cbind, lapply(map.ryan, function(x) as.matrix(x$rea.y.f))), order.by=time(map.ryan[[1]]$rea.y.f))
y.lens.f = y.lens[years_f, ]
rea.y.lens.f = xts(do.call(cbind, lapply(map.lens, function(x) as.matrix(x$rea.y.f))), order.by=time(map.lens[[1]]$rea.y.f))

png(paste('figures/EmuProjectionRyan_all.png', sep=''), width=pwidth, height=pheight)
seasonal.marginal.summary.plot(y.f, rea.y.f, season, season.names)
dev.off()

png(paste('figures/EmuProjectionLens_all.png', sep=''), width=pwidth, height=pheight)
seasonal.marginal.summary.plot(y.lens.f, rea.y.lens.f, season, season.names)
dev.off()

# Make boxplots of the loss separated by Winter/Summer, Mid/Late periods,
# and quantiles
season.days = list(c(335:365, 1:59), c(60:151), c(152:242), c(243:334))
all.days = 1:365; all.years = 1:180;
map.ryan.all = list(map.ryan.mid, map.ryan)
loss = array(0, dim=c( length(season), length(map.ryan.all), dim(y.xts)[2],length(q_all)))
temploss = array(0, dim=c( length(season), length(map.ryan.all), dim(y.xts)[2],length(q_all)))
qqstats = array(0, dim=c(length(season), length(map.ryan.all),dim(y.xts)[2],3))
for (i in 1:length(season)) {
	print(i)
	for (j in 1:length(period)) {
		idx = time.month %in% season[[i]] & time.year %in% period[[j]]		
		temploss[i,j,,] = do.call(rbind, lapply(seq(nfolds), function(k) {
			y.hat = make.quantile.surfaces(qout.ryan[[k]])
			l = 0
			y.hat.mean = apply(y.hat, 3, function(y.hat.q) mean(c(y.hat.q)[idx]))
			loss.ij = apply(y.hat, 3, function(y.hat.q) {
				l <<- l + 1
				y.hat.sub = c(y.hat.q)[idx]
				y = map.ryan.all[[j]][[k]]$rea.y.f
		    	y.sub = y[.indexmon(y) %in% season[[i]],] # Test on projected out-of-sample data
		    	# y.sub = y.xts[idx, foldid==k] # Test on out-of-sample sctual data
		    	loss.q = exceedence_loss_1(y.sub, y.hat.sub) - q_all[l]
		    	# print(loss.q)
		    	dq2temp(loss.q, y.hat.mean, q_all[l], q_all)
		    })		
			loss.ij
	    }))		
		qqstats[i,j,,] = model_qqstats(map.ryan.all[[j]], y.xts[idx,], season[[i]])
	}
}

# png('figures/test.png', 500, 500)
# marginal.summary.plot(y[.indexmon(y) %in% season[[i]],], y.xts[idx, foldid==k])
# dev.off()
# Loss vs. Quantile
library(reshape2)
dfloss = melt(temploss, varnames=c('Season', 'Period', 'Run', 'Quantile'),
value.name='Loss')
dfloss[,'Quantile'] = as.factor(q_all[dfloss[,'Quantile']])
dfloss[,'Season'] = as.factor(season.names)[dfloss[,'Season']]
dfloss[,'Period'] = factor(c('Mid', 'Late'), levels=c('Mid', 'Late'))[dfloss[,'Period']]
dfloss = dfloss[dfloss[,'Period']=='Late',]
dfloss = dfloss[dfloss[,'Quantile'] %in% c(0.01, .1, .5, .9, 0.99), ]
pdf('figures/boxplot.pdf', width=8, height=6)
ggplot(dfloss) + 
geom_boxplot(aes(x=Quantile, y=Loss, fill=Season)) +
# geom_boxplot(aes(x=Quantile, y=Loss, fill=Period)) +
scale_y_continuous(name=expression(paste(hat(tau) - tau, " [",degree*C,"]", sep=''))) + 
coord_flip() + 
# facet_wrap(~Season, ncol=2, nrow=2) +
theme(axis.title.x=element_text(size=25,face="bold"),
	axis.title.y=element_text(size=15))
dev.off()

# Median & Max diff.
qnames = c("Mean Diff.", "Median Diff.", "Max Diff.")
dfqqstats = melt(qqstats, varnames=c('Season', 'Period', 'Run', 'Stat'),
value.name='Loss')
dfqqstats[,'Stat'] = factor(qnames[dfqqstats[,'Stat']], levels=c("Median Diff.", "Max Diff."))
dfqqstats[,'Season'] = as.factor(season.names)[dfqqstats[,'Season']]
dfqqstats[,'Period'] = factor(c('Mid', 'Late'), levels=c('Mid', 'Late'))[dfqqstats[,'Period']]
dfqqstats = dfqqstats[dfqqstats[,'Stat'] %in% c("Median Diff.", "Max Diff."), ]
pdf('figures/QQbox.pdf')
ggplot(dfqqstats) + 
geom_boxplot(aes(x=Season, y=Loss), varwidth=TRUE) +
coord_flip() +
facet_wrap(~Period + Stat, nrow=2, ncol=2, scales="free_y", as.table=TRUE) 
dev.off()

# just compare histograms
# comp.hist(as.numeric(y.te[years_f, ]), as.numeric(out.map[[k]]$rea.y.f), ylab='Temperature (Deg C)', main="Daily Temperature Shift")
png(paste('figures/', 'QuantileDiff', mlat, '_', mlon,'.png', sep=''), width=800, height=500)
k=4
p1 = plot.quantile.change(qout.ryan[[k]],logplot=F, main=list(label="Sriver", vjust=3))
p2 = plot.quantile.change(qout.lens[[k]],logplot=F, main=list(label="LENS", vjust=3))
grid.arrange(p1,p2, ncol=2)
dev.off()

pdf(paste('figures/', 'QuantileDouble', mlat, '_', mlon,'.pdf', sep=''), width=12, height=6)
plot.double.transform(qout.ryan, qout.lens, logplot=T)
dev.off()

# Get tauhat segmented by year and day of year
nbins.days=1; nbins.years=4
out.byyear = lapply(seq(nfolds), function(k) {apply(make.quantile.surfaces(qout.ryan[[k]]), 3, function(y.hat) {
 	exceedence_loss(y.xts[, foldid.lens==k], y.hat, nbins.days=nbins.days, nbins.years=nbins.years)
 })})
nbins.days=4; nbins.years=1
out.byseason = lapply(seq(nfolds), function(k) {apply(make.quantile.surfaces(qout.ryan[[k]]), 3, function(y.hat) {
 	exceedence_loss(y.xts[, foldid.lens==k], y.hat, nbins.days=nbins.days, nbins.years=nbins.years)
 })})

q_all = qout.ryan[[1]]$q_all
pdf("figures/ExceedenceDiffYear.pdf", width=6, height=6)
tauhat = apply(simplify2array(out.byyear), c(1,2), mean)
image.plot(x=q_all, y=seq(dim(tauhat)[1]), z=(apply(tauhat, 1, function(x) x-q_all)), 
	zlim=c(-.0035,.0035), col=rwb.palette(50), xlab="Quantile", ylab="Time (1920-2099)")
dev.off()

# Each panel as a season
pdf("figures/ExceedenceDiffSeason.pdf", width=6, height=6)
tauhat = apply(simplify2array(out.byseason), c(1,2), mean)
image.plot(x=q_all, y=seq(dim(tauhat)[1]), z=(apply(tauhat, 1, function(x) x-q_all)), 
	zlim=c(-.01,.01), col=rwb.palette(50), xlab="Quantile", ylab="Season")
dev.off()

# Plot absolute sums of quantile bias for all quantiles for each projection
nbins.days=4; nbins.years=4
indiv = do.call(rbind, lapply(seq(nfolds), function(k) apply(make.quantile.surfaces(qout.ryan[[k]]), 3, function(y.hat) {
 	exceedence_loss_i(y.xts[, foldid==k], y.hat, nbins.days, nbins.years)
	})))
exceed.ryan = array(indiv, dim=c(16, 50, 21))
# taudiff = apply(exceed.ryan, c(1,2), function(ex) abs(ex - q_all))
taudiff = apply(exceed.ryan, c(1,2), function(ex) ex - q_all)
taudiffmean = array(t(apply(taudiff, c(1,2), mean)), dim=c(nbins.days,nbins.years,21))
loss.byyear = apply(taudiffmean, c(1,3), mean)
pdf("figures/ExceedenceDiffSeason2.pdf")
image.plot(x=q_all, z=t(loss.byyear), col=rwb.palette(50)
	, zlim=c(-0.01, 0.01)
	)
dev.off()
loss.byseason = apply(taudiffmean, c(2,3), mean)
pdf("figures/ExceedenceDiffYear2.pdf")
image.plot(x=q_all, z=t(loss.byseason), col=rwb.palette(50)
	,zlim=c(-.0035,.0035)
	)
dev.off()

png('figures/DetailedComparison.png', width=860, height=860)
k=4
par(mfrow=c(2,2), oma=c(0,1,0,0), xpd=NA)
plot.quantile.fit(y.lens[,foldid.lens==k], qout.lens[[k]], tau=c(.1, .9), 
	plotit=TRUE, normit=F, by.year=T, main="LENS", ylab="Temperature (deg C)")
plot.quantile.fit(y.xts[,foldid==k], qout.ryan[[k]], tau=c(.1, .9), 
	plotit=TRUE, normit=F, by.year=T, main="Sriver", ylab="")
plot.quantile.fit(y.lens[,foldid.lens==k], qout.lens[[k]], tau=c(.1, .9), 
	plotit=TRUE, normit=F, by.year=F, main="LENS", ylab="Temperature (deg C)")
plot.quantile.fit(y.xts[,foldid==k], qout.ryan[[k]], tau=c(.1, .9), 
	plotit=TRUE, normit=F, by.year=F, main="Sriver", ylab="")
dev.off()


colMeans(do.call(rbind, tauhat))
sd.tau = apply(do.call(rbind, tauhat), 2, sd)
bias.tau = qout[[1]]$q_all - colMeans(do.call(rbind, tauhat))
abs(bias.tau) > sd.tau


# Compare initial distributions
pdf('figures/MarginalComparison.pdf')
comp.histlike(list(as.numeric(y.lens[years_i,]), 
	as.numeric(y.xts[years_i,]),
	as.numeric(rea.y)), breaks=60, ylab="Density", xlab="Temperature (deg C)", lwd=2)
legend('topleft', c('LENS', 'Sriver', 'Reanalysis'), col=c(1,2,3), lwd=2)
dev.off()


# Plot the fit of the whole time span of the training set
yqhat = make.quantile.surfaces(qout.lens[[4]])
plot.surface(yqhat[,,qout.lens[[1]]$q_all==.5])
y.lens.k = y.lens[,foldid.lens==2]
doy = as.numeric(strftime(time(y.lens.k), '%j'))
year = as.numeric(strftime(time(y.lens.k), '%Y')) - 1920 + 1
co2 = read.table("RCP85_MIDYR_CONC.DAT", skip=38, header=TRUE)
b2 = co2[co2$YEARS>=year.range[1] & co2$YEARS<=year.range[2], 2]
plot(year, y.lens.k[,1])
lines(year, lm(y.lens.k[,1] ~ b2[year])$fit, col=2, lwd=2)

# plot just the initial years, target time span, with a fitted surface
x = getPredictors(10, -1, 3, get.volc=T)
plot.rea.fit(map.lens[[2]], y.lens[years_i,foldid.lens==2], 18, normit=T)
y.rea.norm = xts(matrix(map.lens[[2]]$rea.y.norm, 365*38, p.lens/nfolds), order.by=tseq_f)
doy = as.numeric(strftime(time(y.rea.norm), '%j'))
plot(rep(doy, 8), as.numeric(map.lens[[2]]$rea.y.norm))

# comp.hist(as.numeric(y.lens), as.numeric(rea.y), breaks=50)
# rea_fit_i = rq(as.numeric(rea.y)~rea.x-1, tau=q_norm)
# yqhat_rea = x.te %*% rea_fit_i$coef
# m = yqhat_rea[, 2]
# r = yqhat_rea[, 3] - yqhat_rea[, 1]
# plot(strftime(time(rea.y), '%j'), as.numeric(rea.y))
# lines(1:365, m[1:365], lwd=2, col=2)
# lines(1:365, yqhat_rea[1:365, 1], lwd=2, col=2)
# y_rea_norm_i = (y.te.i[, which.run] - m) / r
# rea_fit_f = rq(as.numeric(y.te[years_f,which.run])~x.te-1, tau=q_norm)
# yqhat_rea_f = x.te %*% rea_fit_f$coef
# m = yqhat_rea_f[, 2]
# r = yqhat_rea_f[, 3] - yqhat_rea[, 1]
# y_rea_norm_f = (y.te[years_f,which.run] - m) / r
# qqplot(as.numeric(y_rea_norm_f), as.numeric(rea.y.f))
# comp.hist(as.numeric(y_rea_norm_f), as.numeric(transform))

# Fit a GPD to the tails
# It turns out that the GPD does not change very much as a function
# of long term change. Therefore, we will skip it.
library(ismev)
# persp(dlow[seq(1,365, by=5), seq(1,250, by=5)], theta=25, phi=25)
lowtail = (y_norm - rep(c(yqhat_norm[,,q==.01]), 50) )
hightail = (y_norm - rep(c(yqhat_norm[,,q==.99]), 50) )
timestep = 1:length(y_norm)
year = (timestep-1) %/% 365 %% 250 + 1
par(mfrow=c(1,2))
doy.low = timestep[lowtail<0] %% 365; my.low = lowtail[lowtail<0]
idx = sample(1:length(my), 10000)
plot(doy.low[idx], my[idx], main='Low tail')
doy.high = timestep[hightail>0] %% 365; my.high = hightail[hightail>0]
myear = year[hightail>0]
plot(doy.high[idx], my[idx], main='High tail')
# outlow = gpd.fit(-lowtail[lowtail<0], ydat = as.matrix(year[lowtail<0]), sigl=1, threshold=0)
# outlow = gpd.fit(-lowtail[lowtail<0], ydat = as.matrix(year[lowtail<0]), shl=1, threshold=0)
# outlow = gpd.fit(-lowtail[lowtail<0], ydat = as.matrix(doy.low), threshold=0)
# outhigh = gpd.fit(my.high, ydat = as.matrix(myear), sigl=1, threshold=0)
# outhigh = gpd.fit(my.high, ydat = (as.matrix(doy.high)), shl=1, sigl=1, threshold=0)
outhigh = gpd.fit(my.high, threshold=0)
# gpd.diag(outlow)
gpd.diag(outhigh)


# Plot qq-plot of the reanalysis tails using the model fitted GPD parameters
library(TLMoments)
high_rea = y_rea_norm - yqhat_norm[dayofyear,year_i,q==.99]
low_rea = (y_rea_norm - yqhat_norm[dayofyear,year_i,q==.01])
high_rea = high_rea[high_rea>0]
low_rea = low_rea[low_rea<0]
par(mfrow=c(2,2))
hist(low_rea[low_rea<0])
hist(high_rea[high_rea>0])
plot(, ylab = "Empirical", xlab = "Model")
qqplot2(as.numeric(-low_rea), outlow)
qqplot2(as.numeric(high_rea), outhigh)
outhigh = gpd.fit(high_rea, threshold=0)
gpd.diag(outhigh)


