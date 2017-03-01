setwd("~/Extremes/data")
source("~/Extremes/r_code/helper.R")
library(splines)
library(pbs)
q  = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975)
n_pixels = 234
nq = length(q)
np = 35
B = 20
which_pixel = 175
n_years = 250
days_per_year = 365
day_of_year = 1:365


start_idx = B * np * nq * (which_pixel - 1) + 1
end_idx = start_idx + B * np * nq
a = read.table("outBoot1")
my_pixel = a[start_idx:end_idx,]

my_beta = array(my_pixel, dim=c(np, nq, B))

x = rep(day_of_year, n_years)
t = c(t(matrix(rep(1:n_years, days_per_year), ncol=days_per_year)))
x.basis = as.matrix(pbs(x, df=6))
t.basis = ns(t, df=4)
X = model.matrix(~x.basis*t.basis)

bad = c(apply(my_beta, c(2,3), function(b) {
	all(b==0)
}))
if (any(bad)) printf("There was some coefs with all zeros")

d_years = 10
x = c(seq(1,365, 5))
y = c(seq(1,250, d_years))
Xidx = c(outer(x, (y-1) * 365, "+"))
Xsubsample = X[Xidx,] 
rm(t, x, x.basis, t.basis, X)
Yhigh = apply(my_beta, 2, function(a) {	
	z = apply(Xsubsample %*% a, 2, function(b) {
		c(matrix(b, ncol=n_years/d_years))
	})
	apply(z, 1, quantile, 0.95)
})
Ylow = apply(my_beta, 2, function(a) {
	z = apply(Xsubsample %*% a, 2, function(b) {
		c(matrix(b, ncol=n_years/d_years))
	})
	apply(z, 1, quantile, 0.05)
})
Ypoint = apply(my_beta, 2, function(a) {
	c(matrix(Xsubsample %*% a[,1], ncol=n_years/d_years))
})

# Take differences between lower quantiles(LQ), IQR and
# upper quantile(UQ), then LQ / IQR and UQ / IQR
dbeta = aperm(apply(my_beta, c(1,3), function(a) {
	c(a[2] - a[1], a[4] - a[3], a[6] - a[5], 
		(a[2] - a[1]) / (a[4] - a[3]),
		(a[6] - a[5]) / (a[4] - a[3]))
}), c(2,1,3))
dYhigh = apply(dbeta, 2, function(a) {	
	z = apply(Xsubsample %*% a, 2, function(b) {
		c(matrix(b, ncol=n_years/d_years))
	})
	apply(z, 1, quantile, 0.95)
})
dYlow = apply(dbeta, 2, function(a) {
	z = apply(Xsubsample %*% a, 2, function(b) {
		c(matrix(b, ncol=n_years/d_years))
	})
	apply(z, 1, quantile, 0.05)
})
dYpoint = apply(dbeta, 2, function(a) {
	c(matrix(Xsubsample %*% a[,1], ncol=n_years/d_years))
})

qi = 2
y.low = matrix(dYlow[,qi], ncol = n_years / d_years)
y.high = matrix(dYhigh[,qi], ncol = n_years / d_years)
y.point = matrix(dYpoint[,qi], ncol=n_years / d_years)

slices = seq(0, 250, 40) / d_years + 1
n.colors = length(slices)
col.map = cm.colors(n.colors)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
col.map = jet.colors(n.colors)
col = 1
for (yi in slices) {
	z.point = y.point[,yi]
	z.low = y.low[,yi]
	z.high = y.high[,yi]
	if (yi == 1) {
		plot(x, z.point, col=col.map[col], ylim=range(c(y.point, y.low, y.high)), type="l")	
	} else {
		lines(x, z.point, col=col.map[col], lty=1)
	}
	lines(z.low, col=col.map[col], lty=2)
	lines(z.high, col=col.map[col], lty=2)
	col = col + 1
}

