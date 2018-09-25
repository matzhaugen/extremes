library(ncdf4)
library(splines)
library(Hmisc)
library(pbs)
library(quantreg)
library(e1071)  
library(maps)
library(ggplot2)
library(doMC)
library(ggthemes)
library(xts)
library(fields)
library(latex2exp)
library(Rcpp)

load("data/latlon.Rdata")
load("data/latlon_big.Rdata")
# a <- read.table("data/outBoot100_nb") # uncomment to run 
load("data/hist.Rdata")

rwb.palette=colorRampPalette(c('darkblue','white','red2'),interpolate='spline')
wr.palette=colorRampPalette(c('white','red2'),interpolate='spline')

# a_mn <- read.table("data/outBoot100_mn")
# a_mx <- read.table("data/outBoot100_mx")
# The data is stored as a set of coefficients filling in the following order:
# p, q, bootstrap, pixel
# For large spatial window
lats = unique(latlon_big[,1])
lons = unique(latlon_big[,2])
mla = matrix(rep(lats, length(lons)), length(lats), length(lons))
mlo = matrix(rep(lons, length(lats)), length(lons), length(lats))
latlon_big = cbind(c(t(mla)), c(mlo))

# lats = unique(latlon[,1])
# lons = unique(latlon[,2])
# For small spatial window
la = round(restricted_lat)
lo = round(restricted_lon) + 360
# la = restricted_lat
# lo = restricted_lon + 360
mla = matrix(rep(la, length(lo)), length(la), length(lo))
mlo = matrix(rep(lo, length(la)), length(lo), length(la))
latlon_small = cbind(c(t(mla)), c(mlo))

lat_idx=1
lon_idx=1
pos_idx=1
n_lats = length(lats)
n_lons = length(lons)
KtoC = 273.15
q  = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
nq = length(q)

# mat_lat = matrix(rep(lats, n_lons), ncol=n_lons)
# mat_lon = matrix(rep(lons, n_lats), byrow=T, ncol=n_lons)
# latlon = cbind(c(t(mat_lat)), c(t(mat_lon)))

B = 100
n_years = 250
days_per_year = 365
day_of_year = 1:365
d_years = 40
d_days_per_year = 5

load('data/subsamples.Rdata')
# x = rep(day_of_year, n_years)
# t = c(t(matrix(rep(1:n_years, days_per_year), ncol=days_per_year)))
# x.int.basis = as.matrix(pbs(x, df=3))
# x.main.basis = as.matrix(pbs(x, df=15))
# t.basis = ns(t, df=4)
# X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis)

# coef = matrix(0, dim(X)[2], length(q))
# # bad = c(apply(my_beta, c(2,3), function(b) {
# #   all(b==0)
# # }))
# # if (any(bad)) printf("There was some coefs with all zeros")

# x_sub = c(seq(1,days_per_year, d_days_per_year))
# y_sub = c(seq(1,n_years, d_years), n_years)
# Xidx = c(outer(x_sub, (y_sub-1) * days_per_year, "+"))
# Xsubsample = X[Xidx,] 
np = dim(Xsubsample)[2]
# rm(t,x, t.basis, x.main.basis, x.int.basis)
# invisible(gc(verbose = FALSE))
breaks = seq(-70,50, .5) # for the histogram plot
# save(x_sub, y_sub, Xidx, Xsubsample, file='subsamples.Rdata')

get.latlon <- function(data_type='avg') {
    if (data_type != 'big') {
    latlon = latlon_small
  } else {
    latlon = latlon_big
  }
  latlon
}

get.n_pixels <- function(data_type='avg') {
  latlon = get.latlon(data_type)
  n_lats = length(unique(latlon[,1])) 
  n_lons = length(unique(latlon[,2])) 
  n_pixels = n_lats * n_lons
  n_pixels
}

# cppFunction('NumericVector rcpp_clip( NumericVector x, double a, double b){
#     return clamp( a, x, b ) ;
# }')

plot.winter.summer.avg <- function(data=rep("avg", 3), limits,  normalize=TRUE, relnorm=FALSE) {
  if (data[1] != 'big') {
    latlon = latlon_small
  } else {
    latlon = latlon_big
  }
  lats = unique(latlon[,1])
  lons = unique(latlon[,2])
  n_lats = length(lats) 
  n_lons = length(lons)
  n_pixels = n_lats * n_lons
  for (i in 1:length(day)) {
  # for (i in 1) {
    for (j in 1:length(q1)) {
      out_ddq = getddq(q1[j], q2[j], day[i], 250, data[j], n_pixels=n_pixels, normalize=normalize, relnorm=relnorm)
      if (missing(limits)) limits = out_ddq[[1]]
      zlims = make_zlims(limits, 2)
      par(new=!(i==1 && j==1), plt = placement[[i]][[j]], oma=c(1,1,4,1), xpd=FALSE)    
      
      # d_iqr = getddq(0.25, 0.75, day[i], year, data[j])
      point_over_sd = out_ddq[[1]] / out_ddq[[4]]
      sig_pts = abs(point_over_sd) > n_stds
      
      # my.filled.contour3( lons, lats, 
      #   z=matrix(out_ddq[[1]], ncol=n_lats), levels=zlims)
      col.info = get.color.info(limits, q=1.5)
      upper_limit = max(zlims)
      lower_limit = min(zlims)
      zvals = rcpp_clip(out_ddq[[1]], lower_limit, upper_limit)
      image(lons, lats, matrix(zvals, ncol=n_lats), 
        breaks=col.info$breaks, col=col.info$col, xaxt='n', yaxt='n', xlab="", ylab="")
      par(xpd=FALSE)
      map('world2', add=TRUE)
      par(xpd=NA)
      if (i == 1) title(q_string[j], cex.main=labelsize, line=1)
      
      points(latlon[!sig_pts,2], latlon[!sig_pts,1], pch=3, cex=1, lwd=1, col=alpha('gray', .5))
      text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=lettersize)
      if (j==1) mtext(season[i], 2, line=1, ylbias=0, cex=labelsize)
      if (i==1){
        par(new = "TRUE", plt = legend.placement[[j]])
        plot.window(xlim = range(zlims), ylim = c(0, 1), xaxs = "i", 
              yaxs = "i")
          rect(zlims[-length(zlims)], 0, zlims[-1L], 1, col = rgb.palette(length(zlims) - 1))
          axis(1, padj=-0.7)
      }
      par(xpd=FALSE)
    }
  }
}


###
# Plot 3 specific locations
### 
plot.3locations <- function(with.rea=TRUE) {
  title = c('a', 'b', 'c')
  season = c('Winter', 'Spring', 'Summer', 'Fall')
  p = list(list()); length(p) = length(season)
  my_colors = c(alpha('green', .4), alpha('red', .4), alpha('blue', .4))
  labelsize = 3
  lettersize = 3
  axissize = 1.6
  par(mfrow=c(2,3), oma=c(0,3,0,2), mar=c(3,1,3,1))
  for (j in c(1,3)){
    for (i in 1:length(my_lon)) {
      # nbreaks=length(h_i[[j]][[i]]$mids)
      dx = h_i[[j]][[i]]$mids[2] - h_i[[j]][[i]]$mids[1]
      lwd= c(2,2.2,1.8)[i]*0.9
      print(dx)
      distr_i = data.frame(Temperature=h_i[[j]][[i]]$mids,
                 Density=h_i[[j]][[i]]$density,
                 Type='Initial')
      distr_f = data.frame(Temperature=h_f[[j]][[i]]$mids,
                 Density=h_f[[j]][[i]]$density,
                 Type='Final')
      distr_rea = data.frame(Temperature=h_ra[[j]][[i]]$mids,
                 Density=h_ra[[j]][[i]]$density,
                 Type='Rea')
      distr = rbind(distr_i, distr_f, distr_rea)
      plot(distr_i$Temp, distr_i$Density, type='h', lwd=lwd, lend=1, 
        col=my_colors[1], ylim=range(distr$Density), 
        ylab="", xlab='Temperature', yaxt='n', cex.axis=axissize)
      
      if (with.rea){
        lines(distr_rea$Temp, distr_rea$Density, type='h', lwd=lwd, lend=1, 
          col=my_colors[3])
      } else {
        lines(distr_f$Temp, distr_f$Density, type='h', lwd=lwd, lend=1, 
        col=my_colors[2])
      }
      
      # base = ggplot(distr, aes(x=Temperature, y=Density, fill=Type))
      # p = base + geom_bar(stat="identity", alpha=0.5,
      #   position="identity", show.legend=TRUE) +
      #  theme_bw() + theme(axis.title.y = element_blank())
      
      #  # annotate('text', label=title, 
       #    y=max(h$density), 
       #    x = min(h$mids), 
       #    size=10)  
      if (with.rea) {
          if (i==1 && j==1) legend('topright', 
            c('Initial', 'Reanalysis'), col=my_colors[c(1,3)], text.col=my_colors[c(1,3)], lty=0, lwd=2, cex=1.5, bty='n')
        } else {
          if (i==1 && j==1) legend('topright', 
            c('Initial', 'Final'), col=my_colors[1:2],text.col=my_colors[1:2], lty=0, lwd=2, cex=1.5, bty='n')
        }
      text(distr_i$Temp[4], max(distr$Dens)*.9, title[i], cex=lettersize)
      if (i==1) title(ylab=season[j], xpd=NA, line=1, cex.lab=labelsize)
      # if (i==1) mtext(season[j], 2, line=1, ylbias=0, cex=2, srt=90)
    }
  }
}



plot_dist_quantiles <- function(pos_idx, ylim) { 
  cexaxis = 1.2
  # Minimums
  q_min = getQuantile(pos_idx,  0.05, 'avg')
  q_iqr = getQuantile(pos_idx, 0.5, 'avg')
  q_max = getQuantile(pos_idx, 0.95, 'avg')
  
  mnmx = c(unlist(q_min[1:3]), unlist(q_max[1:3]))
  ylims.iqr = c(unlist(q_iqr[1:3]))
  if (!missing(ylim)) {
    mnmx = ylim
    ylims.iqr = ylim
  }
  plot_q_or_dq(q_min, ylim=range(mnmx), cex.axis=cexaxis)
  # IQR
  plot_q_or_dq(q_iqr, ylim=range(ylims.iqr), cex.axis=cexaxis)
  # Maximums
  plot_q_or_dq(q_max, ylim=range(mnmx), cex.axis=cexaxis)
}


plot_dist_diff <- function(pos_idx, ylim, div_by_iqr=TRUE, div_by_initial=FALSE) { 
  cexaxis = 1.2
  # Minimums
  dq_min = getdq(pos_idx, 0.025, 0.05, 'avg', div_by_iqr, div_by_initial)
  dq_iqr = getdq(pos_idx, 0.25, 0.75, 'avg', div_by_iqr=FALSE, div_by_initial)
  dq_max = getdq(pos_idx, 0.95, 0.975, 'avg', div_by_iqr, div_by_initial)
  
  mnmx = c(unlist(dq_min[1:3]), unlist(dq_max[1:3]))
  ylims.iqr = c(unlist(dq_iqr[1:3]))
  if (!missing(ylim)) {
    mnmx = ylim
    ylims.iqr = ylim
  }
  plot_q_or_dq(dq_min, ylim=range(mnmx), cex.axis=cexaxis)
  # IQR
  plot_q_or_dq(dq_iqr, ylim=range(ylims.iqr), cex.axis=cexaxis)
  # Maximums
  plot_q_or_dq(dq_max, ylim=range(mnmx), cex.axis=cexaxis)
}
# Map of skewness difference between final and initial
plot.map <- function(x, y, z, colormap, zlim, lon_min=0, lon_max=0, lat_min=0, lat_max=0, ylab="", legend=TRUE, ...) {
  par(xpd=NA)
  if (legend) {
    image.plot(x, y, z=matrix(z, ncol=length(y)), col=colormap,
      zlim=c(-zlim, zlim), ylab=ylab, xlab="", legend.cex=2,
      ...)  
  } else {
    image(x, y, z=matrix(z, ncol=length(y)), col=colormap,
      zlim=c(-zlim, zlim), ylab=ylab, xlab="",
      ...)
  } 
  par(xpd=FALSE)
  map('world2', add=T)
  # map(xlim=c(lon_min, lon_max),
  #     ylim=c(lat_min, lat_max), add=TRUE, xaxt='n')
  # # map.axes()
}

plot.contour <- function(x, y, z, zlim,  ...) {
  
  contour(x, y, z=matrix(z, ncol=length(y)),
    zlim=c(-zlim, zlim), ...) 
  map('world2', add=TRUE)
}

# ifd - (1,3) initial(1), final(2) or difference(3)
# mss - mean(1), std(2), skew(3)
analyse_moments <- function(moments, ifd=1, mss=3) {
  mom_med = unlist(lapply(moments, function(mom) mom$mom[1,mss]))
  mom_lower = unlist(lapply(moments, function(mom) mom[2,mss]))
  mom_upper = unlist(lapply(moments, function(mom) mom[3,mss]))
  sig_pts = !(mom_lower < 0 & mom_upper > 0)
  cbind(mom_med, sig_pts)
}
# mss - mean(1), std(2), skew(3)
plot_x_season <- function(x_winter, x_summer, mss=1, max_x=2.5, max_dx=1.5, plot_sig_pts=FALSE) {
  
  colormap = colorRampPalette(c("blue", "white", "red"))(100)
  cex.lab = 3
  cex.main = 3
  m_mar = c(3,2,1,1)
  x.i = do.call(rbind, lapply(x_winter, function(el) el$moments.i[,mss]))
  x.f = do.call(rbind, lapply(x_winter, function(el) el$moments.f[,mss]))
  x.d = do.call(rbind, lapply(x_winter, function(el) el$moments.d[,mss]))
  sig.pts.d = !apply(x.d, 1, function(row) {!(row[2] < 0 & row[3] > 0)})
  x.i.summer = do.call(rbind, lapply(x_summer, function(el) el$moments.i[,mss]))
  x.f.summer = do.call(rbind, lapply(x_summer, function(el) el$moments.f[,mss]))
  x.d.summer = do.call(rbind, lapply(x_summer, function(el) el$moments.d[,mss]))
  sig.pts.d.summer = !apply(x.d.summer, 1, function(row) {!(row[2] < 0 & row[3] > 0)})
  par(mfrow=c(2,3), oma=c(0,4,3,2), mar=c(1,1,4,2))
  # Winter
  plot.map(unique(latlon[,2]), unique(latlon[,1]), x.i[,1], colormap, max_x, 
    lon_min, lon_max, lat_min, lat_max, main="1850", ylab="", cex.lab=cex.lab, 
    xaxt='n', yaxt='n', legend=TRUE, cex.main=cex.main, horizontal=TRUE, legend.mar=1)
  mtext("Winter",2, line=1, ylbias=0, cex=cex.lab)
  text(my_lon, my_lat, labels=letters[1:length(my_lat)], cex=2)
  plot.map(unique(latlon[,2]), unique(latlon[,1]), x.f[,1], colormap, max_x, 
    lon_min, lon_max, lat_min, lat_max, main="2100", legend=TRUE, 
    xaxt='n', yaxt='n', cex.main=cex.main, horizontal=TRUE, legend.mar=1)
  text(my_lon, my_lat, labels=letters[1:length(my_lat)], cex=2)
  plot.map(unique(latlon[,2]), unique(latlon[,1]), x.d[,1], colormap, max_dx, 
    lon_min, lon_max, lat_min, lat_max, main="Change", legend=TRUE, 
    xaxt='n', yaxt='n', cex.main=cex.main, horizontal=TRUE, legend.mar=1)
  text(my_lon, my_lat, labels=letters[1:length(my_lat)], cex=2)
  if (plot_sig_pts) points(latlon[sig.pts.d,2], latlon[sig.pts.d,1], pch=3, col='gray')
  # Summer
  par(mar=c(4,1,2,2))
  plot.map(unique(latlon[,2]), unique(latlon[,1]), x.i.summer[,1], cex.lab=cex.lab, colormap, max_x, 
    lon_min, lon_max, lat_min, lat_max, main="", ylab="", xaxt='n', yaxt='n', legend=FALSE, cex.main=cex.main)
  text(my_lon, my_lat, labels=letters[1:length(my_lat)], cex=2)
  mtext("Summer",2, line=1, ylbias=0, cex=cex.lab)
  plot.map(unique(latlon[,2]), unique(latlon[,1]), x.f.summer[,1], colormap, max_x, 
    lon_min, lon_max, lat_min, lat_max, main="", legend=FALSE, xaxt='n', yaxt='n', cex.main=cex.main)
  text(my_lon, my_lat, labels=letters[1:length(my_lat)], cex=2)
  plot.map(unique(latlon[,2]), unique(latlon[,1]), x.d.summer[,1], colormap, max_dx, 
    lon_min, lon_max, lat_min, lat_max, main="", legend=FALSE, xaxt='n', yaxt='n', cex.main=cex.main)
  text(my_lon, my_lat, labels=letters[1:length(my_lat)], cex=2)
  if (plot_sig_pts) points(latlon[sig.pts.d.summer,2], latlon[sig.pts.d.summer,1], pch=3, col='gray')
}

rgb.palette=colorRampPalette(c('darkblue','white','red2'),interpolate='spline')

filled.legend <- function(x = seq(0, 1, length.out = nrow(z)), 
    y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes, 
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, ...) 
{
  # modification of filled.contour by Carey McGilliard and Bridget Ferris
  # designed to just plot the legend
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
  #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #  on.exit(par(par.orig))
  #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  #  par(las = las)
  #  mar <- mar.orig
  #  mar[4L] <- mar[2L]
  #  mar[2L] <- 1
  #  par(mar = mar)
   # plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
        yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
        if (axes) 
            axis(4, ...)
    }
    else key.axes
    box()
}
    #
#    if (!missing(key.title)) 
#        key.title
#    mar <- mar.orig
#    mar[4L] <- 1
#    par(mar = mar)
#    plot.new()
#    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
#    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
#        stop("no proper 'z' matrix specified")
#    if (!is.double(z)) 
#        storage.mode(z) <- "double"
#    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
#        col = col))
#    if (missing(plot.axes)) {
#        if (axes) {
#            title(main = "", xlab = "", ylab = "")
#            Axis(x, side = 1)
#            Axis(y, side = 2)
#        }
#    }
#    else plot.axes
#    if (frame.plot) 
#        box()
#    if (missing(plot.title)) 
#        title(...)
#    else plot.title
#    invisible()
#}


filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page

  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
 # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
 # on.exit(par(par.orig))
 # w <- (3 + mar.orig[2]) * par("csi") * 2.54
 # par(las = las)
 # mar <- mar.orig
 plot.new()
 # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

plot.reanalysis.temp <- function(h, title) {
	reanalysis = data.frame(Temperature=h$mids, Density=h$density)
	par(mar=c(0,0,3,0))
	base = ggplot(reanalysis, aes(x=Temperature, y=Density))
	p = base + geom_bar(stat="identity", alpha=0.5,
		position="identity", show.legend=TRUE) +
	 theme_bw() + theme(axis.title.y = element_blank()) + 
	 annotate('text', label=title, 
	 		y=max(h$density), 
	 		x = min(h$mids), 
	 		size=10)
}

outer_concat <- function(la, lo) {
  mla = matrix(rep(la, length(lo)), length(la), length(lo))
  mlo = matrix(rep(lo, length(la)), length(lo), length(la))
  latlon = cbind(c(t(mla)), c(mlo))
  latlon
}

# this plots subplots in larger figure with image as the base plotting function. Just specify the main window,
# and the legend will automatically be added.
plot.mini.panel <- function(x, y, z, 
  my_lon, my_lat, main.window, sig_pts, add=FALSE, main='', ylab='',
  labelsize=1.7, lettersize=1.5, colbar_fontsize=1) {
  xrange = diff(range(x))
  yrange = diff(range(y))
  par(new=add,
    plt = main.window,   # using plt instead of mfcol (compare                                  
      oma=c(1,1,2,1))

  col.info = get.color.info(z)
  breaks = col.info$breaks
  image(x, y, z, breaks=breaks, col=col.info$col, xaxt='n', yaxt='n', xlab="", ylab="")
  if (!missing(sig_pts)){
    latlon = outer_concat(y, x)
    points(latlon[!c(sig_pts),2], latlon[!c(sig_pts),1], pch=3, cex=1, lwd=1, 
      col=alpha('gray', .95))
  }
  if (any(c(z) > max(breaks))) {print(z[c(z) > max(breaks)]); print("values out of color range! oops")}
  par(xpd=FALSE)
  map('world2', add=T)
  text(my_lon, my_lat, labels=letters[1:length(my_lat)],  cex=lettersize)
  
  text(x=x[1] - xrange*0.09, y=y[1] + yrange*0.42, adj=2, label=ylab, srt=90, 3,
    cex=labelsize, xpd=NA)
  mtext(main, 3, line=1, cex=labelsize)
  par(xpd=NA)
  par(plt = c(main.window[2] + 0.01,
              main.window[2] + 0.02, 
              main.window[3],
              main.window[4]), cex.axis = 1,las=1)
  plot.window(xlim = c(0, 1), ylim = range(breaks), xaxs = "i", 
       yaxs = "i")
  rect(0, breaks[-length(breaks)], 1, breaks[-1L], col = col.info$col)
  axis(side=4, hadj=.4)
  box()
}

get.color.info <- function(z, q = 1.5) {
  
  maxabs =  max(abs(z))
  if (range(z)[1] > 0) { 
    col = wr.palette(28)
    # breaks=seq(min(z), maxabs, length=29)
    breaks = make_zlims(z, q)
    col.palette = wr.palette
  } else {
    col = rwb.palette(28)
    # breaks=seq(-maxabs, maxabs, length=29)
    breaks = make_zlims(z, q)
    col.palette = rwb.palette
  } 
  # print(breaks)
  list(col=col, breaks=breaks, col.palette=col.palette)
}

make_zlims <- function(v, q=1, symmetric=TRUE, include.zero=FALSE) {
  al = 0.95
  ah = 1.05
	if (any(v < 0)) {
    if (symmetric) {
      maxabs = max(abs(v))*ah
      res = c(-(abs(seq(-maxabs^(1/q), 0, length=15)[-15])^q), 
        0, (seq(0, maxabs^(1/q), length=15)[-1])^q)  
    } else {
      res = c(-(abs(seq(-abs(min(v))^(1/q), 0, length=15)[-15])^q), 
        0, (seq(0, max(v)^(1/q), length=15)[-1])^q)  
    }
		
	} else {
    if (!include.zero) {
      res = seq(min(v)*al, (max(v)*ah)^(1/q), length=29)^(q)
    } else {
      res = seq(0, (max(v)*ah)^(1/q), length=29)^(q)
    }
	}
	res
}

my.filled.contour3 <- function(x, y, z, levels,ylab='', xlab='', main='', 
color.palette=rgb.palette, ...) {
  filled.contour3(x, y, z, color.palette=color.palette, 
  levels=levels,
  plot.title=title(main=main, xlab=xlab, ylab=ylab, line=1, cex.lab=2, cex.main=2),
    plot.axes={map('world2', add=TRUE)},
    ...)
}

my.filled.contour <- function(x, y, z, levels, ...) {
	filled.contour(x, y, z, color.palette=rgb.palette, 
	levels=levels,
	plot.title=title(main="", xlab='', ylab=''),
    plot.axes={axis(1); axis(2);map('world2', add=TRUE);grid()})
}

hist_stats <- function(out, alpha=0.025) {
	med = apply(out, 2, median)
	lower = apply(out, 2, quantile, alpha)
	upper = apply(out, 2, quantile, 1-alpha)
	rbind(med, lower, upper)
}

exeedence_var <- function(y, yhat, breaks=30, check.seasonality=T) {
	day_of_year = 1:365
	n = length(yhat) %/% (days_per_year * upper_year)
	res = y - yhat
	posres = res[res > 0]
	# Check either in the seasonality variable or in the long term change variable
	if (check.seasonality) {
		postime = rep(day_of_year, n*upper_year)[res > 0]
	} else {
		t = rep(c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year))), n_files)
		postime = t[res > 0]	
	}
	
	h = hist(postime, breaks=breaks, plot=F)
	ti = h$counts[-length(h$counts)]
	T = sd(ti / sum(ti))
	T
}

printf <- function(...) invisible(print(sprintf(...)))

# x appears in blue and y appears in red.
comp.hist <- function(x,y, breaks=20, ...) {
	h.x = hist(x, breaks=breaks, plot=FALSE)
	h.y = hist(y, breaks=breaks, plot=FALSE)
	dots = list(...)
	dx = min(abs(h.x$breaks[2]-h.x$breaks[1]),
	         abs(h.y$breaks[2]-h.y$breaks[1]))

	ylim = range(c(h.x$density,h.y$density))
	u = range(c(h.x$breaks,h.y$breaks))
	breaks = seq(u[1],u[2]+dx,dx)
	if (max(breaks)<max(ylim)) breaks = c(breaks,max(breaks)+dx)
		if ("plot" %in% names(dots) && dots$plot == FALSE) {
			h.x = hist(x,breaks=breaks, plot=FALSE)	
		} else if ("xlim" %in% names(dots)) {
			h.x = hist(x,breaks=breaks,density=15,ylim=ylim,col="black",freq=FALSE, ...)
		} else {
			h.x = hist(x,breaks=breaks,density=15,ylim=ylim,col="black",xlim=u,freq=FALSE, ...)
		}

	# abline(v=median(x), col="blue")
	h.y = hist(y,breaks=breaks,density=15,col="red",add=TRUE,freq=FALSE,angle=-45)
	list(h.x=h.x, h.y=h.y)
}


print.memory.usage <- function() {
	size = 0
	for (o in objects(envir=.GlobalEnv)) {
		size = size + object.size(get(o))
		message(o); print(object.size(get(o)), units='auto')
	}
	print(paste("Total memory usage:", size / 1000 / 1000, "Mb"))
}

plot.location.scale <- function(data) {
	par(mfrow=c(2,1))
	plot(data$mu)
	plot(data$scale)
}


fit.2d.location.scale <- function(t, x, y, q1, q2, days_per_year=365, n_years=250) {	
	n_days = (days_per_year * n_years)
	sample = seq(1, n_days, 1)
	rq.fit = rq(y ~ x*t, tau=c(0.25, 0.75, q1, q2), method="pfn")

	mu1 = matrix(rq.fit$fitted.values[sample,3], ncol=n_years)
	mu2 = matrix(rq.fit$fitted.values[sample,4], ncol=n_years)
	coef1 = rq.fit$coef[,3]
	coef2 = rq.fit$coef[,4]

	scale = mu2 - mu1
	list(mu=mu1, scale=scale, 
		coef_low=rq.fit$coef[,1], coef_high=rq.fit$coef[,2],
		coef1=coef1, coef2=coef2)
}

fit.location.scale <- function(x, y, q1, q2, days_per_year=365) {	
	
	rq.fit1 = rq('y ~ .', tau=q1, method="pfn", data=data.frame(x=x, y=y))

	rq.fit2 = rq('y ~ .', tau=q2, method="pfn", data=data.frame(x=x, y=y))
	mu1 = rq.fit1$fitted.values[1:days_per_year]
	mu2 = rq.fit2$fitted.values[1:days_per_year]
	scale = mu2 - mu1
	data.frame(mu=mu1, scale=scale)
}

get.ncdf.timefirst <- function(filename, lat_start, lat_count, lon_start, lon_count, t_start, t_count, varname='TREFHT') {
	ncdata <- nc_open(filename)	
		data = ncvar_get(ncdata, 
		varname, 
		start=c(t_start, lon_start, lat_start), 
		count=c(t_count, lon_count, lat_count))
	nc_close(ncdata)
	data
}
get.ncdf.light <- function(filename, lat_start, lat_count, lon_start, lon_count, t_start, t_count) {
	ncdata <- nc_open(filename)	
		data = ncvar_get(ncdata, 
		'TREFHT', 
		start=c(lon_start, lat_start, t_start), 
		count=c(lon_count, lat_count, t_count))
	nc_close(ncdata)
	data
}

get.ncdf <- function(filename, lat_min, lat_max, lon_min, lon_max) {
	ncdata <- nc_open(filename)	
	lons = ncvar_get(ncdata, 'lon') - 180
	lats = ncvar_get(ncdata, 'lat')
	lon_idx = which((lons < lon_max) & (lons > lon_min))
	lat_idx = which((lats < lat_max) & (lats > lat_min))
	data = ncvar_get(ncdata, 
		'TREFHT', 
		start=c(lon_idx[1], lat_idx[1], 1), 
		count=c(length(lon_idx), length(lat_idx), -1))
	nc_close(ncdata)
	data
}


# get.color.info <- function(z) {
#   q = 1.5
#   maxabs =  max(abs(z))
#   if (range(z)[1] > 0) { 
#     col = wr.palette(28)
#     # breaks=seq(min(z), maxabs, length=29)
#     breaks = make_zlims(z, q)
#     col.palette = wr.palette
#   } else {
#     col = rwb.palette(28)
#     # breaks=seq(-maxabs, maxabs, length=29)
#     breaks = make_zlims(z, q)
#     col.palette = rwb.palette
#   } 
#   # print(breaks)
#   list(col=col, breaks=breaks, col.palette=col.palette)
# }

plot_ddq <- function(latlon, i_lat, i_lon, i_q1, i_q2, i_day, year, i_data, main="") {
  n_stds = 4 # number of standard deviations from initial distribution to be significant 
  
  lats = unique(latlon[,1])
  lons = unique(latlon[,2])
  n_lats = length(lats)
  n_lons = length(lons)
  n_pixels = n_lats * n_lons
  out_ddq = getddq(i_q1, i_q2, i_day, year, i_data, n_pixels=n_pixels)
  zlims = make_zlims(out_ddq[[1]], 2)
  point_over_sd = out_ddq[[1]] / out_ddq[[4]]
  sig_pts = abs(point_over_sd) > n_stds
  col.info = get.color.info(out_ddq[[1]])
  point_over_sd = out_ddq[[1]] / out_ddq[[4]]
  par(new=FALSE, plt = c(0.05,0.82,0.10,0.99), oma=c(1,1,4,1), xpd=FALSE)   
  image(lons, lats, matrix(out_ddq[[1]], ncol=n_lats), 
      breaks=col.info$breaks, col=col.info$col, xaxt='n', yaxt='n', xlab="", ylab="")
  par(xpd=FALSE)
  map('world2', add=T)
  par(xpd=NA)
  points(as.numeric(i_lon), as.numeric(i_lat), col="black", cex=3, lwd=3, pch=2, bg=2)
  text(i_lon, i_lat, labels=letters[1:length(i_lat)], pos=4, cex=2)
  

  points(latlon[!sig_pts,2], latlon[!sig_pts,1], pch=3, cex=2, lwd=1, 
    col=alpha('gray', .95))
  par(new = "TRUE", plt = c(0.83,0.85,0.10,0.99), las=1)
  plot.window(xlim = c(0, 1), ylim = range(col.info$breaks), xaxs = "i", 
       yaxs = "i")
  rect(0, col.info$breaks[-length(col.info$breaks)], 1, col.info$breaks[-1L], col = col.info$col)
  axis(side=4, hadj=.4)
}

get_skew <- function(counts, values) {
  sample_mean = 1/sum(counts)*sum(counts * values)
  sample_std = sqrt(1/sum(counts)*sum(counts * (values - sample_mean)^2))
  sample_skew = 1/sum(counts)*sum(counts * ((values - sample_mean)/sample_std)^3)
  c(sample_mean, sample_std, sample_skew)
}
hist_moments <- function(mat) {
    out = t(apply(mat, 1, function(counts) {
        mi = min(which(counts != 0))
        ma = max(which(counts != 0))
        avg_vals_sub = avg_vals[mi:ma]  
        get_skew(counts[mi:ma], avg_vals_sub)
      }))
}
hist_stats <- function(out, alpha=0.025) {
  med = apply(out, 2, median)
  lower = apply(out, 2, quantile, alpha)
  upper = apply(out, 2, quantile, 1-alpha)
  rbind(med, lower, upper)
}

moments <- function(day) {
  time_idx = floor(as.numeric(day) / 10) + 1
  avg_vals = (breaks[2:length(breaks)] + breaks[1:(length(breaks)-1)]) / 2
  alpha = 0.025

  moments.i = lapply(1:n_pixels, function(pos_idx) {
      output = hist_moments(hist.data[[pos_idx]]$counts.i[[time_idx]])
      hist_stats(output)
    })
  moments.f = lapply(1:n_pixels, function(pos_idx) {
      output = hist_moments(hist.data[[pos_idx]]$counts.f[[time_idx]])
      hist_stats(output)
    })
  moments.d = lapply(1:n_pixels, function(pos_idx) {
      output.i = hist_moments(hist.data[[pos_idx]]$counts.i[[time_idx]])
      output.f = hist_moments(hist.data[[pos_idx]]$counts.f[[time_idx]])
      hist_stats(output.f - output.i)
    })

    # yval.i = hist.data[[pos_idx]]$counts.i[[time_idx]]
    # yval.f = hist.data[[pos_idx]]$counts.f[[time_idx]]
    
    # moments_f = get_skew(yval.f[mi:ma], avg_vals_sub)
  list(moments.i, moments.f, moments.d)
}

before_after_hist <- function(day, pos_idx, title="") {
  ymax = 1300
  
  my_colors = c(alpha('green', .8), alpha('red', .8), alpha('blue', .8))
  time_idx = floor(as.numeric(day) / 10) + 1
  data = hist.data[[pos_idx]][[time_idx]]
  density.i = data$hist_stats_i[1,] / sum(data$hist_stats_i[1,])
  density.f = data$hist_stats_f[1,] / sum(data$hist_stats_f[1,])
  n.i = length(density.i)
  n.f = length(density.f)
  nbreaks_i = length(data$breaks_i)
  nbreaks_f = length(data$breaks_f)
  dx = data$breaks_i[2] - data$breaks_i[1]
  # lwd = 100 / diff(range(c(data$breaks_i, data$breaks_f)))
  lwd = 8
  # plot(data$breaks_i[1:n.i], density.i, type='h', lwd=lwd, lend=1, 
 #        col=my_colors[1], xlim=range(c(data$breaks_i, data$breaks_f)), 
 #        ylim = range(c(density.i, density.f)),
 #        ylab="", xlab='Temperature', yaxt='n', cex.axis=1.2)
    histlike(data$breaks_i, density.i, col=my_colors[1], 
      xlim=range(c(data$breaks_i, data$breaks_f)), 
        ylim=range(c(density.i, density.f)),
        ylab="", xlab='Temperature', yaxt='n', cex.axis=1.2)
  histlike(data$breaks_f, density.f, col=my_colors[2], 
        add=TRUE)
    # lines(data$breaks_f[1:n.f], density.f, type='h', lwd=lwd, lend=1, 
     #  col=my_colors[2])
    legend('topright', 
    c('Initial', 'Final'), col=my_colors[1:2], lty=1, lwd=2, cex=1.5)
}

histlike <- function(breaks, values, col=1, xlim=range(breaks), 
  ylim=c(0, max(values)), add=FALSE, ...) {
  if (!add) {
    plot(NA, ylim=ylim, xlim=xlim, ...)
  } else {

  }
  nb = length(breaks)
  nv = length(values)
  segments(x0=breaks[1:(nb-1)], x1=breaks[2:nb], y0=values, y1=values, col=col)
  segments(x0=breaks[2:(nb-1)], x1=breaks[2:(nb-1)], y0=values[1:(nv-1)], y1=values[2:nv], col=col)
  segments(x0=breaks[c(1,nb)], x1=breaks[c(1,nb)], y0=c(0,0), y1=values[c(1,nv)], col=col)
  abline(h=0)
}

printf <- function(...) invisible(print(sprintf(...)))

#plot incremental mean increase as a function of incremental standard deviation increase
rise_over_runs <- function(q1, q2, day, pos_idx, data_type="avg") {
  q1_idx = c(outer(c(1:np), c(nq*np*(1:B))), "+") + (which(q1 == q) - 1) * np + nq*np*B*(pos_idx - 1)
  q2_idx = c(outer(c(1:np), c(nq*np*(1:B))), "+") + (which(q2 == q) - 1) * np + nq*np*B*(pos_idx - 1)

  if (data_type == "avg") {
    b_q1 = a[q1_idx,] 
    b_q2 = a[q2_idx,] 
  } else if (data_type == "min") {
    b_q1 = a_mn[q1_idx,]  
    b_q2 = a_mn[q2_idx,]  
  } else if (data_type == "max") {
    b_q1 = a_mx[q1_idx,]  
    b_q2 = a_mx[q2_idx,]  
  } else {
    print("data_type not selected")
  }
  betaq1 = matrix(y_q1, np, B)
  betaq2 = matrix(y_q2, np, B)
  beta_dq = betaq2 - betaq1
  const_d_idx = seq(1, n_years*days_per_year, days_per_year) + (n_years - 1)*day
  X_const_d = X[const_d_idx, ]

  y_q1 = X_const_d %*% betaq1
  dy = X_const_d %*% beta_dq

  list(y_q1=y_q1, dy=dy)
}

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='', ...) {
    scale = (length(lut)-1)/(max-min)

    # dev.new(width=1.55, height=4)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title, ...)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}



getddq <- function(q1, q2, day, year, data_type, n_pixels, normalize=TRUE, relnorm=FALSE) {
  first_year_idx = which(x_sub == day) 
  last_year_idx = which(x_sub == day) + (which(y_sub==year) - 1) * length(x_sub)

  q1_idx = c(outer(c(1:np), c(nq*np * (1:(B*n_pixels) - 1)), "+")) + (which(q1 == q) - 1) * np
  q2_idx = c(outer(c(1:np), c(nq*np * (1:(B*n_pixels) - 1)), "+")) + (which(q2 == q) - 1) * np
  low_idx = c(outer(c(1:np), c(nq*np * (1:(B*n_pixels) - 1)), "+")) + (which(0.25 == q) - 1) * np
  high_idx = c(outer(c(1:np), c(nq*np * (1:(B*n_pixels) - 1)), "+")) + (which(0.75 == q) - 1) * np
  
  betaq1 = array(getData(data_type, q1_idx), dim=c(np, B, n_pixels))
  betaq2 = array(getData(data_type, q2_idx), dim=c(np, B, n_pixels))
  beta_low = array(getData(data_type, low_idx), dim=c(np, B, n_pixels))
  beta_high = array(getData(data_type, high_idx), dim=c(np, B, n_pixels))
  
  beta_dq = betaq2 - betaq1
  beta_iqr = beta_high - beta_low
  
  # The following extracts two 33-length vectors
  # we will multiply these by the beta matrix representing all pixels
  x_first = Xsubsample[first_year_idx, ]
  x_last = Xsubsample[last_year_idx, ]  
  first_y = apply(beta_dq, c(2,3), function(b) x_first %*% b)
  first_y_iqr = apply(beta_iqr, c(2,3), function(b) x_first %*% b)
  last_y = apply(beta_dq, c(2,3), function(b) x_last %*% b)
  last_y_iqr = apply(beta_iqr, c(2,3), function(b) x_last %*% b)
  if (relnorm) {
      d_y = (last_y  - first_y) / first_y
  } else {
    if (!normalize || (q1 == 0.25 && q2 == 0.75)) {
      d_y = last_y - first_y
    } else {
      d_y = last_y / last_y_iqr - first_y / first_y_iqr
    }
  }
  y_i = apply(first_y, 2, quantile, 0.5)
  y_f = apply(last_y, 2, quantile, 0.5)
  d_ypoint = apply(d_y, 2, quantile, 0.5, na.rm=TRUE)
  d_yhigh = apply(d_y, 2, quantile, 0.95, na.rm=TRUE)
  d_ylow = apply(d_y, 2, quantile, 0.05, na.rm=TRUE)
  sd_dy = apply(d_y, 2, sd, na.rm=TRUE)

  list(d_ypoint=d_ypoint, d_yhigh=d_yhigh, d_ylow=d_ylow, sd_dy=sd_dy, 
    y_i=y_i, y_f=y_f)
}

getcrossings <- function(quantiles) {
  apply(quantiles, 1, function(x) which(diff(x)<0))
}

getq <- function(pos_idx, data_type) {
  start_idx = B * np * nq * (pos_idx - 1) + 1
  end_idx = start_idx + B * np * nq
  # my_pixel = updateData(data_type, start_idx - 1, end_idx)
  my_pixel = getData(data_type, start_idx:end_idx)
  
  my_beta = array(my_pixel, dim=c(np, nq, B))
  quantiles = array(apply(my_beta, 3, function(bi) {
    bi[1, ] = bi[1, ] - KtoC
    Xsubsample %*% bi
  }), dim=c(dim(Xsubsample)[1], nq, B))
  # if (div_by_iqr) {
  #   biqr1 = my_beta[,0.25==q,,drop=T]
  #   biqr1 = biqr1[,!apply(biqr1, 2, function(a) all(a==0))]
  #   biqr2 = my_beta[,0.75==q,,drop=T]
  #   biqr2 = biqr2[,!apply(biqr2, 2, function(a) all(a==0))]
  #   qiqr =  Xsubsample %*% (biqr2 - biqr1)
  #   quantiles = quantiles / qiqr
  # }

  yhigh = apply(quantiles, c(1,2), quantile, 0.95)
  ylow = apply(quantiles, c(1,2), quantile, 0.05)
  ypoint = apply(quantiles, c(1,2), quantile, 0.5)
  
  list(point=ypoint, high=yhigh, low=ylow)
}

nearest_idx <- function(v, v0) {
  which(abs(v-v0) == min(abs(v-v0)))
}


getdq <- function(pos_idx, q1, q2, data_type, div_by_iqr=TRUE, div_by_initial=FALSE) {
  start_idx = B * np * nq * (pos_idx - 1) + 1
  end_idx = start_idx + B * np * nq
  # my_pixel = updateData(data_type, start_idx - 1, end_idx)
  my_pixel = getData(data_type, start_idx:end_idx)
  
  my_beta = array(my_pixel, dim=c(np, nq, B))
  b1 = my_beta[,q1==q,,drop=T]
  b1 = b1[,!apply(b1, 2, function(a) all(a==0))]
  b2 = my_beta[,q2==q,,drop=T]
  b2 = b2[,!apply(b2, 2, function(a) all(a==0))]
  db = b2 - b1
  dq = Xsubsample %*% db

  if (div_by_iqr) {
    biqr1 = my_beta[,0.25==q,,drop=T]
    biqr1 = biqr1[,!apply(biqr1, 2, function(a) all(a==0))]
    biqr2 = my_beta[,0.75==q,,drop=T]
    biqr2 = biqr2[,!apply(biqr2, 2, function(a) all(a==0))]
    qiqr =  Xsubsample %*% (biqr2 - biqr1)
    dq = dq / qiqr
  }
  if (div_by_initial) {
    dqi = dq[rep(1:(365/d_days_per_year), length(y_sub)), ]
    dq = dq / dqi
  }

  yhigh = apply(dq, 1, quantile, 0.95)
  ylow = apply(dq, 1, quantile, 0.05)
  ypoint = apply(dq, 1, quantile, 0.5)
  
  list(point=ypoint, high=yhigh, low=ylow)
}

# Gets the differences between initial and final time of any quantile
getdt <- function(q1, day, data_type) {
  n_pixels = get.n_pixels(data_type)
  year=241
  first_year_idx = which(x_sub == day) 
  last_year_idx = which(x_sub == day) + (which(y_sub==year) - 1) * length(x_sub)
  q1_idx = c(outer(c(1:np), c(nq*np * (1:(B*n_pixels) - 1)), "+")) + (which(q1 == q) - 1) * np
  
  y_q1 = getData(data_type, q1_idx)
  bi = array(y_q1, dim=c(np, B, n_pixels))
  
  # The following extracts two 33-length vectors
  # we will multiply these by the beta matrix representing all pixels
  x_first = Xsubsample[first_year_idx, ]
  x_last = Xsubsample[last_year_idx, ]  
  first_y = apply(bi, c(2,3), function(b) x_first %*% b)
  last_y = apply(bi, c(2,3), function(b) x_last %*% b)
  
  d_y = last_y - first_y
  y_i = apply(first_y, 2, function(sample) {
        quantile(sample[sample!=0], 0.5)})
  y_f = apply(first_y, 2, function(sample) {
        quantile(sample[sample!=0], 0.5)})
  d_ypoint = d_y[1,]
  d_yhigh = apply(d_y, 2, function(sample) {
        quantile(sample[sample!=0], 0.95)})
  d_ylow = apply(d_y, 2, function(sample) {
        quantile(sample[sample!=0], 0.05)})
  sd_dy = apply(d_y, 2, function(sample) {
        sd(sample[sample!=0])})

  list(d_ypoint=d_ypoint, d_yhigh=d_yhigh, d_ylow=d_ylow, sd_dy=sd_dy, 
    y_i=y_i, y_f=y_f)
}

make_zlims <- function(v, q=1, symmetric=TRUE, include.zero=FALSE) {
  al = 0.95
  ah = 1.05
  if (any(v < 0)) {
    if (symmetric) {
      maxabs = max(abs(v))*ah
      res = c(-(abs(seq(-maxabs^(1/q), 0, length=15)[-15])^q), 
        0, (seq(0, maxabs^(1/q), length=15)[-1])^q)  
    } else {
      res = c(-(abs(seq(-abs(min(v))^(1/q), 0, length=15)[-15])^q), 
        0, (seq(0, max(v)^(1/q), length=15)[-1])^q)  
    }
    
  } else {
    if (!include.zero) {
      res = seq(min(v)*al, (max(v)*ah)^(1/q), length=29)^(q)
    } else {
      res = seq(0, (max(v)*ah)^(1/q), length=29)^(q)
    }
  }
  res
}


getQuantile <- function(pos_idx, qi, data_type) {
  start_idx = B * np * nq * (pos_idx - 1) + 1
  end_idx = start_idx + B * np * nq
  # my_pixel = a[start_idx:end_idx,]
  # my_pixel = updateData(data_type, start_idx - 1, end_idx)
  my_pixel = getData(data_type, start_idx:end_idx)
  my_beta = array(my_pixel, dim=c(np, nq, B))

  bi = my_beta[,qi==q,,drop=T]
  
  if (any(bi==0)) {
    print("there were zero coefficients. These runs will be ignored")
    idx = apply(bi, 2, function(a) {
      all(a==0)
    })
    print(paste("Removing", sum(idx), "run(s)"))
    # print(idx)
    bi = bi[,!idx]
  }

  bi[1, ] = bi[1, ] - KtoC

  yhigh = t(apply(Xsubsample %*% bi, 1, quantile, 0.95))
  ymed = t(apply(Xsubsample %*% bi, 1, quantile, 0.5))
  ylow = t(apply(Xsubsample %*% bi, 1, quantile, 0.05))
  ypoint = Xsubsample %*% bi[,1]
  

  list(point=ymed, high=yhigh, low=ylow)
}

updateData <- function(type, start_idx, end_idx) {
  nrows = end_idx - start_idx
  if (type=="avg") {
    a <- read.table("../data/outBoot100_nb", skip=start_idx, nrows=nrows)
  } else if (type=="min") {
    a <- read.table("../data/outBoot100_mn", skip=start_idx, nrows=nrows)
  } else if (type=="max") {
    a <- read.table("../data/outBoot100_mx", skip=start_idx, nrows=nrows)
  }
  a[,1]
}

appendList <- function (x, val) 
{
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    for (v in names(val)) {
        x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
            appendList(x[[v]], val[[v]])
        else c(x[[v]], val[[v]])
    }
    x
}

plot_q_or_dq <- function(data, ...) {
  # qidx = which(q == qi)
  # y.low = matrix(data$low[,qidx], ncol = n_years / d_years)
  # y.high = matrix(data$high[,qidx], ncol = n_years / d_years)
  # y.point = matrix(data$point[,qidx], ncol=n_years / d_years)
  y.low = matrix(data$low, ncol = length(y_sub))
  y.high = matrix(data$high, ncol = length(y_sub))
  y.point = matrix(data$point, ncol= length(y_sub))

  dots = list(...)
  if (!('ylim' %in% names(dots))) {
    dots$ylim = range(c(c(y.low), c(y.high)))
  }

  slices = seq(0, n_years, 40) / d_years + 1
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
    if (col == 1) {
      plotargs = appendList(list(x=x_sub, y=z.point, col=col.map[col], xlab="Day of the year", 
              ylab="Degrees (C)",
              type="l"), dots)  
      do.call('plot', plotargs)
    } else {
      lines(x_sub, z.point, col=col.map[col], lty=1)
    }
    lines(x_sub, z.low, col=col.map[col], lty=2)
    lines(x_sub, z.high, col=col.map[col], lty=2)
    col = col + 1
  }
}

getData <- function(data_type, idx) {
  if (data_type == "avg") {
    my_pixel = a[idx,]  
  } else if (data_type == "big") {
    my_pixel = a_big[idx,]  
  } else if (data_type == "min") {
    my_pixel = a_mn[idx,]
  } else if (data_type == "max") {
    my_pixel = a_mx[idx,]
  } else {
    print("data_type not selected")
  }
  my_pixel
}

getPos_idx <- function(lat, lon, lats=la, lons=lo) {
  n_lons = length(lo)
  lat_idx = nearest_idx(lats, lat)
  lon_idx = nearest_idx(lons, lon)
  pos_idx = (lat_idx - 1) * n_lons + lon_idx   
}


