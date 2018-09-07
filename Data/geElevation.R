tab= d.calib
known.coord = c(lon.calib, lat.calib)
plot=TRUE
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  tab <- geolight.convert(tab[,1], tab[,2], tab[,3])  
  
  sun  <- solar(tab[,1])
  z    <- refracted(zenith(sun, known.coord[1], known.coord[2]))
  plot(z)
  
  inc = 0
  repeat {
  twl_t   <- twilight(tab[,1], known.coord[1], known.coord[2], rise = tab[,2], zenith = max(z)+inc)
  twl_dev <- ifelse(tab$Rise, as.numeric(difftime(tab[,1], twl_t, units = "mins")),
                    as.numeric(difftime(twl_t, tab[,1], units = "mins")))
  if(all(twl_dev>=0)) {
    break
  } else {
    inc <- inc+0.01
  }
  }
  z0 <- max(z)+inc
  
  seq <- seq(0, max(twl_dev), length = 100)

  fitml_ng <- suppressWarnings(fitdistr(twl_dev, "gamma"))
  lns      <- dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2])
  
  diffz <- as.data.frame(cbind(min = apply(cbind(tab[,1], twilight(tab[,1], known.coord[1], known.coord[2], rise = tab[,2], zenith = z0)), 1, function(x) abs(x[1]-x[2]))/60, z = z))
  mod  <- lm(z~min, data = diffz)
  mod2 <- lm(min~z, data = diffz)
  
  a1.0 <- seq[which.max(lns)]
  a1.1 <- 90-predict(mod, newdata = data.frame(min = a1.0))
  
  if(plot) {
  opar <- par(mar = c(10, 4, 1, 1))
  hist(twl_dev, freq = F, breaks = 26, main = "Twilight Model", xlab = "twilight error (min)")
  lines(seq, lns, col = "firebrick", lwd = 3, lty = 2)
  
  points(predict(mod2, newdata=data.frame(z = z0)), 0, pch = 21, cex = 5, bg = "white", lwd = 2)
  text(predict(mod2, newdata=data.frame(z = z0)), 0, "0")

  points(a1.0, 0, pch = 21, cex = 5, bg = "white", lwd = 2)
  text(a1.0, 0, "1")
  
  axis(1, at = seq(0, max(twl_dev), 6), labels = round(90-predict(mod, newdata = data.frame(min = seq(0, max(twl_dev), 6))),1), line = 5)
  mtext("sun elevation angle (degrees)", 1, line = 8)
  
  legend("topright", paste(c("0. Sun elevation angle (zero)", "1. Sun elevation angle (median)", "log-mean", "log-sd"), round(c(90-z0, a1.1, fitml_ng$estimate[1], fitml_ng$estimate[2]),3)), bty = "n")
  
  par(opar)
  }

  data.frame(a1 = a1.1, e0 = 90-z0, log.mean =  fitml_ng$estimate[1], log.sd =  fitml_ng$estimate[2])
}
