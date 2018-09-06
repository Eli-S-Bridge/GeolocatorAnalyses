site <- cL$site
degElevation <- gE[1]
distThreshold = 600 
ll.radius.degrees = 30
fixed = NULL
alpha = c(0, 15)
plot = TRUE



  gl.loglik <- function(crds, Twilight, Rise, degElevation, parms, twilight = NULL) {
    t.tw <- twilight(Twilight, lon = crds[1], lat = crds[2], 
                     rise = Rise, zenith = 94, 
                     iters = 6)
    
    diff.sr <- as.numeric(difftime(Twilight[Rise], t.tw[Rise], units = "mins"))
    diff.ss <- as.numeric(difftime(t.tw[!Rise], Twilight[!Rise], units = "mins"))
    
    if(is.null(twilight)) {
      ll <- sum(dlnorm(c(diff.sr, diff.ss), parms[1], parms[2], log = F), na.rm = T)
    } else {
      if(twilight=="sr"){
        ll <- sum(dlnorm(diff.sr, parms[1], parms[2], log = F), na.rm = T)
      }
      if(twilight=="ss") {
        ll <- sum(dlnorm(diff.ss, parms[1], parms[2], log = F), na.rm = T)
      }}
    
    return(ifelse(!is.infinite(abs(ll)), ll, -10000))
  }

  
  
mergeSites2 <- function(tFirst, tSecond, type, twl, site, degElevation,
                        distThreshold = 250, ll.radius.degrees = 30, fixed = NULL, alpha = c(0, 15), plot = TRUE) {
  
  mycl <- parallel::makeCluster(parallel::detectCores()-1)
  tmp  <- parallel::clusterSetRNGStream(mycl)

  tab <- twl.gl
  
  if(is.null(fixed)) fixed <- matrix(FALSE, ncol = 2, nrow = nrow(tab))
  fixed.ind <- apply(fixed, 1, function(x) any(x))
  
  site0 <- site
  tw <- data.frame(datetime = .POSIXct(c(tab$tFirst, tab$tSecond), "GMT"), 
                   type     = c(tab$type, ifelse(tab$type == 1, 2, 1)), site = site, fixed = fixed.ind)
  tw <- tw[!duplicated(tw$datetime), ]
  tw <- tw[order(tw[, 1]), ]
  tw <- tw[1:nrow(tab),]
  tw$Rise <- ifelse(tw$type==1, TRUE, FALSE)
    site.old <- tw$site  
  
  crds0 <- coord(tab, degElevation = degElevation, note = F)
  tab$lon <- crds0[,1]
  tab$lat <- crds0[,2]
  tw$lon  <- crds0[,1]
  tw$lat  <- crds0[,2]
  
  lonlim <- range(crds0[, 1], na.rm = T)
  lon.seq <- seq(lonlim[1] - 1, lonlim[2] + 1, by = 1)
  latlim <- range(crds0[, 2], na.rm = T)
  lat.seq <- seq(latlim[1] - 1, latlim[2] + 1, by = 1)
  
  tmp  <- parallel::clusterEvalQ(mycl, library("GeoLight"))   
  
  mod <- function(x) {
    
  lons <- seq(-180, 180, by = 2.5)
  lats <- seq(-75, 75, by = 2.5)
    
  crdsT <- expand.grid(lons, lats)
  
  ll.sr    <- parallel::parRapply(mycl, crdsT, FUN = gl.loglik, Twilight = x$datetime, Rise = ifelse(x$type==1, TRUE, FALSE),
                               degElevation = gE[1], parms = parms,  twilight = "sr")
  ll.ss    <- parallel::parRapply(mycl, crdsT, FUN = gl.loglik, Twilight = x$datetime, Rise = ifelse(x$type==1, TRUE, FALSE),
                                  degElevation = gE[1], parms = parms,  twilight = "ss")
  ll <- apply(cbind(ll.sr/max(ll.sr, na.rm = T), ll.ss/max(ll.sr, na.rm = T)), 1, function(x) ifelse(any(x<0.01), NA, sum(x)))
  ll <- ll/max(ll, na.rm = T)
  
  crdsT <- crdsT[!is.na(ll),]
  
  r0     <- raster(xmn =-180, xmx = 180, ymn = -75, ymx = 75, nrow = length(lats), ncol = length(lons))
  r      <- rasterize(crdsT[!is.na(ll),], r0, field = ll[!is.na(ll)])
  plot(r)

  r0     <- raster(xmn = min(crdsT[,1])-3, xmx = max(crdsT[,1])+3, ymn = min(crdsT[,2])-3, ymx = max(crdsT[,2])+3, res = 0.5)
  crdsT  <- coordinates(r0)
  ll.sr  <- parallel::parRapply(mycl, crdsT, FUN = gl.loglik, Twilight = x$datetime, Rise = ifelse(x$type==1, TRUE, FALSE),
                               degElevation = gE[1], parms = parms, twilight = "sr")
  ll.ss    <- parallel::parRapply(mycl, crdsT, FUN = gl.loglik, Twilight = x$datetime, Rise = ifelse(x$type==1, TRUE, FALSE),
                               degElevation = gE[1], parms = parms, twilight = "ss")
  
  ll <- apply(cbind(ll.sr/max(ll.sr, na.rm = T), ll.ss/max(ll.sr, na.rm = T)), 1, function(x) ifelse(any(x<0.001), NA, sum(x)))
  centre <- crdsT[which.max(ll),]
  
  
  r   <- rasterize(crdsT[!is.na(ll),], r0, field = ll[!is.na(ll)])
  r[] <- r[]/max(r[], na.rm = T)
  
  plot(r)
  contour(r, add = T, levels = c(1, 0.495))
  
  crdsRange <- coordinates(r)[!is.na(r[]) & r[]>0.495,]
  
  matrix(c(centre[1], centre[2], min(crdsRange[,1]), max(crdsRange[,1]), 
           min(crdsRange[,2]), max(crdsRange[,2])), ncol = 6)
    
  }
  
  tw$merge = FALSE
  xTab <- split(tw, f = tw$site)
    x0 <- xTab[[1]]
  xTab <- xTab[-1]
    
  sm <- matrix(ncol = 6, nrow = max(site.old))
  repeat{
    
    for(i in 1:(length(xTab)-1)) {
    
      if(all(!xTab[[i]]$merge)) {
      
      cat(sprintf('\n comparing site %d and %d', i, i+1)) 
      out  <- do.call("rbind", lapply(xTab[c(i, i+1)], function(x) mod(x)))
      
      sm[i,] <- out[1,]
      if(i == (length(xTab)-1)) sm[i+1,] <- out[2,]
      
      # plot(out[,1], out[,2], xlim = range(out[,c(1,3,4)]), ylim = range(out[,c(2,5,6)]))
      # plot(wrld_simpl, add = T)
      # arrows(out[,1], out[,5], out[,1], out[,6], length = 0)
      # arrows(out[,3], out[,2], out[,4], out[,2], length = 0)
      
      if(all(!xTab[[i]]$fixed)) {
      
      dist <- fields:::rdist.earth(out[,1:2])[2,1]
      
      if(dist<distThreshold & 
         (out[2,3] < out[1,1] & out[2,4] > out[1,1]) & (out[2,1] > out[1,3] & out[2,1] < out[1,4]) &
         (out[2,5] < out[1,2] & out[2,6] > out[1,2]) & (out[2,2] > out[1,5] & out[2,2] < out[1,6])) {
        
      tw$site[min(which(tw$site==i)):max(which(tw$site==(i+1)))] <- i
      tw$site[tw$site>0 & tw$site>i] <-  tw$site[tw$site>0 & tw$site>i]-1
        
      xTab <- split(tw, f = tw$site)
        x0 <- xTab[[1]]
      xTab <- xTab[-1]
      
      cat(".... merge (back to site 1)")
      
      break
      
      } else {
        tw$merge[tw$site==i] <- TRUE
        cat(".... no action")
      }   
      
      }
      }
    }
    if(i==(length(xTab)-1)) {
      xTab <- split(tw, f = tw$site)
      break
    }
  }
  
  parallel::stopCluster(mycl) 
  
  sm  <- cbind(1:sum(!is.na(sm[,1])), sm[!is.na(sm[,1]),])
  
  out <- do.call("rbind", xTab)[,c(1,5,3,4)]
    out <- out[order(out$datetime),]       

  
  if (plot) {
    site   <- out$site
    hours0 <- as.numeric(format(out[, 1], "%H")) + as.numeric(format(out[, 1], "%M"))/60
    crd0   <- sm[match(out$site, sm[,1]), 2:3]
    # crd0[is.na(crd0[,1]),] <- crds0[is.na(crd0[,1]),]
    hours1 <- twilight(out$datetime, rise = out$Rise, zenith = 90 - degElevation, 
                       lon = crd0[,1], lat = crd0[,2])
    hours1 <- as.numeric(format(hours1, "%H")) + as.numeric(format(hours1, "%M"))/60
    hours2 <- twilight(out$datetime, rise = out$Rise, zenith = 90 - degElevation, 
                       lon = crd0[,1], lat = crd0[,2])
    
    hours2 <- as.numeric(format(hours2, "%H")) + as.numeric(format(hours2,"%M"))/60
    hours3 <- twilight(out$datetime, rise = out$Rise, zenith = 90 - degElevation, 
                       lon = crd0[,1], lat = crd0[,2])
    hours3 <- as.numeric(format(hours3, "%H")) + as.numeric(format(hours3,"%M"))/60
    
    
    for (t in c(TRUE, FALSE)) {
      cor <- rep(NA, 24)
      for (i in 0:23) {
        cor[i + 1] <- max(abs((c(hours0[out$Rise==t][1], 
                                 hours0[out$Rise == t]) + i)%%24 - 
                                (c(hours0[out$Rise == t], 
                                   hours0[out$Rise == t][length(hours0)]) +i)%%24), 
                          na.rm = T)
      }
      hours0[out$Rise == t] <- (hours0[out$Rise == t] + (which.min(round(cor,2))) - 1)%%24
      hours1[out$Rise == t] <- (hours1[out$Rise == t] + (which.min(round(cor,2))) - 1)%%24
      hours2[out$Rise == t] <- (hours2[out$Rise == t] + (which.min(round(cor,2))) - 1)%%24
      hours3[out$Rise == t] <- (hours3[out$Rise == t] + (which.min(round(cor,2))) - 1)%%24
    }
    
    opar <- par(mfrow = c(5, 1), oma = c(5, 0, 0, 0), mar = c(1.5,5, 1, 1))
    mig1 <- site.old
    mig1[mig1 > 0] <- 1
    mig2 <- out$site
    mig2[mig2 > 0] <- 1
    plot(out[, 1], ifelse(mig2 > 0, 1, 0), type = "l", yaxt = "n", 
         ylab = NA, ylim = c(0, 1.5), col = "firebrick", lwd = 2, 
         xaxt = "n")
    lines(out[, 1], ifelse(mig1 > 0, 1, 0), type = "l", lty = 2)
    rect(out[site > 0 & !duplicated(site), 1], 1.1, 
         out[site > 0 & !duplicated(site, fromLast = T), 1], 1.4, col = "grey90", lwd = 0)
    rect(out[site > 0 & !duplicated(site), 1], 1.1, 
         out[site > 0 & !duplicated(site, fromLast = T), 1], 1.4, 
         col = ifelse(apply(fixed[site>0,], 1, function(x) any(x))[!duplicated(site[site>0])], "red", "transparent"), 
         density = 60)
    axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
         labels = FALSE)
    
    plot(out[out$Rise, 1], hours1[out$Rise], type = "l", 
         lwd = 2, col = "firebrick", ylab = "Sunrise (red)", 
         xlim = range(out[, 1]), ylim = range(hours0[out[, 2] == 1]), xaxt = "n")
    lines(out[out$Rise, 1], hours2[out$Rise], type = "l", 
          lwd = 1, lty = 2)
    lines(out[out[, 2] == 1, 1], hours3[out[, 2] == 1], type = "l", 
          lwd = 1, lty = 2)
    points(out[out[, 2] == 1 & !fixed.ind, 1], hours0[out[, 2] == 1 & !fixed.ind], cex = 0.5, 
           pch = 21, col = "black", bg = "firebrick", lwd = 0.5)
    axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
         labels = FALSE)
    plot(out[!out$Rise, 1], hours1[!out$Rise], type = "l", 
         lwd = 2, col = "cornflowerblue", ylab = "Sunset (blue)", 
         xlim = range(out[, 1]), ylim = range(hours0[!out$Rise]), xaxt = "n")
    lines(out[!out$Rise, 1], hours2[!out$Rise], type = "l", 
          lwd = 1, lty = 2)
    lines(out[!out$Rise, 1], hours3[!out$Rise], type = "l", 
          lwd = 1, lty = 2)
    points(out[!out$Rise & !fixed.ind, 1], hours0[!out$Rise & !fixed.ind], cex = 0.5, 
           pch = 21, col = "black", bg = "cornflowerblue", lwd = 0.5)
    axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
         labels = FALSE)
    
    plot(out[, 1], ifelse(fixed.ind, NA, crds0[, 1]), type = "o", pch = 16, cex = 0.5, 
         xaxt = "n", ylab = "Longitude", cex.lab = 1.7, xlab = "")
    abline(v = c(out[site.old > 0 & !duplicated(site.old), 1], 
                 out[site.old > 0 & !duplicated(site.old, fromLast = T), 1]), lty = 2)
    abline(v = c(out[out$site  > 0 & !duplicated(out$site), 1], 
                 out[out$site  > 0 & !duplicated(out$site, fromLast = T), 1]), lwd = 1.5, col = "firebrick")
    axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
                 labels = FALSE)
    plot(out[, 1], ifelse(fixed.ind, NA, crds0[, 2]), type = "o", pch = 16, cex = 0.5, 
         xaxt = "n", ylab = "Latitude", cex.lab = 1.7, xlab = "")
    abline(v = c(out[site.old > 0 & !duplicated(site.old), 1], 
                 out[site.old > 0 & !duplicated(site.old, fromLast = T), 1]), lty = 2)
    abline(v = c(out[out$site > 0 & !duplicated(out$site), 1], 
                 out[out$site > 0 & !duplicated(out$site, fromLast = T), 1]), lwd = 1.5, 
           col = "firebrick")
    axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
         labels = format(seq(out[1, 1], out[nrow(out), 1], length = 10), "%d-%b"))
    mtext("Date", 1, outer = T, line = 1.6, cex = 1.2)
    par(opar)
  }
  sm        <- as.data.frame(sm)
  names(sm) <- c("site", "Lon", "Lat", "Lon.upper", "Lat.upper", 
                  "Lon.lower", "Lat.lower")
  
  diff <- c(apply(cbind(out$datetime[-nrow(out)], 
                        out$datetime[-1]), 1, function(x) c(x[2] - x[1])/60/60), 0)
  
  out.gl <- data.frame(tFirst = as.POSIXct("1900-01-01 01:01", "GMT"), tSecond = as.POSIXct("1900-01-01 01:01", "GMT"),  type = 0, site = 0, fixed = 0, diff.max = 0)
  rw <- 1
  for (k in 1:(nrow(out) - 1)) {
    if (as.numeric(difftime(out$datetime[k], out$datetime[k + 1])) < 24 & out$datetime[k] != out$datetime[k + 1]) {
      out.gl[rw, 1] <- out$datetime[k]
      out.gl[rw, 2] <- out$datetime[k + 1]
      out.gl[rw, 3] <- ifelse(out$Rise[k], 1, 2)
      out.gl[rw, 4] <- out$site[k]
      out.gl[rw, 5] <- out$fixed[k]
      out.gl[rw, 6] <- max(diff[k:(k + 1)])
      rw <- rw + 1
    }
  }
  out.gl <- subset(out.gl, diff.max < 23)
  
  list(twl = out.gl[,c(1,2)], site = out.gl$site, summary = sm)
}
