
library(GeoLocTools)
setupGeolocation()

Species <- "PuMa"
ID = "PuMa01"
wd <- "Data"

lat.calib <-  33.9
lon.calib <- -96.8


raw <- readLig(paste0(wd, "/RawData/", Species, "/", ID, ".lig"))
# raw <- read.csv(paste0(wd, "/RawData/", Species, "/", ID, ".csv"))
#   names(raw)[1:2] <- c("Date", "Light")
#   raw$Date <- as.POSIXct(raw$Date, format = "%Y-%m-%dT%T", tz = "GMT")
raw$Light <- log(raw$Light)

with(raw[2000:5000,], plot(Date, Light, type = "o", pch = 16))  
abline(h = 1.75)

offset <- 15
twl <- preprocessLight(raw, threshold = 1.75, offset = offset, lmax =  5, gr.Device = "x11")
twl <- subset(twl, !is.na(Twilight))

# write.csv(twl, paste0(wd, "/Results/", Species, "/", ID, "_twl.csv"), row.names = FALSE)
twl <- read.csv(paste0(wd, "/Results/", Species, "/", ID, "_twl.csv"))
  twl$Twilight <- as.POSIXct(twl$Twilight, tz = "GMT")

  
  # twl0 <- read.csv("C:/Users/SLi/Dropbox/Science/Projects/RedKnot_tracking/FLightR/examples/tree_swallow_BAS_tag_example/749.csv")
  # raw  <- data.frame(Date = as.POSIXct(twl0$datetime, format = "%Y-%m-%dT%T", tz = "GMT"), Light = twl0$light)
  
  # twl  <- data.frame(Twilight = as.POSIXct(twl0$datetime, format = "%Y-%m-%dT%T", tz = "GMT"), Rise = twl0$twilight)[twl0$twilight!=0,]
  # twl$Rise  <- ifelse(twl$Rise==1, TRUE, FALSE)
    
  
  
twl <- twl[!twl$Deleted,]
# twl <- twilightAdjust(twl)

### Calibration
tm <- seq(min(raw$Date), max(raw$Date), by = "day")
rise <- rep(c(TRUE, FALSE), length(tm))

c.dat <- data.frame(Twilight = twilight(rep(tm, each = 2), lon = lon.calib, lat = lat.calib, 
                                        rise = rise, zenith = 96), Rise = rise)


lightImage(tagdata = raw, offset = offset, zlim = c(0, 6))
tsimagePoints(c.dat$Twilight, offset = 15, pch = 16, cex = 0.25,
              col = "red")

tm.calib <- as.POSIXct(c(min(raw$Date)+5*24*60*60, min(raw$Date)+50*24*60*60))
abline(v = tm.calib, col = "red")


###
twl.gl  <- export2GeoLight(twl) 

d.calib <- subset(twl.gl, tFirst>=tm.calib[1] & tSecond<=tm.calib[2])

gE      <- getElevation(twl = d.calib, known.coord = c(lon.calib, lat.calib), method = "gamma")

crds <- GeoLight::coord(twl.gl, degElevation = gE[1])
tripMap(crds, xlim = range(crds[,1], na.rm = T), ylim = range(crds[,2], na.rm = T))
points(lon.calib, lat.calib, pch = 21, cex = 1.5, bg = "white")



### Movement analysis
cL <- changeLight(twl = twl.gl, quantile = 0.78, days = 1)
mS <- mergeSites2(twl = twl.gl, site = cL$site, degElevation = gE[2]-1, distThreshold = 450, alpha = gE[3:4], method = "gamma", mask = NULL)


Seasonal_palette <- grDevices::colorRampPalette(grDevices::hsv(1 - ((1:365) + (365/4))%%365/365, s = 0.8, v = 0.8), 
                                                space = "Lab")


day  <- as.POSIXlt(aggregate(mS$twl$tFirst[mS$site>0], by = list(mS$site[mS$site>0]), FUN = median)$x, origin = "1970-01-01", tz  ="GMT")$yday
stp  <- as.numeric(aggregate(mS$twl$tFirst[mS$site>0], by = list(mS$site[mS$site>0]), FUN = function(x) difftime(x[length(x)],x[1], units = "days"))$x)
cexf <- approxfun(range(stp), c(3, 9), rule = 3)

plot(NA, xlim = range(mS$summary[,c(2,4:7)]), ylim = range(mS$summary[,c(3,8:10)]), type = "n", bty = "n", mgp = c(3,2,1), xlab = "", ylab = "", las = 1)
plot(wrld_simpl, add = T, col = "grey90", border = "grey90")

lines(mS$summary[,2], mS$summary[,3], lwd = 2, lty = 2)
arrows(mS$summary[,2], mS$summary[,8], mS$summary[,2], mS$summary[,11], lwd = 0.8, length = 0, col = adjustcolor("black", alpha.f = 0.6))
arrows(mS$summary[,2], mS$summary[,9], mS$summary[,2], mS$summary[,10], lwd = 2, length = 0)

arrows(mS$summary[,4], mS$summary[,3], mS$summary[,7], mS$summary[,3], lwd = 0.8, length = 0, col = adjustcolor("black", alpha.f = 0.6))
arrows(mS$summary[,5], mS$summary[,3], mS$summary[,6], mS$summary[,3], lwd = 2, length = 0)


arrows(mS$summary[,4], mS$summary[,3], mS$summary[,5], mS$summary[,3], lwd = 1.2, length = 0)
points(mS$summary[,2], mS$summary[,3], pch = 21, cex = c(1, cexf(stp[-c(1, length(stp))]), 1), bg = Seasonal_palette(365)[day], lwd = 2)
text(mS$summary[,2], mS$summary[,3], 1:nrow(mS$summary), col = c("transparent", rep(1, nrow(mS$summary)-2), "transparent"))
points(lon.calib, lat.calib, pch  = 21, bg = "white", lwd = 2, cex = 2)
mapplots::add.pie(x = -90, y = -50, z = rep(1, 12), radius = 10, col = Seasonal_palette(12), init.angle = day[1])



