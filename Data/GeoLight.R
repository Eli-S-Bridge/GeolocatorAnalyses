
Species <- "OnOn"
ID = "OnOn52"

lat.calib <- 48.2
lon.calib <- 99.7


raw <- read.csv(paste0(wd, "/RawData/", Species, "/", ID, ".csv"))
  names(raw) <- c("Date", "Light")
  raw$Date  <- as.POSIXct(raw$Date, tz = "GMT")
  raw$Light <- log(raw$Light)

with(raw[8000:10000,], plot(Date, Light, type = "o", pch = 16))  
abline(h = 1.5)

offset <- 8
twl <- preprocessLight(raw, threshold = 2, offset = offset, lmax = 6)
twl <- subset(twl, !is.na(Twilight))

# write.csv(twl, paste0(wd, "/Results/", Species, "/", ID, "_twl.csv"), row.names = FALSE)
twl <- read.csv(paste0(wd, "/Results/", Species, "/", ID, "_twl.csv"))
  twl$Twilight <- as.POSIXct(twl$Twilight, tz = "GMT")

twl <- twl[!twl$Deleted,]
  
  
twl <- twilightAdjust(twl, interval = 600)

### Calibration
raw <- subset(raw, Date>=min(twl$Twilight) & Date<=max(twl$Twilight))

tm <- seq(min(raw$Date), max(raw$Date), by = "day")
rise <- rep(c(TRUE, FALSE), length(tm))

c.dat <- data.frame(Twilight = twilight(rep(tm, each = 2), lon = lon.calib, lat = lat.calib, 
                                        rise = rise, zenith = 96), Rise = rise)


lightImage(tagdata = raw, offset = offset, zlim = c(0, 6))
tsimagePoints(c.dat$Twilight, offset = offset, pch = 16, cex = 0.25,
              col = "orange")

tm.calib <- as.POSIXct(c(min(raw$Date), min(raw$Date)+40*24*60*60))
abline(v = tm.calib, col = "red")


###
twl.gl  <- export2GeoLight(twl) 

d.calib <- subset(twl.gl, tFirst>=tm.calib[1] & tSecond<=tm.calib[2])

gE <- getElevation(twl = d.calib, known.coord = c(lon.calib, lat.calib), lnorm.pars = TRUE)

crds <- coord(twl.gl, degElevation = gE[1])
tripMap(crds, xlim = range(crds[,1], na.rm = T), ylim = range(crds[,2], na.rm = T))
points(lon.calib, lat.calib, pch = 21, cex = 1.5, bg = "white")



### Movement analysis
cL <- changeLight(twl = twl.gl, quantile = 0.85, days = 2)
mS <- mergeSites(twl = twl.gl, site = cL$site, degElevation = gE[1])





