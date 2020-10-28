
################################################################################
### WARNING!!! Right now this is our best way to do this, we are currently working 
###            on another solution but this is the best we can do at the moment
################################################################################

library(dplyr)
library(sp)
library(crawl)
library(gdistance)

### The fix_path function currently in crawl doesn't work very well use the old
### version available here:
devtools::source_gist("ca62cbdbb773a828d4a41faa2c2746cd")


data("harborSeal")

harborSeal %>% dplyr::filter(!is.na(latitude)) %>% 
  as.data.frame() -> harborSeal_sp

sp::coordinates(harborSeal_sp) <- ~longitude + latitude
sp::proj4string(harborSeal_sp) <- CRS("+init=epsg:4326")
harborSeal_sp <- sp::spTransform(harborSeal_sp, CRS("+init=epsg:3338")) %>% 
  as("data.frame")

harborSeal %>%  
  dplyr::select( -longitude, -latitude) %>% 
  left_join(harborSeal_sp) -> harborSeal

levels(harborSeal$Argos_loc_class) = c("3","2","1","0","A","B")

#######################################################
# fit crawl model
######################################################

initial = list(
  a=c(harborSeal$longitude[1],0,
      harborSeal$latitude[1],0),
  P=diag(c(10000^2,5400^2,10000^2,5400^2))
)

fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0)

constr=list(
  lower=c(rep(log(1500),3), rep(-Inf,2)),
  upper=rep(Inf,5)
)

fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal, coord=c("longitude","latitude"), Time.name="Time", 
  initial.state=initial, fixPar=fixPar, theta=c(rep(log(5000),3),log(3*3600), 0),
  constr=constr,
  control=list(maxit=2000, trace=1, REPORT=1)#,
  #initialSANN=list(maxit=200, trace=1, REPORT=1)
)

predTimes <- seq(min(harborSeal$Time), max(harborSeal$Time), by = 0.5)

# Created predicted path
hs_pred <- crawl::crwPredict(fit1, predTime=predTimes)
coordinates(hs_pred) <- ~mu.x+mu.y
proj4string(hs_pred) <- CRS("+init=epsg:3338")

########################################################
# Obtain coastline using ptolemy
# If not installed:
# devtools::install_github("jmlondon/ptolemy")
######################################################

library(ptolemy)
library(raster)
library(gdistance)
# Get Alaska polygon
ak <- ptolemy::alaska()

## This is a pretty coarse raster, just for demonstration purposes
land <- raster(
  ext =extend(extent(bbox(hs_pred)),100000),
  resolution = 5000,
  crs = CRS(proj4string(hs_pred))
) %>% rasterize(ak,., getCover=TRUE)/100

plot(land)
lines(hs_pred$mu.x, hs_pred$mu.y)

#Create transition matrix
trans = transition(1-land, prod, directions = 16)

# Move path based on shortest distance around land
new_hs_pred = fix_path(hs_pred, hs_pred$Time, land, trans)

plot(as(new_hs_pred, "SpatialLines"), col="red")
plot(ak, col=gray(0.5), add=TRUE)

