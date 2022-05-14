library(crawl)
library(sf)
library(mapview); mapviewOptions(fgb=FALSE)
library(tidyverse)

data(northernFurSeal)

## Project data with sf##
northernFurSeal <- st_as_sf(northernFurSeal, coords=c("long","lat"), crs=4326)
northernFurSeal <- northernFurSeal %>% st_transform(
  crs="+proj=aea +lat_1=30 +lat_2=70 +lat_0=52 +lon_0=-170 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  )


fixPar <- c(log(250), log(500), log(1500), rep(NA,3), NA, rep(NA,2))
displayPar( mov.model=~1, err.model=list(x=~loc_class-1), drift = TRUE,
            data=northernFurSeal,fixPar=fixPar)
constr <- list(
  lower=c(rep(log(1500),2), -Inf, rep(-Inf,3)),
  upper=rep(Inf,6)
)

ln.prior <- function(theta){-abs(theta[4]+3)/0.5}

set.seed(123)
fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~loc_class-1), drift=TRUE,
  data=northernFurSeal, Time.name="GMT", time.scale="hour",
  fixPar=fixPar, constr=constr,# prior=ln.prior,
  theta = c(rep(log(1500)+1,2), log(3000), 0, 0,0),
  method="L-BFGS-B",
  control=list(trace=1, REPORT=1)
)
fit1
##Make hourly location predictions
predObj <- crwPredict(fit1, predTime="6 hours", return.type = "flat") 
plt <- mapview(predObj %>% crw_as_sf(ftype="LINESTRING"), cex=1, color='red', map.types="Esri.WorldImagery")
print(plt)

##Create simulation object with 100 parameter draws
set.seed(123)
simObj <- crwSimulator(fit1, predTime="6 hour", method="quadrature", quad.ask=F)


## Sample 20 tracks from posterior predictive distribution
samp_df <- tibble(rep=1:20) %>% mutate(
  samp = map(rep, ~crwPostIS(simObj, fullPost = TRUE) %>% crw_as_sf("LINESTRING"))
) 
plt <- plt + mapview(samp_df$samp, color = 'red', alpha=0.3)
print(plt)
      