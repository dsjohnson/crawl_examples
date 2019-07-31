library(crawl)
library(sf)
library(mapview)
library(tidyverse)

data("beardedSeals")

## Project data with sf##
beardedSeals <- st_as_sf(beardedSeals, coords=c("longitude","latitude"), crs=4326) %>% 
  st_transform(3338)

# Process ellipse information for crawl error structure
beardedSeals <- crawl::argosDiag2Cov(
  Major=beardedSeals$error_semimajor_axis, 
  Minor=beardedSeals$error_semiminor_axis,
  Orientation=beardedSeals$error_ellipse_orientation
) %>% bind_cols(beardedSeals, .)

# Break into individual data sets
beardedSeals <- beardedSeals %>% group_by(deployid) %>% nest()

# Fit models
beardedSeals <- beardedSeals %>% mutate(
  crw_fit = map(data,
                ~{
                  crwMLE(
                    err.model = list(x=~0+ln.sd.x, y=~0+ln.sd.y, rho=~error.corr),
                    data=.x, Time.name = "date_time",
                    fixPar = c(1,1,NA,NA), theta=c(8,-1),
                    attempts = 10
                  )
                }
  )
)

## Fit to just one animal if you don't use tidyverse functions
dat <- beardedSeals$data[[1]]
head(dat)
fit1 <- crwMLE(
  err.model = list(x=~0+ln.sd.x, y=~0+ln.sd.y, rho=~error.corr),
  data=dat, Time.name = "date_time",
  fixPar = c(1,1,NA,NA), theta=c(8,-1),
  attempts = 10
)

fit1

# predict hourly
beardedSeals <- beardedSeals %>% mutate(
  crw_predict = map(crw_fit, ~{
    crwPredict(object.crwFit=.x, predTime = "1 hour") %>% 
      crw_as_sf(ftype = "POINT") # make it an sf object
  }
  )
)

# or for just the first animal:
pred <- crwPredict(object.crwFit=fit1, predTime = "1 hour")
head(pred)

# Make a picture
beardedSeals %>% select(deployid, crw_predict) %>% unnest() %>% st_as_sf(crs=3338) %>% 
  mapview(.,zcol="deployid", burst=TRUE)

