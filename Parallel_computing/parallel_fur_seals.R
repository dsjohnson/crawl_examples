library(tidyverse)
library(sf)
library(crawl)
library(mapview)
library(furrr)


### Obtain fur seal data  
download.file("https://ndownloader.figshare.com/files/10364766", "nfs_pups_2005.csv")

### Read in pup telemetry data
pup_frame = read_csv("nfs_pups_2005.csv", 
                     col_types=cols(
                       dbid = col_integer(),
                       site = col_character(),
                       GMT = col_datetime(format = "%m/%d/%y %H:%M"),
                       Habitat = col_character(),
                       DateDeploy = col_date(format = "%m/%d/%y"),
                       loc_class = col_factor(levels = c('3', '2', '1', '0', 'A', 'B','Z')),
                       lat = col_double(),
                       long = col_double(),
                       sex = col_factor(levels = c('F','M'))
                     )
) %>% 
  dplyr::filter(loc_class!='Z', loc_class!="B") %>% 
  droplevels() %>% 
  dplyr::select(dbid, site, sex, GMT, long, lat, loc_class) 

### Convert to sf object and project
pup_frame = pup_frame %>% st_as_sf(coords=c("long","lat")) %>% 
  st_set_crs(4326) %>% st_transform(3832)

### Split data by individual
pup_frame = pup_frame %>% group_by(dbid, site, sex) %>% nest()


### Fit crawl model in parallel over animals
### Will be based on parallel R sessions 
plan("multisession")

pup_frame %>% #slice(90:97) %>% 
  mutate(
    # Fit a ctcrw model and make some predictions
    fit_pred = future_map(data, 
                         .f=~{
                           # fit ctcrw model
                           fit = crwMLE(
                             err.model=list(x=~loc_class-1), 
                             drift=TRUE,
                             data = .x,
                             Time.name="GMT",
                             fixPar=c(log(250), log(500), log(1500), NA, NA, NA,4, NA, NA),
                             theta=c(log(2000), log(2000), 5, 0, 0),
                             constr=list(
                               lower=c(rep(log(1500),2), rep(-Inf,3)), upper=rep(Inf,5)
                             ), method='L-BFGS-B', attempts = 5)
                           # Make some predictions at 6 hour intervals
                           pred = crwPredict(fit, predTime = "6 hour") %>% 
                             crw_as_tibble() %>% filter(locType=="p")
                           list(fit=fit, pred=pred)
                         }
    ),
    # Just bookeeping to split fit and pred apart
    fit = map(fit_pred, ~.x$fit),
    pred = map(fit_pred, ~.x$pred)
  ) %>% 
  select(-fit_pred) -> pup_frame

### Make some plots
# get shoreline 
np = nPacMaps::npac()

pred_data = pup_frame %>% select(-fit) %>% unnest(pred) %>% 
  st_as_sf(coords=c("mu.x","mu.y")) %>% st_set_crs(3832)

ggplot() + geom_sf(data=pred_data) + geom_sf(data=np)



