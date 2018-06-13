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

pup_frame %>% 
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
                               lower=c(rep(log(1500),2), rep(-Inf,3)), upper=Inf
                             ), method='L-BFGS-B', attempts = 5)
                           # Make some predictions at 6 hour intervals
                           pred = crwPredict(fit, predTime = "1 day") %>% 
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
np = nPacMaps::npac() # install with-- devtools::install_github("jmlondon/nPacMaps")

pred_data = pup_frame %>% select(-fit) %>% unnest(pred)
bb = c(range(pred_data$mu.x), range(pred_data$mu.y)) + c(-1, -1, 1, 1)*100000


### Plot by sex
ggplot(data=pred_data) + geom_point(aes(x=mu.x, y=mu.y, color=sex, group=dbid), alpha=0.2) + geom_sf(data=np, fill=1, color=1) +
  coord_sf(xlim =bb[1:2], ylim=bb[3:4]) + xlab("Longitude") + ylab("Latitude")
ggsave(file="nfs_sex.png", width = 8, height=6)

### Plot by site
ggplot(data=pred_data) + geom_point(aes(x=mu.x, y=mu.y, color=site, group=dbid), alpha=0.2) + geom_sf(data=np, fill=1, color=1) +
  coord_sf(xlim =bb[1:2], ylim=bb[3:4]) + xlab("Longitude") + ylab("Latitude")
ggsave(file="nfs_site.png", width = 8, height=6)

