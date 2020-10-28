library(tidyverse)
library(sf)
library(crawl)
library(mapview)
library(foreach)
library(doFuture)

registerDoFuture()


### Obtain fur seal data  
download.file("https://ndownloader.figshare.com/files/10364766", "nfs_pups_2005.csv")

### Read in pup telemetry data
pup_frame <- read_csv("nfs_pups_2005.csv", 
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
pup_frame <- pup_frame %>% st_as_sf(coords=c("long","lat")) %>% 
  st_set_crs(4326) %>% st_transform(3832)

### Split data by individual
pup_frame <- pup_frame %>% group_by(dbid, site, sex) %>% nest()


### Fit crawl model in parallel over animals
### Will be based on parallel R sessions 
plan("multisession", workers=6)

pup_frame$fit_pred <- foreach(i=1:nrow(pup_frame)) %dopar% {
      # fit ctcrw model
      set.seed(123)
      .x <- pup_frame$data[[i]]
      if(nrow(.x)<20) return(NA)
      fit <- crwMLE(
        err.model=list(x=~loc_class-1), 
        data = .x,
        Time.name="GMT",
        fixPar=c(log(250), log(500), log(1500), NA, NA, NA,NA),
        theta=c(log(2000), log(2000), 8, -3),
        constr=list(
          lower=c(rep(log(1500),2), rep(-Inf,2)), upper=Inf
        ), method='L-BFGS-B', attempts = 5)
      # Make some predictions at 6 hour intervals
      pred = crwPredict(fit, predTime = "1 day") %>% 
        crw_as_tibble() %>% filter(locType=="p")
      return(list(fit=fit, pred=pred))
}

plan("sequential")

pup_frame <- pup_frame %>% filter(!is.na(fit_pred)) %>% mutate(
# Just bookeeping to split fit and pred apart
    fit = map(fit_pred, ~{.x$fit}),
    pred = map(fit_pred, ~{.x$pred})
  ) %>% 
  select(-fit_pred) -> pup_frame



### Make some plots
# get shoreline 
np <- ptolemy::npac() # install with-- devtools::install_github("jmlondon/ptolemy")

pred_data <- pup_frame %>% filter(!is.na(pred)) %>% select(-fit) %>% unnest(cols=pred)
bb = c(range(pred_data$mu.x), range(pred_data$mu.y)) + c(-1, -1, 1, 1)*100000


### Plot by sex
ggplot(data=pred_data) + geom_point(aes(x=mu.x, y=mu.y, color=sex, group=dbid), alpha=0.2) + geom_sf(data=np, fill=1, color=1) +
  coord_sf(xlim =bb[1:2], ylim=bb[3:4]) + xlab("Longitude") + ylab("Latitude")
# ggsave(file="nfs_sex.png", width = 8, height=6)

### Plot by site
ggplot(data=pred_data) + geom_point(aes(x=mu.x, y=mu.y, color=site, group=dbid), alpha=0.2) + geom_sf(data=np, fill=1, color=1) +
  coord_sf(xlim =bb[1:2], ylim=bb[3:4]) + xlab("Longitude") + ylab("Latitude")
# ggsave(file="nfs_site.png", width = 8, height=6)

