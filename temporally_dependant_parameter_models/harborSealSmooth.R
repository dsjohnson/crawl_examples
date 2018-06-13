library(crawl)
library(ggplot2)
library(splines)
library(sp)
library(dplyr)
data(harborSeal)
head(harborSeal)
harborSeal$Argos_loc_class = factor(harborSeal$Argos_loc_class, levels=c("3","2","1","0","A","B"))
ak = nPacMaps::alaska()

## Project data ##

toProj = harborSeal[!is.na(harborSeal$latitude),c("Time","latitude","longitude")]
coordinates(toProj) = ~longitude+latitude
proj4string(toProj) <- CRS("+proj=longlat")
toProj <- spTransform(toProj, CRS("+init=epsg:3338"))
toProj = as.data.frame(toProj)
colnames(toProj)[2:3] = c("x","y")
harborSeal = merge(toProj, harborSeal, by="Time", all=TRUE)
harborSeal = harborSeal[order(harborSeal$Time),]


##Fit model as given in Johnson et al. (2008) Ecology 89:1208-1215
## Start values for theta come from the estimates in Johnson et al. (2008)

### Show all the parameters to provide start values and define a prior...
df=50
# 2m/s =2*60*60 m/h 
theta.start = c(rep(log(2000),3),log(2*60*60),rep(0,df),rep(0,df+1))

fixPar = c(log(250), log(500), log(1500), rep(NA,2*df+8-3), 0)

displayPar( mov.model=~bs(harborSeal$Time, df=df), err.model=list(x=~Argos_loc_class-1),data=harborSeal, 
            activity=~I(1-DryTime),fixPar=fixPar, theta=theta.start
)

constr=list(lower=c(rep(log(1500),3),rep(-Inf,2*df+5-3)), upper=rep(Inf,2*df+5))
omega_beta = 10
omega_sigma = 10

ln_prior = function(par){
  sum(dnorm(par[5:(df+4)], 0, omega_sigma, log=T)) + # normal prior for sigma coefs
    sum(dnorm(par[(df+6):(2*df+5)], 0, omega_beta, log=T)) # normal prior for beta_coefs
}

set.seed(321)
fit1 <- crwMLE(
  mov.model=~bs(harborSeal$Time, df=df), 
  err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal, coord=c("x","y"), Time.name="Time", 
  fixPar=fixPar, constr=constr, method="L-BFGS-B",
  theta = theta.start,
  prior=ln_prior,
  control=list(maxit=2000, trace=1, REPORT=1))

print(fit1)

predTimes = floor(seq(min(harborSeal$Time), max(harborSeal$Time), by=6))[-1]

pred1 = crwPredict(fit1, return.type = "flat") %>% crw_as_tibble()

require(ggplot2)
p1=ggplot(aes(x=mu.x, y=mu.y), data=pred1) + geom_path(col="red") + geom_point(aes(x=x, y=y), col="blue") + coord_fixed()
p2=ggplot(aes(x=Time, y=mu.x), data=pred1) + geom_ribbon(aes(ymin=mu.x-2*se.mu.x, ymax=mu.x+2*se.mu.x), fill="green", alpha=0.5)  + 
  geom_path(, col="red") + geom_point(aes(x=Time, y=x), col="blue", size=1)
p3=ggplot(aes(x=Time, y=mu.y), data=pred1) + geom_ribbon(aes(ymin=mu.y-2*se.mu.y, ymax=mu.y+2*se.mu.y), fill="green", alpha=0.5)  + 
  geom_path(, col="red") + geom_point(aes(x=Time, y=y), col="blue", size=1)
print(p1)
print(p2)
print(p3)

displayPar( mov.model=~bs(harborSeal$Time, df=df), err.model=list(x=~Argos_loc_class-1),data=harborSeal, 
            activity=~I(1-DryTime),fixPar=fit1$par)

# time series of velocity correlation
vel.cor = exp(-exp(fit1$mov.mf%*%fit1$par[(df+8):(2*df+8)]))
#time series of velosity variance parameter
vel.sigma = exp(fit1$mov.mf%*%fit1$par[7:(df+7)])
# plot of correlation
plot(fit1$data$Time, vel.cor, type='l')
# plot of corrleation v. sigma
plot(vel.cor, vel.sigma, type='l')


### Experimental clustering for calculating behavior states
library(cluster)
mov_par = data.frame(ln_beta=fit1$mov.mf%*%fit1$par[(df+8):(2*df+8)], ln_sigma=fit1$mov.mf%*%fit1$par[7:(df+7)])
d = dist(mov_par, method = "euclidean") # distance matrix
clust = hclust(d, method="ward") 
states = cutree(clust, 2)
pred1$states = factor(states)
mov_par$states = factor(states)

###
bb = c(range(pred1$mu.x),range(pred1$mu.y)) + c(-1,1,-1,1)*50000
ggplot(data=pred1) + geom_sf(data=ak, col=1, fill=1) + 
  geom_point(aes(x=mu.x, y=mu.y, col=states), alpha=0.3) +
  coord_sf(xlim=bb[1:2], ylim=bb[3:4]) + xlab("Longitude") + ylab("Latitude")
# ggsave(file="hs_smooth_map.png", dpi=300, width = 6, height=6, units="in")
ggplot(data=mov_par) + geom_point(aes(x=pred1$Time, y=exp(-exp(ln_beta)), col=states)) + 
  geom_path(aes(x=pred1$Time, y=exp(-exp(ln_beta)))) + xlab("Time") + ylab("OU vel corr")
# ggsave(file="hs_smooth_beta.png", dpi=300, width=4, height=1.75)
ggplot(data=mov_par) + geom_point(aes(x=pred1$Time, y=exp(ln_sigma), col=states)) + 
  geom_path(aes(x=pred1$Time, y=exp(ln_sigma))) + xlab("Time") + ylab("OU vel var")
# ggsave(file="hs_smooth_sigma.png", dpi=300, width=4, height=1.75)
