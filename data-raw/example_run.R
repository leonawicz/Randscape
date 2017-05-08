library(parallel)
library(raster)
library(dplyr)
library(purrr)
rasterOptions(chunksize=10^12,maxmemory=10^11)

agg <- 1 # aggregate from 1 sq. km to agg sq. km
use_files <- TRUE # flammability maps year by year via raster files from disk vs. brick in memory
setwd("/workspace/UA/mfleonawicz/projects/randscape/workspaces")
source("../code/functions.R")
load(paste0("example_inputs_", agg, "km.RData"))

if(use_files){
  b.flam <- list.files(paste0("../data/flam", agg, "km"), full=TRUE)
  n <- length(b.flam)
  yrs <- as.numeric(substr(basename(b.flam), 11, 14))
  r.flam <- raster(b.flam[1])
} else {
  n <- nlayers(b.flam)
  yrs <- as.numeric(substring(names(b.flam), 11))
  r.flam <- subset(b.flam, 1)
}

n.sim <- 4
n.parsim <- 16
r.veg.list <- mclapply(1:n.parsim, function(i, ...) lapply(1:n.sim, function(x) r.veg), n.sim, r.veg)
r.age.list <- mclapply(1:n.parsim, function(i, ...) lapply(1:n.sim, function(x) r.age), n.sim, r.age)
r.spruce.list <- mclapply(1:n.parsim, setSpruceTypes, n=n.sim, r=r.spruce, slope=r.slope, aspect=r.site, mc.cores=n.parsim)
inputs <- map(seq_along(r.veg.list), ~list(Veg=r.veg.list[[.x]], Age=r.age.list[[.x]], Spruce=r.spruce.list[[.x]], iter=.x))
n.strikes <- 500
ignit <- 5 # Ignition factor
sens <- 10 # Sensitivity factor
set.seed(856)
options(expressions=5e4)

#system.time({
#  Rprof("/workspace/UA/mfleonawicz/profiling.txt")
#  results <- simulate(inputs[[1]], n.strikes=n.strikes, ignit=ignit, sens=sens, b.flam=b.flam, years=yrs,
#    r.burn=r.burn, prob=fire.prob, tr.br=tr.br, ignore.veg=0, keep.maps=FALSE, verbose=FALSE)
#  Rprof()
#})
#summaryRprof("/workspace/UA/mfleonawicz/profiling.txt")

system.time({
  results <- mclapply(inputs, simulate, n.strikes=n.strikes, ignit=ignit, sens=sens, b.flam=b.flam, years=yrs,
    r.burn=r.burn, prob=fire.prob, tr.br=tr.br, ignore.veg=0, keep.maps=FALSE, verbose=FALSE, mc.cores=n.parsim)
})

fire <- map(results, ~.x$Fire) %>% bind_rows()
veg <- map(results, ~.x$Veg) %>% bind_rows()

ba.obs <- data.frame(Year=1950:2013, BA=as.integer(cellStats(fire.obs, sum)), Replicate=0L) %>% tbl_df
ba <- filter(fire, Year >= 1950) %>% group_by(Replicate, Year) %>% summarise(BA=sum(!is.na(Cell))*agg^2)
ba <- bind_rows(ba.obs, ba)
saveRDS(ba, "../data/outputs/ba.rds")

SimBurnProbByYear[[z]] <- Reduce("+", lapply(results, "[[", 1))/(n.sim*n.parsim) #results$SimBurnProb

# Analyze and save results
flam.out <- as.numeric(cellStats(do.call(brick, SimBurnProbByYear), mean))
flam.in <- as.numeric(cellStats(b.flam, mean))
cor(flam.in, flam.out)^2

keep.obj <- c("fire.prob", "b.flam", "r.flam.age", "SimBurnProbByYear", "r.age", "Sim1AgeByYear", "r.veg")
save(list=keep.obj, file="example_outputs.RData")
