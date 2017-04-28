library(raster)
rasterOptions(chunksize=10^12,maxmemory=10^11)
setwd("/workspace/UA/mfleonawicz/projects/randscape/workspaces")
#setwd("C:/github/Randscape/workspaces")

# Load functions
source("../code/functions.R")

agg <- TRUE

Mode <- function(x, na.rm){
  if(na.rm) x <- x[!is.na(x)]
  if(!length(x)) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

yrs <- 1900:2013
flamDir <- "../../Flammability/data"
# Objects related to flammability map
flam.files <- list.files(file.path(flamDir, "gbmFlammability/samples_based/historical/CRU32/3m100n_cavmDistTrunc_loop_Lmap"), pattern="\\.tif$", full=TRUE)
flam.files <- flam.files[match(yrs, as.numeric(gsub("[a-z._]", "", basename(flam.files))))]
b.flam <- readAll(stack(flam.files))#, quick=TRUE))

r.veg <- readAll(raster(file.path(flamDir, "alf2005.cavm.merged.030212.tif")))

# r.age <- readAll(raster("../data/Age_0_1900.tif")) weird map
r.age <- r.veg
r.age[!is.na(r.age) & r.age!=0] <- 10 # set all non-NA and non-zero veg cell ages to 10
r.age[r.age==0] <- NA

# What to do with these?
r.site <- readAll(raster("../data/centaksite.tif"))
r.slope <- readAll(raster("../data/centakslope.tif"))
r.tden <- readAll(raster("../data/centaktreedensity.tif"))
projection(r.site) <- projection(b.flam)

if(agg){
  b.flam <- aggregate(b.flam, 10)
  r.veg <- aggregate(r.veg, 10, fun=Mode, na.rm=T)
  r.age <- aggregate(r.age, 10)
  r.site <- aggregate(r.site, 10)
  r.slope <- aggregate(r.slope, 10, fun=Mode, na.rm=T)
  r.tden <- aggregate(r.tden, 10, fun=Mode, na.rm=T)
}

b.flam <- mask(b.flam, r.veg)
r.flam <- subset(b.flam, 1)
r.veg <- mask(r.veg, r.flam)
r.age <- mask(r.age, r.flam)
#r.age[r.age > 200] <- 200 # Testing only: truncate ages
#r.age[r.age==0] <- 1 # Testing only: truncate ages
#r.flam[r.veg==0] <- 0 # Why are there postive flammablity probabilies in the flammability maps where veg ID is 0? # this has been moved to the z loop below

#r.spruce <- readAll(raster("/workspace/Shared/Users/mfleonawicz/big_scratch_backup_101614/mfleonawicz/Alf_Files_20121129/Spinup300Year_32Reps/originalVeg/Veg_0_1900.tif"))
r.spruce <- r.veg
r.spruce[!(r.spruce %in% c(2,3,4))] <- NA
nc <- ncol(r.flam)
tr.br <- list(rep(c(-nc, nc), each=3) + c(-1:1))

# Create empty raster of burned cells
r.burn <- r.flam
r.burn[!is.na(r.burn)] <- 0

# Settings
fire.prob <- list(
  "Alpine tundra"=c(0.75, 4, 0.3),
  "Black spruce"=c(0.85, 4, 0.08),
  "White spruce"=c(0.7, 8, 0.08),
  Deciduous=c(0.1, 25, 10), # Not age dependent?
  Shrub=c(0.2, 25, 1.4),
  Graminoid=c(0.3, 25, 1.4),
  Wetland=c(0.1, 25, 1.4)
)

# Test age multipler for flammability
r.flam.age <- r.flam
for(i in 1:length(fire.prob)){
  pars <- fire.prob[[i]]
  k <- pars[1]
  a <- pars[2]
  b <- pars[3]
  ind <- r.veg[]==i
  r.flam.age[ind] <- (k/(1+exp(a-b*r.age[ind])))*r.flam[ind]
}

file <- "example_inputs.RData"
save(b.flam, fire.prob, r.age, r.burn, r.flam, r.site, r.slope, r.spruce, r.tden, r.veg, tr.br, yrs, file=file)
