# Still need to incorporate:
# lightning density maps
# historical observed fire perimeters
# weaken spread by vegetation class
# ignition/sensitivity by vegetation class
# site? slope? tree density? what are these and how should they affect burning
# how to deal with deciduous?
# vegetation transition rules

#' Set spruce types
#'
#' This function establishes black vs. white spruce trajectories on the landscape.
#'
#' Existing black and white spruce cells in a raster vegetation map layer remain as they are.
#' Deciduous cells are converted randomly to black or white spruce as a function of each cell's slope and aspect.
#' This function should only be called by \code{setSpruceTypesByRep}, not directly.
#'
#' @param k integer, a dummy variable (iterator).
#' @param r raster layer, vegetation IDs.
#' @param slope raster layer, slope.
#' @param aspect raster layer, aspect.
#' @param bs.probs.slope vector of black spruce probabilities associated with aspect bins.
#' White spruce probabilities are their reciprocal.
#' @param bs.probs.aspect vector of black spruce probabilities associated with aspect bins.
#' White spruce probabilities are their reciprocal.
#' @param bwd.ids integer vector of ID codes for black spruce, white spruce and deciduous, respectively.
#'
#' @return A raster layer of spruce trajectories.
.setSpruceTypes <- function(k, r, slope, aspect, bs.probs.slope, bs.probs.aspect, bwd.ids){
	# First aspect
	ind <- which(r[]==bwd.ids[3] & !is.na(aspect[]))
	bins <- seq(-2, 360, length=9)
	names(bs.probs.aspect) <- 1:length(bs.probs.aspect)
	aspect <- cut(aspect, breaks=bins)
	aspect.v <- aspect[ind]
	bs.ind <- bs.probs.aspect[aspect.v] > runif(length(ind))
	r[ind[bs.ind]] <- bwd.ids[1]
	# Then slope
	ind <- which(r[]==bwd.ids[3] & !is.na(slope[]))
	names(bs.probs.slope) <- 1 + 0:1
	slope.v <- slope[ind] + 1
	bs.ind <- bs.probs.slope[slope.v] > runif(length(ind))
	r[ind[bs.ind]] <- bwd.ids[1]
	r[r==bwd.ids[3]] <- bwd.ids[2]
	r
}

#' Set spruce types
#'
#' This function establishes black vs. white spruce trajectories on the landscape.
#'
#' Existing black and white spruce cells in a raster vegetation map layer remain as they are.
#' Deciduous cells are converted randomly to black or white spruce as a function of each cell's slope and aspect.
#'
#' @param k integer, a dummy variable (iterator).
#' @param i integer vector (inner iterator).
#' @param r raster layer, vegetation IDs.
#' @param slope raster layer, slope.
#' @param aspect raster layer, aspect.
#' @param bs.probs.slope vector of black spruce probabilities associated with aspect bins.
#' Defaults provided. White spruce probabilities are their reciprocal.
#' @param bs.probs.aspect vector of black spruce probabilities associated with aspect bins.
#' Defaults provided. White spruce probabilities are their reciprocal.
#' @param bwd.ids integer vector of ID codes for black spruce, white spruce and deciduous, respectively.
#'
#' @return A length-\code{i} list of raster layers of spruce trajectories.
#' @export
#'
#' @examples
#' # not run
setSpruceTypesByRep <- function(k, i, r, slope, aspect,
  bs.probs.slope=c(0.8, 0.95, 0.95, 0.8, 0.2, 0.05, 0.05, 0.2),
  bs.probs.aspect=c(0.9, 0.1), bwd.ids=c(2, 3, 4)){
	lapply(i, .setSpruceTypes, r=r, slope=slope, aspect=aspect,
	       bs.probs.slope=bs.probs.slope, bs.probs.aspect=bs.probs.aspect, bwd.ids=bwd.ids)
}

#' Sequential fire spread.
#'
#' Spread fire on the landscape one fire at a time.
#'
#' This function spreads fire out from their origin cell sequentially
#' so that individual fire IDs can be associated with bunred cells.
#' Fire spread is recursive.
#'
#' @param x a vector of ignition locations (grid cell indices).
#' @param x.hold storage for \code{x} when function is called recursively.
#' @param xid integer vector of fire IDs to be associated with each unique fire origin cell in \code{x}.
#' @param xid.hold storage for \code{x} when function is called recursively.
#' @param sens numeric constant, the spread sensitivity of neighboring cells.
#' @param v a vector of vegetation flammabilities of all grid cells.
#' @param tr.br a vector of column number breaks assisting with recursive fire spread within the raster domain.
#' @param j recursion iterator.
#' @param weaken boolean, whether sensitivity should weaken with each round of burning.
#' @param ... arguments passed from \code{spreadByRep}.
#'
#' @return a table containing a column of burned cell indices and a column of fire IDs.
#' @export
#'
#' @examples
#' # not run
spread <- function(x, x.hold=x, xid=1:length(x), xid.hold=xid, sens, v, tr.br, j=1, weaken=TRUE, ...){
  if(!length(x)) return(x)
  #if(sens > 1) stop("sens cannot exceed 1.")
  print(paste("Recursion level:", j))
  if(weaken) sens <- sens*(1-exp(-50/j))
  j <- j+1
  smax <- 1/sens

  .prep_cells_ids <- function(x, tr.br){
    n <- length(x)
    x <- c(x-1, x+1, mapply("+", x, MoreArgs=tr.br)) # boundary cells
    xid <- c(rep(1:n, 2), rep(1:n, each=6)) # corresponding fire IDs
    x.ind <- which(!duplicated(x)) # remove overlap (randomized input order)
    xid <- xid[x.ind]
    x <- x[x.ind]
    xid <- xid[order(x)] # Sort by cell index
    x <- sort(x)
    list(x, xid)
  }

  .rm_burned_na <- function(x, x.hold, v){
    ind <- !(x[[1]] %in% x.hold | is.na(v[x[[1]]])) # remove cells previously burned or NA-valued
    x[[1]] <- x[[1]][ind]
    x[[2]] <- x[[2]][ind]
    x
  }

  .spread_to_cells <- function(x, v, smax){
    ind <- runif(length(x[[1]]), 0, smax) < v[x[[1]]]
    x[[1]] <- x[[1]][ind]
    x[[2]] <- x[[2]][ind]
    x
  }

  x.list <- .prep_cells_ids(x, tr.br)
  x.list <- .rm_burned_na(x.list, x.hold, v)
  x.list <- .spread_to_cells(x.list, v, smax)
  if(length(x.list[[1]])){
    x.hold <- c(x.hold, x.list[[1]])
    xid.hold <- c(xid.hold, x.list[[2]])
    Recall(x=x.list[[1]], x.hold=x.hold, xid=x.list[[2]], xid.hold=xid.hold, sens=sens, v=v, tr.br=tr.br, j=j)
  } else {
    ord <- order(xid.hold)
    x <- dplyr::tbl_df(data.frame(x.hold, xid.hold)) %>% dplyr::slice(ord)
    names(x) <- c("Cell", "FID")
    x
  }
}

#' Mediate vegetation flammability by vegetation age.
#'
#' Scale climate-mediated vegetation flammability input raster
#' using vegetation age-mediated vegetation flammability.
#'
#' This function takes a raster layer representing vegetation flammability
#' already mediated by climate and further scales the flammability values
#' in the raster grid cells based on vegetation age. The effect of vegetation
#' age on flammability is also vegetation class-dependent. There is a different
#' set of parameters, \code{k}, \code{a} and \code{b} for each vegetation class, hence
#' why \code{prob} is a list of length-\code{3} vectors.
#' The coefficient applied to each cell in \code{r.flam} is \eqn{(k/(1+exp(a[i]-b[i]*r.age[idx[i]]))} where
#' \code{i} refers to each vegetation class and \code{idx} refers to the cell indices with that vegetation.
#'
#' @param r.flam a raster layer of climate-mediated vegetation flammability.
#' @param r.veg a raster layer of vegetation cover type class IDs.
#' @param r.age a raster layer of vegetation age.
#' @param prob a list of vectors of parameters influencing the logistic curve that
#' represents the probability of fire based on vegetation age.
#' @param ignore.veg integer vector of vegetation ID codes to be ignored.
#'
#' @return a raster layer of vegetation flammabilities.
#' @export
#'
#' @examples
#' # not run
flamByAge <- function(r.flam, r.veg, r.age, prob, ignore.veg=0){
	n.v <- length(prob)
	for(i in 1:n.v){
		pars <- prob[[i]]
		k <- pars[1]
		a <- pars[2]
		b <- pars[3]
		ind <- r.veg[]==i
		if(!(i %in% ignore.veg))
		  r.flam[ind] <- ( k / (1 + exp(a - b*r.age[ind])) ) * r.flam[ind]
	}
	r.flam
}

#' Mediate vegetation flammability by vegetation age.
#'
#' Scale climate-mediated vegetation flammability input raster
#' using vegetation age-mediated vegetation flammability.
#'
#' This function takes a raster layer representing vegetation flammability
#' already mediated by climate and further scales the flammability values
#' in the raster grid cells based on vegetation age. The effect of vegetation
#' age on flammability is also vegetation class-dependent. There is a different
#' set of parameters, \code{k}, \code{a} and \code{b} for each vegetation class, hence
#' why \code{prob} is a list of length-\code{3} vectors.
#' The coefficient applied to each cell in \code{r.flam} is \eqn{(k/(1+exp(a[i]-b[i]*r.age[idx[i]]))} where
#' \code{i} refers to each vegetation class and \code{idx} refers to the cell indices with that vegetation.
#'
#' @param i integer, a dummy variable (iterator).
#' @param r.flam a length-\code{i} list of raster layers of climate-mediated vegetation flammability.
#' @param r.veg a length-\code{i} list of raster layers of vegetation cover type class IDs.
#' @param r.age a length-\code{i} list of raster layer of vegetation age.
#' @param prob a list of vectors of parameters influencing the logistic curve that
#' represents the probability of fire based on vegetation age. The length is equal to the number of
#' vegetation classes under consideration. If an ID code is not present in \code{r.veg} and is not ignored
#' via \code{ignore.veg}, an error is thrown.
#' @param ignore.veg integer vector of vegetation ID codes to be ignored.
#'
#' @return a raster layer of vegetation flammabilities.
#' @export
#'
#' @examples
#' # not run
flamByAgeByRep <- function(i, r.flam, r.veg, r.aage, prob, ignore.veg=0){
  flamByAge(r.flam=r.flam, r.veg=r.veg[[i]], r.age=r.age[[i]], prob=prob, ignore.veg=ignore.veg)
}



#' Title
#'
#' @param r.age a raster layer of vegetation age.
#' @param cells.burned a table containing burned grid cell indices.
#'
#' @return a raster layer of updated vegetation age.
#' @export
#'
#' @examples
#' # not run
updateAge <- function(r.age, cells.burned){
  r.age <- r.age + 1
  if(length(cells.burned)) r.age[cells.burned$Cell] <- 0
  r.age
}

#' Update vegetation
#'
#' Update vegetation following burning, succession and colonization.
#'
#' This function currently only implements updating vegetation classes of grid cells
#' in resonse to fire and general ecological succession. Colonization is not yet implemented.
#' While there are several hardcoded vegetation IDs and logisitic curve parameters,
#' they are provided as defaults to function arguments so they can still be changed if needed.
#' However, the various forms of vegetation transition from one vegetation class to another
#' are assumed fixed and coded into the function. They involve black spruce, white spruce,
#' deciduous, shrub tundra and graminoid tundra.
#'
#' @param r.veg a raster layer of vegetation IDs.
#' @param r.spruce a raster layer of spruce trajectory IDs.
#' @param r.age a raster layer of vegetation age.
#' @param cells.burned a table containing burned cell indices.
#' @param trans.succeed a list mapping the direction of succession from one vegetation
#' class to possible different classes. Not currently in use.
#' @param trans.colonize a list mapping the direction of colonization from one vegetation
#' class to possible different classes. Not currently implemented.
#' @param trans.burn a list mapping the direction of fire transition from one vegetation
#' class to possible different classes.
#' @param bs.pars a vector of parameters \code{k}, \code{a} and \code{b} describing the logistic
#' curve representing age-based ecological succession to black spruce from
#' black spruce-trajectory deciduous using the equation \eqn{(k/(1+exp(a-b*age))}.
#' @param ws.pars the same as \code{bs.pars} but succeeding to white spruce from
#' white spruce-trajectory deciduous.
#' @param gram.pars the same as \code{bs.pars} but succeeding to shrub tundra from graminoid tundra.
#' @param bwdsg integer vector, ID codes for black spruce, white spruce, deciduous, shrub tundra
#' and graminoid tundra, respectively. These must correspond to the \code{trans*} arguments
#' and vice versa.
#'
#' @return a raster layer of updated vegetation classes.
#' @export
#'
#' @examples
#' # not run
updateVeg <- function(r.veg, r.spruce, r.age, cells.burned,
                           trans.succeed=list('1'=1, '2'=2, '3'=3, '4'=c(2,3), '5'=5, '6'=5, '7'=7),
                           trans.colonize=list('1'=1, '2'=2, '3'=3, '4'=4, '5'=3, '6'=3, '7'=7),
                           trans.burn=list('1'=1, '2'=4, '3'=4, '4'=4, '5'=6, '6'=6, '7'=7),
                           bs.pars=c(1, 6, 0.2),
                           ws.pars=c(1, 10, 0.25),
                           gram.pars=c(1, 15, 0.25),
                           bwdsg=c(2, 3, 4, 5, 6)){
  i <- cells.burned$Cell
  # Burn transition
  v0 <- r.veg[i]
  r.veg[i] <- trans.burn[[v0]]
  # Succession
  v1 <- r.veg[-i]
  spr <- r.spruce[-i]
  a1 <- r.age[-i]
  succession <- function(age, p) p[1]/(1+exp(p[2]-p[3]*age)) > runif(length(age))
  bs.ind <- v1==bwdsg[3] & spr==bwdsg[1]
  ws.ind <- v1==bwdsg[3] & spr==bwdsg[2]
  gram.ind <- v1==bwdsg[5]
  v1[bs.ind & succession(age=a1[bs.ind], p=bs.pars)] <- bwdsg[1]
  v1[ws.ind & succession(age=a1[ws.ind], p=ws.pars)] <- bwdsg[2]
  v1[gram.ind & succession(age=a1[gram.ind], p=gram.pars)] <- bwdsg[4]
  r.veg[-i] <- v1
  # Colonization
  # Establish conditions for colonization
  r.veg
}

#' Landscape lightning strikes
#'
#' Strike the landscape randomly.
#'
#' This function is only used internally. It generates \code{n} lightning strikes where \code{n}
#' sampled from a normal distribution centered on \code{n.strikes} with a standard deviation of five.
#' See the source code for \code{simulate} for context.
#'
#' @param n.sim number of samples of strike count drawn from a eqn\{~Normal(n.strikes, 25)} distribution.
#' @param n.strikes mean number of strikes to center strikes sample around.
#' @param index integer vector indicating grid cells in raster layers which contain data values
#' (i.e., not \code{NA}-valued).
#' @param ignit numeric, ignition probability.
#'
#' @return a list of two vectors: a vector of struck cells and a vector of strike intensities.
#'
#' @examples
#' # not run
strike <- function(n.sim, n.strikes, index, ignit){
	n <- round(rnorm(n.sim, n.strikes, 5)) # Number of lightning strikes
	cells <- purrr::map(1:n.sim, ~sample(index, size=n[.x])) # Location of strikes
	intensity <- purrr::map(1:n.sim, ~runif(n[.x], 0, 1/ignit)) # Strike intensity
	return(list(cells=cells, intensity=intensity))
}

#' Ignite cells by lightning strikes on landscape
#'
#' Ignite landscape cells based on their lightning sensitivty and vegetation flammability.
#'
#' Ignition is based on a simple boolean check of whether lightning sensitivity, \code{s}, exceeds
#' vegetation flammability, \code{r}, for each location in \code{x}. The subset of ignited locations is returned.
#'
#' @param x a vector of strike locations.
#' @param s a vector of strike intensity.
#' @param r a vector of vegetation flammability.
#'
#' @return a vector of ignited locations.
#' @export
#'
#' @examples
#' # not run
ignite <- function(x, s, r) x[s < r[x]] # x=location, s=lightning intensity, r=vegetation flammability

# Multiple simulation replicates
simulate <- function(par.iter, prob, Maps=F, simMaps, r.spruce.type, r.age, r.veg, b.flam, verbose=TRUE, ...){
	r.spruce.type <- r.spruce.type[[par.iter]]
	r.age <- r.age[[par.iter]]
	r.veg <- r.veg[[par.iter]]
	veg0.ind <- which(r.veg[[1]][]==0)
	ignore.veg <- list(...)$ignore.veg
	n.strikes <- list(...)$n.strikes
	ignit <- list(...)$ignit
	sens <- list(...)$sens
	index <- which(!is.na(subset(b.flam, 1)[])) # all layers have same data vs NA cell indices
	n.yrs <- nlayers(b.flam)
	n.sim <- length(r.age)
	if(missing(simMaps)) simMaps <- 1:n.sim
	burned.list <- n.fires.list <- fs.list <- tba.list <- vector("list", n.yrs)
	r.age.list <- r.veg.list <- vector("list", n.yrs)

	for(z in 1:n.yrs){
	print(z)
		r.flam <- subset(b.flam, z)
		r.flam[veg0.ind] <- 0 # Why are there postive flammablity probabilies in the flammability maps where veg ID is 0?
		strikes <- strike(n.sim, n.strikes, index, ignit)
		if(verbose) cat("Completed lightning strikes on landscape.\n")
		r.flam <- purrr::map(1:n.sim, ~flamByAge(r.flam, r.veg[[.x]], r.age[[.x]], prob=prob, ignore.veg=ignore.veg))
		if(verbose) cat("Completed adjustment of climate- and vegetation-mediated flammability by vegetation age.\n")
		ig.pts <- purrr::map(1:n.sim, ~ignite(x=strikes$cells[[.x]], s=strikes$intensity[[.x]], r=r.flam[[.x]]))
		if(!length(ig.pts[[1]])) next

		if(verbose) cat("Completed ignitions on landscape.\n")
		v.flam <- rapply(r.flam, getValues, how="replace")
		if(verbose) cat("Obtained list of vegetation flammability vectors.\n")

		n.fires <- lapply(ig.pts, length)  # store the total number of fires by sim
		fs <- cells.burned <- tba <- vector("list", n.sim)
		for(i in 1:n.sim){
		  v.flam.tmp <- v.flam[[i]]
		  cells.burned.tmp <- vector("list", n.fires[[i]])
		  for(j in 1:n.fires[[i]]){
		    cells.burned.tmp[[j]] <- spread(x=ig.pts[[i]][j], v=v.flam.tmp, ...)
		    v.flam.tmp[cells.burned.tmp[[j]]] <- 0
		  }
		  fs[[i]] <- sapply(cells.burned.tmp, length) # store the fire sizes by sim
		  tba[[i]] <- sum(fs[[i]]) # store the total burn area by sim
		  cells.burned[[i]] <- dplyr::bind_rows(cells.burned.tmp) # store the cells burned by sim
		}

		if(verbose) cat("Completed fire spread.\n")
		r.age <- purrr::map(1:n.sim, ~updateAge(a=r.age[[.x]], i=cells.burned[[.x]]))
		if(verbose) cat("Updated vegetation ages.\n")
		r.veg <- purrr::map(1:n.sim, ~updateVeg(r.veg[[.x]], r.spruce.type[[.x]], r.age[[.x]], cells.burned[[.x]]))
		if(verbose) cat("Updated vegetation types.\n")
		cells.burned.list[[z]] <- cells.burned
		n.fires.list[[z]] <- n.fires
		fs.list[[z]] <- fs
		tba.list[[z]] <- tba
		if(Maps){ r.age.list[[z]] <- r.age[simMaps]; r.veg.list[[z]] <- r.veg[simMaps] }
	}

	list(
	  CellsBurned=cells.burned.list, N=n.fires.list, FS=fs.list, TBA=tba.list,
    AgeByRep=r.age.list[simMaps], VegByRep=r.veg.list[simMaps]
  )
}
