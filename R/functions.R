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
	lapply(i, .setSpruceTypes, r=r, slope=slope, aspect=aspect, bs.probs=bs.probs, bwd.ids=bwd.ids)
}

#' Simultaneous fire spread.
#'
#' Spread all fires on the landscape simultaneously.
#'
#' This function spreads fires out from their origin cells simultaneously.
#' This is algorithmically and computationally a bit simpler than \code{.spreadSequential}
#' but it does not allow tracking of individual fires on the landcaspe by fire IDs.
#' Fire spread is recursive. This function should only be called by \code{spreadByRep},
#' not directly.
#'
#' @param x a vector of ignition locations (grid cell indices).
#' @param x.hold storage for \code{x} when function is called recursively.
#' @param sens numeric constant, the spread sensitivity of neighboring cells.
#' @param v a vector of vegetation flammabilities of all grid cells.
#' @param tr.br a vector of column number breaks assisting with recursive fire spread within the raster domain.
#' @param j recursion iterator.
#' @param weaken boolean, whether sensitivity should weaken with each round of burning.
#' @param ... arguments passed from \code{spreadByRep}.
#'
#' @return
#'
#' @examples
#' # not run
.spreadSimultaneous <- function(x, x.hold=x, sens, v, tr.br, j=1, weaken=TRUE, ...){
	if(!length(x)) return(x)
	#if(sens > 1) stop("sens cannot exceed 1.")
	print(paste("Recursion level:", j))
	if(weaken) sens <- sens*(1-exp(-50/j))
	j <- j+1
	smax <- 1/sens
	b <- sort(unique(c(x-1, x+1, mapply("+", x, MoreArgs=tr.br))))
	b <- b[!(b %in% x.hold | is.na(v[b]))]
	x.new <- b[runif(length(b), 0, smax) < v[b]]
	if(length(x.new)){
		x.hold <- c(x.hold, x.new)
		Recall(x=x.new, x.hold=x.hold, sens=sens, v=v, tr.br=tr.br, j=j)
	} else x.hold
}

#' Sequential fire spread.
#'
#' Spread fire on the landscape one fire at a time.
#'
#' This function spreads fire out from their origin cell sequentially
#' so that individual fire IDs can be associated with bunred cells.
#' Fire spread is recursive. This function should only be called by \code{spreadByRep},
#' not directly.
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
#' @return a matrix containing a column of burned cell indices and a column of fire IDs.
#'
#' @examples
#' # not run
.spreadSequential <- function(x, x.hold=x, xid=1:length(x), xid.hold=xid, sens, v, tr.br, j=1, weaken=TRUE, ...){
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
    x <- cbind(x.hold, xid.hold)[ord,]
    colnames(x) <- c("Cell", "FID")
    x
  }
}

#' Spread fire to surrounding grid cells
#'
#' Fire is spread recursively to surrounding grid cells from starting cells.
#'
#' If \code{sequential=TRUE} (default), each fire is ignited on the landscape
#' and runs its course before the next is ignited.
#' Otherwise all fires burn the landscape simultaneously.
#'
#' @param k integer, a dummy variable (iterator).
#' @param x a list of vectors of ignition points (fire origin grid cells).
#' @param v a list of vectors of vegetation flammability (all grid cells).
#' @param sequential boolean, whether fires burn the landscape one a time or all at once.
#' @param ... arguments passed to fire spread functions.
#'
#' @return if \code{sequential=TRUE}, a matrix containing a column of burned cell indices
#'  and a column of fire IDs. Otherwise, a vector of burned cell indices.
#' @export
#'
#' @examples
#' # not run
spreadByRep <- function(k, x, v, sequential=TRUE, ...){
  if(sequential)
    .spreadSequential(x=x[[k]], v=v[[k]], ...) else
      .spreadSimultaneous(x=x[[k]], v=v[[k]], ...)
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
		if(!(i %in% ignore.veg)) r.flam[ind] <- (k/(1+exp(a-b*r.age[ind])))*r.flam[ind]
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
#' @param k integer, a dummy variable (iterator).
#' @param a a list of raster layers of vegetation age.
#' @param i a list of vectors or matrices containing burned grid cell indices.
#'
#' @return a raster layer of updated vegetation age for iteration \code{k}.
#' @export
#'
#' @examples
#' # not run
updateAgeByRep <- function(k, a, i){
  updateAge <- function(a, i){
    a <- a + 1
    if(length(i)){
      if(is.matrix(i)) a[i["Cell"]] <- 0 else a[i] <- 0
    }
    a
  }
  updateAge(a=a[[k]], i=i[[k]])
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
#' @param k integer, a dummy variable (iterator).
#' @param v a list of raster layers of vegetation IDs.
#' @param spruce.type a list of raster layers of spruce trajectory IDs.
#' @param a a list of raster layers of vegetation age.
#' @param i a list of vectors or matrices containing burned cell indices.
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
#' @return a list of raster layers of updated vegetation classes.
#' @export
#'
#' @examples
#' # not run
updateVegByRep <- function(k, v, spruce.type, a, i,
                           trans.succeed=list('1'=1, '2'=2, '3'=3, '4'=c(2,3), '5'=5, '6'=5, '7'=7),
                           trans.colonize=list('1'=1, '2'=2, '3'=3, '4'=4, '5'=3, '6'=3, '7'=7),
                           trans.burn=list('1'=1, '2'=4, '3'=4, '4'=4, '5'=6, '6'=6, '7'=7),
                           bs.pars=c(1, 6, 0.2),
                           ws.pars=c(1, 10, 0.25),
                           gram.pars=c(1, 15, 0.25),
                           bwdsg=c(2, 3, 4, 5, 6)){

  updateVeg <- function(v, spruce.type, a, i,
                        trans.succeed=list('1'=1, '2'=2, '3'=3, '4'=c(2,3), '5'=5, '6'=5, '7'=7),
                        trans.colonize=list('1'=1, '2'=2, '3'=3, '4'=4, '5'=3, '6'=3, '7'=7),
                        trans.burn=list('1'=1, '2'=4, '3'=4, '4'=4, '5'=6, '6'=6, '7'=7),
                        bs.pars=c(1, 6, 0.2),
                        ws.pars=c(1, 10, 0.25),
                        gram.pars=c(1, 15, 0.25),
                        bwdsg=c(2, 3, 4, 5, 6)){
    if(is.matrix(i)) i <- i["Cell"]
    # Burn transition
    v0 <- v[i]
    v[i] <- trans.burn[[v0]]
    # Succession
    v1 <- v[-i]
    st1 <- spruce.type[-i]
    a1 <- a[-i]
    succession <- function(age, p) p[1]/(1+exp(p[2]-p[3]*age)) > runif(length(age))
    bs.ind <- v1==bwdsg[3] & st1==bwdsg[1]
    ws.ind <- v1==bwdsg[3] & st1==bwdsg[2]
    gram.ind <- v1==bwdsg[5]
    v1[bs.ind & succession(age=a1[bs.ind], p=bs.pars)] <- bwdsg[1]
    v1[ws.ind & succession(age=a1[ws.ind], p=ws.pars)] <- bwdsg[2]
    v1[gram.ind & succession(age=a1[gram.ind], p=gram.pars)] <- bwdsg[4]
    v[-i] <- v1
    # Colonization
    # Establish conditions for colonization
    v
  }

	updateVeg(v=v[[k]], spruce.type[[k]], a[[k]], i[[k]],
    trans.succeed=trans.succeed, trans.colonize=trans.colonize, trans.burn=trans.burn,
    bs.pars=bs.pars, ws.pars=ws.pars, gram.pars=gram.pars, bwdsg=bwdsg)
}

#' Landscape lightning strikes
#'
#' Strike the landscape randomly.
#'
#' This function is only used internally. It generates \code{n} lightning strikes where \code{n}
#' sampled from a normal distribution centered on \code{n.strikes} with a standard deviattion of five.
#' See the source code for \code{simulate} for context.
#'
#' @param ... arguments passed from \code{simulate}.
#'
#' @return a list of two vectors: a vector of struck cells and a vector of strike intensities.
#'
#' @examples
#' # not run
strike <- function(...){
	dots <- list(...)
	index=dots$index
	n.strikes=dots$n.strikes
	n <- round(rnorm(n.sim, n.strikes, 5)) # Number of lightning strikes
	cells <- lapply(1:n.sim, function(i, x, n) sample(x, size=n[i]), x=index, n=n) # Location of strikes
	intensity <- lapply(1:n.sim, function(i, n, ig) runif(n[i], 0, 1/ig), n=n, ig=ignit) # Strike intensity
	return(list(cells=cells, intensity=intensity))
}

#' Ignite cells by lightning strikes on landscape
#'
#' Ignite landscape cells based on their lightning sensitivty and vegetation flammability.
#'
#' Ignition is based on a simple boolean check of whether lightning sensitivity, \code{s}, exceeds
#' vegetation flammability, \code{r}, for each location in \code{x}. The subset of ignited locations is returned.
#'
#' @param i integer, a dummy variable (iterator).
#' @param x a list of vectors of strike locations.
#' @param s a list of vectors of strike intensity.
#' @param r a list of vectors of vegetation flammability.
#'
#' @return a vector of ignited locations.
#' @export
#'
#' @examples
#' # not run
ignite <- function(i, x, s, r) x[[i]][s[[i]] < r[[i]][x[[i]]]] # x=location, s=lightning intensity, r=vegetation flammability

# Multiple simulation replicates
simulate <- function(par.iter, n.sim, prob, Maps=F, simMaps=1:n.sim, r.spruce.type, r.age, r.veg, ...){
	r.spruce.type <- r.spruce.type[[par.iter]]
	r.age <- r.age[[par.iter]]
	r.veg <- r.veg[[par.iter]]
	ignore.veg <- list(...)$ignore.veg
	strikes <- strike(...)
	#r.spruce.type <- list(...)$r.spruce.type
	print("Completed lightning strikes on landscape.")
	r.flam <- lapply(1:n.sim, flamByAgeByRep, r.f=r.flam, r.v=r.veg, r.a=r.age, p=prob, iv=ignore.veg)
	print("Completed adjustment of climate- and vegetation-mediated flammability by vegetation age.")
	ig.pts <- lapply(1:n.sim, ignite, x=strikes$cells, s=strikes$intensity, r=r.flam) # cells where strike intensity overcomes vegetation flammability
	print("Completed ignitions on landscape.")
	v.flam <- rapply(r.flam, getValues, how="replace")
	print("Obtained list of vegetation flammability vectors.")
	cells.burned <- lapply(1:n.sim, spreadByRep, x=ig.pts, v=v.flam, ...) # Spread fire from ignition points
	print("Completed fire spread.")
	r.age <- lapply(1:n.sim, updateAgeByRep, a=r.age, i=cells.burned)
	print("Updated vegetation ages.")
	r.veg <- lapply(1:n.sim, updateVegByRep, v=r.veg, spruce.type=r.spruce.type, a=r.age, i=cells.burned)
	print("Updated vegetation types.")
	# Compile burn probability across simulations
	burn.table <- table(do.call(c, cells.burned))
	#burn.table <- (burn.table - min(burn.table))/max(burn.table - min(burn.table))
	r.burn[as.integer(names(burn.table))] <- burn.table
	n.fires <- lapply(ig.pts, length)
	if(Maps) list(SimBurnProb=r.burn, AgeByRep=r.age[simMaps], VegByRep=r.veg[simMaps]) else list(SimBurnProb=r.burn)
}

# Multiple simulation replicates
simulate2 <- function(par.iter, n.sim, prob, Maps=F, simMaps=1:n.sim, r.spruce.type, r.age, r.veg, b.flam, ...){
	r.spruce.type <- r.spruce.type[[par.iter]]
	r.age <- r.age[[par.iter]]
	r.veg <- r.veg[[par.iter]]
	ignore.veg <- list(...)$ignore.veg
	strikes <- strike(...)
	#r.spruce.type <- list(...)$r.spruce.type
	print("Completed lightning strikes on landscape.")
	r.flam <- lapply(1:n.sim, flamByAgeByRep, r.f=r.flam, r.v=r.veg, r.a=r.age, p=prob, iv=ignore.veg)
	print("Completed adjustment of climate- and vegetation-mediated flammability by vegetation age.")
	ig.pts <- lapply(1:n.sim, ignite, x=strikes$cells, s=strikes$intensity, r=r.flam) # cells where strike intensity overcomes vegetation flammability
	print("Completed ignitions on landscape.")
	v.flam <- rapply(r.flam, getValues, how="replace")
	print("Obtained list of vegetation flammability vectors.")

	#cells.burned <- lapply(1:n.sim, spreadByRep, x=ig.pts, v=v.flam, ...) # Spread fire from ignition points
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
		cells.burned[[i]] <- do.call(c, cells.burned.tmp) # store the cells burned by sim
	}

	print("Completed fire spread.")
	r.age <- lapply(1:n.sim, updateAgeByRep, a=r.age, i=cells.burned)
	print("Updated vegetation ages.")
	r.veg <- lapply(1:n.sim, updateVegByRep, v=r.veg, spruce.type=r.spruce.type, a=r.age, i=cells.burned)
	print("Updated vegetation types.")
	# Compile burn probability across simulations
	#burn.table <- table(do.call(c, cells.burned))
	####burn.table <- (burn.table - min(burn.table))/max(burn.table - min(burn.table))
	#r.burn[as.integer(names(burn.table))] <- burn.table
	if(Maps) list(CellsBurned=cells.burned, N=n.fires, FS=fs, TBA=tba, AgeByRep=r.age[simMaps], VegByRep=r.veg[simMaps]) else list(CellsBurned=cells.burned, N=n.fires, FS=fs, TBA=tba)
}

# Multiple simulation replicates
simulate3 <- function(par.iter, n.sim, prob, Maps=F, simMaps=1:n.sim, r.spruce.type, r.age, r.veg, b.flam, ...){
	r.spruce.type <- r.spruce.type[[par.iter]]
	r.age <- r.age[[par.iter]]
	r.veg <- r.veg[[par.iter]]
	veg0.ind <- which(r.veg[[1]][]==0)
	ignore.veg <- list(...)$ignore.veg
	n.yrs <- nlayers(b.flam)
	cells.burned.list <- n.fires.list <- fs.list <- tba.list <- vector("list", n.yrs)
	if(Maps) r.age.list <- r.veg.list <- vector("list", n.yrs)
	for(z in 1:n.yrs){
	print(z)
		r.flam <- subset(b.flam, z)
		r.flam[veg0.ind] <- 0 # Why are there postive flammablity probabilies in the flammability maps where veg ID is 0?

		strikes <- strike(...)
		#r.spruce.type <- list(...)$r.spruce.type
		print("Completed lightning strikes on landscape.")
		r.flam <- lapply(1:n.sim, flamByAgeByRep, r.f=r.flam, r.v=r.veg, r.a=r.age, p=prob, iv=ignore.veg)
		print("Completed adjustment of climate- and vegetation-mediated flammability by vegetation age.")
		ig.pts <- lapply(1:n.sim, ignite, x=strikes$cells, s=strikes$intensity, r=r.flam) # cells where strike intensity overcomes vegetation flammability
		print("Completed ignitions on landscape.")
		v.flam <- rapply(r.flam, getValues, how="replace")
		print("Obtained list of vegetation flammability vectors.")

		#cells.burned <- lapply(1:n.sim, spreadByRep, x=ig.pts, v=v.flam, ...) # Spread fire from ignition points
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
			cells.burned[[i]] <- do.call(c, cells.burned.tmp) # store the cells burned by sim
		}

		print("Completed fire spread.")
		r.age <- lapply(1:n.sim, updateAgeByRep, a=r.age, i=cells.burned)
		print("Updated vegetation ages.")
		r.veg <- lapply(1:n.sim, updateVegByRep, v=r.veg, spruce.type=r.spruce.type, a=r.age, i=cells.burned)
		print("Updated vegetation types.")
		# Compile burn probability across simulations
		#burn.table <- table(do.call(c, cells.burned))
		####burn.table <- (burn.table - min(burn.table))/max(burn.table - min(burn.table))
		#r.burn[as.integer(names(burn.table))] <- burn.table
		cells.burned.list[[z]] <- cells.burned
		n.fires.list[[z]] <- n.fires
		fs.list[[z]] <- fs
		tba.list[[z]] <- tba
		if(Maps){ r.age.list[[z]] <- r.age[simMaps]; r.veg.list[[z]] <-r.veg[simMaps] }
	}

	if(Maps) list(CellsBurned=cells.burned.list, N=n.fires.list, FS=fs.list, TBA=tba.list, AgeByRep=r.age.list[simMaps], VegByRep=r.veg.list[simMaps]) else list(CellsBurned=cells.burned.list, N=n.fires.list, FS=fs.list, TBA=tba.list)
}

# Multiple simulation replicates
simulate4 <- function(par.iter, n.sim, prob, Maps=F, simMaps=1:n.sim, r.spruce.type, r.age, r.veg, b.flam, ...){
	r.spruce.type <- r.spruce.type[[par.iter]]
	r.age <- r.age[[par.iter]]
	r.veg <- r.veg[[par.iter]]
	veg0.ind <- which(r.veg[[1]][]==0)
	ignore.veg <- list(...)$ignore.veg
	n.yrs <- nlayers(b.flam)
	burned.list <- n.fires.list <- fs.list <- tba.list <- vector("list", n.yrs)
	if(Maps) r.age.list <- r.veg.list <- vector("list", n.yrs)
	for(z in 1:n.yrs){
	print(z)
		r.flam <- subset(b.flam, z)
		r.flam[veg0.ind] <- 0 # Why are there postive flammablity probabilies in the flammability maps where veg ID is 0?
		strikes <- strike(...)
		#r.spruce.type <- list(...)$r.spruce.type
		print("Completed lightning strikes on landscape.")
		r.flam <- lapply(1:n.sim, flamByAgeByRep, r.f=r.flam, r.v=r.veg, r.a=r.age, p=prob, iv=ignore.veg)
		print("Completed adjustment of climate- and vegetation-mediated flammability by vegetation age.")
		ig.pts <- lapply(1:n.sim, ignite, x=strikes$cells, s=strikes$intensity, r=r.flam) # cells where strike intensity overcomes vegetation flammability
		if(!length(ig.pts[[1]])) next

		print("Completed ignitions on landscape.")
		v.flam <- rapply(r.flam, getValues, how="replace")
		print("Obtained list of vegetation flammability vectors.")

		n.fires <- lapply(ig.pts, length)  # store the total number of fires by sim
		burned <- lapply(1:n.sim, spreadByRep2, x=ig.pts, v=v.flam, ...) # Spread fire from ignition points
		cells.burned <- lapply(burned, function(x) x[,1])
		#print(head(burned[[1]]))
		#print(head(burned[[2]]))
		fs <- lapply(1:length(burned), function(i, a) as.numeric(tapply(a[[i]][,2], a[[i]][,2], length)), a=burned)
		tba <- lapply(fs, sum)

		print("Completed fire spread.")
		r.age <- lapply(1:n.sim, updateAgeByRep, a=r.age, i=cells.burned)
		print("Updated vegetation ages.")
		r.veg <- lapply(1:n.sim, updateVegByRep, v=r.veg, spruce.type=r.spruce.type, a=r.age, i=cells.burned)
		print("Updated vegetation types.")
		# Compile burn probability across simulations
		#burn.table <- table(do.call(c, cells.burned))
		####burn.table <- (burn.table - min(burn.table))/max(burn.table - min(burn.table))
		#r.burn[as.integer(names(burn.table))] <- burn.table
		burned.list[[z]] <- burned
		n.fires.list[[z]] <- n.fires
		fs.list[[z]] <- fs
		tba.list[[z]] <- tba
		if(Maps){ r.age.list[[z]] <- r.age[simMaps]; r.veg.list[[z]] <-r.veg[simMaps] }
	}

	if(Maps) list(Burn=burned.list, N=n.fires.list, FS=fs.list, TBA=tba.list, AgeByRep=r.age.list[simMaps], VegByRep=r.veg.list[simMaps]) else list(Burn=burned.list, N=n.fires.list, FS=fs.list, TBA=tba.list)
}
