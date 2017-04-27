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

# Recursive fire spread function
spread <- function(x, x.hold=x, sens, v, tr.br, j=1, weaken=TRUE, ...){
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

spreadByRep <- function(k, x, v, ...) spread(x=x[[k]], v=v[[k]], ...)

# Recursive fire spread function
prep_cells_ids <- function(x, tr.br){
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

rm_burned_na <- function(x, x.hold, v){
	ind <- !(x[[1]] %in% x.hold | is.na(v[x[[1]]])) # remove cells previously burned or NA-valued
	x[[1]] <- x[[1]][ind]
	x[[2]] <- x[[2]][ind]
	x
}

spread_to_cells <- function(x, v, smax){
	ind <- runif(length(x[[1]]), 0, smax) < v[x[[1]]]
	x[[1]] <- x[[1]][ind]
	x[[2]] <- x[[2]][ind]
	x
}

spread2 <- function(x, x.hold=x, xid=1:length(x), xid.hold=xid, sens, v, tr.br, j=1, weaken=TRUE, ...){
	if(!length(x)) return(x)
	#if(sens > 1) stop("sens cannot exceed 1.")
	print(paste("Recursion level:", j))
	if(weaken) sens <- sens*(1-exp(-50/j))
	j <- j+1
	smax <- 1/sens
	x.list <- prep_cells_ids(x, tr.br)
	x.list <- rm_burned_na(x.list, x.hold, v)
	x.list <- spread_to_cells(x.list, v, smax)
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

spreadByRep2 <- function(k, x, v, ...) spread2(x=x[[k]], v=v[[k]], ...)

# Scale climate-mediated vegetation flammability by age-mediated vegetation flammability
flamByAge <- function(r.flam, r.veg, r.age, prob, ignore.veg=0){
	#if(!all(c(hasArg(r.veg), hasArg(r.age), hasArg(prob)))) return(r.flam)
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

flamByAgeByRep <- function(k, r.f, r.v, r.a, p, iv=0) flamByAge(r.flam=r.f, r.veg=r.v[[k]], r.age=r.a[[k]], prob=p, ignore.veg=iv)

updateAge <- function(a, i){
	a <- a + 1
	if(length(i)) a[i] <- 0
	a
}

updateAgeByRep <- function(k, a, i)	updateAge(a=a[[k]], i=i[[k]])

trans.succeed <- c('1'=1, '2'=2, '3'=3, '4'=c(2,3), '5'=5, '6'=5, '7'=7)
trans.colonize <- c('1'=1, '2'=2, '3'=3, '4'=4, '5'=3, '6'=3, '7'=7)
trans.burn <- c('1'=1, '2'=4, '3'=4, '4'=4, '5'=6, '6'=6, '7'=7)

updateVeg <- function(v, spruce.type, a, i){
	# Burn transition
	v0 <- v[i]
	v[i] <- trans.burn[v0]
	# Succession
	v1 <- v[-i]
	st1 <- spruce.type[-i]
	a1 <- a[-i]
	bs.pars <- c(1, 6, 0.2)
	ws.pars <- c(1, 10, 0.25)
	gram.pars <- c(1, 15, 0.25)
	succession <- function(age, p) p[1]/(1+exp(p[2]-p[3]*age)) > runif(length(age))
	bs.ind <- v1==4 & st1==2
	ws.ind <- v1==4 & st1==3
	gram.ind <- v1==6
	v1[bs.ind & succession(age=a1[bs.ind], p=bs.pars)] <- 2
	v1[ws.ind & succession(age=a1[ws.ind], p=ws.pars)] <- 3
	v1[gram.ind & succession(age=a1[gram.ind], p=gram.pars)] <- 5
	v[-i] <- v1
	# Colonization
	# Establish conditions for colonization
	v
}

updateVegByRep <- function(k, v, spruce.type, a, i){
	updateVeg(v=v[[k]], spruce.type[[k]], a[[k]], i[[k]])
}

strike <- function(...){
	dots <- list(...)
	index=dots$index
	n.strikes=dots$n.strikes
	n <- round(rnorm(n.sim, n.strikes, 5)) # Number of lightning strikes
	cells <- lapply(1:n.sim, function(i, x, n) sample(x, size=n[i]), x=index, n=n) # Location of strikes
	intensity <- lapply(1:n.sim, function(i, n, ig) runif(n[i], 0, 1/ig), n=n, ig=ignit) # Strike intensity
	return(list(cells=cells, intensity=intensity))
}

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
