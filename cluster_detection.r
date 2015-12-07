#Author: Vishakh Gopu
#An attempt to combine connected components in a Self organizing map 
#If they represent the same cluster. 
require(popsom)

#Get the unique centroids
get_centroids <- function(map, coords){
	xdim <- map$xdim
	ydim <- map$ydim
	xlist <- c()
	ylist <- c()
	xcoords <- coords$x
	ycoords <- coords$y

	for(ix in 1:xdim){
		for(iy in 1:ydim) {
			cx <- xcoords[ix, iy]
			cy <- ycoords[ix, iy]
			if(!(cx %in% xlist) || !(cy %in% ylist)){
				xlist <- c(xlist,cx)
				ylist <- c(ylist,cy)
			}
			if(!(cy %in% ylist) || !(cx %in% xlist) ){
				xlist <- c(xlist,cx)
				ylist <- c(ylist,cy)
			}
		}
	}
		list(xvals=xlist, yvals=ylist)
}	

#Get average distance from centroid by cluster
distance_from_centroids <- function(map, centroids, id_centroids, heat){
	xdim <- map$xdim
	ydim <- map$ydim
	xcoords <- id_centroids$xvals
	ycoords <- id_centroids$yvals
	within <- c()
	for (i in 1:length(xcoords)){
		cx <- xcoords[i]
		cy <- ycoords[i]
		distance <- cluster_spread(cx, cy, heat, centroids, map)
		within <- c(within, distance)
	}
	within
}

#Function to calculate the average distance in
#one cluster given one centroid
cluster_spread <- function(x, y, heat, centroids, map){
	centroidx <- x
	centroidy <- y
	sum <- 0
	elements <- 0
	xdim <- map$xdim
	ydim <- map$ydim
	centroid_weight <- heat[centroidx, centroidy]
	for(xi in 1:xdim){
		for(yi in 1:ydim){
			cx <- centroids$x[xi, yi]
			cy <- centroids$y[xi, yi]
			if(cx == centroidx && cy == centroidy){
				cweight <- heat[xi,yi]
				sum <- sum+abs(cweight-centroid_weight)
				elements <- elements+1
			}
		}
	}
	
	average <- sum/elements
	average
}

#The average pairwise distance between clusters
distance_between_clusters <- function(map, coords, centroids, umat){
	cluster_elements <- list_clusters(map, coords, centroids, umat)
	cluster_elements <- sapply(cluster_elements,'[',
							   seq(max(sapply(cluster_elements, length))))
	
	columns<-ncol(cluster_elements)
	cluster_elements <- matrix(unlist(cluster_elements), 
		                       ncol = ncol(cluster_elements), byrow = FALSE)
	cluster_elements <- apply(combn(ncol(cluster_elements), 2), 2, function(x)
							  abs(cluster_elements[, x[1]] - cluster_elements[, x[2]]))
	mean <- colMeans(cluster_elements, na.rm=TRUE)
	index <- 1
	mat <- matrix(data=NA, nrow=columns, ncol=columns)
	for(xi in 1:(columns-1)){
		for (yi in xi:(columns-1)){
			mat[xi, yi+1] <- mean[index]
			mat[yi+1, xi] <- mean[index]
		    index <- index+1
		}
	}
	mat
}

#Get the clusters as a list of lists
list_clusters <- function(map, coords, centroids, umat){
	cent_x <- centroids$xvals
	cent_y <- centroids$yvals
	componentx <- coords$x
	componenty <- coords$y
	cluster_list <- list()
	for(i in 1:length(cent_x)){
		cx <- cent_x[i]
		cy <- cent_y[i]
		cluster_list[i] <- list_from_centroid(cx, cy, coords, umat)
	}
 cluster_list
}

#Get all cluster elements associated to one centroid
list_from_centroid <- function(x, y, components, heat){
	centroidx <- x
	centroidy <- y
	sum <- 0
	xdim <- map$xdim
	ydim <- map$ydim
	centroid_weight <- heat[centroidx, centroidy]
	cluster_list <- c()
	for(xi in 1:xdim){
		for(yi in 1:ydim){
			cx <- components$x[xi,yi]
			cy <- components$y[xi,yi]

			if(cx == centroidx && cy == centroidy){
				cweight <- heat[xi,yi]
				cluster_list <- c(cluster_list, cweight)
			}
		}
	}
	list(cluster_list)
}

#Boolean matrix representing which clusters should be combined
combine_decision <- function(within_cluster_dist, distance_between_clusters){
	inter_cluster <- distance_between_clusters
	centroid_dist <- within_cluster_dist
	dim <- dim(inter_cluster)[1]
	to_combine <- matrix(data=FALSE, nrow=dim, ncol=dim)
	
	for(xi in 1:dim){
		for(yi in 1:dim){
			cdist <- inter_cluster[xi,yi]
			if(! is.na(cdist)){
				rx <- centroid_dist[xi]*(.1)
				ry <- centroid_dist[yi]*(.1)
				if( cdist < (centroid_dist[xi]+rx) || cdist < (centroid_dist[yi]+ry)){
					to_combine[xi,yi] <- TRUE
				}
			}
		}
	}
	to_combine
}

#Changes every instance of a centroid to one that it should be combined with
swap_centroids <- function(x1, y1, x2, y2, components){
	xdim <- map$xdim
	ydim <- map$ydim
	compn_x <- components$x
	compn_y <- components$y
	for(xi in 1:xdim){
		for(yi in 1:ydim){
			if(compn_x[xi] == x1 && compn_y[yi] == y1){
				compn_x[xi] <- x2
				compn_y[yi] <- y2
			}
		}

	}
	list(xcoords=compn_x, ycoords=compn_y)
}


#Combine centroids based on matrix of booleans
new_centroid <- function(bools, heat, components, centroids, map){
	xdim <- dim(bools)[1]
	print(xdim)
	ydim <- dim(bools)[2]
	centroids_x <- centroids$xvals
	centroids_y <- centroids$yvals
	components <- components
	for(xi in 1:xdim){
		for(yi in 1:ydim){
			if(bools[xi,yi] == TRUE){
				x1 <- centroids_x[xi]
				y1 <- centroids_y[xi]
				
				x2 <- centroids_x[yi]
				y2 <- centroids_y[yi]

				components <- swap_centroids(x1, y1, x2, y2, components)
				
			}
		} 
	}
	components
}	
################## Modifications of POPSOM to evaluate clustering #############	

#Modifying the heatmap method to take the combined components
plot_heat_mod <- function(map, heat, components, explicit=FALSE, comp=TRUE) {
	components <- components
	labels <- map$labels
	if (is.null(labels))
		stop("plot.heat: no labels available")

	x <- map$xdim
	y <- map$ydim
	nobs <- nrow(map$data)
	count <- array(data=0,dim=c(x,y))

	# need to make sure the map doesn't have a dimension of 1
	if (x > 1 && y > 1) {
		# bin the heat values into 100 bins used for the 100 heat colors below
		heat.v <- as.vector(heat)
		heat.v <- cut(heat.v, breaks=100, labels=FALSE)
		heat <- array(data=heat.v, dim=c(x,y))
	}

	# set up the graphics window
	par.v <- map.graphics.set()
	plot.new()
	plot.window(xlim=c(0,x),ylim=c(0,y))
	box()
	
	title(xlab="x",ylab="y")
	
	xticks <- seq(0.5,x-0.5,1)
	yticks <- seq(0.5,y-0.5,1)
	xlabels <- seq(1,x,1)
	ylabels <- seq(1,y,1)
	axis(1,at=xticks,labels=xlabels)
	axis(3,at=xticks,labels=xlabels)
	axis(2,at=yticks,labels=ylabels)
	axis(4,at=yticks,labels=ylabels)
		
	# plot heat
	colors<- heat.colors(100)
	
	for (ix in 1:x) {
		for (iy in 1:y) {
			rect(ix-1,iy-1,ix,iy,col=colors[100 - heat[ix,iy] + 1],border=NA)
		}
	}
	
	# put the connected component lines on the map
	if (comp) {

		# compute the connected components
		coords <- components

		for(ix in 1:x){
			for (iy in 1:y) {
				cx <- coords$xcoords[ix,iy]
				cy <- coords$ycoords[ix,iy]
				points(c(ix,cx)-.5,c(iy,cy)-.5,type="l",col="grey")
			}
		}
	}

	# put the labels on the map
	# count the labels in each map cell
	for(i in 1:nobs){
		ix <- map$visual$x[i]
		iy <- map$visual$y[i]
		count[ix+1,iy+1] <- count[ix+1,iy+1]+1
	}
	
	for(i in 1:nobs){
		ix <- map$visual$x[i]
		iy <- map$visual$y[i]
		# we only print one label per cell
		if (count[ix+1,iy+1] > 0) {
			count[ix+1,iy+1] <- 0
			ix <- ix + .5
			iy <- iy + .5
			l <- labels[i,1]
			text(ix,iy,labels=l)
		}
	}

	map.graphics.reset(par.v)
}
#Modification of the starbust  function to plot modified conected components
plot_starburst_mod <- function(map,umat,components,explicit=FALSE,smoothing=2) {

	if (class(map) != "map")
		stop("map.starburst: first argument is not a map object.")

	umat <- umat
	plot_heat_mod(map,umat,components,explicit=explicit,comp=TRUE)
}

####################  TOP LEVEL for testing ###########################
data <- read.csv("/Users/Vishakh/Desktop/Projects/Cluster_detection/iris.csv", header=TRUE)
labels <- data[,5]
data <- data[0:4]
map <- map.build(data,labels=labels,xdim=25, ydim=20, alpha=.6, train=100000)
png(filename="/Users/Vishakh/Desktop/Projects/Cluster_detection/old_starburst.png")
#Plot Starburst without modification
map.starburst(map)
dev.off()
umat <- compute.umat(map, smoothing=2)
coords <- compute.internal.nodes(map, umat, explicit=FALSE)
#Get unique centroids
centroids<-get_centroids(map, coords)
#Get distance from centroid to cluster elements for all centroids
within_cluster_dist <- distance_from_centroids(map, coords, centroids,umat)
#Get average pairwise distance between clusters
between_cluster_dist <- distance_between_clusters(map, coords, centroids, umat)
#Get a boolean matrix of whether two components should be combined
combine_cluster_bools <- combine_decision(within_cluster_dist, between_cluster_dist)
#Create the modified connected components grid
new_centroid <- new_centroid(combine_cluster_bools,heat,coords,centroids,map)
png(filename="/Users/Vishakh/Desktop/Projects/Cluster_detection/new_starburst.png")
#Plot modified starburst
plot_heat_mod(map,umat,new_centroid)
dev.off()