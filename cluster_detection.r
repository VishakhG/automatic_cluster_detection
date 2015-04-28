require(popsom)


get_centroids <- function(map, centroids){
	xdim <- map$xdim
	ydim <- map$ydim
	xlist <- c()
	ylist <- c()
	xcoords <- coords$x
	ycoords <- coords$y

	for(ix in 1:xdim){
		for(iy in 1:ydim) {
			cx <- coords$xcoords[ix, iy]
			cy <- coords$ycoords[ix, iy]
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

cluster_spread <- function(x, y, heat, centroids, map){
	centroidx <- x
	centroidy <- y
	sum <- 0
	elements <- 0
	xdim <- map$xdim
	ydim <- map$ydim
	centroid_weight <- heat[centroidx,centroidy]
	for(xi in 1:xdim){
		for(yi in 1:ydim){
			cx <- centroids$x[xi,yi]
			cy <- centroids$y[xi,yi]

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

distance_between_clusters<-function(map, coords, centroids, umat){
	cluster_elements<-list_clusters(map,coords,centroids,umat)
	print(length(cluster_elements))

}
list_clusters<-function(map,coords,centroids,umat){
	cent_x <- centroids$xvals
	cent_y <- centroids$yvals
	componentx <- coords$x
	componenty <- coords$y
	cluster_list <- list()
	for(i in 1:length(centroids)){
		cx <- cent_x[i]
		cy <- cent_y[i]
		cluster_list[i] <- list_from_centroid(cx,cy,coords,umat)
	}
 cluster_list
}

list_from_centroid <- function(x, y, components, heat){
	centroidx <- x
	centroidy <- y
	sum <- 0
	xdim <- map$xdim
	ydim <- map$ydim
	centroid_weight <- heat[centroidx,centroidy]
	cluster_list<-c()
	for(xi in 1:xdim){
		for(yi in 1:ydim){
			cx <- components$x[xi,yi]
			cy <- components$y[xi,yi]

			if(cx == centroidx && cy == centroidy){
				cweight <- heat[xi,yi]
				cluster_list <- c(cluster_list,cweight)
			}
		}
	}
	
	list(cluster_list)
}


data <- read.csv("iris.csv", header=TRUE)
labels <- data[,3]
data <- data[0:3]
map <- map.build(data,xdim=10, ydim=5, alpha=.6, train=100)
umat <- compute.umat(map, smoothing=2)
coords <- compute.internal.nodes(map,umat, explicit=FALSE)
#plot(coords$x,coords$y)
#Get unique centroids
centroids<-get_centroids(map, coords)
print(centroids)
#get distance from centroid to cluster elements
#for each cluster (centroid)
within_cluster_dist <- distance_from_centroids(map, coords, centroids,umat)
between_cluster_dist <- distance_between_clusters(map, coords, centroids, umat)



	
