require(ape)

dbcv <- function(partition, dataset, dist_method = "euclidean") {

	partition <- as.integer(partition)
	
	partitionOrg <- partition
	nObjsOrg     <- length(partition)       

	dims      <- dim(dataset)
	nFeatures <- dims[2]
	distance  <- as.matrix(dist(dataset, method = dist_method)^2)

	if (dims[1] != nObjsOrg) {
		stop('The number of objects has to the be the same from partition.')
	}

	# assign all singletons to noise class
	partition[which(partition %in% as.integer(names(which(table(partition)==1))))] <- -1

	# remove noisy objects from data
	noiseIdx  <- which(partition %in% -1)
	if (length(noiseIdx) != 0) {
		dataset   <- dataset[-noiseIdx,]
		partition <- partition[-noiseIdx]
		distance  <- distance[-noiseIdx,-noiseIdx]
		colnames(distance) <- seq(1,length(partition))
		rownames(distance) <- seq(1,length(partition))
	}

	clusters  <- sort(unique(partition))	# cluster's ids
	nClusters <- length(clusters)	        # number of unique clusters
	nObjs     <- length(partition)			  # number of objects in dataset


	if (nObjs <= 1 || is.null(clusters) || nClusters <= 1) {
	  stop('You need more than one object/cluster, other than noise.')
	}
	
	compactness <- c()
	internalVertices <- c()
	
	ducore <- rep(0,nObjs)
	
	for (i in 1:nClusters) {
	  objCluster <- which(partition == clusters[i])
	  nCluster   <- length(objCluster)
	  
	  mrd <- mmrd(distance[objCluster,objCluster],nFeatures,nCluster)
	  ducore[objCluster] <- mrd$ducore
	  
	  clusterMST <- spantree(mrd$objDist)
	  
	  edgeList <- c()
	  
	  for (e in 1:length(clusterMST$dist)) {
	    edgeList <- rbind(edgeList,c(as.integer(clusterMST$labels[e]),clusterMST$dist[e],as.integer(clusterMST$labels[clusterMST$kid[e]])))
	  }
	  
	  vertexDegree <- table(c(edgeList[,1],edgeList[,3]))

	  internalVertices[[i]] <- as.integer(names(which(vertexDegree > 1)))
	  
	  internalEdges    <- intersect(which(edgeList[,1]%in%internalVertices[[i]]),which(edgeList[,3]%in%internalVertices[[i]]))
	  
	  if (is.null(internalEdges)) {
	    compactness[i] <- max(edgeList[,2])
	  } else {
	    compactness[i] <- max(edgeList[internalEdges,2])
	  }
	  
	  if (is.null(internalVertices[i])) {
	    internalVertices[[i]] <- objCluster
	  }
	  
	}
	
	separationPoint <- matrix(0,nrow = nObjs, ncol = nObjs)
	for (i in 1:(nObjs-1)) {
	  for (j in i:nObjs) {
	    separationPoint[i,j] <- max(c(distance[i,j],ducore[i],ducore[j]))
	    separationPoint[j,i] <- separationPoint[i,j]
    }
  }
	
	cv <- 0
	
	clusterSeparation <- rep(Inf,nClusters)
	
	for (i in 1:nClusters) {
	  otherClusters <- setdiff(clusters,clusters[i])
	  
	  separation <- rep(0,length(otherClusters))
	  
	  for (j in 1:length(otherClusters)) {
  	  separation[j] = min(min(separationPoint[internalVertices[[i]],internalVertices[[otherClusters[j]]]])) 
	  }
	  clusterSeparation[i] = min(separation)
    dbcvcl = (clusterSeparation[i] - compactness[i]) / max(compactness[i],clusterSeparation[i]);
    cv <- cv + (dbcvcl * sum(partition == clusters[i]))
	}
	
	cv <- cv / nObjsOrg
	
	cv
}

mmrd <- function(objDist, dim, minPts) {
  
  kNNDist <- objDist^(-dim)
  kNNDist[is.infinite(kNNDist)] <- 0
  
  ducore <- colSums(kNNDist)
  ducore <- ducore / (minPts - 1)   
  ducore <- (1/ducore)^(1/(1*dim))
  
  ducore[is.infinite(ducore)] <- 0
  
  for (i in 1:minPts) {
    for (j in 1:minPts) {
      objDist[i,j] <- max(ducore[i],ducore[j],objDist[i,j])
      objDist[j,i] <- objDist[i,j]
    }
  }
  
  mrd <- c()
  mrd$ducore  <- ducore
  mrd$objDist <- objDist
  
  mrd
}
