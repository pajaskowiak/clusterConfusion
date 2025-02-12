dbcv <- function(data, partition, noiseLabel = -1) {
  clusters <- unique(partition)
  dist <- as.matrix(dist(data))^2
  
  for (i in seq_along(clusters)) {
    if (sum(partition == clusters[i]) == 1) {
      partition[partition == clusters[i]] <- noiseLabel
      clusters[i] <- noiseLabels
    }
  }
  
  clusters <- setdiff(clusters, noiseLabel)
  
  if (length(clusters) == 0 || length(clusters) == 1) {
    return(0)
  }
  
  data <- data[partition != noiseLabel, ]
  dist <- dist[partition != noiseLabel, partition != noiseLabel]
  poriginal <- partition
  partition <- partition[partition != noiseLabel]
  
  nclusters <- length(clusters)
  nobjects <- nrow(data)
  nfeatures <- ncol(data)
  
  d_ucore_cl <- rep(0, nobjects)
  compcl <- rep(0, nclusters)
  int_edges <- vector("list", nclusters)
  int_node_data <- vector("list", nclusters)
  
  for (i in seq_along(clusters)) {
    objcl <- which(partition == clusters[i])
    nuobjcl <- length(objcl)
    
    mr <- matrix_mutual_reachability_distance(nuobjcl, dist[objcl, objcl], nfeatures)
    d_ucore_cl[objcl] <- mr$d_ucore
    G <- list(no_vertices = nuobjcl, MST_edges = matrix(0, nrow = nuobjcl - 1, ncol = 3), 
              MST_degrees = rep(0, nuobjcl), MST_parent = rep(0, nuobjcl))
    
    mst_results <- MST_Edges(G, 1, mr$G_edges_weights)
    Edges <- mst_results$Edg
    Degrees <- mst_results$Degr
    
    int_node <- which(Degrees != 1)
    int_edg1 <- which(Edges[, 1] %in% int_node)
    int_edg2 <- which(Edges[, 2] %in% int_node)
    int_edges[[i]] <- intersect(int_edg1, int_edg2)
    
    if (length(int_edges[[i]]) > 0) {
      compcl[i] <- max(Edges[int_edges[[i]], 3])
    } else {
      compcl[i] <- max(Edges[, 3])
    }
    
    int_node_data[[i]] <- objcl[int_node]
    if (length(int_node_data[[i]]) == 0) {
      int_node_data[[i]] <- objcl
    }
  }
  
  sep_point <- matrix(0, nrow = nobjects, ncol = nobjects)
  for (i in 1:(nobjects - 1)) {
    for (j in i:nobjects) {
      sep_point[i, j] <- max(c(dist[i, j], d_ucore_cl[i], d_ucore_cl[j]))
      sep_point[j, i] <- sep_point[i, j]
    }
  }
  
  valid <- 0
  sepcl <- rep(Inf, nclusters)
  for (i in seq_along(clusters)) {
    other_cls <- setdiff(clusters, clusters[i])
    sep <- sapply(other_cls, function(cls) {
      min(sep_point[int_node_data[[i]], int_node_data[[which(clusters == cls)]]])
    })
    sepcl[i] <- min(sep)
    dbcvcl <- (sepcl[i] - compcl[i]) / max(compcl[i], sepcl[i])
    valid <- valid + (dbcvcl * sum(partition == clusters[i]))
  }
  
  valid <- valid / length(poriginal)
  return(valid)
}

matrix_mutual_reachability_distance <- function(MinPts, G_edges_weights, d) {
  No <- nrow(G_edges_weights)
  
  K_NN_Dist <- G_edges_weights^(-1 * d)
  K_NN_Dist[K_NN_Dist == Inf] <- 0
  
  d_ucore <- colSums(K_NN_Dist)
  d_ucore <- d_ucore / (No - 1)
  d_ucore <- (1 / d_ucore)^(1 / (1 * d))
  d_ucore[d_ucore == Inf] <- 0
  
  for (i in 1:No) {
    for (j in 1:MinPts) {
      G_edges_weights[i, j] <- max(c(d_ucore[i], d_ucore[j], G_edges_weights[i, j]))
      G_edges_weights[j, i] <- G_edges_weights[i, j]
    }
  }
  
  return(list(d_ucore = d_ucore, G_edges_weights = G_edges_weights))
}

MST_Edges <- function(G, start, G_edges_weights) {
  intree <- rep(0, G$no_vertices)
  d <- rep(Inf, G$no_vertices)
  G$MST_parent <- 1:G$no_vertices
  
  d[start] <- 0
  v <- start
  counter <- 0
  
  while (counter < (G$no_vertices - 1)) {
    intree[v] <- 1
    dist <- Inf
    
    for (w in 1:G$no_vertices) {
      if (w != v && intree[w] == 0) {
        weight <- G_edges_weights[v, w]
        if (d[w] > weight) {
          d[w] <- weight
          G$MST_parent[w] <- v
        }
        if (dist > d[w]) {
          dist <- d[w]
          next_v <- w
        }
      }
    }
    
    counter <- counter + 1
    G$MST_edges[counter, ] <- c(G$MST_parent[next_v], next_v, G_edges_weights[G$MST_parent[next_v], next_v])
    G$MST_degrees[G$MST_parent[next_v]] <- G$MST_degrees[G$MST_parent[next_v]] + 1
    G$MST_degrees[next_v] <- G$MST_degrees[next_v] + 1
    v <- next_v
  }
  
  Edg <- G$MST_edges
  Degr <- G$MST_degrees
  
  return(list(Edg = Edg, Degr = Degr))
}


