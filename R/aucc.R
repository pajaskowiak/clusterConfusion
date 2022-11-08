#' @title aucc
#'
#' @description Computes the Area Under the (ROC) Curve for Clustering (AUCC)
#'
#' @param partition A hard/crisp partition
#' @param dataset The dataset from which the partition was obtained
#' @param distance Distances between object pairs, obtained with the dist function
#' @param distanceMethod The distance method employed to compute distances from the dataset (see ?dist)
#'
#' @return The AUCC of the partition (clustering solution)
#'
#' @details
#'
#' Function computes the Area Under the (ROC) Curve for Clustering (AUCC) of a clustering solution.
#'
#' You need to provide a partitioning of the data (clustering solution) and one of the following:
#' (i)   the dataset from which the partition is originated;
#' (ii)  the distance matrix from which the partition was originated.
#'
#' In case of (i), you may also specify the distance function that will be employed to
#' compute the distances internally by aucc, using the distanceMethod parameter. If you do not
#' specify any string, it defaults to the Euclidean Distance. Make sure to pass strings accepted
#' by the function dist (package stats).
#'
#' In case of (ii), make sure that the distance matrix was originated/converted with the function
#' dist from the package stats, using the parameters diag = FALSE, upper = FALSE (these are the
#' defaults for dist). Other than that, you may choose any distance / dissimilarity you find
#' appropriate for your data. If you have to perform multiple computations of aucc for the same
#' dataset, with different partitions, prefer this approach as computing the distance matrix
#' is computationally expensive.
#'
#' Make sure to validate your partition with the same distance you used to generate it.
#'
#' @references
#'
#' Jaskowiak, P.A., Costa, I.G. & Campello, R.J.G.B. The area under the ROC curve as a measure of clustering quality. Data Mining and Knowledge Discovery 36, 1219–1245 (2022). https://doi.org/10.1007/s10618-022-00829-0
#'
#' @examples
#'
#' \dontrun{
#' library(cluster) #for ruspini dataset and k-means
#'
#' # Computing AUCC from the dataset itself, partitions with k = 2
#'
#' clustering <- kmeans(ruspini,2)$cluster
#' aucc(clustering, ruspini)
#'
#' # Computing AUCC from a precomputed distance matrix and selecting best k
#' myDistance <- dist(ruspini)
#'
#' auccValues <- c()
#'
#' #varying k in a reasonable range...
#' for (k in 2:ceiling(sqrt(dim(ruspini)[1]))) {
#'    clustering <- kmeans(ruspini,k)$cluster
#'    auccValues <- append(auccValues, aucc(clustering, distance = myDistance))
#' }
#'
#' sprintf('Suggested number of clusters is: %d', which.max(auccValues)+1)
#' }
#'
#' @export
#' @importFrom precrec evalmod
#' @importFrom stats dist
#' @importFrom rdist rdist
#'

aucc <- function(partition, dataset = NULL, distance = NULL, distanceMethod = 'euclidean') {

  if (missing(dataset) && missing(distance)) {
    stop('You need to specify a distance matrix or a dataset.')
  }

  if (missing(partition)) {
    stop('You need to specify a hard partition - clustering solution.')
  }

  if (!is.null(distance) && !is.null(dataset)) {
    stop('You can only specify a dataset or a distance, not both.')
  }

  if (!missing(distance) && !is.null(distance) && class(distance) != 'dist') {
    stop('The distance has to be computed with the dist function (stats).
    If you want other distances or dissimilarities no available,
    check its documentation in in order to use the desired one.
    Please, see: ?dist')

  }

  if (!is.null(dataset) && is.null(distance)) {
    if (dim(dataset)[1] != length(partition)) {
      stop('The number of objects has to the be the same from partition.')
    }
    distance <- c(dist(dataset, diag = FALSE, upper = FALSE, method = distanceMethod))
  } else {
    if (is.null(dataset) && !is.null(distance)) {
      distance <- c(distance)
      if (length(partition)*(length(partition)-1)/2 != length(distance)) {
        stop('Distance and partitions sizes don\'t match. If n is the number of objects, the required length of distance is n(n-1)/2.')
      }
    }
  }

  pred     <- c(rdist(partition,metric = 'hamming'))
  distance <- (distance - min(distance)) / (max(distance) - min(distance))

  as.numeric(as.data.frame(evalmod(scores = 1-distance, labels = 1-pred,mode='aucroc'))['aucs'])

}
