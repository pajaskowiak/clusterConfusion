#' @title aucc
#'
#' @description Computes the Area Under the (ROC) Curve for Clustering (AUCC).
#'
#' @param partition A hard/crisp partition (array of integers specifying cluster memberships).
#' @param dataset The dataset from which the partition was obtained.
#' @param distance The distances between object pairs. This should be obtained with the dist function (please, see ?dist).
#' @param distanceMethod The distance method that should be employed to compute distances from the dataset (please, check ?dist for more information on options). This parameter is used only if distance is not provided, given that they should be computed internally by the aucc function from the dataset.
#' @param returnRates If TRUE, TPR and FPR are returned alongside AUCC in an object of class "aucc" (useful for plot purposes). Default is FALSE, in which case the function directly outputs a numeric value that corresponds to the AUCC of the partition.
#'
#' @returns If returnRates=FALSE (default), the function returns a numerical value, which corresponds to the AUCC of the partition (clustering solution). If returnRates=TRUE, the function returns an object of class "aucc", with the following components/fields:
#'
#' * aucc: the AUCC value of the partition.
#' * tpr: the corresponding True Positive Rate (TPR) values.
#' * fpr: the corresponding False Positive Rate (FPR) values.
#'
#' @details
#'
#' Function computes the Area Under the (ROC) Curve for Clustering (AUCC) of a partition.
#'
#' You need to provide a partitioning of the data (clustering solution) and one of the following:
#'
#' (i)   the dataset from which the partition was originated;
#'
#' (ii)  the distance matrix for the data from which the partition was originated.
#'
#' In case of (i), you may also specify the distance function that will be employed to compute the distances internally by aucc, using the distanceMethod parameter. If you do not specify any string, it defaults to the Euclidean Distance. Make sure to pass strings accepted by the dist function (package stats), see ?dist.
#'
#' In case of (ii), make sure that the distance matrix was originated/converted with the function
#' dist from the package stats, using the parameters diag = FALSE, upper = FALSE (these are the
#' defaults for dist). Other than that, you may choose any distance / dissimilarity you find
#' appropriate for your data. If you have to perform multiple computations of aucc for the same
#' dataset, with different partitions, prefer this approach as computing the distance matrix
#' repeatedly can be computationally expensive.
#'
#'
#' @note
#' * Make sure to validate your partition(s) with the same distance you used to generate it.
#'
#' @references
#'
#' Jaskowiak, P.A., Costa, I.G. & Campello, R.J.G.B. The area under the ROC curve as a measure of clustering quality. Data Mining and Knowledge Discovery 36, 1219–1245 (2022). [https://doi.org/10.1007/s10618-022-00829-0](https://doi.org/10.1007/s10618-022-00829-0).
#'
#' @examples
#' # Example 1 (basic usage, returnRates=FALSE)
#'
#' # Computing AUCC from the dataset itself, partition with k = 2
#' clustering <- kmeans(ruspini,2)$cluster
#' print(aucc(clustering, dataset = ruspini))
#'
#' # Computing AUCC from a precomputed distance matrix and selecting best k
#' ruspiniDistances <- dist(ruspini)
#'
#' auccValues <- c()
#'
#' #varying k in a reasonable range...
#' maxK <- ceiling(sqrt(dim(ruspini)[1]))
#'
#' for (k in 2:maxK) {
#'    clustering <- kmeans(ruspini,k)$cluster
#'    auccValues <- append(auccValues, aucc(clustering, distance = ruspiniDistances))
#' }
#'
#' # Printing the suggested number of clusters and plotting the results
#' sprintf('Suggested number of clusters is: %d', which.max(auccValues)+1)
#' plot(2:maxK,auccValues)
#'
#' # Example 2 (getting full report, returnRates=TRUE)
#'
#' clustering <- kmeans(ruspini,2)$cluster
#' evaluation <- aucc(clustering, dataset = ruspini, returnRates=TRUE)
#'
#' print(evaluation$aucc)
#'
#' #calling package function to plot the corresponding ROC Curve
#' rocPlot(evaluation)
#'
#' @seealso [rocPlot()] [clusteringROCs()]
#'
#' @export
#' @importFrom precrec evalmod
#' @importFrom stats dist
#' @importFrom rdist rdist
#'

aucc <- function(partition, dataset = NULL, distance = NULL, distanceMethod = 'euclidean', returnRates = FALSE) {

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

  pred     <- c(rdist::rdist(partition,metric = 'hamming'))
  distance <- (distance - min(distance)) / (max(distance) - min(distance))

  if (!returnRates) {
    r <- as.numeric(as.data.frame(precrec::evalmod(scores = 1 - distance, labels = 1 - pred,mode='aucroc'))['aucs'])
  } else {
    basicEval <- precrec::evalmod(scores = 1-distance, labels = 1-pred,mode='basic')
    r         <- c()
    r$tpr     <- basicEval$sn[[1]]$y
    r$fpr     <- 1-basicEval$sp[[1]]$y
    r$aucc    <- as.numeric(as.data.frame(precrec::evalmod(scores = 1-distance, labels = 1-pred,mode='aucroc'))['aucs'])
    class(r) <- 'aucc'
  }

  r

}
