#' @title clusteringROCs
#'
#' @description Given a list of partitions computes their AUCC and corresponding ROC Curves.
#'
#' @param partitions A list of partitions.
#' @param dataset The dataset from which the partitions were obtained.
#' @param distance The distances between object pairs. This should be obtained with the dist function (please, see ?dist).
#' @param distanceMethod The distance method that should be employed to compute distances from the dataset (please, check ?dist for more information on options). This parameter is used only if distance is not provided, given that they should be computed internally by the clusteringROCs function from the dataset.
#' @param showPlot Defines whether the plot should be displayed by the function or not.
#' @param decoratePlot Whether the plot should be decorated with the best AUCC values or not.
#' @param alphaPlot If TRUE, the intensity of the ROC Curves will be proportional to their AUCC. Otherwise, all ROC Curves will be plotted with the same intensity (alpha/transparency).
#' @param rLine Whether to plot or not the random line that accounts for a random clustering solution.
#' @param sample Defines if all FPR and TPR should be used to generate the plot (slower), or just a subsample (faster).
#' @param sampleSize If sample=TRUE, defines the number of FPR and TPR pairs that will be used to generate the plot.
#'
#' @returns Returns a 'clusteringROCs' object containing three fields:
#'
#' * auccs: The corresponding values of AUCC for each one of the partitions in the list provided as input.
#' * rocplot: A ROC Plot, with one curve for each one of the partitions provided in the input list.
#' * auccplot: A Plot of the k values x AUCC values.
#'
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
#' Jaskowiak, P.A., Costa, I.G. & Campello, R.J.G.B. The area under the ROC curve as a measure of clustering quality. Data Mining and Knowledge Discovery 36, 1219–1245 (2022). [https://doi.org/10.1007/s10618-022-00829-0](https://doi.org/10.1007/s10618-022-00829-0).
#'
#' @examples
#'
#' set.seed(666)
#'
#' partitions <- list()
#'
#' #varying k in a reasonable range...
#' for (k in 2:ceiling(sqrt(dim(ruspini)[1]))) {
#'    partitions <- append(partitions,list(kmeans(ruspini,k)$cluster))
#' }
#' res <- clusteringROCs(partitions,dataset = ruspini,show = FALSE, decorate = TRUE, alphaPlot = FALSE)
#'
#' require(patchwork) #for combining the plots.
#' res$auccplot + res$rocplot
#'
#' # Now we just add more and more partitions to the list, with same number of k
#' # Same
#'
#' for (r in 1:10) {
#'    for (k in 2:ceiling(sqrt(dim(ruspini)[1]))) {
#'       partitions <- append(partitions,list(kmeans(ruspini,k)$cluster))
#'    }
#' }
#'
#' res <- clusteringROCs(partitions,dataset = ruspini,show = FALSE, decorate = TRUE, alphaPlot = FALSE)
#'
#' require(patchwork) #for combining the plots.
#' res$auccplot + res$rocplot

#'
#'
#' @export
#' @importFrom precrec evalmod
#' @importFrom stats dist
#' @importFrom rdist rdist
#' @importFrom reshape2 melt
#' @importFrom grDevices colors
#' @import ggplot2
#' @import ggthemes
#' @import grid
#'

clusteringROCs <- function(partitions, dataset = NULL, distance = NULL, distanceMethod = 'euclidean', showPlot = TRUE, decoratePlot = TRUE, alphaPlot = TRUE, rLine = TRUE, sample = TRUE, sampleSize = 5000) {

  curves <- FP <- TP <- NULL

  if (missing(dataset) && missing(distance)) {
    stop('You need to specify a distance matrix or a dataset.')
  }

  if (missing(partitions)) {
    stop('You need to specify a list of hard partitions - clustering solutions.')
  } else {
    if (class(partitions) != 'list') {
      stop('This function is supposed to be used with a list of partitions. For single partition evaluations use aucc.')
    }
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

  if (length(unique(do.call(rbind,lapply(partitions, length)))) != 1) {
    stop('All partitions should have the same size.')
  }

  if (!is.null(dataset) && is.null(distance)) {
    if (class(partitions) == 'list' && dim(dataset)[1] != length(partitions[[1]])) {
      stop('The number of objects has to the be the same from partition.')
    }
    distance <- c(dist(dataset, diag = FALSE, upper = FALSE, method = distanceMethod))
  } else {
    if (is.null(dataset) && !is.null(distance)) {
      distance <- c(distance)
      if (length(partitions[[1]])*(length(partitions[[1]])-1)/2 != length(distance)) {
        stop('Distance and partitions sizes don\'t match. If n is the number of objects, the required length of distance is n(n-1)/2.')
      }
    }
  }

  auccs <- unlist(lapply(partitions, aucc, dataset=dataset))

  auccsOriginalOrder <- auccs
  kOriginalOrder    <- unlist(lapply(lapply(partitions, unique),length))

  #plot higher auc curver later on, so they show better
  plotOrder  <- order(auccs)
  partitions <- partitions[plotOrder]
  auccs       <- auccs[plotOrder]

  pairwisePartitions <- lapply(lapply(partitions, rdist::rdist, metric = 'hamming'),c)
  distance           <- (distance - min(distance)) / (max(distance) - min(distance))

  #Calculations based on precrec
  performances <- lapply(pairwisePartitions, function(x) precrec::evalmod(scores=1-distance, labels=1-x, mode='basic'))
  df <- lapply(performances, function(x) data.frame(1-x$sp[[1]]$y,x$sn[[1]]$y))

  #Calculations based on ROCR (legacy)
  #predictions  <- lapply(pairwisePartitions, ROCR::prediction, predictions = distance)
  #performances <- lapply(predictions, performance, measure = "tpr", x.measure = "fpr")
  #df <- lapply(performances, function(x) data.frame(FalsePositive=x@x.values[[1]],TruePositive=x@y.values[[1]]))

  if (sample && (sampleSize < dim(df[[1]])[1])) {
    subsample <- c(1,sample(1:dim(df[[1]])[1],sampleSize),dim(df[[1]])[1])
  }

  for (i in 1:length(df)) {
    if (sample && (sampleSize < dim(df[[1]])[1])) {
      df[[i]] <- df[[i]][subsample,]
    }
    colnames(df[[i]]) <- c('FP','TP')
  }

  normalize <- function(x){(x-min(x))/(max(x)-min(x))}

  if (alphaPlot) {
    nauccs <- normalize(auccs)
  } else {
    nauccs <- rep(1,length(auccs))
  }

  melted <- reshape2::melt(df,id.vars = c('TP','FP'))
  colnames(melted) <- c('TP','FP','curves')
  melted$alpha <- nauccs[melted$curves]

  k <- unlist(lapply(lapply(partitions, unique),length))

  kbest <- k[which.max(auccs)]

  curveColors <- ifelse(k==kbest&auccs==max(auccs), 'gold', ifelse(k>kbest, 'red', 'blue'))
  melted$colors <- curveColors[melted$curves]

  colorLines <- c("gold" = "gold", "blue" = "blue", "red" = "red")
  colorLegend <- c("gold" = "k*", "blue" = "k < k*", "red" = "k > k*")

  if (length(unique(melted$colors)) < 3) {
    if (!any(curveColors=='red')) {
      colorLines <- c("gold" = "gold", "blue" = "blue")
      colorLegend <- c("gold" = "k*", "blue" = "k < k*")
    } else {
      if (!any(curveColors=='blue')) {
        colorLines <- c("gold" = "gold", "red" = "red")
        colorLegend <- c("gold" = "k*", "red" = "k > k*")
      }
    }
  }

  pltROCs <- ggplot2::ggplot(melted, aes(curves,x=FP, y=TP, col=colors,alpha = alpha)) +
    geom_line(size = 1,aes(group=curves,x=FP, y=TP)) +
    xlab('FPR') +
    ylab('TPR') +
    scale_color_manual(values=colorLines,labels=colorLegend) +
    theme_minimal() +
    scale_alpha(guide='none')

  linePlot <- data.frame(cbind(auccsOriginalOrder,kOriginalOrder))

  if (length(kOriginalOrder) == length(unique(kOriginalOrder))) {
    pltAUCCs <- ggplot2::ggplot(linePlot, aes(linePlot,x=kOriginalOrder, y=auccsOriginalOrder)) +
      geom_point(alpha = auccsOriginalOrder) +
      geom_line() +
      #geom_smooth(method = 'loess', formula = y~x) +
      geom_point(data=linePlot[which.max(linePlot[,1]),],
                 aes(x=kOriginalOrder, y=auccsOriginalOrder),
                 size=3, colour='gold') +
      geom_point(data=linePlot[which.max(linePlot[,1]),],
                 aes(x=kOriginalOrder, y=auccsOriginalOrder),
                 size=5, pch=21, colour='gold') +
      xlab('k') +
      ylab('AUCC') +
      scale_x_continuous(breaks=sort(unique(kOriginalOrder))) +
      theme_minimal()
  } else {
    pltAUCCs <- ggplot2::ggplot(linePlot, aes(linePlot,x=kOriginalOrder, y=auccsOriginalOrder)) +
      geom_point(col=ifelse(kOriginalOrder<kbest,'blue','red')) +
      #geom_line() +
      geom_smooth(method = 'loess', formula = y~x, span = 0.3, col = 'darkgreen') +
      geom_point(data=linePlot[which.max(linePlot[,1]),],
                 aes(x=kOriginalOrder, y=auccsOriginalOrder),
                 size=3, colour='gold') +
      geom_point(data=linePlot[which.max(linePlot[,1]),],
                 aes(x=kOriginalOrder, y=auccsOriginalOrder),
                 size=5, pch=21, colour='gold') +
      xlab('k values') +
      ylab('AUCC') +
      scale_x_continuous(breaks=sort(unique(kOriginalOrder))) +
      theme_minimal()
  }

  if (decoratePlot) {
    grobROCs <- grid::grobTree(grid::textGrob(paste("Best AUCC: ", formatC(max(auccs), digits = 3,format = "f"),sep = ""),
                              x=0.40,  y=0.10, hjust=0,
                              gp=grid::gpar(col="black", fontsize=13)))

    pltROCs <- pltROCs  + annotation_custom(grobROCs)

    pltAUCCs <- pltAUCCs + geom_text(aes(label=ifelse(auccsOriginalOrder==max(auccsOriginalOrder),formatC(auccsOriginalOrder,digits = 3, format = 'f'),'')),hjust=0,vjust=0)

  }

  if (showPlot) {
    print(pltROCs)
  }

  r <- c()

  r$rocplot   <- pltROCs
  r$auccplot  <- pltAUCCs
  r$auccs     <- auccsOriginalOrder

  class(r) <- 'clusteringROCs'

  r
}
