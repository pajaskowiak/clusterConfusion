#' @title rocPlot
#'
#' @description Plots the ROC Curve of a clustering solution.
#'
#' @param auccResult An object of class "aucc", containing the evaluation, obtained from the aucc function with returnRates=TRUE.
#' @param showPlot Determines if the plot should be displayed or just generated and returned for later use.
#' @param rLine Whether to plot or not the random line that accounts for a random clustering solution.
#' @param sample Defines if all FPR and TPR should be used to generate the plot (slower), or just a subsample (faster).
#' @param sampleSize If sample=TRUE, defines the number of FPR and TPR pairs that will be used to generate the plot.
#' @param color The color of the curve.
#' @param size The size/width of the line/curve.
#' @param lineType The type of line used in the curve (solid, dashed, ...).
#' @param addToPlot If NULL just plots the curve. If you pass a previous plot obtained with the function, you can overlay another curve (or as many as you want, with subsequent calls).
#'
#' @return A 'ggplot' object with the corresponding ROC Curve. Note that you can manipulate and customize the plot with subsequent ggplot function calls.
#'
#' @details
#'
#' Plots the ROC curve for a clustering solution. Note that this function does no evaluate the clustering result by itself, it only plots the ROC Curve. In order to use this function you first need to call the [aucc()] with the parameter returnRates=TRUE.
#'
#' Function calls are flexible. You may choose whether to show the plot or just get it for later use. You may choose whether to add the random line (the expected "curve" from a random clustering solution, which has an AUCC of 0.5) from (0,0) to (1,1). Other than that, you can select the curve color and its line width. Given that the function returns a 'ggplot' object, you can customize the plot to your taste after the function call too.
#'
#' It is worth mentioning that when called with returnRates=TRUE, The [aucc()] function will return the AUCC value alongside _all_ FPR and TPR values. In case of large datasets, the number of FPR and TPR pairs may translate into high computational costs for obtaining the ROC Plot (with the function taking too long to plotting the Curve). In order to reduce the wating time/computational cost, the sample parameter is used by default (sample=TRUE). This means that only a subsample of FPR and TPR values will be used in the ROC Plot. Note that this has no implication at all in the AUCC value itself, it only affects the ROC Plot aesthetics. You may define different sample sizes (sampleSize) and see the one that fits you best (balance between time and plot quality). For quick plots, set a small sampleSize. For high fidelity plots, set sample=FALSE and, if you have a large dataset, be patient (or maybe, go for a cup of coffee).
#'
#'Finally, regarding the addToPlot argument. By default it is set as NULL, which indicates that a new ggplot will be generated and returned by the function call. You may, however, add subsequent ROC Curves to the first plot, for other clustering solutions, for instance. In this case, instead of NULL, just pass the plot to which you want to add another curve (see the example section for more).
#'
#' @references
#'
#' Jaskowiak, P.A., Costa, I.G. & Campello, R.J.G.B. The area under the ROC curve as a measure of clustering quality. Data Mining and Knowledge Discovery 36, 1219–1245 (2022).  [https://doi.org/10.1007/s10618-022-00829-0](https://doi.org/10.1007/s10618-022-00829-0).
#'
#' @examples
#'
#' require(ggplot2)
#'
#' #Just picking some nice colors
#' colorsScatter <- c("#0066FF","#FF0099")
#' colorROCS     <- heat.colors(9)
#'
#' data     <- c()
#' distance <- c()
#'
#' #Just place all the Gaussian data in a list so we can iterate
#' #All these datasets are loaded within the package
#'
#' data[[1]] <- gData1
#' data[[2]] <- gData2
#' data[[3]] <- gData3
#' data[[4]] <- gData4
#' data[[5]] <- gData5
#' data[[6]] <- gData6
#'
#' rPlot <- NULL     #We will store the ROC Plots here
#' dPlot <- list()   #And the data scatter plots here
#'
#' #So now we iterate over the six datasets...
#' for (d in 1:6) {
#'
#'   #Store the results of the AUCC considering the ground truth and the distance (aucc object)
#'   r <- aucc(partition = as.integer(data[[d]][,3]), dataset = data[[d]][,-3], returnRates = TRUE)
#'
#'   #For the first dataset, we make a new plot (addToPlot=NULL). For the following plots
#'   #we call the same function, but we add the new curves (addToPlot = rPlot). We choose
#'   #not to show each plot (showPlot=FALSE). We'll compose all togheter in the end
#'
#'   if (d==1) {
#'     rPlot <- rocPlot(r, addToPlot = NULL, col=colorROCS[d],
#'                      lineType = d, sampleSize = 3000, showPlot = FALSE)
#'   } else {
#'     rPlot <- rocPlot(r, addToPlot = rPlot, col=colorROCS[d],
#'                      lineType = d, showPlot = FALSE)
#'   }
#'
#'   #We make a separate scatter plot for each one of the datasets
#'   dPlot[[d]] <- ggplot(data[[d]],aes(x=x,y=y,color=class)) +
#'     geom_point() +
#'     theme_minimal() +
#'     scale_color_manual(values=colorsScatter) +
#'     theme(legend.position="none") +
#'     xlab(element_blank()) +
#'     ylab(element_blank())
#'
#' }
#'
#' #We now put everything together. This is the same plot presented in Fig. 5 of the reference paper.
#'
#' require(patchwork) #we need this library in order to compose the ggplots in the following line
#' print((((dPlot[[1]] | dPlot[[2]] | dPlot[[3]])/
#'          (dPlot[[4]] | dPlot[[5]] | dPlot[[6]]))|rPlot) + plot_layout(widths = c(4,2.5)))
#'
#'
#' @export
#' @import ggplot2
#'
rocPlot <- function(auccResult, showPlot = TRUE, rLine = TRUE, sample = TRUE, sampleSize = 500, color = 'black', size = 1, lineType = 1, addToPlot = NULL) {

  FPR <- TPR <- NULL

  if (class(auccResult) != 'aucc') {
    stop('Object auccResults has to come from function aucc, with returnRates = TRUE')
  }

  if (sample && sampleSize < length(auccResult$fpr)) {
    subsample <- c(1,sample(1:length(auccResult$fpr),sampleSize),length(auccResult$fpr))
    rDF <- data.frame(cbind(auccResult$fpr[subsample],auccResult$tpr[subsample]))
  } else {
    rDF <- data.frame(cbind(auccResult$fpr,auccResult$tpr))
  }

  colnames(rDF) <- c('FPR','TPR')

  if (is.null(addToPlot)) {
    rPlot <- ggplot2::ggplot(rDF,aes(x=FPR,y=TPR)) +
      ggplot2::geom_line(size = size, col = color, linetype = lineType) +
      ggplot2::theme_minimal()
  } else {
    rPlot <- addToPlot +
      ggplot2::geom_line(data=rDF,ggplot2::aes(x=FPR,y=TPR), col = color, size = size, linetype = lineType) +
      ggplot2::theme_minimal()
  }

  if (rLine) {
    rPlot <- rPlot +
      ggplot2::geom_segment(ggplot2::aes(x = 0, y = 0, xend = 1, yend = 1), color = 'darkgrey', linetype=2)
  }

  if (showPlot) {
    print(rPlot)
  }

  rPlot
}
