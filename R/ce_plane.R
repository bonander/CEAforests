#' @title Plot an individual or population-level cost-effectiveness plane using a CEA forest.
#' @description \code{ce_plane} Plots the results from CEA forests, optionally with uncertainty estimates from an NMB forest.
#'
#' @param forest A trained CEA forest.
#' @param WTP The willingness to pay for an incremental increase in the outcome. Is set to zero if NULL.
#' @param alpha The desired confidence level, defaults to 0.05. Only used when nmb is supplied.
#' @param type A string (either "individual" or "population") specifying which type of cost-effectiveness plane to plot. Defaults to individual.
#' @param R Number of bootstrap replicates (only used for population-level planes).
#' @param subset A specified subpopulation to compute the plane for (only used for population-level planes). Defaults to NULL.
#' @param WTP The willingness to pay per one unit increase in the outcome. Only used for population-level plots. Re-estimate the forest to change WTP for individual-level plots. Uses the WTP from the forest object by default.
#' @return A ggplot2 object.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @export

ce_plane = function(forest, alpha=0.05, type="individual", R=999, subset=NULL, WTP=NULL) {

  if (isTRUE(any(class(forest) %in% c("CEAforests")))==FALSE) {stop("Unrecognized class. Please supply a CEAforests object.")}

  obj = forest

  if (isTRUE(type=="individual")) {
  WTP = obj[["WTP"]]
  preds = predict(obj)
  dy = preds$predicted.delta_y
  lowerx = dy-sqrt(preds$variance.delta_y)*qnorm(1-alpha/2)
  upperx = dy+sqrt(preds$variance.delta_y)*qnorm(1-alpha/2)
  dc = preds$predicted.delta_cost
  lowery = dc-sqrt(preds$variance.delta_cost)*qnorm(1-alpha/2)
  uppery = dc+sqrt(preds$variance.delta_cost)*qnorm(1-alpha/2)
  df = as.data.frame(cbind(dc, dy))
  p = ggplot2::ggplot(data=df,ggplot2::aes(y=dc, x=dy))
  p = p + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))


    nmb.point = preds$predicted.nmb

    nmb.lower = nmb.point-sqrt(preds$variance.nmb)*qnorm(1-alpha/2)
    nmb.upper = nmb.point+sqrt(preds$variance.nmb)*qnorm(1-alpha/2)

    certain=as.numeric(rep(0, length(nmb.point)))
    certain[nmb.point>0 & nmb.lower>0 & nmb.upper>0] = 1
    certain[nmb.point<0 & nmb.lower<0 & nmb.upper<0] = -1

    if (isTRUE(length(unique(certain))==3)){
    certain = factor(certain, labels=c(paste("Not cost-effective\n(n = ", sum(as.numeric(certain==-1)),", ",  round(sum(as.numeric(certain==-1))/length(certain)*100, 1),"%)", sep=""),
                                       paste("Undetermined\n(n = ", sum(as.numeric(certain==0)),", ",  round(sum(as.numeric(certain==0))/length(certain)*100, 1),"%)", sep=""),
                                       paste("Cost-effective\n(n = ", sum(as.numeric(certain==1)),", ",  round(sum(as.numeric(certain==1))/length(certain)*100, 1),"%)", sep="")))
    scalecol = c("tomato2", "beige", "palegreen3")
    } else if (isTRUE(identical(base::sort(unique(certain)), c(-1,0)))) {
      certain = factor(certain, labels=c(paste("Not cost-effective\n(n = ", sum(as.numeric(certain==-1)),", ",  round(sum(as.numeric(certain==-1))/length(certain)*100, 1),"%)", sep=""),
                                         paste("Undetermined\n(n = ", sum(as.numeric(certain==0)),", ",  round(sum(as.numeric(certain==0))/length(certain)*100, 1),"%)", sep="")))
      scalecol = c("tomato2", "beige")
    } else if (isTRUE(identical(base::sort(unique(certain)), c(0,1)))) {
      certain = factor(certain, labels=c(paste("Cost-effective\n(n = ", sum(as.numeric(certain==1)),", ",  round(sum(as.numeric(certain==1))/length(certain)*100, 1),"%)", sep=""),
                                         paste("Undetermined\n(n = ", sum(as.numeric(certain==0)),", ",  round(sum(as.numeric(certain==0))/length(certain)*100, 1),"%)", sep="")))
      scalecol = c("beige", "palegreen3")
    } else if (isTRUE(identical(unique(certain), c(0)))) {
      certain = factor(certain, labels=c(paste("Undetermined\n(n = ", sum(as.numeric(certain==0)),", ",  round(sum(as.numeric(certain==0))/length(certain)*100, 1),"%)", sep="")))
      scalecol = c("beige")
    } else if (isTRUE(identical(unique(certain), c(-1)))) {
      certain = factor(certain, labels=c(paste("Not cost-effective\n(n = ", sum(as.numeric(certain==-1)),", ",  round(sum(as.numeric(certain==-1))/length(certain)*100, 1),"%)", sep="")))
      scalecol = c("tomato2")
    } else if (isTRUE(identical(unique(certain), c(1)))) {
      certain = factor(certain, labels=c(paste("Cost-effective\n(n = ", sum(as.numeric(certain==1)),", ",  round(sum(as.numeric(certain==1))/length(certain)*100, 1),"%)", sep="")))
      scalecol = c("palegreen3")
    }
    p = p + ggplot2::geom_point(ggplot2::aes(y=dc, x=dy, color=certain)) +
      ggplot2::scale_color_manual(values=scalecol) +
      ggplot2::labs(color=paste("Certainty (WTP = ", WTP, ", Certainty level: ", round((1-alpha)*100,2),"%):", sep="")) +
      ggplot2::theme(legend.position="top") +
      ggplot2::geom_abline(intercept=0, slope=ifelse(is.null(WTP), 0, WTP)) +
      ggplot2::geom_hline(yintercept=0) +
      ggplot2::geom_vline(xintercept=0) +
      ggplot2::ylab(expression(paste(Delta,"Cost")[i])) +
      ggplot2::xlab(expression(paste(Delta,"Y")[i]))

  } else if (isTRUE(type=="population")) {
    if (is.null(WTP)) {WTP = obj[["WTP"]]} else {WTP=WTP}
    dy=cate.prepare(forest[["outcome.forest"]], subset=subset)
    dc=cate.prepare(forest[["cost.forest"]], subset=subset)
    boot.estimates = boot.dr_scores(dc, dy, R=R, WTP=WTP, alpha=alpha)[[2]]
    dc = boot.estimates[,2]
    dy = boot.estimates[,1]
    df = as.data.frame(cbind(dy,dc))
    below_threshold=mean(as.numeric((dy*WTP-dc)>0))

    p = ggplot2::ggplot(data=df,ggplot2::aes(y=dc, x=dy)) +
      ggplot2::geom_point(color="gray") +
      ggplot2::geom_abline(intercept=0, slope=ifelse(is.null(WTP), 0, WTP)) +
      ggplot2::geom_hline(yintercept=0) +
      ggplot2::geom_vline(xintercept=0) +
      ggplot2::ylab(expression(paste(Delta,"Cost"))) +
      ggplot2::xlab(expression(paste(Delta,"Y"))) +
      ggplot2::ggtitle(paste("Cost-effectiveness plane (WTP: ", WTP,"; ", "% Cost-effective: ", round(below_threshold*100,1),")", sep=""))

    p = p + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    p = p + ggplot2::stat_ellipse(linetype=2, level=1-alpha)

  } else (stop("Unrecognized plot type. Select either individual or population."))

  return(p)
}
