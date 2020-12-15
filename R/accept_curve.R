#' @title (Sub-)population acceptability curves.
#' @description \code{accept_curve} Provides acceptability curves for the entire sample or for a specified subpopulation.
#'
#' @param forest A trained CEA forest.
#' @param from Lowest WTP to consider.
#' @param to Highest WTP to consider.
#' @param length.out Number of increments from lowest to highest WTP.
#' @param subset A specified subpopulation. See grf::average_treatment_effect for details.
#' @param compare A logical statement. If TRUE, the subpopulation curve is compared to the sample average curve. Defaults to FALSE.
#' @param robust.se Whether or not robust (sandwich) standard errors should be used to compute the acceptance probability. Defaults to FALSE.
#'
#' @return A matrix with acceptance probabilities given the data, subpopulation and different choices of WTP.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @import grf
#' @import stats
#' @import graphics
#' @export

accept_curve = function(forest, from, to, length.out, subset=NULL, compare=FALSE, robust.se=FALSE, labels=NULL) {


  if (isTRUE(any(class(forest) %in% c("CEAforests")))==FALSE) {stop("Unrecognized class. Please supply a CEAforests object.")}

  if (("instrumental_forest" %in% class(forest[["outcome.forest"]]))) {
  fobj = forest[["outcome.forest"]]
  compliance.forest.y <- grf::causal_forest(fobj$X.orig, Y = fobj$W.orig,
                                            W = fobj$Z.orig, Y.hat = fobj$W.hat, W.hat = fobj$Z.hat)
  compliance.score.y <- predict(compliance.forest.y)$predictions
  }
  if (("instrumental_forest" %in% class(forest[["cost.forest"]]))) {
    fobj = forest[["cost.forest"]]
    compliance.forest.c <- grf::causal_forest(fobj$X.orig, Y = fobj$W.orig,
                                              W = fobj$Z.orig, Y.hat = fobj$W.hat, W.hat = fobj$Z.hat)
    compliance.score.c <- predict(compliance.forest.c)$predictions
  }
  if (("instrumental_forest" %in% class(forest[["cost.forest"]]))) {
    compliance.scores = cbind(compliance.score.y, compliance.score.c)
  } else {compliance.scores = NULL}


  WTP.seq = base::seq(from=from, to=to, length.out=length.out)

  reslist = list()
  for (i in 1:length(WTP.seq)) {
    WTP = WTP.seq[i]
    reslist[i]=avg_effects(forest, WTP=WTP, subset=subset, robust.se=robust.se, compliance.scores=compliance.scores, icer.ci=FALSE)[3,5]
  }
  res1 = do.call("rbind", reslist)
  if (isTRUE(compare==TRUE && !is.null(subset)==TRUE)) {
    reslist2 = reslist3 = list()
    for (i in 1:length(WTP.seq)) {
      WTP = WTP.seq[i]
      reslist2[i]=avg_effects(forest, WTP=WTP, subset=NULL, robust.se=robust.se, compliance.scores=compliance.scores, icer.ci=FALSE)[3,5]
      reslist3[i]=avg_effects(forest, WTP=WTP, subset=!subset, robust.se=robust.se, compliance.scores=compliance.scores, icer.ci=FALSE)[3,5]
    }
    res2 = do.call("rbind", reslist2)
    res3 = do.call("rbind", reslist3)
    ret=cbind(WTP=WTP.seq, accept.prob=res1,accept.prob.average=res2,accept.prob.not.in.subset=res3)
    colnames(ret) = c("WTP", "accept.prob", "accept.prob.average","accept.prob.not.in.subset")
    class(ret) = c("accept_compare","accept_curve")
    return(plot.accept_curve(ret, labels=labels))
  } else {ret=cbind(WTP=WTP.seq, accept.prob=res1)
          colnames(ret)=c("WTP", "accept.prob")
          class(ret) = c("accept_main", "accept_curve")
          return(plot(ret))}
}

#' @title Plotting function for (sub)-population acceptability curves.
#' @description \code{plot.accept_curve} Plots acceptability curves using ggplot2.
#'
#' @param x An acceptability curve object.
#' @param labels A character vector of labels for the selected subset, the entire sample and the unselected subset, defaults to c("Selected subgroup", "Full sample", "Not in selected subgroup").
#' @param ... Additional options (currently ignored).
#' @keywords internal
#' @return A ggplot.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @import ggplot2
#' @method plot accept_curve
#' @export

plot.accept_curve = function(x, labels=NULL, ...) {
  curve = x
  if (isTRUE(base::identical(class(curve), c("accept_main","accept_curve")))) {
    df = as.data.frame(curve[,c(1,2)])

    p = ggplot(df, aes(y=accept.prob, x=WTP)) + theme_bw() +
      xlab("Willingness to pay") + ylab("Pr(Cost-Effective)") +
      geom_line(size=1) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size=15))

  }
  if (isTRUE(base::identical(class(curve), c("accept_compare","accept_curve")))) {
    df1 = as.data.frame(curve[,c(1,2)])
    df2 = as.data.frame(curve[,c(1,3)])
    df3 = as.data.frame(curve[,c(1,4)])
    if (is.null(labels)) {
    df1$Group = "Selected subgroup"
    df2$Group = "Full sample"
    df3$Group = "Not in selected subgroup"
    } else {

      df1$Group = labels[1]
      df2$Group = labels[2]
      df3$Group = labels[3]

    }
    colnames(df1) = colnames(df2) = colnames(df3) = c("WTP", "accept.prob", "Group")
    df = rbind(df1, df2, df3)
    p = ggplot(df, aes(y=accept.prob, x=WTP, linetype=Group)) + theme_bw() +
      xlab("Willingness to pay") + ylab("Pr(Cost-Effective)") +
      geom_line(size=1) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size=15),
            legend.title=element_blank(),
            legend.position="bottom") + scale_linetype_manual(values=c("dashed", "dotted", "solid"))
  }
  return(p)
}
