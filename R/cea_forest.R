#' @title Estimate causal forests for outcomes, costs and net monetary benefits.
#' @description \code{cea_forest} Runs causal forests for outcomes, costs and net monetary benefits given a specified willingness to pay (a wrapper for grf::causal_forest).
#'
#' @param Y The outcome vector.
#' @param C The cost vector.
#' @param X The covariate matrix.
#' @param W The treatment vector.
#' @param Z An instrumental variable. (Optional)
#' @param WTP Willingness to pay per one-unit increase in the outcome. Defaults to 1.
#' @param W.hat Pre-fitted propensity scores for treatment (W). If NULL, the algorithm fits a regression forest to estimate W.hat.
#' @param tune.parameters Which hyperparameters to tune. Defaults to "all". See grf::causal_forest for other options. Option "none" uses default settings for all parameters.
#' @param num.trees The number of trees in each forest. Defaults to 5000. Can (and probably should) be set to a higher number to reduce Monte Carlo errors.
#' @param ... Other options to be passed to grf::causal_forest() or grf::instrumenal_forest() if instrument is supplied.
#'
#' @references Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random forests. The Annals of Statistics, 47(2), 1148-1178.
#'
#' @return Returns a list containing two causal forest objects (one for the outcome and one for costs). If an instrument is supplied, the code returns two instrumental forest objects.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @import grf
#' @export
cea_forest = function(Y, C, X, W, Z=NULL, WTP=NULL, W.hat=NULL, tune.parameters="all", num.trees=5000, ...) {
  if (is.null(WTP)) {
    message("No willingness to pay (WTP) per one-unit increase in Y supplied. Setting WTP to 1.")
    WTP = 1
  }
  if (length(unique(W))>2) {
    message("You seem to have supplied a non-binary treatment variable. The causal forest algorithm will still run, but other functions in the CEAforest package may not work as intended. Use with caution.")
  }
  if (length(unique(Z))>2) {
    message("You seem to have supplied a non-binary instrument. The instrumental forest algorithm will still run, but other functions in the CEAforest package may not work as intended. Use with caution.")
  }
  if (base::exists("clusters")==TRUE) {
    message("You seem to have supplied a cluster variable. The forest algorithm will still run, but other functions in the CEAforest package have not yet been extended to clustered data. Use with caution.")
  }
  if (is.null(W.hat)) {
    #Unless custom propensity scores are provided, pre-fit a regression forest for W to speed up algorithm.
    w_forest = grf::regression_forest(Y=W, X=X, tune.parameters=tune.parameters, num.trees=num.trees, ...)
    W.hat = predict(w_forest)$predictions
  }
  if (is.null(Z)) {
  y_forest = grf::causal_forest(X=X, Y=Y, W=W, W.hat=W.hat, tune.parameters=tune.parameters, num.trees=num.trees, ...)
  c_forest = grf::causal_forest(X=X, Y=C, W=W, W.hat=W.hat, tune.parameters=tune.parameters, num.trees=num.trees, ...)
  nmb_forest = grf::causal_forest(X=X, Y=Y*WTP-C, W=W, W.hat=W.hat, tune.parameters=tune.parameters, num.trees=num.trees, ...)
  forest = list()
  forest[["outcome.forest"]] = y_forest
  forest[["cost.forest"]] = c_forest
  forest[["nmb.forest"]] = nmb_forest
  forest[["WTP"]] = WTP
  class(forest) = c("cea_forest", "CEAforests")
  } else {
  y_forest = grf::instrumental_forest(X=X, Y=Y, W=W, Z=Z, W.hat=W.hat, tune.parameters=tune.parameters, num.trees=num.trees, ...)
  c_forest = grf::instrumental_forest(X=X, Y=C, W=W, Z=Z, W.hat=W.hat, tune.parameters=tune.parameters, num.trees=num.trees, ...)
  nmb_forest = grf::instrumental_forest(X=X, Y=Y*WTP-C, W=W, Z=Z, W.hat=W.hat, tune.parameters=tune.parameters, num.trees=num.trees, ...)
  forest = list()
  forest[["outcome.forest"]] = y_forest
  forest[["cost.forest"]] = c_forest
  forest[["nmb.forest"]] = nmb_forest
  forest[["WTP"]] = WTP
  class(forest) = c("cea_forest_instrumental", "CEAforests")
  }
  return(forest)
}

#' @title Predict with a CEA forest.
#' @description \code{predict.CEAforests} Gets estimates of conditional incremental outcomes and costs given x using a cea_forest object (a wrapper for grf::predict.causal_forest).
#'
#' @param object The trained CEA forest.
#' @param ... Other options to be passed to grf::predict.causal_forest() or grf::predict.instrumental_forest(). See grf documentation for additional information.
#'
#' @return A matrix of predictions of conditional average treatment effects for the outcome and costs, along with variance estimates.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @import grf
#' @import stats
#' @export
predict.CEAforests = function(object, ...) {
  obj = object
    yp = predict(obj[["outcome.forest"]], estimate.variance=TRUE, ...)
    predicted.delta_y = yp$predictions; variance.delta_y = yp$variance.estimates
    cp = predict(obj[["cost.forest"]], estimate.variance=TRUE, ...)
    predicted.delta_cost = cp$predictions; variance.delta_cost = cp$variance.estimates
    nmb = predict(obj[["nmb.forest"]], estimate.variance=TRUE, ...)
    predicted.nmb = nmb$predictions; variance.nmb = nmb$variance.estimates
    return(as.data.frame(cbind(predicted.delta_y, variance.delta_y,
                               predicted.delta_cost, variance.delta_cost,
                               predicted.nmb, variance.nmb)))
}

#' @title Plotting function for CEA forests.
#' @description Provides histograms, scatter plots or (partial) effects plots to assess heterogeneity with respect to a covariate after a CEAforest.
#' @param forest A trained CEA forest.
#' @param which.y A string or string vector naming which outcomes to plot (any combination of "outcome", "cost" and "nmb" is acceptable). Defaults to "all", which produces three plots in one-row, three-column grid.
#' @param which.x A column number or string naming a single variable from the X matrix in the CEAforests object. If null, the function outputs histograms of the out-of-bag estimates.
#' @param conditional Whether or not to keep all other variables in the X matrix constant at their mean in bivariate plots (defaults to FALSE). Ignored if which.x is null.
#' @param smooth Whether or not to plot a semi-parametric smooth function fit to doubly robust scores for tau(x) (via the mgcv package) instead of out-of-bag estimates (unconditional) or non-parametric predictions (conditional). Defaults to FALSE.
#' @param ci.level The desired confidence level for confidence intervals (used when applicable). Defaults to 0.95.
#' @param x.range A two-element numeric vector that controls the range of the x-axis. Defaults to min(which.x) and max(which.x).
#' @param length.out The length of the sequence of X values to be plotted (used when conditional=TRUE). Defaults to 100.
#' @param xlab Label for the X axis. Defaults to the column name for the focal X variable, or to "X" if no column name is available.
#' @param labels Labels to be passed to cowplot::plot_grid when more than one outcome type (which.y) is supplied. Defaults to "AUTO", which labels the plots using letters.
#' @param ... Additional arguments to be passed to the gam function in mgcv. Ignored if smooth=FALSE.
#'
#' @return A ggplot via ggplot2 or grid of ggplots via the cowplot package.
#' @examples
#' \dontrun{
#' To be added...
#' }
#' @import stats
#' @import grf
#' @import ggplot2
#' @export
plot.CEAforests = function(forest, which.y="all", which.x=NULL, conditional=FALSE, smooth=FALSE, ci.level=0.95, x.range=NULL, length.out=100, xlab="X", labels="AUTO", ...) {

  preds = predict(forest)
  dr.taus = as.data.frame(debias_effects(forest))


  if (isTRUE(any(c("cost", "outcome", "nmb", "all") %in% which.y))==FALSE) {
    stop("which.y appears to  be misspecified. See help file for valid options.")
  }

  if (isTRUE(which.y=="all")) {which.y=c("outcome", "cost", "nmb")}

    if (is.null(which.x)) {#Plot histogram(s)

      if (isTRUE(any(c("outcome") %in% which.y))) {
      yplot = ggplot(preds, ggplot2::aes(x=predicted.delta_y, ..count../sum(..count..))) +
                    ggplot2::geom_histogram(binwidth=2*stats::IQR(preds$predicted.delta_y)/(length(preds$predicted.delta_y)^(1/3)), boundary = 0, color="black", fill="gray50") +
                    ggplot2::ylab("Density") + ggplot2::xlab(expression(paste(Delta,"Y"))) +
                    ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                             panel.grid.minor = ggplot2::element_blank(),
                                             axis.line = ggplot2::element_line(colour = "black"),
                                             text = ggplot2::element_text(size=12))
      }
      if (isTRUE(any(c("cost") %in% which.y))) {
      cplot = ggplot(preds, ggplot2::aes(x=predicted.delta_cost, ..count../sum(..count..))) +
                    ggplot2::geom_histogram(binwidth=2*stats::IQR(preds$predicted.delta_cost)/(length(preds$predicted.delta_cost)^(1/3)), boundary = 0, color="black", fill="gray50") +
                    ggplot2::ylab("Density") + ggplot2::xlab(expression(paste(Delta,"Cost"))) +
                    ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                             panel.grid.minor = ggplot2::element_blank(),
                                             axis.line = ggplot2::element_line(colour = "black"),
                                             text = ggplot2::element_text(size=12))
      }
      if (isTRUE(any(c("nmb") %in% which.y))) {
      nplot = ggplot(preds, ggplot2::aes(x=predicted.nmb, ..count../sum(..count..))) +
                    ggplot2::geom_histogram(binwidth=2*stats::IQR(preds$predicted.nmb)/(length(preds$predicted.nmb)^(1/3)), boundary = 0, color="black", fill="gray50") +
                    ggplot2::ylab("Density") + ggplot2::xlab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) +
                    ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                             panel.grid.minor = ggplot2::element_blank(),
                                             axis.line = ggplot2::element_line(colour = "black"),
                                             text = ggplot2::element_text(size=12))
      }

      if (isTRUE(setequal(which.y, c("outcome", "cost", "nmb")))) {
        p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
      else if (isTRUE("all" %in% which.y)) {
        p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
      else if (isTRUE(setequal(which.y, c("outcome", "cost")))) {
        p = cowplot::plot_grid(yplot, cplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
      else if (isTRUE(setequal(which.y, c("outcome", "nmb")))) {
        p = cowplot::plot_grid(yplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
      else if (isTRUE(setequal(which.y, c("cost", "nmb")))) {
        p = cowplot::plot_grid(cplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
      else if (isTRUE(setequal(which.y, c("outcome")))) {
        p = yplot}
      else if (isTRUE(setequal(which.y, c("cost")))) {
        p = cplot}
      else if (isTRUE(setequal(which.y, c("nmb")))) {
        p = nplot}

    } else { # Bivariate plots

      #Prepare inputs
      alpha = 1-ci.level
      Xmat = forest[["outcome.forest"]]$X.orig
      x = Xmat[,which.x]
      if (isTRUE(all(class(which.x)=="character" & xlab=="X"))) {xlab = which.x}
      if (isTRUE(class(which.x)=="character")) {which.x = which(colnames(Xmat)==which.x)}
      if (is.null(x.range)) {x.range=c(min(x), max(x))}

      if (isTRUE(conditional)) {#Keep all other variables in X-matrix constant at their mean

        if (isTRUE(smooth)) {#Use mgcv to produce a spline plot with CIs
          if (isTRUE("mgcv" %in% rownames(installed.packages())==FALSE)) {
            stop("The mgcv package must be installed to plot smooth functions.")
          }

          xdf = as.data.frame(Xmat)
          xvars = colnames(Xmat)
          Xmeans <- apply(xdf, 2, mean)

          is.bin = apply(xdf, 2, function(x) isTRUE(length(unique(x))==2)) #Check for binary covariates before fitting model
          is.factor = apply(xdf, 2, function(x) isTRUE(class(x)=="factor" | class(x)=="character"))

          focal.x = colnames(xdf)[which.x]
          avars = colnames(xdf)[is.bin|is.factor]
          cvars = colnames(xdf)[!(is.bin|is.factor)]

          splineterms = paste("s(", cvars,", k=-1)", sep="")
          b = paste(c(avars, splineterms), collapse="+")

          #This should be conditioned on continuous variables; otherwise predict for each category in the fVar.

          Xmeans <- apply(Xmat, 2, mean)
          if (isTRUE(focal.x %in% cvars)) {
          X.test = matrix(rep(Xmeans, length.out), length.out, ncol(Xmat), byrow=T)
          X.test[,which.x] = seq(x.range[1], x.range[2], length.out=length.out)
          } else {
            X.test = matrix(rep(Xmeans, length(unique(Xmat[,which.x]))), length(unique(Xmat[,which.x])), ncol(Xmat), byrow=T)
            X.test[,which.x] = unique(Xmat[,which.x])
          }
          colnames(X.test) = colnames(xdf)
          X.test = as.data.frame(X.test)
          x.pred = X.test[,which.x]

          #Outcomes

          if (isTRUE(any(c("outcome") %in% which.y))) {

          y.score = dr.taus$debiased.delta_outcome
          pdf = as.data.frame(cbind(y.score, xdf))
          form = as.formula(paste("y.score~",b,sep=""))
          sgam = mgcv::gam(form, data=pdf, ...)

          gam.preds = predict(sgam, newdata=X.test, se.fit=TRUE)
          tau.fit = gam.preds$fit
          tau.se = gam.preds$se.fit
          lower.tau = tau.fit-tau.se*qnorm(1-alpha/2)
          upper.tau = tau.fit+tau.se*qnorm(1-alpha/2)
          df.y = as.data.frame(cbind(tau.fit, lower.tau, upper.tau, x.pred))

          if (isTRUE(focal.x %in% cvars)) {#If continuous, plot spline

          yplot = ggplot2::ggplot(data=df.y,ggplot2::aes(y=tau.fit, x=x.pred)) +
            geom_line(size=1) + geom_ribbon(aes(ymin=lower.tau, ymax=upper.tau), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
            ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                 panel.grid.minor = ggplot2::element_blank(),
                                                 axis.line = ggplot2::element_line(colour = "black"),
                                                 text = ggplot2::element_text(size=12)) +
            ggplot2::ylab(expression(paste(Delta,"Y"))) + ggplot2::xlab(xlab)

          } else { #Else plot point with SE
            yplot = ggplot2::ggplot(data=df.y,ggplot2::aes(y=tau.fit, x=as.factor(x.pred))) +
              geom_pointrange(aes(ymin=lower.tau, ymax=upper.tau), size=1) + geom_hline(yintercept=0, linetype="dashed") +
              ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                   panel.grid.minor = ggplot2::element_blank(),
                                                   axis.line = ggplot2::element_line(colour = "black"),
                                                   text = ggplot2::element_text(size=12)) +
              ggplot2::ylab(expression(paste(Delta,"Y"))) + ggplot2::xlab(xlab)
          }
          }

          #Costs

          if (isTRUE(any(c("cost") %in% which.y))) {

          y.score = dr.taus$debiased.delta_cost
          pdf = as.data.frame(cbind(y.score, xdf))
          form = as.formula(paste("y.score~",b,sep=""))
          sgam = mgcv::gam(form, data=pdf, ...)

          gam.preds = predict(sgam, newdata=X.test, se.fit=TRUE)
          tau.fit = gam.preds$fit
          tau.se = gam.preds$se.fit
          lower.tau = tau.fit-tau.se*qnorm(1-alpha/2)
          upper.tau = tau.fit+tau.se*qnorm(1-alpha/2)

          df.c = as.data.frame(cbind(tau.fit, lower.tau, upper.tau, x.pred))

          if (isTRUE(focal.x %in% cvars)) {#If continuous, plot spline

            cplot = ggplot2::ggplot(data=df.c,ggplot2::aes(y=tau.fit, x=x.pred)) +
              geom_line(size=1) + geom_ribbon(aes(ymin=lower.tau, ymax=upper.tau), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
              ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                   panel.grid.minor = ggplot2::element_blank(),
                                                   axis.line = ggplot2::element_line(colour = "black"),
                                                   text = ggplot2::element_text(size=12)) +
              ggplot2::ylab(expression(paste(Delta,"Cost"))) + ggplot2::xlab(xlab)

          } else { #Else plot point with SE
            cplot = ggplot2::ggplot(data=df.c,ggplot2::aes(y=tau.fit, x=as.factor(x.pred))) +
              geom_pointrange(aes(ymin=lower.tau, ymax=upper.tau), size=1) + geom_hline(yintercept=0, linetype="dashed") +
              ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                   panel.grid.minor = ggplot2::element_blank(),
                                                   axis.line = ggplot2::element_line(colour = "black"),
                                                   text = ggplot2::element_text(size=12)) +
              ggplot2::ylab(expression(paste(Delta,"Cost"))) + ggplot2::xlab(xlab)
          }
          }

          #NMB

          if (isTRUE(any(c("nmb") %in% which.y))) {

          y.score = dr.taus$debiased.nmb
          pdf = as.data.frame(cbind(y.score, xdf))
          form = as.formula(paste("y.score~",b,sep=""))
          sgam = mgcv::gam(form, data=pdf, ...)

          gam.preds = predict(sgam, newdata=X.test, se.fit=TRUE)
          tau.fit = gam.preds$fit
          tau.se = gam.preds$se.fit
          lower.tau = tau.fit-tau.se*qnorm(1-alpha/2)
          upper.tau = tau.fit+tau.se*qnorm(1-alpha/2)

          df.n = as.data.frame(cbind(tau.fit, lower.tau, upper.tau, x.pred))

          if (isTRUE(focal.x %in% cvars)) {#If continuous, plot spline

            nplot = ggplot2::ggplot(data=df.n,ggplot2::aes(y=tau.fit, x=x.pred)) +
              geom_line(size=1) + geom_ribbon(aes(ymin=lower.tau, ymax=upper.tau), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
              ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                   panel.grid.minor = ggplot2::element_blank(),
                                                   axis.line = ggplot2::element_line(colour = "black"),
                                                   text = ggplot2::element_text(size=12)) +
              ggplot2::ylab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) + ggplot2::xlab(xlab)

          } else { #Else plot point with SE
            nplot = ggplot2::ggplot(data=df.n,ggplot2::aes(y=tau.fit, x=as.factor(x.pred))) +
              geom_pointrange(aes(ymin=lower.tau, ymax=upper.tau), size=1) + geom_hline(yintercept=0, linetype="dashed") +
              ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                   panel.grid.minor = ggplot2::element_blank(),
                                                   axis.line = ggplot2::element_line(colour = "black"),
                                                   text = ggplot2::element_text(size=12)) +
              ggplot2::ylab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) + ggplot2::xlab(xlab)
          }
          }

          if (isTRUE(setequal(which.y, c("outcome", "cost", "nmb")))) {
            p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
          else if (isTRUE("all" %in% which.y)) {
            p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome", "cost")))) {
            p = cowplot::plot_grid(yplot, cplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome", "nmb")))) {
            p = cowplot::plot_grid(yplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("cost", "nmb")))) {
            p = cowplot::plot_grid(cplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome")))) {
            p = yplot}
          else if (isTRUE(setequal(which.y, c("cost")))) {
            p = cplot}
          else if (isTRUE(setequal(which.y, c("nmb")))) {
            p = nplot}

        } else {#Else do not plot the smooth (non-parametric plot grf-style)

          Xmeans <- apply(Xmat, 2, mean)

          is.bin = apply(Xmat, 2, function(x) isTRUE(length(unique(x))==2)) #Check for binary covariates before fitting model
          is.factor = apply(Xmat, 2, function(x) isTRUE(class(x)=="factor" | class(x)=="character"))

          focal.x = colnames(Xmat)[which.x]
          avars = colnames(Xmat)[is.bin|is.factor]
          cvars = colnames(Xmat)[!(is.bin|is.factor)]

          if (isTRUE(focal.x %in% cvars)) {
            X.test = matrix(rep(Xmeans, length.out), length.out, ncol(Xmat), byrow=T)
            X.test[,which.x] = seq(x.range[1], x.range[2], length.out=length.out)
          } else {
            X.test = matrix(rep(Xmeans, length(unique(Xmat[,which.x]))), length(unique(Xmat[,which.x])), ncol(Xmat), byrow=T)
            X.test[,which.x] = unique(Xmat[,which.x])
          }
          x.pred = X.test[,which.x]
          preds_n = predict(forest, newdata=X.test)

          if (isTRUE(any(c("outcome") %in% which.y))) {

          dy = preds_n$predicted.delta_y
          vy = preds_n$variance.delta_y
          lower.y = dy-sqrt(vy)*qnorm(1-alpha/2)
          upper.y = dy+sqrt(vy)*qnorm(1-alpha/2)
          df.y = as.data.frame(cbind(dy, lower.y, upper.y, x.pred))
          if (isTRUE(focal.x %in% cvars)) {

          yplot = ggplot2::ggplot(data=df.y,ggplot2::aes(y=dy, x=x.pred)) +
            geom_line(size=1) + geom_ribbon(aes(ymin=lower.y, ymax=upper.y), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
            ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                 panel.grid.minor = ggplot2::element_blank(),
                                                 axis.line = ggplot2::element_line(colour = "black"),
                                                 text = ggplot2::element_text(size=12)) +
            ggplot2::ylab(expression(paste(Delta,"Y"))) + ggplot2::xlab(xlab)

             } else {

               yplot = ggplot2::ggplot(data=df.y,ggplot2::aes(y=dy, x=as.factor(x.pred))) +
                 geom_pointrange(aes(ymin=lower.y, ymax=upper.y), size=1) + geom_hline(yintercept=0, linetype="dashed") +
                 ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                      panel.grid.minor = ggplot2::element_blank(),
                                                      axis.line = ggplot2::element_line(colour = "black"),
                                                      text = ggplot2::element_text(size=12)) +
                 ggplot2::ylab(expression(paste(Delta,"Y"))) + ggplot2::xlab(xlab)

            }

          }

          if (isTRUE(any(c("cost") %in% which.y))) {

          dc = preds_n$predicted.delta_cost
          vc = preds_n$variance.delta_cost
          lower.c = dc-sqrt(vc)*qnorm(1-alpha/2)
          upper.c = dc+sqrt(vc)*qnorm(1-alpha/2)

          df.c = as.data.frame(cbind(dc, lower.c, upper.c, x.pred))

          if (isTRUE(focal.x %in% cvars)) {

          cplot = ggplot2::ggplot(data=df.c,ggplot2::aes(y=dc, x=x.pred)) +
            geom_line(size=1) + geom_ribbon(aes(ymin=lower.c, ymax=upper.c), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
            ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                 panel.grid.minor = ggplot2::element_blank(),
                                                 axis.line = ggplot2::element_line(colour = "black"),
                                                 text = ggplot2::element_text(size=12)) +
            ggplot2::ylab(expression(paste(Delta,"Cost"))) + ggplot2::xlab(xlab)

          } else {

            cplot = ggplot2::ggplot(data=df.c,ggplot2::aes(y=dc, x=as.factor(x.pred))) +
              geom_pointrange(aes(ymin=lower.c, ymax=upper.c), size=1) + geom_hline(yintercept=0, linetype="dashed") +
              ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                   panel.grid.minor = ggplot2::element_blank(),
                                                   axis.line = ggplot2::element_line(colour = "black"),
                                                   text = ggplot2::element_text(size=12)) +
              ggplot2::ylab(expression(paste(Delta,"Cost"))) + ggplot2::xlab(xlab)
          }


          }

          if (isTRUE(any(c("nmb") %in% which.y))) {
          dn = preds_n$predicted.nmb
          vn = preds_n$variance.nmb
          lower.n = dn-sqrt(vn)*qnorm(1-alpha/2)
          upper.n = dn+sqrt(vn)*qnorm(1-alpha/2)

          df.n = as.data.frame(cbind(dn, lower.n, upper.n, x.pred))

          if (isTRUE(focal.x %in% cvars)) {
          nplot = ggplot2::ggplot(data=df.n,ggplot2::aes(y=dn, x=x.pred)) +
            geom_line(size=1) + geom_ribbon(aes(ymin=lower.n, ymax=upper.n), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
            ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                 panel.grid.minor = ggplot2::element_blank(),
                                                 axis.line = ggplot2::element_line(colour = "black"),
                                                 text = ggplot2::element_text(size=12)) +
            ggplot2::ylab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) + ggplot2::xlab(xlab)

          } else {

            nplot = ggplot2::ggplot(data=df.n,ggplot2::aes(y=dn, x=as.factor(x.pred))) +
              geom_pointrange(aes(ymin=lower.n, ymax=upper.n), size=1) + geom_hline(yintercept=0, linetype="dashed") +
              ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                   panel.grid.minor = ggplot2::element_blank(),
                                                   axis.line = ggplot2::element_line(colour = "black"),
                                                   text = ggplot2::element_text(size=12)) +
              ggplot2::ylab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) + ggplot2::xlab(xlab)

            }
          }
          if (isTRUE(setequal(which.y, c("outcome", "cost", "nmb")))) {
            p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
          else if (isTRUE("all" %in% which.y)) {
            p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome", "cost")))) {
            p = cowplot::plot_grid(yplot, cplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome", "nmb")))) {
            p = cowplot::plot_grid(yplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("cost", "nmb")))) {
            p = cowplot::plot_grid(cplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome")))) {
            p = yplot}
          else if (isTRUE(setequal(which.y, c("cost")))) {
            p = cplot}
          else if (isTRUE(setequal(which.y, c("nmb")))) {
            p = nplot}
        }

      } else {#Else plot unconditional

        if (isTRUE(smooth)) {#Use mgcv to produce a spline plot with CIs

            if (isTRUE("mgcv" %in% rownames(installed.packages())==FALSE)) {
              stop("The mgcv package must be installed to plot smooth functions.")
            }

            xdf = as.data.frame(Xmat)
            xvars = colnames(Xmat)
            Xmeans <- apply(xdf, 2, mean)

            is.bin = apply(xdf, 2, function(x) isTRUE(length(unique(x))==2)) #Check for binary covariates before fitting model
            is.factor = apply(xdf, 2, function(x) isTRUE(class(x)=="factor" | class(x)=="character"))

            focal.x = colnames(xdf)[which.x]
            avars = colnames(xdf)[is.bin|is.factor]
            cvars = colnames(xdf)[!(is.bin|is.factor)]

            if (focal.x %in% cvars) {xb = paste("s(", focal.x,", k=-1)", sep="")} else {xb = focal.x}

            Xmeans <- apply(Xmat, 2, mean) # Not necessary here, can be removed. Leaving for now.
            if (isTRUE(focal.x %in% cvars)) {
              X.test = matrix(rep(Xmeans, length.out), length.out, ncol(Xmat), byrow=T)
              X.test[,which.x] = seq(x.range[1], x.range[2], length.out=length.out)
            } else {
              X.test = matrix(rep(Xmeans, length(unique(Xmat[,which.x]))), length(unique(Xmat[,which.x])), ncol(Xmat), byrow=T)
              X.test[,which.x] = unique(Xmat[,which.x])
            }
            colnames(X.test) = colnames(xdf)
            X.test = as.data.frame(X.test)
            x.pred = X.test[,which.x]


            #Outcomes

            if (isTRUE(any(c("outcome") %in% which.y))) {

            y.score = dr.taus$debiased.delta_outcome
            pdf = as.data.frame(cbind(y.score, xdf))
            form = as.formula(paste("y.score~",xb,sep=""))
            sgam = mgcv::gam(form, data=pdf, ...)

            gam.preds = predict(sgam, newdata=X.test, se.fit=TRUE)
            tau.fit = gam.preds$fit
            tau.se = gam.preds$se.fit
            lower.tau = tau.fit-tau.se*qnorm(1-alpha/2)
            upper.tau = tau.fit+tau.se*qnorm(1-alpha/2)
            df.y = as.data.frame(cbind(tau.fit, lower.tau, upper.tau, x.pred))

            if (isTRUE(focal.x %in% cvars)) {#If continuous, plot spline

              yplot = ggplot2::ggplot(data=df.y,ggplot2::aes(y=tau.fit, x=x.pred)) +
                geom_line(size=1) + geom_ribbon(aes(ymin=lower.tau, ymax=upper.tau), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12)) +
                ggplot2::ylab(expression(paste(Delta,"Y"))) + ggplot2::xlab(xlab)

            } else { #Else plot point with SE
              yplot = ggplot2::ggplot(data=df.y,ggplot2::aes(y=tau.fit, x=as.factor(x.pred))) +
                geom_pointrange(aes(ymin=lower.tau, ymax=upper.tau), size=1) + geom_hline(yintercept=0, linetype="dashed") +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12)) +
                ggplot2::ylab(expression(paste(Delta,"Y"))) + ggplot2::xlab(xlab)
            }
            }

            #Costs

            if (isTRUE(any(c("cost") %in% which.y))) {

            y.score = dr.taus$debiased.delta_cost
            pdf = as.data.frame(cbind(y.score, xdf))
            form = as.formula(paste("y.score~",xb,sep=""))
            sgam = mgcv::gam(form, data=pdf, ...)

            gam.preds = predict(sgam, newdata=X.test, se.fit=TRUE)
            tau.fit = gam.preds$fit
            tau.se = gam.preds$se.fit
            lower.tau = tau.fit-tau.se*qnorm(1-alpha/2)
            upper.tau = tau.fit+tau.se*qnorm(1-alpha/2)

            df.c = as.data.frame(cbind(tau.fit, lower.tau, upper.tau, x.pred))

            if (isTRUE(focal.x %in% cvars)) {#If continuous, plot spline

              cplot = ggplot2::ggplot(data=df.c,ggplot2::aes(y=tau.fit, x=x.pred)) +
                geom_line(size=1) + geom_ribbon(aes(ymin=lower.tau, ymax=upper.tau), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12)) +
                ggplot2::ylab(expression(paste(Delta,"Cost"))) + ggplot2::xlab(xlab)

            } else { #Else plot point with SE
              cplot = ggplot2::ggplot(data=df.c,ggplot2::aes(y=tau.fit, x=as.factor(x.pred))) +
                geom_pointrange(aes(ymin=lower.tau, ymax=upper.tau), size=1) + geom_hline(yintercept=0, linetype="dashed") +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12)) +
                ggplot2::ylab(expression(paste(Delta,"Cost"))) + ggplot2::xlab(xlab)
            }
            }

            #NMB

            if (isTRUE(any(c("nmb") %in% which.y))) {

            y.score = dr.taus$debiased.nmb
            pdf = as.data.frame(cbind(y.score, xdf))
            form = as.formula(paste("y.score~",xb,sep=""))
            sgam = mgcv::gam(form, data=pdf, ...)

            gam.preds = predict(sgam, newdata=X.test, se.fit=TRUE)
            tau.fit = gam.preds$fit
            tau.se = gam.preds$se.fit
            lower.tau = tau.fit-tau.se*qnorm(1-alpha/2)
            upper.tau = tau.fit+tau.se*qnorm(1-alpha/2)

            df.n = as.data.frame(cbind(tau.fit, lower.tau, upper.tau, x.pred))

            if (isTRUE(focal.x %in% cvars)) {#If continuous, plot spline

              nplot = ggplot2::ggplot(data=df.n,ggplot2::aes(y=tau.fit, x=x.pred)) +
                geom_line(size=1) + geom_ribbon(aes(ymin=lower.tau, ymax=upper.tau), alpha=0.3, fill="gray") + geom_hline(yintercept=0, linetype="dashed") +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12)) +
                ggplot2::ylab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) + ggplot2::xlab(xlab)

            } else { #Else plot point with SE
              nplot = ggplot2::ggplot(data=df.n,ggplot2::aes(y=tau.fit, x=as.factor(x.pred))) +
                geom_pointrange(aes(ymin=lower.tau, ymax=upper.tau), size=1) + geom_hline(yintercept=0, linetype="dashed") +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12)) +
                ggplot2::ylab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) + ggplot2::xlab(xlab)
            }
            }

            if (isTRUE(setequal(which.y, c("outcome", "cost", "nmb")))) {
              p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
            else if (isTRUE("all" %in% which.y)) {
              p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
            else if (isTRUE(setequal(which.y, c("outcome", "cost")))) {
              p = cowplot::plot_grid(yplot, cplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
            else if (isTRUE(setequal(which.y, c("outcome", "nmb")))) {
              p = cowplot::plot_grid(yplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
            else if (isTRUE(setequal(which.y, c("cost", "nmb")))) {
              p = cowplot::plot_grid(cplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
            else if (isTRUE(setequal(which.y, c("outcome")))) {
              p = yplot}
            else if (isTRUE(setequal(which.y, c("cost")))) {
              p = cplot}
            else if (isTRUE(setequal(which.y, c("nmb")))) {
              p = nplot}


        } else {#Else do not plot the smooth (scatter plot or box plot using out-of-bag estimates)

          x = Xmat[,which.x]

          is.bin = apply(Xmat, 2, function(x) isTRUE(length(unique(x))==2)) #Check for binary covariates before fitting model
          is.factor = apply(Xmat, 2, function(x) isTRUE(class(x)=="factor" | class(x)=="character"))

          focal.x = colnames(Xmat)[which.x]
          avars = colnames(Xmat)[is.bin|is.factor]
          cvars = colnames(Xmat)[!(is.bin|is.factor)]

          if (isTRUE(focal.x %in% cvars)) { # If continuous focal X, use scatter

          if (isTRUE(any(c("outcome") %in% which.y))) { # Something is wrong with x-axis...
          yplot = ggplot(preds, ggplot2::aes(y=predicted.delta_y, x=as.numeric(as.character(x)))) +
            ggplot2::geom_point(size=1) +
            ggplot2::ylab(expression(paste(Delta,"Y"))) + ggplot2::xlab(xlab) +
            ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                 panel.grid.minor = ggplot2::element_blank(),
                                                 axis.line = ggplot2::element_line(colour = "black"),
                                                 text = ggplot2::element_text(size=12)) + ggplot2::xlim(x.range)
          }
          if (isTRUE(any(c("cost") %in% which.y))) {
          cplot = ggplot(preds, ggplot2::aes(y=predicted.delta_cost, x=as.numeric(as.character(x)))) +
            ggplot2::geom_point(size=1) +
            ggplot2::ylab(expression(paste(Delta,"Cost"))) + ggplot2::xlab(xlab) +
            ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                 panel.grid.minor = ggplot2::element_blank(),
                                                 axis.line = ggplot2::element_line(colour = "black"),
                                                 text = ggplot2::element_text(size=12)) + ggplot2::xlim(x.range)
          }
          if (isTRUE(any(c("nmb") %in% which.y))) {
          nplot = ggplot(preds, ggplot2::aes(y=predicted.nmb, x=as.numeric(as.character(x)))) +
            ggplot2::geom_point(size=1) +
            ggplot2::ylab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) + ggplot2::xlab(xlab) +
            ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                 panel.grid.minor = ggplot2::element_blank(),
                                                 axis.line = ggplot2::element_line(colour = "black"),
                                                 text = ggplot2::element_text(size=12)) + ggplot2::xlim(x.range)
          }
          } else { # Else use box plots for binary/factors in x

            if (isTRUE(any(c("outcome") %in% which.y))) {
              yplot = ggplot(preds, ggplot2::aes(y=predicted.delta_y, x=as.factor(x))) +
                ggplot2::geom_boxplot() +
                ggplot2::ylab(expression(paste(Delta,"Y"))) + ggplot2::xlab(xlab) +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12))
            }
            if (isTRUE(any(c("cost") %in% which.y))) {
              cplot = ggplot(preds, ggplot2::aes(y=predicted.delta_cost, x=as.factor(x))) +
                ggplot2::geom_boxplot() +
                ggplot2::ylab(expression(paste(Delta,"Cost"))) + ggplot2::xlab(xlab) +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12))
            }
            if (isTRUE(any(c("nmb") %in% which.y))) {
              nplot = ggplot(preds, ggplot2::aes(y=predicted.nmb, x=as.factor(x))) +
                ggplot2::geom_boxplot() +
                ggplot2::ylab(paste("Net monetary benefit (WTP = ", forest[["WTP"]], ")", sep="")) + ggplot2::xlab(xlab) +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                     panel.grid.minor = ggplot2::element_blank(),
                                                     axis.line = ggplot2::element_line(colour = "black"),
                                                     text = ggplot2::element_text(size=12))
            }

          }

          if (isTRUE(setequal(which.y, c("outcome", "cost", "nmb")))) {
          p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
          else if (isTRUE("all" %in% which.y)) {
          p = cowplot::plot_grid(yplot, cplot, nplot, align="h", axis="b", nrow=1, ncol=3, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome", "cost")))) {
          p = cowplot::plot_grid(yplot, cplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome", "nmb")))) {
          p = cowplot::plot_grid(yplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("cost", "nmb")))) {
          p = cowplot::plot_grid(cplot, nplot, align="h", axis="b", nrow=1, ncol=2, labels=labels)}
          else if (isTRUE(setequal(which.y, c("outcome")))) {
          p = yplot}
          else if (isTRUE(setequal(which.y, c("cost")))) {
          p = cplot}
          else if (isTRUE(setequal(which.y, c("nmb")))) {
          p = nplot}

        }

      }

    }

  return(p)
}
